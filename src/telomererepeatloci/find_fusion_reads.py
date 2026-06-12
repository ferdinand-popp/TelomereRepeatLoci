#!/usr/bin/env python3

import argparse
import re

import pandas as pd
import pysam

from pipeline.tables import FUSION_READS_COLUMNS, read_tsv, write_tsv


TELOMERE_PATTERN = re.compile(r"TTAGGG|CCCTAA")

# CIGAR op codes
# 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
READ_CONSUME_OPS = {0, 1, 4, 7, 8}  # ops that advance the query position
REF_CONSUME_OPS = {0, 2, 3, 7, 8}  # ops that advance the reference position

WINDOW_EXTENSION = 300


def read_pair_label(read):
    if read.is_read1:
        return "READ1"
    if read.is_read2:
        return "READ2"
    return ""


def alignment_end(start0, cigartuples):
    """Return 0-based exclusive reference end for an alignment."""
    if not cigartuples:
        return start0
    ref_len = sum(length for op, length in cigartuples if op in REF_CONSUME_OPS)
    return start0 + ref_len


def reverse_complement(seq):
    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(trans)[::-1]


def clipped_sequences_from_cigar(seq, cigartuples):
    if not seq or not cigartuples:
        return []

    clips = []
    qpos = 0
    for op, length in cigartuples:
        if op == 4:  # soft clip — bases ARE in query_sequence
            clips.append(seq[qpos : qpos + length])
            qpos += length
        elif op == 5:  # hard clip — bases NOT in query_sequence; skip
            pass
        elif op in READ_CONSUME_OPS:
            qpos += length
    return [c for c in clips if c]


def expected_pos_fusion(cigar):
    if re.match(r"^\d+M.*\d+[HS]$", cigar):
        return "downstream"
    if re.match(r"^\d+[HS].*\d+M$", cigar):
        return "upstream"
    return ""


def parse_sa_tag(sa_tag: str) -> tuple[str, int, str]:
    """
    Parse the first entry of an SA tag.
    Format: chr,pos,strand,cigar,mapq,nm;...
    Returns (chrom, 1-based pos, strand).
    """
    first = sa_tag.split(";")[0]
    fields = first.split(",")
    if len(fields) < 3:
        return "", 0, "+"
    chrom = fields[0]
    try:
        pos = int(fields[1])
    except ValueError:
        pos = 0
    strand = fields[2] if fields[2] in ("+", "-") else "+"
    return chrom, pos, strand


def telomere_counts(seq):
    t = len(re.findall("TTAGGG", seq))
    c = len(re.findall("CCCTAA", seq))
    return t, c


def _has_bad_encoding(value) -> bool:
    if value is None:
        return False
    try:
        text = str(value)
    except Exception:
        return True
    if "\x00" in text:
        return True
    try:
        text.encode("utf-8")
    except UnicodeEncodeError:
        return True
    return False


def _row_has_bad_encoding(row: dict) -> bool:
    return any(_has_bad_encoding(v) for v in row.values())


def _strip_nuls(value: str) -> str:
    if value and "\x00" in value:
        return value.replace("\x00", "")
    return value


# ---------------------------------------------------------------------------
# Primary-sequence retrieval (for supplementary alignments)
# ---------------------------------------------------------------------------


def build_primary_cache(
    bam: pysam.AlignmentFile,
    chrom: str,
    start0: int,
    end0: int,
) -> dict[tuple, pysam.AlignedSegment]:
    """
    Index all non-supplementary, non-secondary reads overlapping a window
    by (query_name, is_read1).  Used to resolve primary sequences for
    supplementary alignments without a fetch-per-read.
    """
    cache: dict[tuple, pysam.AlignedSegment] = {}
    for read in bam.fetch(chrom, start0, end0):
        if read.is_unmapped or read.is_supplementary or read.is_secondary:
            continue
        key = (read.query_name, read.is_read1)
        if key not in cache:
            cache[key] = read
    return cache


def get_primary_sequence(
    bam: pysam.AlignmentFile,
    sa_read: pysam.AlignedSegment,
    primary_chr: str,
    primary_pos: int,  # 1-based, from SA tag
    primary_strand: str,
    primary_cache: dict | None = None,
) -> str:
    """
    Return the full query sequence of the primary alignment that corresponds
    to *sa_read*.

    Strategy (mirrors the R original's samtools grep logic):
      1. Try the pre-built cache for the primary region if provided.
      2. Fall back to a targeted bam.fetch() around primary_pos.
      3. Last resort: return the supplementary read's own stored sequence.

    Strand correction: if the primary is on the opposite strand from the
    supplementary, reverse-complement the sequence (matching R's
    reverseComplement call).
    """
    supp_strand = "-" if sa_read.is_reverse else "+"

    def _strand_correct(seq: str) -> str:
        if primary_strand != supp_strand:
            return reverse_complement(seq)
        return seq

    # 1. Try cache (built from the *primary* chromosome window)
    if primary_cache is not None:
        key = (sa_read.query_name, sa_read.is_read1)
        cached = primary_cache.get(key)
        if cached and cached.query_sequence:
            return _strand_correct(cached.query_sequence)

    # 2. Targeted fetch around the SA-tag coordinate
    if primary_chr and primary_pos > 0:
        fetch_start = max(0, primary_pos - 1)  # 0-based inclusive
        fetch_end = primary_pos  # 0-based exclusive (1 bp window)
        for read in bam.fetch(primary_chr, fetch_start, fetch_end):
            if read.query_name != sa_read.query_name:
                continue
            if read.is_supplementary or read.is_secondary:
                continue
            if sa_read.is_read1 != read.is_read1 or sa_read.is_read2 != read.is_read2:
                continue
            # Verify the read actually starts at primary_pos (matches R's awk '$4 == coord')
            if read.reference_start != fetch_start:
                continue
            seq = read.query_sequence or ""
            if seq:
                return _strand_correct(seq)

    # 3. Fallback
    return _strip_nuls(sa_read.query_sequence or "")


# ---------------------------------------------------------------------------
# Core row builder
# ---------------------------------------------------------------------------


def _build_row(
    window: str,
    read: pysam.AlignedSegment,
    sequence: str,
    chr_primary: str,
    coord_primary_1based: int,
    strand_primary: str,
) -> dict:
    """Build one output row from a processed read."""
    cigar = read.cigarstring or ""
    read_start = read.reference_start  # 0-based
    read_end = alignment_end(read_start, read.cigartuples)  # 0-based exclusive

    sequence = _strip_nuls(sequence)
    clipped_parts = clipped_sequences_from_cigar(sequence, read.cigartuples)
    clipped_sequence = _strip_nuls(", ".join(clipped_parts))

    part_telomere = bool(TELOMERE_PATTERN.search(clipped_sequence))
    t_count, c_count = telomere_counts(clipped_sequence)

    # Store coord_primary_align as 1-based to match R original
    return {
        "window": window,
        "read_name": read.query_name,
        "read_1_2": read_pair_label(read),
        "start": read_start,
        "end": read_end,
        "cigar": cigar,
        "chr_primary_align": chr_primary,
        "coord_primary_align": coord_primary_1based if coord_primary_1based else "",
        "strand_primary_align": strand_primary,
        "sequence": sequence,
        "clipped_sequence": clipped_sequence,
        "part_telomere": str(part_telomere),
        "TTAGGG_count": t_count,
        "CCCTAA_count": c_count,
        "expected_pos_fusion": expected_pos_fusion(cigar),
    }


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------


def find_fusion_reads(candidate_region_file: str, bamfile: str) -> pd.DataFrame:
    candidate_regions = read_tsv(candidate_region_file).to_dict("records")
    bam = pysam.AlignmentFile(bamfile, "rb")
    out_rows: list[dict] = []

    for region in candidate_regions:
        window = region.get("window", "")
        chrom = region.get("chrom", "")
        try:
            chrom_start = int(float(region.get("chromStart", 0)))
            chrom_end = int(float(region.get("chromEnd", 0)))
        except ValueError:
            continue

        # Use distinct names so they are never overwritten inside the loop
        win_start = max(0, chrom_start - WINDOW_EXTENSION - 1)
        win_end = chrom_end + WINDOW_EXTENSION

        # ------------------------------------------------------------------
        # Single pass: collect soft-clipped reads and supplementary
        # alignments together to avoid fetching the BAM region twice.
        # ------------------------------------------------------------------
        # We also cache primary reads encountered in this window so that
        # supplementary alignments whose primary falls in the same region
        # don't need an extra fetch.
        window_primary_cache: dict[tuple, pysam.AlignedSegment] = {}

        soft_clipped_reads: list[pysam.AlignedSegment] = []
        supplementary_reads: list[pysam.AlignedSegment] = []

        for read in bam.fetch(chrom, win_start, win_end):
            if read.is_unmapped:
                continue

            if not read.is_supplementary and not read.is_secondary:
                key = (read.query_name, read.is_read1)
                if key not in window_primary_cache:
                    window_primary_cache[key] = read

            cigar = read.cigarstring or ""

            if not read.is_supplementary and "S" in cigar:
                soft_clipped_reads.append(read)
            elif read.is_supplementary:
                supplementary_reads.append(read)

        # ------------------------------------------------------------------
        # Process soft-clipped reads
        # ------------------------------------------------------------------
        for read in soft_clipped_reads:
            row = _build_row(
                window=window,
                read=read,
                sequence=read.query_sequence or "",
                chr_primary="",
                coord_primary_1based=0,
                strand_primary="",
            )
            if not _row_has_bad_encoding(row):
                out_rows.append(row)

        # ------------------------------------------------------------------
        # Process supplementary alignments
        # ------------------------------------------------------------------
        for read in supplementary_reads:
            try:
                sa_tag = read.get_tag("SA")
            except KeyError:
                continue

            primary_chr, primary_pos, primary_strand = parse_sa_tag(sa_tag)

            # Build a targeted primary cache for reads on a *different*
            # chromosome (not covered by window_primary_cache).
            if primary_chr and primary_chr != chrom and primary_pos > 0:
                remote_cache = build_primary_cache(
                    bam,
                    primary_chr,
                    max(0, primary_pos - 1),
                    primary_pos,
                )
            else:
                remote_cache = window_primary_cache

            sequence = get_primary_sequence(
                bam,
                read,
                primary_chr,
                primary_pos,
                primary_strand,
                primary_cache=remote_cache,
            )

            row = _build_row(
                window=window,
                read=read,
                sequence=sequence,
                chr_primary=primary_chr,
                coord_primary_1based=primary_pos,  # keep 1-based as in R
                strand_primary=primary_strand,
            )
            if not _row_has_bad_encoding(row):
                out_rows.append(row)

    bam.close()
    return pd.DataFrame(out_rows)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract fusion reads from a BAM file around candidate regions."
    )
    parser.add_argument("candidate_region_file")
    parser.add_argument("bamfile")
    parser.add_argument("outfile")
    args = parser.parse_args()

    df = find_fusion_reads(args.candidate_region_file, args.bamfile)
    write_tsv(df, args.outfile, FUSION_READS_COLUMNS)


if __name__ == "__main__":
    main()
