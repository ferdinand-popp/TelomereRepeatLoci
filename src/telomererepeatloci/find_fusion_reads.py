#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd
import pysam

from pipeline.tables import (
    FUSION_READS_COLUMNS,
    read_tsv,
    sanitize_tsv_values,
    write_tsv,
)


TELOMERE_PATTERN = re.compile(r"TTAGGG|CCCTAA")
READ_CONSUME_OPS = {0, 1, 4, 7, 8}
REF_CONSUME_OPS = {0, 2, 3, 7, 8}
WINDOW_EXTENSION = 300
# Number of rows to buffer before flushing to outfile. Candidate regions can
# span very wide/high-depth loci (e.g. after window fusion), and each row
# carries a full read sequence -- buffering the whole result set in memory
# for the entire candidate-region file scales with total reads fetched, not
# with output size, and can exhaust memory at high coverage. The flush check
# must run per-row (not just between regions, e.g. via buffer.extend(generator))
# or a single high-depth/fused region can still dump its entire row set into
# `buffer` in one shot before the size check ever fires.
FLUSH_ROWS = 5000
# Max reads collected per SA-tag primary locus in _primary_reads_at(). A
# repeat-collapsed/high-depth site can have thousands of reads overlapping
# one exact position, but only a handful of specific (read_name, read1/2)
# combinations will ever actually be looked up there -- capping the scan
# bounds memory at exactly the loci that would otherwise blow it up. A read
# not found within the cap falls back to its own (possibly hard-clip-
# truncated) sequence, the same fallback already used for a true miss.
MAX_READS_PER_PRIMARY_LOCUS = 2000
# Max distinct (chrom, pos) primary loci kept in primary_seq_cache across a
# whole find_fusion_reads.py run. The cache never otherwise shrinks, so a
# candidate-region file touching many distinct high-depth loci could grow
# it unboundedly; FIFO eviction bounds total cache memory at the cost of
# occasionally re-fetching an evicted-then-revisited locus.
MAX_CACHED_PRIMARY_LOCI = 5000


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
    """Extract the soft/hard-clipped portions of `seq` per `cigartuples`.

    For a record's own query_sequence, hard-clipped (H) bases are absent from
    `seq` entirely, so H must not consume query positions. But for a
    supplementary alignment whose `sequence` has been substituted with the
    full primary read (see get_primary_sequence), the hard-clipped bases
    *are* present in `seq` -- H must then be treated like a soft clip (both
    consuming query space and itself extractable), matching R's
    cigarRangesAlongQuerySpace(ops=c("S","H"), before.hard.clipping=TRUE)
    against that same full sequence. Detect which case applies by comparing
    `seq`'s length to the cigar-implied length with and without H.
    """
    if not seq or not cigartuples:
        return []

    len_without_h = sum(length for op, length in cigartuples if op in READ_CONSUME_OPS)
    len_with_h = len_without_h + sum(
        length for op, length in cigartuples if op == 5
    )
    hard_clip_present_in_seq = len_with_h != len_without_h and len(seq) == len_with_h

    qpos = 0
    clips = []
    for op, length in cigartuples:
        if op in READ_CONSUME_OPS:
            if op == 4:  # soft clip
                clips.append(seq[qpos : qpos + length])
            qpos += length
        elif op == 5 and hard_clip_present_in_seq:  # hard clip, seq is unclipped
            clips.append(seq[qpos : qpos + length])
            qpos += length
    return [c for c in clips if c]


def expected_pos_fusion(cigar):
    if re.match(r"^\d+M.*\d+[HS]$", cigar):
        return "downstream"
    if re.match(r"^\d+[HS].*\d+M$", cigar):
        return "upstream"
    return ""


def _primary_reads_at(primary_bam, chrom, pos, cache):
    """Return {(read_name, is_read1, is_read2): sequence} for non-supplementary,
    non-secondary reads at a 1bp SA-tag locus, fetching it at most once per
    (chrom, pos) for the lifetime of `cache` instead of once per supplementary
    read -- fusion breakpoints often have many supplementary alignments
    pointing back at the same primary locus."""
    key = (chrom, pos)
    cached = cache.get(key)
    if cached is not None:
        return cached

    seqs = {}
    start0 = max(0, pos - 1)
    # Query exactly the SA-tag primary position in 0-based half-open coordinates.
    end0 = pos
    for read in primary_bam.fetch(chrom, start0, end0):
        if read.is_supplementary or read.is_secondary:
            continue
        seq = read.query_sequence or ""
        if not seq:
            continue
        seqs.setdefault((read.query_name, read.is_read1, read.is_read2), seq)
        if len(seqs) >= MAX_READS_PER_PRIMARY_LOCUS:
            break

    if len(cache) >= MAX_CACHED_PRIMARY_LOCI:
        cache.pop(next(iter(cache)))
    cache[key] = seqs
    return seqs


def get_primary_sequence(
    primary_bam, sa_read, primary_chr, primary_pos, primary_strand, primary_seq_cache
):
    if not primary_chr or primary_pos <= 0:
        return sa_read.query_sequence or ""

    seqs = _primary_reads_at(primary_bam, primary_chr, primary_pos, primary_seq_cache)
    seq = seqs.get((sa_read.query_name, sa_read.is_read1, sa_read.is_read2))
    if not seq:
        return sa_read.query_sequence or ""

    supp_strand = "-" if sa_read.is_reverse else "+"
    if primary_strand != supp_strand:
        seq = reverse_complement(seq)
    return seq


def parse_sa_tag(sa_tag):
    # chr,pos,strand,cigar,mapq,nm;...
    first = sa_tag.split(";")[0]
    fields = first.split(",")
    if len(fields) < 3:
        return "", 0, "+"
    chrom = fields[0]
    try:
        pos = int(fields[1])
    except ValueError:
        pos = 0
    strand = fields[2]
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
    for value in row.values():
        if _has_bad_encoding(value):
            return True
    return False


def _strip_nuls(value: str) -> str:
    if not value:
        return value
    if "\x00" in value:
        return value.replace("\x00", "")
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_region_file")
    parser.add_argument("bamfile")
    parser.add_argument("outfile")
    args = parser.parse_args()

    write_fusion_reads_streaming(args.candidate_region_file, args.bamfile, args.outfile)


def _fusion_rows_for_region(bam, region, primary_seq_cache, primary_bam):
    """Yield soft-clip and supplementary-alignment fusion-read rows for one
    candidate region, without accumulating them anywhere.

    Both row kinds are derived from the same window fetch (a read can produce
    either, or both, e.g. a supplementary alignment that also has a soft clip)
    -- iterating the window once instead of twice halves the read I/O/parsing
    cost per region.

    `primary_bam` must be a *different* AlignmentFile handle than `bam`, used
    only for the supplementary branch's primary-sequence lookup
    (get_primary_sequence). pysam does not support two concurrent fetch()
    iterators on the same handle: a nested bam.fetch() call on `bam` while
    the outer `for read in bam.fetch(...)` below is still iterating silently
    corrupts/resets that outer iterator -- this was happening on every
    supplementary read and could truncate the rest of the region's scan
    (observed dropping a region's yield from ~9500 candidate rows to 9 at a
    high-depth, supplementary-alignment-rich locus).
    """
    window = region.get("window", "")
    chrom = region.get("chrom", "")
    try:
        chrom_start = int(float(region.get("chromStart", 0)))
        chrom_end = int(float(region.get("chromEnd", 0)))
    except ValueError:
        return

    window_start0 = max(0, chrom_start - WINDOW_EXTENSION - 1)
    window_end0 = chrom_end + WINDOW_EXTENSION

    for read in bam.fetch(chrom, window_start0, window_end0):
        if read.is_unmapped:
            continue

        cigar = read.cigarstring or ""
        start0 = read.reference_start
        end0 = alignment_end(start0, read.cigartuples)

        # soft-clipped read
        if "S" in cigar:
            sequence = _strip_nuls(read.query_sequence or "")
            clipped_parts = clipped_sequences_from_cigar(sequence, read.cigartuples)
            clipped_sequence = _strip_nuls(", ".join(clipped_parts))
            part_telomere = bool(TELOMERE_PATTERN.search(clipped_sequence))
            t_count, c_count = telomere_counts(clipped_sequence)
            row = {
                "window": window,
                "read_name": read.query_name,
                "read_1_2": read_pair_label(read),
                "start": start0,
                "end": end0,
                "cigar": cigar,
                "chr_primary_align": "",
                "coord_primary_align": "",
                "strand_primary_align": "",
                "sequence": sequence,
                "clipped_sequence": clipped_sequence,
                "part_telomere": str(part_telomere),
                "TTAGGG_count": t_count,
                "CCCTAA_count": c_count,
                "expected_pos_fusion": expected_pos_fusion(cigar),
            }
            if not _row_has_bad_encoding(row):
                yield row

        # supplementary alignment (hard-clipped candidate)
        if read.is_supplementary:
            try:
                sa_tag = read.get_tag("SA")
            except KeyError:
                continue
            primary_chr, primary_pos, primary_strand = parse_sa_tag(sa_tag)
            primary_pos0 = primary_pos - 1 if primary_pos else 0
            sequence = _strip_nuls(
                get_primary_sequence(
                    primary_bam,
                    read,
                    primary_chr,
                    primary_pos,
                    primary_strand,
                    primary_seq_cache,
                )
            )

            clipped_parts = clipped_sequences_from_cigar(sequence, read.cigartuples)
            clipped_sequence = _strip_nuls(", ".join(clipped_parts))
            part_telomere = bool(TELOMERE_PATTERN.search(clipped_sequence))
            t_count, c_count = telomere_counts(clipped_sequence)
            row = {
                "window": window,
                "read_name": read.query_name,
                "read_1_2": read_pair_label(read),
                "start": start0,
                "end": end0,
                "cigar": cigar,
                "chr_primary_align": primary_chr,
                "coord_primary_align": primary_pos0,
                "strand_primary_align": primary_strand,
                "sequence": sequence,
                "clipped_sequence": clipped_sequence,
                "part_telomere": str(part_telomere),
                "TTAGGG_count": t_count,
                "CCCTAA_count": c_count,
                "expected_pos_fusion": expected_pos_fusion(cigar),
            }
            if not _row_has_bad_encoding(row):
                yield row


def find_fusion_reads(candidate_region_file: str, bamfile: str) -> pd.DataFrame:
    """Return every fusion-read row as one in-memory DataFrame.

    Kept for callers that want the full result set directly (e.g. tests);
    the CLI entry point uses write_fusion_reads_streaming() instead, which
    never holds the whole result set in memory at once.
    """
    candidate_regions = read_tsv(candidate_region_file).to_dict("records")

    bam = pysam.AlignmentFile(bamfile, "rb")
    # Separate handle for primary-sequence lookups -- see _fusion_rows_for_region's
    # docstring on why this can't share `bam`'s fetch() iterator.
    primary_bam = pysam.AlignmentFile(bamfile, "rb")
    primary_seq_cache = {}
    out_rows = []
    try:
        for region in candidate_regions:
            out_rows.extend(
                _fusion_rows_for_region(bam, region, primary_seq_cache, primary_bam)
            )
    finally:
        bam.close()
        primary_bam.close()

    return pd.DataFrame(out_rows)


def _flush_rows(rows, outfile, wrote_header):
    if not rows:
        return wrote_header
    df = sanitize_tsv_values(pd.DataFrame(rows, columns=FUSION_READS_COLUMNS))
    df.to_csv(
        outfile,
        sep="\t",
        index=False,
        encoding="utf-8",
        mode="a" if wrote_header else "w",
        header=not wrote_header,
    )
    return True


def write_fusion_reads_streaming(
    candidate_region_file: str,
    bamfile: str,
    outfile: str,
    flush_rows: int = FLUSH_ROWS,
) -> None:
    """Write fusion-read rows straight to outfile in bounded-size batches
    instead of buffering the whole candidate-region file's worth of reads
    (each carrying a full sequence) in memory before a single write."""
    candidate_regions = read_tsv(candidate_region_file).to_dict("records")

    out_path = Path(outfile)
    if out_path.exists():
        out_path.unlink()

    bam = pysam.AlignmentFile(bamfile, "rb")
    # Separate handle for primary-sequence lookups -- see _fusion_rows_for_region's
    # docstring on why this can't share `bam`'s fetch() iterator.
    primary_bam = pysam.AlignmentFile(bamfile, "rb")
    primary_seq_cache = {}
    buffer = []
    wrote_header = False
    try:
        for region in candidate_regions:
            for row in _fusion_rows_for_region(
                bam, region, primary_seq_cache, primary_bam
            ):
                buffer.append(row)
                if len(buffer) >= flush_rows:
                    wrote_header = _flush_rows(buffer, outfile, wrote_header)
                    buffer = []
        wrote_header = _flush_rows(buffer, outfile, wrote_header)
    finally:
        bam.close()
        primary_bam.close()

    if not wrote_header:
        write_tsv(
            pd.DataFrame(columns=FUSION_READS_COLUMNS), outfile, FUSION_READS_COLUMNS
        )


if __name__ == "__main__":
    main()
