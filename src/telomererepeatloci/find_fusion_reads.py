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
    full primary read (see _supp_sequence), the hard-clipped bases
    *are* present in `seq` -- H must then be treated like a soft clip (both
    consuming query space and itself extractable), matching R's
    cigarRangesAlongQuerySpace(ops=c("S","H"), before.hard.clipping=TRUE)
    against that same full sequence. Detect which case applies by comparing
    `seq`'s length to the cigar-implied length with and without H.
    """
    if not seq or not cigartuples:
        return []

    len_without_h = sum(length for op, length in cigartuples if op in READ_CONSUME_OPS)
    len_with_h = len_without_h + sum(length for op, length in cigartuples if op == 5)
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


def _resolve_primary_sequences(primary_bam, pending_by_locus):
    """For each (chrom, pos) SA-tag primary locus referenced by this region's
    supplementary reads, fetch that locus and capture only the sequences of
    the specific reads pending there (from pending_by_locus), stopping as
    soon as every one of them has been found.

    Unlike scanning a locus generically and capping how much of its pileup
    to keep, the wanted set here is known up front and is bounded by how
    many supplementary reads in this region actually reference that locus
    (typically a handful) -- not by the locus's total read depth. So a read
    is only ever missing from the result because it genuinely isn't at that
    locus, never because an arbitrary cap was hit first. Early-exiting once
    every wanted key is found keeps this cheap on the common case without
    needing any cap at all.
    """
    resolved = {}
    for (chrom, pos), wanted_keys in pending_by_locus.items():
        found = {}
        start0 = max(0, pos - 1)
        # Query exactly the SA-tag primary position in 0-based half-open coordinates.
        end0 = pos
        for read in primary_bam.fetch(chrom, start0, end0):
            if read.is_supplementary or read.is_secondary:
                continue
            key = (read.query_name, read.is_read1, read.is_read2)
            if key in wanted_keys and key not in found:
                seq = read.query_sequence or ""
                if seq:
                    found[key] = seq
                    if len(found) >= len(wanted_keys):
                        break
        resolved[(chrom, pos)] = found
    return resolved


def _supp_sequence(meta, resolved):
    """Return the correctly-stranded primary sequence for one supplementary
    read's metadata (see _fusion_rows_for_region), falling back to the
    read's own (possibly hard-clip-truncated) sequence only when the primary
    read genuinely isn't found at its SA-tag locus."""
    locus = meta["primary_locus"]
    seq = resolved.get(locus, {}).get(meta["key"]) if locus is not None else None
    if not seq:
        return meta["own_sequence"]

    supp_strand = "-" if meta["is_reverse"] else "+"
    if meta["strand_primary_align"] != supp_strand:
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


def _fusion_rows_for_region(bam, region, primary_bam):
    """Yield soft-clip and supplementary-alignment fusion-read rows for one
    candidate region.

    Soft-clip rows are self-contained and yielded immediately while scanning
    the region's window (a read can also be a supplementary alignment at the
    same time, e.g. one with its own soft clip -- both rows get produced from
    this same single pass). Supplementary (hard-clip) rows additionally need
    their SA-tag primary read's full sequence, which requires a second,
    targeted fetch of that primary locus; those are collected as
    `pending_by_locus` while scanning and resolved together, once, after the
    window scan finishes, by _resolve_primary_sequences() -- so multiple
    supplementary reads sharing one primary locus still cost only one fetch
    of it, and a read is never dropped just because that locus happens to be
    deep (see _resolve_primary_sequences's docstring).

    `primary_bam` must be a *different* AlignmentFile handle than `bam`.
    pysam does not support two concurrent fetch() iterators on the same
    handle: a nested bam.fetch() call on `bam` while the outer
    `for read in bam.fetch(...)` below is still iterating silently
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

    pending_by_locus = {}
    supp_meta_list = []

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

        # supplementary alignment (hard-clipped candidate) -- deferred until
        # the primary sequence has been resolved below.
        if read.is_supplementary:
            try:
                sa_tag = read.get_tag("SA")
            except KeyError:
                continue
            primary_chr, primary_pos, primary_strand = parse_sa_tag(sa_tag)
            primary_pos0 = primary_pos - 1 if primary_pos else 0
            key = (read.query_name, read.is_read1, read.is_read2)
            locus = (
                (primary_chr, primary_pos) if primary_chr and primary_pos > 0 else None
            )
            supp_meta_list.append(
                {
                    "window": window,
                    "read_name": read.query_name,
                    "read_1_2": read_pair_label(read),
                    "start": start0,
                    "end": end0,
                    "cigar": cigar,
                    "chr_primary_align": primary_chr,
                    "coord_primary_align": primary_pos0,
                    "strand_primary_align": primary_strand,
                    "expected_pos_fusion": expected_pos_fusion(cigar),
                    "own_sequence": read.query_sequence or "",
                    "is_reverse": read.is_reverse,
                    "cigartuples": read.cigartuples,
                    "key": key,
                    "primary_locus": locus,
                }
            )
            if locus is not None:
                pending_by_locus.setdefault(locus, set()).add(key)

    resolved = _resolve_primary_sequences(primary_bam, pending_by_locus)

    for meta in supp_meta_list:
        sequence = _strip_nuls(_supp_sequence(meta, resolved))
        clipped_parts = clipped_sequences_from_cigar(sequence, meta["cigartuples"])
        clipped_sequence = _strip_nuls(", ".join(clipped_parts))
        part_telomere = bool(TELOMERE_PATTERN.search(clipped_sequence))
        t_count, c_count = telomere_counts(clipped_sequence)
        row = {
            "window": meta["window"],
            "read_name": meta["read_name"],
            "read_1_2": meta["read_1_2"],
            "start": meta["start"],
            "end": meta["end"],
            "cigar": meta["cigar"],
            "chr_primary_align": meta["chr_primary_align"],
            "coord_primary_align": meta["coord_primary_align"],
            "strand_primary_align": meta["strand_primary_align"],
            "sequence": sequence,
            "clipped_sequence": clipped_sequence,
            "part_telomere": str(part_telomere),
            "TTAGGG_count": t_count,
            "CCCTAA_count": c_count,
            "expected_pos_fusion": meta["expected_pos_fusion"],
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
    out_rows = []
    try:
        for region in candidate_regions:
            out_rows.extend(_fusion_rows_for_region(bam, region, primary_bam))
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
    buffer = []
    wrote_header = False
    try:
        for region in candidate_regions:
            for row in _fusion_rows_for_region(bam, region, primary_bam):
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
