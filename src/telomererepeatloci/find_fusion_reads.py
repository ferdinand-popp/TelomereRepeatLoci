#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd
import pysam

from pipeline.tables import FUSION_READS_COLUMNS, read_tsv, sanitize_tsv_values, write_tsv


TELOMERE_PATTERN = re.compile(r"TTAGGG|CCCTAA")
READ_CONSUME_OPS = {0, 1, 4, 7, 8}
REF_CONSUME_OPS = {0, 2, 3, 7, 8}
WINDOW_EXTENSION = 300
# Number of rows to buffer before flushing to outfile. Candidate regions can
# span very wide/high-depth loci (e.g. after window fusion), and each row
# carries a full read sequence -- buffering the whole result set in memory
# for the entire candidate-region file scales with total reads fetched, not
# with output size, and can exhaust memory at high coverage.
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
    if not seq or not cigartuples:
        return []

    qpos = 0
    clips = []
    for op, length in cigartuples:
        if op in READ_CONSUME_OPS:
            if op == 4:  # soft clip
                clips.append(seq[qpos : qpos + length])
            qpos += length
    return [c for c in clips if c]


def expected_pos_fusion(cigar):
    if re.match(r"^\d+M.*\d+[HS]$", cigar):
        return "downstream"
    if re.match(r"^\d+[HS].*\d+M$", cigar):
        return "upstream"
    return ""


def get_primary_sequence(bam, sa_read, primary_chr, primary_pos, primary_strand):
    if not primary_chr or primary_pos <= 0:
        return sa_read.query_sequence or ""
    start0 = max(0, primary_pos - 1)
    # Query exactly the SA-tag primary position in 0-based half-open coordinates.
    end0 = primary_pos
    for read in bam.fetch(primary_chr, start0, end0):
        if read.query_name != sa_read.query_name:
            continue
        if read.is_supplementary or read.is_secondary:
            continue
        if sa_read.is_read1 != read.is_read1 or sa_read.is_read2 != read.is_read2:
            continue

        seq = read.query_sequence or ""
        if not seq:
            continue

        supp_strand = "-" if sa_read.is_reverse else "+"
        if primary_strand != supp_strand:
            seq = reverse_complement(seq)
        return seq
    return sa_read.query_sequence or ""


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


def _fusion_rows_for_region(bam, region):
    """Yield soft-clip and supplementary-alignment fusion-read rows for one
    candidate region, without accumulating them anywhere."""
    window = region.get("window", "")
    chrom = region.get("chrom", "")
    try:
        chrom_start = int(float(region.get("chromStart", 0)))
        chrom_end = int(float(region.get("chromEnd", 0)))
    except ValueError:
        return

    window_start0 = max(0, chrom_start - WINDOW_EXTENSION - 1)
    window_end0 = chrom_end + WINDOW_EXTENSION

    # soft-clipped reads
    for read in bam.fetch(chrom, window_start0, window_end0):
        if read.is_unmapped:
            continue
        cigar = read.cigarstring or ""
        if "S" not in cigar:
            continue

        start0 = read.reference_start
        end0 = alignment_end(start0, read.cigartuples)
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

    # supplementary alignments (hard-clipped candidates)
    for read in bam.fetch(chrom, window_start0, window_end0):
        if read.is_unmapped or not read.is_supplementary:
            continue

        try:
            sa_tag = read.get_tag("SA")
        except KeyError:
            continue
        primary_chr, primary_pos, primary_strand = parse_sa_tag(sa_tag)
        primary_pos0 = primary_pos - 1 if primary_pos else 0
        sequence = _strip_nuls(
            get_primary_sequence(bam, read, primary_chr, primary_pos, primary_strand)
        )

        cigar = read.cigarstring or ""
        start0 = read.reference_start
        end0 = alignment_end(start0, read.cigartuples)
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
    out_rows = []
    try:
        for region in candidate_regions:
            out_rows.extend(_fusion_rows_for_region(bam, region))
    finally:
        bam.close()

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
    buffer = []
    wrote_header = False
    try:
        for region in candidate_regions:
            buffer.extend(_fusion_rows_for_region(bam, region))
            if len(buffer) >= flush_rows:
                wrote_header = _flush_rows(buffer, outfile, wrote_header)
                buffer = []
        wrote_header = _flush_rows(buffer, outfile, wrote_header)
    finally:
        bam.close()

    if not wrote_header:
        write_tsv(pd.DataFrame(columns=FUSION_READS_COLUMNS), outfile, FUSION_READS_COLUMNS)


if __name__ == "__main__":
    main()
