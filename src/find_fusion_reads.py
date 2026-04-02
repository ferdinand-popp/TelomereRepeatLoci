#!/usr/bin/env python3

import argparse
import csv
import re

import pysam


TELOMERE_PATTERN = re.compile(r"TTAGGG|CCCTAA")
READ_CONSUME_OPS = {0, 1, 4, 7, 8}
REF_CONSUME_OPS = {0, 2, 3, 7, 8}
WINDOW_EXTENSION = 300


def read_pair_label(read):
    if read.is_read1:
        return "READ1"
    if read.is_read2:
        return "READ2"
    return ""


def alignment_end(start_1based, cigartuples):
    """Return 1-based inclusive reference end for an alignment."""
    if not cigartuples:
        return start_1based
    ref_len = sum(length for op, length in cigartuples if op in REF_CONSUME_OPS)
    return start_1based + ref_len - 1


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_region_file")
    parser.add_argument("bamfile")
    parser.add_argument("outfile")
    args = parser.parse_args()

    with open(args.candidate_region_file, newline="") as handle:
        candidate_regions = list(csv.DictReader(handle, delimiter="\t"))

    bam = pysam.AlignmentFile(args.bamfile, "rb")
    out_rows = []

    for region in candidate_regions:
        window = region.get("window", "")
        chrom = region.get("chrom", "")
        try:
            chrom_start = int(float(region.get("chromStart", 0)))
            chrom_end = int(float(region.get("chromEnd", 0)))
        except ValueError:
            continue

        start0 = max(0, chrom_start - WINDOW_EXTENSION - 1)
        end0 = chrom_end + WINDOW_EXTENSION

        # soft-clipped reads
        for read in bam.fetch(chrom, start0, end0):
            if read.is_unmapped:
                continue
            cigar = read.cigarstring or ""
            if "S" not in cigar:
                continue

            start_1based = read.reference_start + 1
            end_1based = alignment_end(start_1based, read.cigartuples)
            sequence = read.query_sequence or ""
            clipped_parts = clipped_sequences_from_cigar(sequence, read.cigartuples)
            clipped_sequence = ", ".join(clipped_parts)
            part_telomere = bool(TELOMERE_PATTERN.search(clipped_sequence))
            t_count, c_count = telomere_counts(clipped_sequence)

            out_rows.append(
                {
                    "window": window,
                    "read_name": read.query_name,
                    "read_1_2": read_pair_label(read),
                    "start": start_1based,
                    "end": end_1based,
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
            )

        # supplementary alignments (hard-clipped candidates)
        for read in bam.fetch(chrom, start0, end0):
            if read.is_unmapped or not read.is_supplementary:
                continue

            try:
                sa_tag = read.get_tag("SA")
            except KeyError:
                continue
            primary_chr, primary_pos, primary_strand = parse_sa_tag(sa_tag)
            sequence = get_primary_sequence(
                bam, read, primary_chr, primary_pos, primary_strand
            )

            cigar = read.cigarstring or ""
            start_1based = read.reference_start + 1
            end_1based = alignment_end(start_1based, read.cigartuples)
            clipped_parts = clipped_sequences_from_cigar(sequence, read.cigartuples)
            clipped_sequence = ", ".join(clipped_parts)
            part_telomere = bool(TELOMERE_PATTERN.search(clipped_sequence))
            t_count, c_count = telomere_counts(clipped_sequence)

            out_rows.append(
                {
                    "window": window,
                    "read_name": read.query_name,
                    "read_1_2": read_pair_label(read),
                    "start": start_1based,
                    "end": end_1based,
                    "cigar": cigar,
                    "chr_primary_align": primary_chr,
                    "coord_primary_align": primary_pos,
                    "strand_primary_align": primary_strand,
                    "sequence": sequence,
                    "clipped_sequence": clipped_sequence,
                    "part_telomere": str(part_telomere),
                    "TTAGGG_count": t_count,
                    "CCCTAA_count": c_count,
                    "expected_pos_fusion": expected_pos_fusion(cigar),
                }
            )

    bam.close()

    fieldnames = [
        "window",
        "read_name",
        "read_1_2",
        "start",
        "end",
        "cigar",
        "chr_primary_align",
        "coord_primary_align",
        "strand_primary_align",
        "sequence",
        "clipped_sequence",
        "part_telomere",
        "TTAGGG_count",
        "CCCTAA_count",
        "expected_pos_fusion",
    ]

    with open(args.outfile, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)


if __name__ == "__main__":
    main()
