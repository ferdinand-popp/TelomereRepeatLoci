#!/usr/bin/env python3
"""Standalone diagnostic: check what pysam's bam.fetch() actually returns for a
given region, independent of find_fusion_reads.py or main.py. Used to isolate
whether a read-count gap seen in find_fusion_reads.py's output comes from
fetch() itself or from its per-read row-building logic.

Usage:
    uv run python scripts/diagnose_fetch.py /path/to/tumor.bam CHROM START END

CHROM/START/END are 0-based half-open, matching find_fusion_reads.py's own
window_start0/window_end0 (candidate region chromStart/chromEnd +/- 300bp).
"""

import argparse

import pysam


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam_path")
    parser.add_argument("chrom")
    parser.add_argument("start0", type=int)
    parser.add_argument("end0", type=int)
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam_path, "rb")

    total = 0
    unmapped = 0
    has_cigar = 0
    has_s_in_cigarstring = 0

    for read in bam.fetch(args.chrom, args.start0, args.end0):
        total += 1
        if read.is_unmapped:
            unmapped += 1
            continue
        cigar = read.cigarstring or ""
        if cigar:
            has_cigar += 1
        if "S" in cigar:
            has_s_in_cigarstring += 1

    bam.close()

    print(f"region: {args.chrom}:{args.start0}-{args.end0} (0-based half-open)")
    print(f"total reads from bam.fetch(): {total}")
    print(f"  unmapped: {unmapped}")
    print(f"  with a cigar string: {has_cigar}")
    print(f"  with 'S' in cigarstring: {has_s_in_cigarstring}")


if __name__ == "__main__":
    main()
