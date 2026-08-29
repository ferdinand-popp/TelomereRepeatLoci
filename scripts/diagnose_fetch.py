#!/usr/bin/env python3
"""Standalone diagnostic: check what pysam's bam.fetch() actually returns for a
given region, AND what find_fusion_reads.py's real _fusion_rows_for_region()
yields for the same region -- independent of main.py, candidate-region-file
parsing, or the output-writing path. Used to isolate whether a read-count gap
seen in find_fusion_reads.py's output comes from fetch() itself or from its
per-read row-building logic.

Usage:
    uv run python scripts/diagnose_fetch.py /path/to/tumor.bam CHROM CHROM_START0 CHROM_END0

CHROM_START0/CHROM_END0 are the *candidate region's own* chromStart/chromEnd
(0-based half-open), matching a row from the candidate-region table -- this
script applies the same +/-300bp WINDOW_EXTENSION internally that
find_fusion_reads.py does, so pass the unexpanded region boundaries.
"""

import argparse
import inspect

import pysam

import telomererepeatloci.find_fusion_reads as ffr
from telomererepeatloci.find_fusion_reads import (
    _fusion_rows_for_region,
    _row_has_bad_encoding,
    _strip_nuls,
    alignment_end,
    clipped_sequences_from_cigar,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam_path")
    parser.add_argument("chrom")
    parser.add_argument("chrom_start0", type=int, help="candidate region chromStart (0-based)")
    parser.add_argument("chrom_end0", type=int, help="candidate region chromEnd (0-based, exclusive)")
    args = parser.parse_args()

    print(f"find_fusion_reads module loaded from: {inspect.getsourcefile(ffr)}")
    print(f"WINDOW_EXTENSION = {ffr.WINDOW_EXTENSION}\n")

    window_start0 = max(0, args.chrom_start0 - 300 - 1)
    window_end0 = args.chrom_end0 + 300

    bam = pysam.AlignmentFile(args.bam_path, "rb")

    total = 0
    unmapped = 0
    has_s_in_cigarstring = 0
    is_supplementary = 0
    has_sa_tag = 0

    for read in bam.fetch(args.chrom, window_start0, window_end0):
        total += 1
        if read.is_unmapped:
            unmapped += 1
            continue
        cigar = read.cigarstring or ""
        if "S" in cigar:
            has_s_in_cigarstring += 1
        if read.is_supplementary:
            is_supplementary += 1
            if read.has_tag("SA"):
                has_sa_tag += 1

    bam.close()

    print(f"region searched: {args.chrom}:{window_start0}-{window_end0} (0-based half-open)")
    print("raw bam.fetch() totals:")
    print(f"  total reads: {total}")
    print(f"  unmapped: {unmapped}")
    print(f"  with 'S' in cigarstring: {has_s_in_cigarstring}")
    print(f"  supplementary: {is_supplementary} (with SA tag: {has_sa_tag})")

    print("\nrunning the REAL _fusion_rows_for_region() against this region...")
    bam2 = pysam.AlignmentFile(args.bam_path, "rb")
    bam2_primary = pysam.AlignmentFile(args.bam_path, "rb")
    region = {
        "window": "diagnostic",
        "chrom": args.chrom,
        "chromStart": args.chrom_start0,
        "chromEnd": args.chrom_end0,
    }
    primary_seq_cache = {}

    yielded = 0
    soft_clip_yielded = 0
    supplementary_yielded = 0

    for row in _fusion_rows_for_region(bam2, region, primary_seq_cache, bam2_primary):
        yielded += 1
        if row["chr_primary_align"]:
            supplementary_yielded += 1
        else:
            soft_clip_yielded += 1
    bam2.close()
    bam2_primary.close()

    print(f"rows yielded by _fusion_rows_for_region(): {yielded}")
    print(f"  soft-clip rows: {soft_clip_yielded}")
    print(f"  supplementary rows: {supplementary_yielded}")
    print(
        "\n(compare 'soft-clip rows' above against \"with 'S' in cigarstring\" from "
        "the raw fetch() totals -- a big gap there means rows are being silently "
        "dropped somewhere between the read and the yield, not by fetch() itself.)"
    )

    print("\nstep-by-step reproduction of the soft-clip branch's own logic...")
    bam3 = pysam.AlignmentFile(args.bam_path, "rb")
    reached_s_check = 0
    passed_s_check = 0
    built_row = 0
    failed_encoding = 0
    passed_encoding = 0

    for read in bam3.fetch(args.chrom, window_start0, window_end0):
        if read.is_unmapped:
            continue
        reached_s_check += 1
        cigar = read.cigarstring or ""
        if "S" not in cigar:
            continue
        passed_s_check += 1

        sequence = _strip_nuls(read.query_sequence or "")
        clipped_parts = clipped_sequences_from_cigar(sequence, read.cigartuples)
        clipped_sequence = _strip_nuls(", ".join(clipped_parts))
        start0 = read.reference_start
        end0 = alignment_end(start0, read.cigartuples)
        row = {
            "window": "diagnostic",
            "read_name": read.query_name,
            "read_1_2": "",
            "start": start0,
            "end": end0,
            "cigar": cigar,
            "chr_primary_align": "",
            "coord_primary_align": "",
            "strand_primary_align": "",
            "sequence": sequence,
            "clipped_sequence": clipped_sequence,
            "part_telomere": "False",
            "TTAGGG_count": 0,
            "CCCTAA_count": 0,
            "expected_pos_fusion": "",
        }
        built_row += 1
        if _row_has_bad_encoding(row):
            failed_encoding += 1
        else:
            passed_encoding += 1
    bam3.close()

    print(f"  reached (mapped) reads: {reached_s_check}")
    print(f"  passed 'S' in cigar check: {passed_s_check}")
    print(f"  rows built (should equal the line above): {built_row}")
    print(f"  rows failing _row_has_bad_encoding: {failed_encoding}")
    print(f"  rows passing _row_has_bad_encoding (would be yielded): {passed_encoding}")


if __name__ == "__main__":
    main()
