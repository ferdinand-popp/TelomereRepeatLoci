#!/usr/bin/env python3

import argparse
from collections import Counter

import pandas as pd
import pysam
from telomererepeatloci.pipeline.tables import DISCORDANT_READS_COLUMNS, write_tsv


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract read names and mate positions for discordant telomeric reads."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input BAM with telomeric reads (or any BAM you want to screen).",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output TSV file for discordant reads.",
    )
    return parser.parse_args()


def read_discordant_pairs(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        name_counts = Counter(read.query_name for read in bam.fetch(until_eof=True))

    rows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.mate_is_unmapped:
                continue
            if read.next_reference_id == -1:
                continue
            if name_counts.get(read.query_name, 0) == 2:
                continue
            rows.append(
                {
                    "read_name": read.query_name,
                    "mate_chr": read.next_reference_name,
                    "mate_position": str(read.next_reference_start),
                }
            )
    return rows


def main():
    args = parse_args()
    rows = read_discordant_pairs(args.input)
    df = pd.DataFrame(rows)
    write_tsv(df, args.output, DISCORDANT_READS_COLUMNS)


def run(input_bam: str, output_tsv: str) -> pd.DataFrame:
    rows = read_discordant_pairs(input_bam)
    df = pd.DataFrame(rows)
    write_tsv(df, output_tsv, DISCORDANT_READS_COLUMNS)
    return df


if __name__ == "__main__":
    main()
