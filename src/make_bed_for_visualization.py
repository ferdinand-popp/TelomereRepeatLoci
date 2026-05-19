#!/usr/bin/env python3

import argparse

import pandas as pd

from pipeline.tables import BED_COLUMNS, read_tsv, write_tsv


def build_bed_rows(candidate_rows, pid, flank):
    rows = []
    for row in candidate_rows:
        try:
            insertion_site = int(float(row.get("insertion_site", "")))
        except (TypeError, ValueError):
            continue
        try:
            support = float(row.get("reads_supporting_insertion_pos", 0))
        except (TypeError, ValueError):
            support = 0
        if support <= 2:
            continue

        chrom_start = max(0, insertion_site - flank)
        chrom_end = insertion_site + flank
        rows.append(
            {
                "#chrom": row.get("chrom", ""),
                "chromStart": chrom_start,
                "chromEnd": chrom_end,
                "pos": insertion_site,
                "pid": pid,
            }
        )
    return rows


def write_bed(path, rows):
    df = pd.DataFrame(rows)
    write_tsv(df, path, BED_COLUMNS)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_region_file")
    parser.add_argument("outfile1")
    parser.add_argument("outfile2")
    parser.add_argument("pid")
    args = parser.parse_args()

    build_beds(args.candidate_region_file, args.pid, args.outfile1, args.outfile2)


def build_beds(
    candidate_region_file: str,
    pid: str,
    zoomed_out_path: str,
    zoomed_in_path: str,
):
    candidate_rows = read_tsv(candidate_region_file).to_dict("records")

    write_bed(zoomed_out_path, build_bed_rows(candidate_rows, pid, 500))
    write_bed(zoomed_in_path, build_bed_rows(candidate_rows, pid, 100))


if __name__ == "__main__":
    main()
