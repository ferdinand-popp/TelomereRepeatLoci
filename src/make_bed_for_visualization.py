#!/usr/bin/env python3

import argparse
import csv


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
    fieldnames = ["#chrom", "chromStart", "chromEnd", "pos", "pid"]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_region_file")
    parser.add_argument("outfile1")
    parser.add_argument("outfile2")
    parser.add_argument("pid")
    args = parser.parse_args()

    with open(args.candidate_region_file, newline="") as handle:
        candidate_rows = list(csv.DictReader(handle, delimiter="\t"))

    write_bed(args.outfile1, build_bed_rows(candidate_rows, args.pid, 500))
    write_bed(args.outfile2, build_bed_rows(candidate_rows, args.pid, 100))


if __name__ == "__main__":
    main()
