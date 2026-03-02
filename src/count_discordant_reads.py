#!/usr/bin/env python3

import argparse
import csv
import os
from collections import Counter

MIN_MATE_MAPQ = 30


def read_discordant_counts(path, sample):
    counts = Counter()
    if not path or path == "NULL" or not os.path.exists(path):
        return counts

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                if float(row.get("mate_mapq", 0)) <= MIN_MATE_MAPQ:
                    continue
            except (TypeError, ValueError):
                continue

            chrom = row.get("mate_chr", "")
            strand = row.get("mate_strand", "")
            try:
                pos = int(float(row.get("mate_position", 0)))
            except (TypeError, ValueError):
                continue

            pos_1kb = (pos // 1000) * 1000
            window = f"{chrom}_{pos_1kb}_{strand}"
            counts[window] += 1

    return counts


def parse_window(window):
    chrom, start, strand = window.rsplit("_", 2)
    chrom_start = int(start)
    chrom_end = chrom_start + 1000
    return chrom, chrom_start, chrom_end, strand


def merge_adjacent(rows):
    rows.sort(key=lambda r: (r["chrom"], r["strand"], r["chromStart"]))

    merged = []
    for row in rows:
        if not merged:
            merged.append(dict(row))
            continue

        prev = merged[-1]
        if (
            prev["chrom"] == row["chrom"]
            and prev["strand"] == row["strand"]
            and prev["chromEnd"] == row["chromStart"]
            and prev["tumor_discordant_read_count"] != 0
            and row["tumor_discordant_read_count"] != 0
        ):
            prev["chromEnd"] = row["chromEnd"]
            prev["tumor_discordant_read_count"] += row["tumor_discordant_read_count"]
            prev["control_discordant_read_count"] += row["control_discordant_read_count"]
            if prev["blacklisted"] == "yes" or row["blacklisted"] == "yes":
                prev["blacklisted"] = "yes"
            elif prev["blacklisted"] == "no" or row["blacklisted"] == "no":
                prev["blacklisted"] = "no"
            else:
                prev["blacklisted"] = ""
            prev["window"] = f"{prev['chrom']}_{prev['chromStart']}_{prev['strand']}"
        else:
            merged.append(dict(row))

    return merged


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--discordantReadFileTumor", required=True)
    parser.add_argument("-c", "--discordantReadFileControl", required=True)
    parser.add_argument("-b", "--blacklist_file", required=True)
    parser.add_argument("-o", "--outFile", required=True)
    parser.add_argument("-f", "--function_file", required=False)
    args = parser.parse_args()

    tumor_counts = read_discordant_counts(args.discordantReadFileTumor, "tumor")
    control_counts = read_discordant_counts(args.discordantReadFileControl, "control")

    windows = sorted(set(tumor_counts) | set(control_counts))
    pid = os.path.basename(args.outFile).replace("_discordant_reads_1_kb_windows.tsv", "")

    blacklist = set()
    if args.blacklist_file and os.path.exists(args.blacklist_file):
        with open(args.blacklist_file, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if row.get("window"):
                    blacklist.add(row["window"])

    rows = []
    for window in windows:
        chrom, chrom_start, chrom_end, strand = parse_window(window)
        rows.append(
            {
                "PID": pid,
                "window": window,
                "chrom": chrom,
                "chromStart": chrom_start,
                "chromEnd": chrom_end,
                "strand": strand,
                "tumor_discordant_read_count": int(tumor_counts.get(window, 0)),
                "control_discordant_read_count": int(control_counts.get(window, 0)),
                "blacklisted": "yes" if window in blacklist else ("no" if blacklist else ""),
            }
        )

    rows = merge_adjacent(rows)

    fieldnames = [
        "PID",
        "window",
        "chrom",
        "chromStart",
        "chromEnd",
        "strand",
        "tumor_discordant_read_count",
        "control_discordant_read_count",
        "blacklisted",
    ]

    with open(args.outFile, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
