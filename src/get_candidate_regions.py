#!/usr/bin/env python3

import argparse
import csv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("window_file")
    parser.add_argument("candidate_region_file")
    parser.add_argument("tumor_discordant_read_lower_limit", type=float)
    parser.add_argument("control_discordant_read_upper_limit", type=float)
    parser.add_argument("consider_blacklist")
    parser.add_argument("function_file", nargs="?")
    args = parser.parse_args()

    with open(args.window_file, newline="") as in_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        rows = []
        for row in reader:
            try:
                tumor = float(row.get("tumor_discordant_read_count", 0))
                control = float(row.get("control_discordant_read_count", 0))
            except (TypeError, ValueError):
                continue

            if tumor < args.tumor_discordant_read_lower_limit:
                continue
            if control > args.control_discordant_read_upper_limit:
                continue
            if args.consider_blacklist == "True" and row.get("blacklisted") == "yes":
                continue
            rows.append(row)

        fieldnames = reader.fieldnames or []

    with open(args.candidate_region_file, "w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
