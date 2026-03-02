#!/usr/bin/env python3

import argparse
import csv
import math
import statistics
from collections import Counter, defaultdict


EMPTY_VALUES = {"", "NA", "NaN", "nan", "None", None}


def is_true(value):
    return str(value).strip().lower() in {"true", "t", "1"}


def parse_int(value):
    if value in EMPTY_VALUES:
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def parse_float(value):
    if value in EMPTY_VALUES:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def median(values):
    values = [v for v in values if v is not None]
    if not values:
        return None
    return statistics.median(values)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_region_file")
    parser.add_argument("clipped_reads_file")
    parser.add_argument("discordant_read_file")
    parser.add_argument("outfile")
    parser.add_argument("function_file", nargs="?")
    args = parser.parse_args()

    with open(args.candidate_region_file, newline="") as handle:
        candidate_reader = csv.DictReader(handle, delimiter="\t")
        candidate_rows = list(candidate_reader)
        candidate_fields = candidate_reader.fieldnames or []

    with open(args.clipped_reads_file, newline="") as handle:
        clipped_rows = list(csv.DictReader(handle, delimiter="\t"))

    with open(args.discordant_read_file, newline="") as handle:
        discordant_rows = list(csv.DictReader(handle, delimiter="\t"))

    clipped_by_window = defaultdict(list)
    for row in clipped_rows:
        clipped_by_window[row.get("window", "")].append(row)

    new_fields = [
        "insertion_site",
        "pos_telomeres_from_insertion",
        "reads_supporting_insertion_pos",
        "sum_TTAGGG_count",
        "sum_CCCTAA_count",
        "repeat_forward",
    ]
    output_fields = candidate_fields + [f for f in new_fields if f not in candidate_fields]

    for region in candidate_rows:
        window = region.get("window", "")
        strand = region.get("strand", "")
        chrom = region.get("chrom", "")
        window_start = parse_int(region.get("chromStart"))
        window_end = parse_int(region.get("chromEnd"))

        for field in new_fields:
            region[field] = ""

        clipped_window = clipped_by_window.get(window, [])
        clipped_tel = [r for r in clipped_window if is_true(r.get("part_telomere"))]
        if not clipped_tel:
            continue

        if strand == "+":
            expected_pos = "downstream"
            clipped_col = "end"
        elif strand == "-":
            expected_pos = "upstream"
            clipped_col = "start"
        else:
            continue

        filtered = [r for r in clipped_tel if r.get("expected_pos_fusion") == expected_pos]

        discordant_filtered = []
        for d in discordant_rows:
            if d.get("mate_chr") != chrom:
                continue
            if d.get("mate_strand") != strand:
                continue
            pos = parse_int(d.get("mate_position"))
            if pos is None or window_start is None or window_end is None:
                continue
            if window_start <= pos <= window_end:
                discordant_filtered.append(pos)

        med = median(discordant_filtered)
        if med is None:
            continue
        med += 50

        filtered_pos = []
        for r in filtered:
            start = parse_int(r.get("start"))
            end = parse_int(r.get("end"))
            if strand == "+" and end is not None and end > med:
                filtered_pos.append(r)
            elif strand == "-" and start is not None and start < med:
                filtered_pos.append(r)

        pos_to_cigars = defaultdict(set)
        for r in filtered_pos:
            pos = parse_int(r.get(clipped_col))
            if pos is None:
                continue
            pos_to_cigars[pos].add(r.get("cigar", ""))

        if not pos_to_cigars:
            continue

        best_unique = max(len(cigars) for cigars in pos_to_cigars.values())
        insertion_candidates = [pos for pos, cigars in pos_to_cigars.items() if len(cigars) == best_unique]

        if len(insertion_candidates) != 1:
            continue

        insertion_site = insertion_candidates[0]
        region["insertion_site"] = str(insertion_site)
        region["pos_telomeres_from_insertion"] = expected_pos
        region["reads_supporting_insertion_pos"] = str(best_unique)

        at_insertion = [r for r in filtered_pos if parse_int(r.get(clipped_col)) == insertion_site]
        sum_t = sum(parse_int(r.get("TTAGGG_count")) or 0 for r in at_insertion)
        sum_c = sum(parse_int(r.get("CCCTAA_count")) or 0 for r in at_insertion)

        region["sum_TTAGGG_count"] = str(sum_t)
        region["sum_CCCTAA_count"] = str(sum_c)

        if sum_t > sum_c:
            region["repeat_forward"] = "TTAGGG"
        elif sum_c > sum_t:
            region["repeat_forward"] = "CCCTAA"

    with open(args.outfile, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(candidate_rows)


if __name__ == "__main__":
    main()
