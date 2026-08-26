#!/usr/bin/env python3

import argparse
import bisect
import os
import sys
from collections import defaultdict

import pandas as pd

from pipeline.tables import WINDOWS_COLUMNS, read_tsv, write_tsv

MIN_MATE_MAPQ = 30
# Legacy fixed-grid window size the shipped blacklist files were built
# against (blacklists/*.tsv only carry a "chrom_start_strand" window id, so
# this is needed purely to reconstruct their extent for overlap checks).
WINDOW_SIZE = 1000
# Max gap (bp) between two same-strand discordant mate positions for them to
# be considered the same locus and merged into one candidate region.
MERGE_GAP = 1000
# Floor for a cluster's span, matching the legacy fixed-grid window size. A
# cluster narrower than this (e.g. a handful of reads landing within a few
# bp of each other) is padded outward symmetrically to this width, so sparse
# regions keep at least as much surrounding context for the downstream
# fusion-read search as the old grid gave, instead of being only as wide as
# the exact discordant-read positions happen to span.
MIN_REGION_WIDTH = 1000
# Purely informational: a cluster wider than this gets a warning printed
# (visible in logs) but is never split or truncated.
LARGE_CLUSTER_WARN_BP = 50000

EMPTY_REGION_COLUMNS = [
    "window",
    "chrom",
    "chromStart",
    "chromEnd",
    "strand",
    "tumor_discordant_read_count",
    "control_discordant_read_count",
    "_tumor_read_names",
    "_control_read_names",
]


def load_discordant(path):
    if not path or path == "NULL" or not os.path.exists(path):
        return pd.DataFrame(
            columns=[
                "read_name",
                "mate_chr",
                "mate_position",
                "mate_mapq",
                "mate_strand",
            ]
        )
    df = read_tsv(path)
    for col in ["read_name", "mate_chr", "mate_position", "mate_mapq", "mate_strand"]:
        if col not in df.columns:
            df[col] = ""
    df["mate_mapq"] = pd.to_numeric(df["mate_mapq"], errors="coerce").fillna(0)
    df["mate_position"] = pd.to_numeric(df["mate_position"], errors="coerce")
    df = df[df["mate_mapq"] > MIN_MATE_MAPQ]
    df = df[df["mate_position"].notna()]
    return df


def _finalize_cluster(cluster):
    chrom = cluster["chrom"]
    strand = cluster["strand"]
    start = cluster["start_pos"]
    end = cluster["end_pos"] + 1
    span = end - start
    if span < MIN_REGION_WIDTH:
        pad = MIN_REGION_WIDTH - span
        left_pad = pad // 2
        right_pad = pad - left_pad
        start = max(0, start - left_pad)
        end = end + right_pad
        span = end - start
    if span > LARGE_CLUSTER_WARN_BP:
        print(
            f"Warning: candidate region {chrom}:{start}-{end} ({strand}) spans "
            f"{span}bp of discordant-read support; keeping it as a single "
            "region rather than splitting it.",
            file=sys.stderr,
        )
    return {
        "window": f"{chrom}_{start}_{strand}",
        "chrom": chrom,
        "chromStart": start,
        "chromEnd": end,
        "strand": strand,
        "tumor_discordant_read_count": len(cluster["_tumor_read_names"]),
        "control_discordant_read_count": len(cluster["_control_read_names"]),
        "_tumor_read_names": cluster["_tumor_read_names"],
        "_control_read_names": cluster["_control_read_names"],
    }


def build_regions(tumor_df, control_df, merge_gap=MERGE_GAP):
    """Cluster discordant mate positions directly instead of tiling a fixed,
    overlapping window grid. For each chrom+strand, same-strand mate
    positions within `merge_gap` bp of each other are merged into one
    candidate region (single-linkage, bedtools-merge style), so a region's
    extent reflects the actual read support instead of a grid cell.
    """
    tumor_df = tumor_df.copy()
    tumor_df["_sample"] = "tumor"
    control_df = control_df.copy()
    control_df["_sample"] = "control"
    combined = pd.concat([tumor_df, control_df], ignore_index=True)

    if combined.empty:
        return pd.DataFrame(columns=EMPTY_REGION_COLUMNS)

    combined["mate_position"] = combined["mate_position"].astype(int)
    combined = combined.sort_values(["mate_chr", "mate_strand", "mate_position"])

    regions = []
    cluster = None
    columns = combined[
        ["mate_chr", "mate_strand", "mate_position", "read_name", "_sample"]
    ]
    for chrom, strand, pos, read_name, sample in columns.itertuples(index=False):
        if (
            cluster is not None
            and cluster["chrom"] == chrom
            and cluster["strand"] == strand
            and pos - cluster["end_pos"] <= merge_gap
        ):
            cluster["end_pos"] = pos
        else:
            if cluster is not None:
                regions.append(_finalize_cluster(cluster))
            cluster = {
                "chrom": chrom,
                "strand": strand,
                "start_pos": pos,
                "end_pos": pos,
                "_tumor_read_names": set(),
                "_control_read_names": set(),
            }
        cluster[f"_{sample}_read_names"].add(str(read_name))
    if cluster is not None:
        regions.append(_finalize_cluster(cluster))

    return pd.DataFrame(regions, columns=EMPTY_REGION_COLUMNS)


def _parse_window_id(window_id):
    parts = str(window_id).rsplit("_", 2)
    if len(parts) != 3:
        return None
    chrom, start, strand = parts
    try:
        start = int(start)
    except ValueError:
        return None
    return chrom, start, strand


def _parse_blacklist_intervals(window_ids):
    """Legacy blacklist files carry a fixed-grid window id per row; parse
    each back into a (chrom, strand) -> merged, sorted list of (start, end)
    intervals so dynamically-sized candidate regions can be checked against
    them by genomic overlap instead of exact window-id match.
    """
    raw = defaultdict(list)
    for window_id in window_ids:
        parsed = _parse_window_id(window_id)
        if parsed is None:
            continue
        chrom, start, strand = parsed
        raw[(chrom, strand)].append((start, start + WINDOW_SIZE))

    merged = {}
    for key, spans in raw.items():
        spans.sort()
        merged_spans = []
        for start, end in spans:
            if merged_spans and start <= merged_spans[-1][1]:
                merged_spans[-1] = (merged_spans[-1][0], max(merged_spans[-1][1], end))
            else:
                merged_spans.append((start, end))
        merged[key] = merged_spans
    return merged


def _overlaps_blacklist(chrom, strand, chrom_start, chrom_end, blacklist_intervals):
    spans = blacklist_intervals.get((chrom, strand))
    if not spans:
        return False
    starts = [s for s, _ in spans]
    idx = bisect.bisect_right(starts, chrom_end) - 1
    for i in (idx, idx - 1):
        if 0 <= i < len(spans):
            start, end = spans[i]
            if start < chrom_end and chrom_start < end:
                return True
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--discordantReadFileTumor", required=True)
    parser.add_argument("-c", "--discordantReadFileControl", required=True)
    parser.add_argument("-b", "--blacklist_file", required=True)
    parser.add_argument("-o", "--outFile", required=True)
    parser.add_argument("-f", "--function_file", required=False)
    args = parser.parse_args()

    df = compute_windows(
        args.discordantReadFileTumor,
        args.discordantReadFileControl,
        args.blacklist_file,
        args.outFile,
    )
    write_tsv(df, args.outFile, WINDOWS_COLUMNS)


def compute_windows(
    tumor_path: str,
    control_path: str,
    blacklist_path: str,
    output_path: str,
) -> pd.DataFrame:
    tumor_df = load_discordant(tumor_path)
    control_df = load_discordant(control_path)

    pid = os.path.basename(output_path).replace(
        "_discordant_reads_1_kb_windows.tsv", ""
    )

    merged_counts = build_regions(tumor_df, control_df)
    merged_counts["PID"] = pid

    merged_counts["blacklisted"] = ""
    if blacklist_path and os.path.exists(blacklist_path) and not merged_counts.empty:
        blacklist_df = read_tsv(blacklist_path)
        if "window" in blacklist_df.columns:
            blacklist_intervals = _parse_blacklist_intervals(
                blacklist_df["window"].tolist()
            )
            blacklisted_flags = []
            for row in merged_counts.itertuples(index=False):
                overlaps = _overlaps_blacklist(
                    row.chrom,
                    row.strand,
                    row.chromStart,
                    row.chromEnd,
                    blacklist_intervals,
                )
                blacklisted_flags.append("yes" if overlaps else "no")
            merged_counts["blacklisted"] = blacklisted_flags

    if not merged_counts.empty:
        merged_counts = merged_counts.sort_values(["chrom", "chromStart", "strand"])
    return merged_counts


if __name__ == "__main__":
    main()
