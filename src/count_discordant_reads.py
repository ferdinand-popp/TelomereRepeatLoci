#!/usr/bin/env python3

import argparse
import os

import pandas as pd

from pipeline.tables import WINDOWS_COLUMNS, read_tsv, write_tsv

MIN_MATE_MAPQ = 30
WINDOW_SIZE = 1000
WINDOW_STEP = 500


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


def build_windows(df):
    if df.empty:
        return pd.DataFrame(
            columns=["window", "chrom", "chromStart", "chromEnd", "strand"]
        )

    windows = []
    for (chrom, strand), group in df.groupby(["mate_chr", "mate_strand"], dropna=False):
        positions = group["mate_position"].astype(int).tolist()
        for pos in positions:
            start = (pos // WINDOW_STEP) * WINDOW_STEP
            for offset in [0, -WINDOW_STEP]:
                win_start = start + offset
                if win_start < 0:
                    continue
                win_end = win_start + WINDOW_SIZE
                if win_start <= pos < win_end:
                    windows.append(
                        {
                            "window": f"{chrom}_{win_start}_{strand}",
                            "chrom": chrom,
                            "chromStart": win_start,
                            "chromEnd": win_end,
                            "strand": strand,
                        }
                    )
    return pd.DataFrame(windows).drop_duplicates()


def count_windows(df, windows, name_column):
    if windows.empty:
        windows = windows.copy()
        windows["count"] = 0
        windows[name_column] = [set() for _ in range(len(windows))]
        return windows
    df = df.copy()
    df["mate_position"] = df["mate_position"].astype(int)

    windows = windows.copy()
    windows["count"] = 0
    windows[name_column] = [set() for _ in range(len(windows))]
    for idx, row in windows.iterrows():
        chrom = row["chrom"]
        strand = row["strand"]
        start = row["chromStart"]
        end = row["chromEnd"]
        mask = (
            (df["mate_chr"] == chrom)
            & (df["mate_strand"] == strand)
            & (df["mate_position"] >= start)
            & (df["mate_position"] < end)
        )
        read_names = set(df.loc[mask, "read_name"].astype(str).tolist())
        windows.at[idx, "count"] = int(len(read_names))
        windows.at[idx, name_column] = read_names
    return windows


def merge_overlapping(rows):
    rows = rows.sort_values(["chrom", "strand", "chromStart"]).to_dict("records")
    merged = []
    for row in rows:
        if not merged:
            merged.append(dict(row))
            continue
        prev = merged[-1]
        overlaps = prev["chrom"] == row["chrom"] and prev["strand"] == row["strand"]
        overlaps = overlaps and prev["chromEnd"] >= row["chromStart"]
        tumor_counts = (
            prev["tumor_discordant_read_count"] != 0
            and row["tumor_discordant_read_count"] != 0
        )
        if overlaps and tumor_counts:
            prev["chromEnd"] = max(prev["chromEnd"], row["chromEnd"])
            prev_tumor = prev.get("_tumor_read_names", set())
            row_tumor = row.get("_tumor_read_names", set())
            prev_control = prev.get("_control_read_names", set())
            row_control = row.get("_control_read_names", set())

            prev["_tumor_read_names"] = prev_tumor | row_tumor
            prev["_control_read_names"] = prev_control | row_control
            prev["tumor_discordant_read_count"] = len(prev["_tumor_read_names"])
            prev["control_discordant_read_count"] = len(prev["_control_read_names"])
            if prev["blacklisted"] == "yes" or row["blacklisted"] == "yes":
                prev["blacklisted"] = "yes"
            elif prev["blacklisted"] == "no" or row["blacklisted"] == "no":
                prev["blacklisted"] = "no"
            else:
                prev["blacklisted"] = ""
            prev["window"] = f"{prev['chrom']}_{prev['chromStart']}_{prev['strand']}"
        else:
            merged.append(dict(row))
    return pd.DataFrame(merged)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--discordantReadFileTumor", required=True)
    parser.add_argument("-c", "--discordantReadFileControl", required=True)
    parser.add_argument("-b", "--blacklist_file", required=True)
    parser.add_argument("-o", "--outFile", required=True)
    parser.add_argument("-f", "--function_file", required=False)
    args = parser.parse_args()

    tumor_df = load_discordant(args.discordantReadFileTumor)
    control_df = load_discordant(args.discordantReadFileControl)

    pid = os.path.basename(args.outFile).replace(
        "_discordant_reads_1_kb_windows.tsv", ""
    )

    windows = build_windows(pd.concat([tumor_df, control_df], ignore_index=True))
    tumor_counts = count_windows(tumor_df, windows, "_tumor_read_names").rename(
        columns={"count": "tumor_discordant_read_count"}
    )
    control_counts = count_windows(control_df, windows, "_control_read_names").rename(
        columns={"count": "control_discordant_read_count"}
    )

    merged_counts = windows.merge(
        tumor_counts[["window", "tumor_discordant_read_count", "_tumor_read_names"]],
        on="window",
        how="left",
    ).merge(
        control_counts[
            ["window", "control_discordant_read_count", "_control_read_names"]
        ],
        on="window",
        how="left",
    )
    merged_counts["tumor_discordant_read_count"] = (
        merged_counts["tumor_discordant_read_count"].fillna(0).astype(int)
    )
    merged_counts["control_discordant_read_count"] = (
        merged_counts["control_discordant_read_count"].fillna(0).astype(int)
    )
    merged_counts["_tumor_read_names"] = merged_counts["_tumor_read_names"].apply(
        lambda x: x if isinstance(x, set) else set()
    )
    merged_counts["_control_read_names"] = merged_counts["_control_read_names"].apply(
        lambda x: x if isinstance(x, set) else set()
    )
    merged_counts["PID"] = pid

    merged_counts["blacklisted"] = ""
    if args.blacklist_file and os.path.exists(args.blacklist_file):
        blacklist_df = read_tsv(args.blacklist_file)
        if "window" in blacklist_df.columns:
            blacklist = set(blacklist_df["window"].tolist())
            merged_counts["blacklisted"] = merged_counts["window"].apply(
                lambda w: "yes" if w in blacklist else "no"
            )

    merged_counts = merge_overlapping(merged_counts)
    if "_tumor_read_names" in merged_counts.columns:
        merged_counts = merged_counts.drop(columns=["_tumor_read_names"])
    if "_control_read_names" in merged_counts.columns:
        merged_counts = merged_counts.drop(columns=["_control_read_names"])
    if not merged_counts.empty:
        merged_counts = merged_counts.sort_values(["chrom", "chromStart", "strand"])
    write_tsv(merged_counts, args.outFile, WINDOWS_COLUMNS)


if __name__ == "__main__":
    main()
