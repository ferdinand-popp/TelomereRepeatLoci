#!/usr/bin/env python3

import argparse
import os

import pandas as pd

from pipeline.tables import WINDOWS_COLUMNS, read_tsv, write_tsv

MIN_MATE_MAPQ = 30
WINDOW_SIZE = 1000


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
            win_start = (pos // WINDOW_SIZE) * WINDOW_SIZE
            windows.append(
                {
                    "window": f"{chrom}_{win_start}_{strand}",
                    "chrom": chrom,
                    "chromStart": win_start,
                    "chromEnd": win_start + WINDOW_SIZE,
                    "strand": strand,
                }
            )
    return pd.DataFrame(windows).drop_duplicates()


def merge_adjacent_tumor_windows(df: pd.DataFrame) -> pd.DataFrame:
    """Merge consecutive same-chrom/strand 1kb windows into one region whenever
    both windows have nonzero tumor discordant-read support, so a locus whose
    support straddles a window boundary isn't undercounted below the
    downstream tumor/control thresholds. Mirrors count_discordant_reads.R's
    adjacent-window merge, which ran before thresholding rather than after.
    """
    if df.empty:
        return df

    df = df.sort_values(["chrom", "strand", "chromStart"]).reset_index(drop=True)

    merged_rows = []
    current = df.iloc[0].to_dict()
    for _, row in df.iloc[1:].iterrows():
        row = row.to_dict()
        adjacent = (
            row["chrom"] == current["chrom"]
            and row["strand"] == current["strand"]
            and row["chromStart"] == current["chromEnd"]
        )
        if (
            adjacent
            and current["tumor_discordant_read_count"] != 0
            and row["tumor_discordant_read_count"] != 0
        ):
            current["chromEnd"] = row["chromEnd"]
            current["_tumor_read_names"] = (
                current["_tumor_read_names"] | row["_tumor_read_names"]
            )
            current["_control_read_names"] = (
                current["_control_read_names"] | row["_control_read_names"]
            )
            current["tumor_discordant_read_count"] = len(current["_tumor_read_names"])
            current["control_discordant_read_count"] = len(
                current["_control_read_names"]
            )
            continue
        merged_rows.append(current)
        current = row
    merged_rows.append(current)
    return pd.DataFrame(merged_rows)


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
    merged_counts = merge_adjacent_tumor_windows(merged_counts)
    merged_counts["PID"] = pid

    merged_counts["blacklisted"] = ""
    if blacklist_path and os.path.exists(blacklist_path):
        blacklist_df = read_tsv(blacklist_path)
        if "window" in blacklist_df.columns:
            blacklist = set(blacklist_df["window"].tolist())
            merged_counts["blacklisted"] = merged_counts["window"].apply(
                lambda w: "yes" if w in blacklist else "no"
            )

    if not merged_counts.empty:
        merged_counts = merged_counts.sort_values(["chrom", "chromStart", "strand"])
    return merged_counts


if __name__ == "__main__":
    main()
