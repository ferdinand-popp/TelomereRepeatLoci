#!/usr/bin/env python3

import argparse

import ast

import pandas as pd

from pipeline.tables import read_tsv, write_tsv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("window_file")
    parser.add_argument("candidate_region_file")
    parser.add_argument("tumor_discordant_read_lower_limit", type=float)
    parser.add_argument("control_discordant_read_upper_limit", type=float)
    parser.add_argument("consider_blacklist")
    args = parser.parse_args()

    df = filter_candidates(
        read_tsv(args.window_file),
        args.tumor_discordant_read_lower_limit,
        args.control_discordant_read_upper_limit,
        args.consider_blacklist,
    )
    write_tsv(df, args.candidate_region_file, list(df.columns))


def filter_candidates(
    df: pd.DataFrame,
    tumor_discordant_read_lower_limit: float,
    control_discordant_read_upper_limit: float,
    consider_blacklist: str,
) -> pd.DataFrame:
    df = df.copy()
    if "tumor_discordant_read_count" not in df.columns:
        df["tumor_discordant_read_count"] = "0"
    if "control_discordant_read_count" not in df.columns:
        df["control_discordant_read_count"] = "0"

    df["tumor_discordant_read_count"] = pd.to_numeric(
        df["tumor_discordant_read_count"], errors="coerce"
    ).fillna(0)
    df["control_discordant_read_count"] = pd.to_numeric(
        df["control_discordant_read_count"], errors="coerce"
    ).fillna(0)

    filtered = df[
        (df["tumor_discordant_read_count"] >= tumor_discordant_read_lower_limit)
        & (df["control_discordant_read_count"] <= control_discordant_read_upper_limit)
    ]
    if consider_blacklist == "True" and "blacklisted" in filtered.columns:
        filtered = filtered[filtered["blacklisted"] != "yes"]
    filtered = fuse_overlapping_candidates(filtered)
    return drop_read_name_columns(filtered)


def drop_read_name_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop(
        columns=[
            col
            for col in df.columns
            if col in {"_tumor_read_names", "_control_read_names"}
        ],
        errors="ignore",
    )


def _normalize_read_names(values) -> set:
    if isinstance(values, set):
        return values
    if values is None:
        return set()
    if isinstance(values, list):
        return set(str(value) for value in values if value)
    if isinstance(values, str) and values:
        try:
            parsed = ast.literal_eval(values)
        except (ValueError, SyntaxError):
            parsed = None
        if isinstance(parsed, (set, list, tuple)):
            return {str(value) for value in parsed if value}
        if values in {"set()", "{}", "[]", "()"}:
            return set()
        if "," in values:
            return {part.strip() for part in values.split(",") if part.strip()}
        return {values}
    return set()


def _overlap_ratio(left: set, right: set) -> float:
    if not left or not right:
        return 0.0
    shared = len(left & right)
    return shared / float(min(len(left), len(right)))


def fuse_overlapping_candidates(
    df: pd.DataFrame,
    overlap_ratio_threshold: float = 0.5,
) -> pd.DataFrame:
    if df.empty:
        return df

    df = df.copy()
    df["chromStart"] = pd.to_numeric(df["chromStart"], errors="coerce")
    df["chromEnd"] = pd.to_numeric(df["chromEnd"], errors="coerce")
    df["_tumor_read_names"] = df["_tumor_read_names"].apply(_normalize_read_names)
    df["_control_read_names"] = df["_control_read_names"].apply(_normalize_read_names)
    df = df.sort_values(["chrom", "strand", "chromStart", "chromEnd"]).reset_index(
        drop=True
    )

    fused_rows = []
    current = df.iloc[0].to_dict()
    for _, row in df.iloc[1:].iterrows():
        row = row.to_dict()
        same_block = (
            row["chrom"] == current["chrom"]
            and row["strand"] == current["strand"]
            and row["chromStart"] <= current["chromEnd"]
        )
        if same_block:
            tumor_overlap = _overlap_ratio(
                current["_tumor_read_names"], row["_tumor_read_names"]
            )
            control_overlap = _overlap_ratio(
                current["_control_read_names"], row["_control_read_names"]
            )
            should_fuse = tumor_overlap >= overlap_ratio_threshold
            if not should_fuse and control_overlap >= overlap_ratio_threshold:
                should_fuse = True
            if should_fuse:
                current["chromStart"] = min(
                    current["chromStart"], row["chromStart"]
                )
                current["chromEnd"] = max(current["chromEnd"], row["chromEnd"])
                current["_tumor_read_names"] = (
                    current["_tumor_read_names"] | row["_tumor_read_names"]
                )
                current["_control_read_names"] = (
                    current["_control_read_names"] | row["_control_read_names"]
                )
                current["tumor_discordant_read_count"] = int(
                    len(current["_tumor_read_names"])
                )
                current["control_discordant_read_count"] = int(
                    len(current["_control_read_names"])
                )
                continue
        fused_rows.append(current)
        current = row
    fused_rows.append(current)
    return pd.DataFrame(fused_rows)


if __name__ == "__main__":
    main()
