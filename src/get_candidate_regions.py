#!/usr/bin/env python3

import argparse

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
    return filtered


if __name__ == "__main__":
    main()
