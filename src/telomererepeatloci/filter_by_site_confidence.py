#!/usr/bin/env python3

import argparse

import pandas as pd

from pipeline.tables import read_tsv, write_tsv


EMPTY_VALUES = {"", "NA", "NaN", "nan", "None", None}
DEFAULT_MAX_TUMOR_NOISE_RATIO = 0.8
DEFAULT_CONTROL_MAX_SEQ_DISTANCE = 2
DEFAULT_CONTROL_MAX_TELO_CLIPPED_AT_SITE = 2


def parse_float(value):
    if value in EMPTY_VALUES:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_int(value):
    if value in EMPTY_VALUES:
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def passes_confidence_filters(
    row: dict,
    max_tumor_noise_ratio: float,
    control_max_seq_distance: int,
    control_max_telo_clipped_at_site: int,
) -> bool:
    """Diagnostic columns from assess_site_confidence.py -> keep/drop decision.

    Blank/missing diagnostics (no insertion_site, no control data, etc.) never
    cause a drop on their own -- only a computed value that actually exceeds a
    threshold does.
    """
    noise_ratio = parse_float(row.get("tumor_noise_ratio"))
    if noise_ratio is not None and noise_ratio > max_tumor_noise_ratio:
        return False

    control_count = parse_int(row.get("control_telo_clipped_at_insertion_site"))
    if control_count is None or control_count <= 0:
        return True

    control_distance = parse_int(row.get("control_min_seq_distance_to_tumor"))
    if control_distance is not None and control_distance <= control_max_seq_distance:
        return False
    if control_count > control_max_telo_clipped_at_site:
        return False
    return True


def filter_regions(
    df: pd.DataFrame,
    max_tumor_noise_ratio: float = DEFAULT_MAX_TUMOR_NOISE_RATIO,
    control_max_seq_distance: int = DEFAULT_CONTROL_MAX_SEQ_DISTANCE,
    control_max_telo_clipped_at_site: int = DEFAULT_CONTROL_MAX_TELO_CLIPPED_AT_SITE,
) -> pd.DataFrame:
    rows = df.to_dict("records")
    kept = [
        row
        for row in rows
        if passes_confidence_filters(
            row,
            max_tumor_noise_ratio,
            control_max_seq_distance,
            control_max_telo_clipped_at_site,
        )
    ]
    print(f"Confidence filter: kept {len(kept)} of {len(rows)} regions.")
    return pd.DataFrame(kept, columns=list(df.columns))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_regions_confidence_file")
    parser.add_argument("outfile")
    parser.add_argument(
        "--max-tumor-noise-ratio",
        type=float,
        default=DEFAULT_MAX_TUMOR_NOISE_RATIO,
    )
    parser.add_argument(
        "--control-max-seq-distance",
        type=int,
        default=DEFAULT_CONTROL_MAX_SEQ_DISTANCE,
    )
    parser.add_argument(
        "--control-max-telo-clipped-at-site",
        type=int,
        default=DEFAULT_CONTROL_MAX_TELO_CLIPPED_AT_SITE,
    )
    args = parser.parse_args()

    df = filter_regions(
        read_tsv(args.candidate_regions_confidence_file),
        args.max_tumor_noise_ratio,
        args.control_max_seq_distance,
        args.control_max_telo_clipped_at_site,
    )
    write_tsv(df, args.outfile, list(df.columns))


if __name__ == "__main__":
    main()
