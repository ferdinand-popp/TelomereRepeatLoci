#!/usr/bin/env python3

import argparse

import pandas as pd

from pipeline.tables import read_tsv, write_tsv


EMPTY_VALUES = {"", "NA", "NaN", "nan", "None", None}
DEFAULT_MAX_TUMOR_NOISE_RATIO = 0.8
DEFAULT_CONTROL_MAX_SEQ_DISTANCE = 2
DEFAULT_CONTROL_MAX_TELO_CLIPPED_AT_SITE = 2
DEFAULT_MIN_INSERTION_SUPPORT = 2.0


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
    min_insertion_support: float,
) -> bool:
    """Diagnostic columns from assess_site_confidence.py -> keep/drop decision.

    A region with no predicted insertion_site is dropped outright -- it has no
    locus to review or plot, so it can't be "kept" in any meaningful sense.
    For regions that do have an insertion_site, blank/missing diagnostics (no
    control data, etc.) never cause a drop on their own -- only a computed
    value that actually exceeds a threshold does.

    Regions below min_insertion_support are dropped here too, matching the
    make_bed_for_visualization.py threshold, so this table never contains
    rows that silently never make it into a plot.
    """
    if row.get("insertion_site") in EMPTY_VALUES:
        return False

    support = parse_int(row.get("reads_supporting_insertion_pos")) or 0
    if support < min_insertion_support:
        return False

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
    min_insertion_support: float = DEFAULT_MIN_INSERTION_SUPPORT,
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
            min_insertion_support,
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
    parser.add_argument(
        "--min-insertion-support",
        type=float,
        default=DEFAULT_MIN_INSERTION_SUPPORT,
        help=(
            "Minimum reads_supporting_insertion_pos required to keep a region. "
            "Should match --plot-min-support so this table never contains "
            "regions that make_bed_for_visualization.py would silently drop. "
            "Default: 2."
        ),
    )
    args = parser.parse_args()

    df = filter_regions(
        read_tsv(args.candidate_regions_confidence_file),
        args.max_tumor_noise_ratio,
        args.control_max_seq_distance,
        args.control_max_telo_clipped_at_site,
        args.min_insertion_support,
    )
    write_tsv(df, args.outfile, list(df.columns))


if __name__ == "__main__":
    main()
