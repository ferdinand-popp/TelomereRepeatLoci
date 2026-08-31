#!/usr/bin/env python3

import argparse

import pandas as pd

from pipeline.tables import read_tsv, write_tsv


EMPTY_VALUES = {"", "NA", "NaN", "nan", "None", None}
DEFAULT_MAX_TUMOR_NOISE_RATIO = 0.8
DEFAULT_CONTROL_MAX_SEQ_DISTANCE = 2
DEFAULT_CONTROL_MAX_TELO_CLIPPED_AT_SITE = 3
DEFAULT_MIN_INSERTION_SUPPORT = 2.0
DEFAULT_MIN_CONTROL_READS_AT_SITE = 3
DEFAULT_MAX_READS_AT_SITE = 1600


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
    min_control_reads_at_site: int = DEFAULT_MIN_CONTROL_READS_AT_SITE,
    max_reads_at_site: int = DEFAULT_MAX_READS_AT_SITE,
) -> bool:
    """Diagnostic columns from assess_site_confidence.py -> keep/drop decision.

    A region with no predicted insertion_site is dropped outright -- it has no
    locus to review or plot, so it can't be "kept" in any meaningful sense.
    For regions that do have an insertion_site, blank/missing diagnostics (no
    control data, control BAM not given, etc.) never cause a drop on their
    own -- only a computed value that actually exceeds a threshold does.

    The one exception is a "control looks clean" verdict (zero telomeric
    clipped reads at the site): that's only meaningful if control actually
    had enough read depth at the site to show a signal. A populated but thin
    control_all_reads_at_site (e.g. 0-2 reads) makes the clean verdict
    uninformative, so it's treated as a drop rather than a pass. A *missing*
    control_all_reads_at_site (no control BAM) still doesn't cause a drop.

    A pileup of thousands of reads at one site (collapsed-repeat/mapping
    artifact loci) is dropped outright: visualize_telomere_insertions.py's
    plot_region silently skips plotting once either side's read count in the
    region reaches 3000, so a candidate that deep would never render anyway.

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

    all_reads_at_site = parse_int(row.get("all_reads_at_site"))
    if all_reads_at_site is not None and all_reads_at_site > max_reads_at_site:
        return False

    control_all_reads_at_site = parse_int(row.get("control_all_reads_at_site"))
    if (
        control_all_reads_at_site is not None
        and control_all_reads_at_site > max_reads_at_site
    ):
        return False

    control_count = parse_int(row.get("control_telo_clipped_at_insertion_site"))
    if control_count is None or control_count <= 0:
        if (
            control_all_reads_at_site is not None
            and control_all_reads_at_site < min_control_reads_at_site
        ):
            return False
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
    min_control_reads_at_site: int = DEFAULT_MIN_CONTROL_READS_AT_SITE,
    max_reads_at_site: int = DEFAULT_MAX_READS_AT_SITE,
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
            min_control_reads_at_site,
            max_reads_at_site,
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
    parser.add_argument(
        "--min-control-reads-at-site",
        type=int,
        default=DEFAULT_MIN_CONTROL_READS_AT_SITE,
        help=(
            "Minimum control_all_reads_at_site required to trust a 'no control "
            "telomeric clips at the site' verdict as evidence of tumor "
            "specificity. Below this, control coverage is too thin for the "
            "absence of a clip to mean anything, so the region is dropped. "
            "Does not apply when control_all_reads_at_site is missing "
            "entirely (e.g. no control BAM). Default: 3."
        ),
    )
    parser.add_argument(
        "--max-reads-at-site",
        type=int,
        default=DEFAULT_MAX_READS_AT_SITE,
        help=(
            "Maximum all_reads_at_site or control_all_reads_at_site allowed "
            "before a region is dropped as a collapsed-repeat/mapping-"
            "artifact pileup. visualize_telomere_insertions.py silently "
            "skips plotting once either side's read count in the region "
            "reaches 3000, so candidates this deep would never render "
            "anyway. Default: 1600."
        ),
    )
    args = parser.parse_args()

    df = filter_regions(
        read_tsv(args.candidate_regions_confidence_file),
        args.max_tumor_noise_ratio,
        args.control_max_seq_distance,
        args.control_max_telo_clipped_at_site,
        args.min_insertion_support,
        args.min_control_reads_at_site,
        args.max_reads_at_site,
    )
    write_tsv(df, args.outfile, list(df.columns))


if __name__ == "__main__":
    main()
