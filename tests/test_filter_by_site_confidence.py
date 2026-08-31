"""Regression tests for the opt-in confidence-based filter.

A region with no predicted insertion_site is dropped outright -- it has no
locus to review or plot. For regions that do have one, blank/missing
diagnostic values (e.g. no control BAM) must never cause a drop -- only a
value that actually exceeds its threshold should. This mirrors how the rest
of the pipeline treats "" as "no signal", not "worst case".
"""

from __future__ import annotations

import pandas as pd

from telomererepeatloci.filter_by_site_confidence import (
    filter_regions,
    passes_confidence_filters,
)


def _row(**overrides):
    row = {
        "window": "1_1000_+",
        "insertion_site": "12345",
        "reads_supporting_insertion_pos": "2",
        "tumor_noise_ratio": "",
        "control_telo_clipped_at_insertion_site": "",
        "control_min_seq_distance_to_tumor": "",
    }
    row.update(overrides)
    return row


def test_keeps_region_with_missing_diagnostics():
    row = _row()
    assert passes_confidence_filters(row, 0.8, 2, 2, 2) is True


def test_drops_region_with_no_insertion_site():
    row = _row(insertion_site="")
    assert passes_confidence_filters(row, 0.8, 2, 2, 2) is False


def test_drops_region_below_min_insertion_support():
    # Support of 1 would never clear make_bed_for_visualization.py's default
    # --min-support of 2, so it can never be plotted -- drop it here too.
    row = _row(reads_supporting_insertion_pos="1")
    assert passes_confidence_filters(row, 0.8, 2, 2, 2) is False


def test_keeps_region_at_min_insertion_support_threshold():
    row = _row(reads_supporting_insertion_pos="2")
    assert passes_confidence_filters(row, 0.8, 2, 2, 2) is True


def test_drops_region_with_missing_insertion_support():
    row = _row(reads_supporting_insertion_pos="")
    assert passes_confidence_filters(row, 0.8, 2, 2, 2) is False


def test_drops_region_with_high_tumor_noise_ratio():
    row = _row(tumor_noise_ratio="0.9")
    assert (
        passes_confidence_filters(
            row,
            max_tumor_noise_ratio=0.8,
            control_max_seq_distance=2,
            control_max_telo_clipped_at_site=2,
            min_insertion_support=2,
        )
        is False
    )


def test_keeps_region_at_or_below_tumor_noise_ratio_threshold():
    row = _row(tumor_noise_ratio="0.8")
    assert passes_confidence_filters(row, 0.8, 2, 2, 2) is True


def test_drops_region_where_control_sequence_matches_tumor():
    # A single matching control read (low Hamming distance) is enough.
    row = _row(
        control_telo_clipped_at_insertion_site="1",
        control_min_seq_distance_to_tumor="0",
    )
    assert (
        passes_confidence_filters(
            row,
            0.8,
            control_max_seq_distance=2,
            control_max_telo_clipped_at_site=2,
            min_insertion_support=2,
        )
        is False
    )


def test_keeps_region_with_one_non_matching_control_read():
    # Present in control, but far from the tumor sequence and under the count
    # threshold -- not enough to call it germline on its own.
    row = _row(
        control_telo_clipped_at_insertion_site="1",
        control_min_seq_distance_to_tumor="6",
    )
    assert (
        passes_confidence_filters(
            row,
            0.8,
            control_max_seq_distance=2,
            control_max_telo_clipped_at_site=2,
            min_insertion_support=2,
        )
        is True
    )


def test_drops_region_with_too_many_control_telo_clipped_reads_even_without_match():
    row = _row(
        control_telo_clipped_at_insertion_site="3",
        control_min_seq_distance_to_tumor="6",
    )
    assert (
        passes_confidence_filters(
            row,
            0.8,
            control_max_seq_distance=2,
            control_max_telo_clipped_at_site=2,
            min_insertion_support=2,
        )
        is False
    )


def test_drops_clean_control_with_no_reads_at_site():
    # control_all_reads_at_site=0 is a populated value, not missing data --
    # a "clean" verdict from zero coverage tells you nothing.
    row = _row(control_all_reads_at_site="0")
    assert (
        passes_confidence_filters(
            row,
            max_tumor_noise_ratio=0.8,
            control_max_seq_distance=2,
            control_max_telo_clipped_at_site=2,
            min_insertion_support=2,
            min_control_reads_at_site=3,
        )
        is False
    )


def test_keeps_clean_control_at_min_control_reads_threshold():
    row = _row(control_all_reads_at_site="3")
    assert (
        passes_confidence_filters(
            row,
            max_tumor_noise_ratio=0.8,
            control_max_seq_distance=2,
            control_max_telo_clipped_at_site=2,
            min_insertion_support=2,
            min_control_reads_at_site=3,
        )
        is True
    )


def test_keeps_clean_control_with_missing_reads_at_site():
    # No control BAM at all -- blank stays non-disqualifying.
    row = _row(control_all_reads_at_site="")
    assert (
        passes_confidence_filters(
            row,
            max_tumor_noise_ratio=0.8,
            control_max_seq_distance=2,
            control_max_telo_clipped_at_site=2,
            min_insertion_support=2,
            min_control_reads_at_site=3,
        )
        is True
    )


def test_drops_region_with_tumor_reads_above_max_reads_at_site():
    row = _row(all_reads_at_site="1601")
    assert (
        passes_confidence_filters(
            row,
            max_tumor_noise_ratio=0.8,
            control_max_seq_distance=2,
            control_max_telo_clipped_at_site=2,
            min_insertion_support=2,
            max_reads_at_site=1600,
        )
        is False
    )


def test_keeps_region_at_max_reads_at_site_threshold():
    row = _row(all_reads_at_site="1600")
    assert (
        passes_confidence_filters(
            row,
            max_tumor_noise_ratio=0.8,
            control_max_seq_distance=2,
            control_max_telo_clipped_at_site=2,
            min_insertion_support=2,
            max_reads_at_site=1600,
        )
        is True
    )


def test_drops_region_with_control_reads_above_max_reads_at_site():
    row = _row(control_all_reads_at_site="1601")
    assert (
        passes_confidence_filters(
            row,
            max_tumor_noise_ratio=0.8,
            control_max_seq_distance=2,
            control_max_telo_clipped_at_site=2,
            min_insertion_support=2,
            max_reads_at_site=1600,
        )
        is False
    )


def test_filter_regions_reports_counts_and_preserves_columns(capsys):
    df = pd.DataFrame(
        [
            _row(window="a", tumor_noise_ratio="0.9"),
            _row(window="b", tumor_noise_ratio="0.1"),
        ]
    )

    result = filter_regions(df, 0.8, 2, 2, 2)

    assert list(result["window"]) == ["b"]
    assert list(result.columns) == list(df.columns)
    assert "kept 1 of 2" in capsys.readouterr().out
