"""Regression tests for the opt-in confidence-based filter.

Blank/missing diagnostic values (e.g. no control BAM, or no insertion_site)
must never cause a drop -- only a value that actually exceeds its threshold
should. This mirrors how the rest of the pipeline treats "" as "no signal",
not "worst case".
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
        "tumor_noise_ratio": "",
        "control_telo_clipped_at_insertion_site": "",
        "control_min_seq_distance_to_tumor": "",
    }
    row.update(overrides)
    return row


def test_keeps_region_with_missing_diagnostics():
    row = _row()
    assert passes_confidence_filters(row, 0.8, 2, 2) is True


def test_drops_region_with_high_tumor_noise_ratio():
    row = _row(tumor_noise_ratio="0.9")
    assert (
        passes_confidence_filters(
            row,
            max_tumor_noise_ratio=0.8,
            control_max_seq_distance=2,
            control_max_telo_clipped_at_site=2,
        )
        is False
    )


def test_keeps_region_at_or_below_tumor_noise_ratio_threshold():
    row = _row(tumor_noise_ratio="0.8")
    assert passes_confidence_filters(row, 0.8, 2, 2) is True


def test_drops_region_where_control_sequence_matches_tumor():
    # A single matching control read (low Hamming distance) is enough.
    row = _row(
        control_telo_clipped_at_insertion_site="1",
        control_min_seq_distance_to_tumor="0",
    )
    assert (
        passes_confidence_filters(
            row, 0.8, control_max_seq_distance=2, control_max_telo_clipped_at_site=2
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
            row, 0.8, control_max_seq_distance=2, control_max_telo_clipped_at_site=2
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
            row, 0.8, control_max_seq_distance=2, control_max_telo_clipped_at_site=2
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

    result = filter_regions(df, 0.8, 2, 2)

    assert list(result["window"]) == ["b"]
    assert list(result.columns) == list(df.columns)
    assert "kept 1 of 2" in capsys.readouterr().out
