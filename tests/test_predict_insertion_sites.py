"""Regression tests for the insertion-site support metric.

`reads_supporting_insertion_pos` is used both to pick the winning breakpoint
position among competing candidates and, downstream, to decide whether a
region is even plotted (`make_bed_for_visualization.py`'s `min_support`
filter). It must reflect the number of distinct supporting reads — a plot
of the same region draws one bar per read, so a metric that undercounts
reads sharing an identical CIGAR string will disagree with what is visible
in the plot and can pick the wrong breakpoint entirely.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from telomererepeatloci.predict_insertion_sites import predict_insertions

CANDIDATE_COLUMNS = ["window", "chrom", "chromStart", "chromEnd", "strand"]
CLIPPED_COLUMNS = [
    "window",
    "read_name",
    "start",
    "end",
    "cigar",
    "part_telomere",
    "expected_pos_fusion",
    "TTAGGG_count",
    "CCCTAA_count",
]
DISCORDANT_COLUMNS = ["mate_chr", "mate_strand", "mate_position"]


def _write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def _run(tmp_path, candidate_rows, clipped_rows, discordant_rows):
    candidate_file = tmp_path / "candidates.tsv"
    clipped_file = tmp_path / "clipped.tsv"
    discordant_file = tmp_path / "discordant.tsv"
    _write_tsv(candidate_file, candidate_rows, CANDIDATE_COLUMNS)
    _write_tsv(clipped_file, clipped_rows, CLIPPED_COLUMNS)
    _write_tsv(discordant_file, discordant_rows, DISCORDANT_COLUMNS)

    output_df, _ = predict_insertions(
        str(candidate_file), str(clipped_file), str(discordant_file)
    )
    return output_df.iloc[0].to_dict()


def _clipped_read(read_name, end, cigar="40M10S"):
    return {
        "window": "1_1000_+",
        "read_name": read_name,
        "start": end - 40,
        "end": end,
        "cigar": cigar,
        "part_telomere": "True",
        "expected_pos_fusion": "downstream",
        "TTAGGG_count": 1,
        "CCCTAA_count": 0,
    }


def _discordant_read(position):
    return {"mate_chr": "1", "mate_strand": "+", "mate_position": position}


@pytest.fixture
def candidate_region():
    return [
        {
            "window": "1_1000_+",
            "chrom": "1",
            "chromStart": 1000,
            "chromEnd": 2000,
            "strand": "+",
        }
    ]


def test_reads_supporting_insertion_pos_counts_reads_not_cigars(
    tmp_path, candidate_region
):
    # 4 physically distinct reads happen to share the exact same CIGAR at the
    # breakpoint. A metric based on distinct CIGAR strings would report 1
    # here; the correct read-support count is 4.
    clipped_rows = [_clipped_read(f"read{i}", end=1600) for i in range(4)]
    discordant_rows = [_discordant_read(1300), _discordant_read(1300)]

    result = _run(tmp_path, candidate_region, clipped_rows, discordant_rows)

    assert result["insertion_site"] == "1600"
    assert result["reads_supporting_insertion_pos"] == "4"
    assert result["sum_TTAGGG_count"] == "4"


def test_prefers_position_with_more_supporting_reads_over_more_distinct_cigars(
    tmp_path, candidate_region
):
    # Position 1600 has 5 independent supporting reads that all happen to
    # share one CIGAR. Position 1700 has only 2 supporting reads, but they
    # happen to have different CIGARs. The winning site must be the one
    # with more actual supporting reads (1600), not the one with more
    # distinct CIGAR strings (1700) — otherwise the table would report a
    # lower-coverage breakpoint with a misleadingly low support count.
    well_supported = [_clipped_read(f"read{i}", end=1600) for i in range(5)]
    weakly_supported = [
        _clipped_read("readA", end=1700, cigar="40M9S"),
        _clipped_read("readB", end=1700, cigar="41M8S"),
    ]
    clipped_rows = well_supported + weakly_supported
    discordant_rows = [_discordant_read(1300), _discordant_read(1300)]

    result = _run(tmp_path, candidate_region, clipped_rows, discordant_rows)

    assert result["insertion_site"] == "1600"
    assert result["reads_supporting_insertion_pos"] == "5"
