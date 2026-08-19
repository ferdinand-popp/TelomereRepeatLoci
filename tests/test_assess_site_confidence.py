"""Regression tests for the site-level confidence/false-positive metrics.

These check the two subtleties that came up while designing this step:
1. A read's clip coordinate depends on its OWN clip direction (leading vs
   trailing clip), not on the region's predicted strand — reads_clipped_at_site
   must key off each row's own `expected_pos_fusion`.
2. Comparing clipped sequences between tumor and control must align on the
   base immediately at the breakpoint for BOTH directions: for a trailing
   clip (downstream) that's the first base of the clip; for a leading clip
   (upstream) the clip precedes the alignment, so the breakpoint-adjacent
   base is the LAST base of the raw clip substring, and it must be reversed
   before comparison — mirroring the orientation convention already used by
   get_consensus.py.
"""

from __future__ import annotations

import pysam
import pytest

from telomererepeatloci.assess_site_confidence import (
    count_reads_covering_site,
    junction_oriented_clip,
    min_hamming_distance,
    reads_clipped_at_site,
)


def _row(read_name, start, end, part_telomere, expected_pos_fusion, clipped_sequence):
    return {
        "read_name": read_name,
        "start": start,
        "end": end,
        "part_telomere": part_telomere,
        "expected_pos_fusion": expected_pos_fusion,
        "clipped_sequence": clipped_sequence,
    }


def test_reads_clipped_at_site_uses_each_reads_own_direction():
    # readA is clipped downstream (trailing clip) with its breakpoint at
    # reference position 1600 -> matches via its "end" column.
    # readB is clipped upstream (leading clip) with its breakpoint at 1600
    # too -> matches via its "start" column, not "end".
    rows = [
        _row(
            "readA",
            start=1560,
            end=1600,
            part_telomere="True",
            expected_pos_fusion="downstream",
            clipped_sequence="TTAGGGTTAGGG",
        ),
        _row(
            "readB",
            start=1600,
            end=1640,
            part_telomere="True",
            expected_pos_fusion="upstream",
            clipped_sequence="TTAGGGTTAGGG",
        ),
        # readC's trailing clip lands elsewhere -- must not match.
        _row(
            "readC",
            start=1560,
            end=1700,
            part_telomere="True",
            expected_pos_fusion="downstream",
            clipped_sequence="TTAGGGTTAGGG",
        ),
    ]

    matches = reads_clipped_at_site(rows, insertion_site=1600)

    assert set(matches) == {"readA", "readB"}


def test_reads_clipped_at_site_can_restrict_to_one_direction():
    rows = [
        _row(
            "readA",
            start=1560,
            end=1600,
            part_telomere="True",
            expected_pos_fusion="downstream",
            clipped_sequence="TTAGGGTTAGGG",
        ),
        _row(
            "readB",
            start=1600,
            end=1640,
            part_telomere="True",
            expected_pos_fusion="upstream",
            clipped_sequence="TTAGGGTTAGGG",
        ),
    ]

    matches = reads_clipped_at_site(
        rows, insertion_site=1600, require_direction="downstream"
    )

    assert set(matches) == {"readA"}


def test_junction_oriented_clip_downstream_takes_first_base_as_is():
    row = _row(
        "readA",
        start=1560,
        end=1600,
        part_telomere="True",
        expected_pos_fusion="downstream",
        clipped_sequence="ACGTAA",
    )
    assert junction_oriented_clip(row) == "ACGTAA"


def test_junction_oriented_clip_upstream_reverses_so_breakpoint_base_is_first():
    # Leading clip "AACGTA" precedes the alignment; the base immediately at
    # the breakpoint is the LAST character ("A"), so the oriented sequence
    # must start with that last character, i.e. be reversed.
    row = _row(
        "readB",
        start=1600,
        end=1640,
        part_telomere="True",
        expected_pos_fusion="upstream",
        clipped_sequence="AACGTA",
    )
    assert junction_oriented_clip(row) == "ATGCAA"


def test_min_hamming_distance_respects_orientation_not_raw_string_equality():
    # Both reads have the identical breakpoint-adjacent 12-mer once properly
    # oriented, even though their raw clipped_sequence strings differ because
    # one is a leading clip and one is a trailing clip.
    tumor_rows = {
        "t1": _row(
            "t1",
            start=1560,
            end=1600,
            part_telomere="True",
            expected_pos_fusion="downstream",
            clipped_sequence="TTAGGGCCCTAA",
        ),
    }
    control_rows = {
        "c1": _row(
            "c1",
            start=1600,
            end=1640,
            part_telomere="True",
            expected_pos_fusion="upstream",
            clipped_sequence="AATCCCGGGATT",
        ),
    }

    assert min_hamming_distance(tumor_rows, control_rows) == 0


def test_min_hamming_distance_none_when_either_side_empty():
    assert (
        min_hamming_distance(
            {}, {"c1": _row("c1", 0, 0, "True", "downstream", "A" * 12)}
        )
        is None
    )
    assert (
        min_hamming_distance(
            {"t1": _row("t1", 0, 0, "True", "downstream", "A" * 12)}, {}
        )
        is None
    )


def test_min_hamming_distance_skips_reads_shorter_than_12bp():
    tumor_rows = {
        "t1": _row("t1", 0, 0, "True", "downstream", "ACGT"),
    }
    control_rows = {
        "c1": _row("c1", 0, 0, "True", "downstream", "ACGT"),
    }
    assert min_hamming_distance(tumor_rows, control_rows) is None


@pytest.fixture
def synthetic_bam(tmp_path):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 5000, "SN": "chr1"}]}
    bam_path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:

        def make_read(name, pos, length, flag=0):
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "A" * length
            read.flag = flag
            read.reference_id = 0
            read.reference_start = pos
            read.mapping_quality = 60
            read.cigartuples = [(0, length)]
            return read

        # Covers site 1000 (999-1049).
        bam.write(make_read("covering1", 999, 50))
        # Also covers site 1000, but flagged as a duplicate -> must be excluded.
        bam.write(make_read("dup", 999, 50, flag=1024))
        # Does not cover site 1000 (ends at 1000, half-open -> exclusive).
        bam.write(make_read("adjacent", 950, 50))
        # Covers site 1000 as well.
        bam.write(make_read("covering2", 980, 40))

    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))
    return bam_path


def test_count_reads_covering_site_excludes_duplicates_and_non_covering_reads(
    synthetic_bam,
):
    with pysam.AlignmentFile(str(synthetic_bam), "rb") as bam:
        count = count_reads_covering_site(bam, "chr1", 1000, window=100)
    assert count == 2


def test_count_reads_covering_site_resolves_chr_prefix_alias(synthetic_bam):
    with pysam.AlignmentFile(str(synthetic_bam), "rb") as bam:
        # Table stores the bare chromosome name; BAM header uses "chr1".
        count = count_reads_covering_site(bam, "1", 1000, window=100)
    assert count == 2
