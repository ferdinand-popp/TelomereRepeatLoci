"""Regression tests for find_fusion_reads.py after merging the soft-clip and
supplementary-alignment passes into a single bam.fetch() per region, and
caching primary-alignment lookups by locus instead of re-fetching per
supplementary read.

The caching in particular is correctness-sensitive: two different reads whose
primary alignments happen to sit at the exact same locus must still each get
their OWN primary sequence back, not a swapped/cross-contaminated one.
"""

from __future__ import annotations

import pandas as pd
import pysam
import pytest

from telomererepeatloci.find_fusion_reads import find_fusion_reads


def _make_bam(path):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 10000}]}

    def _segment(name, pos, cigar, seq, flag):
        seg = pysam.AlignedSegment()
        seg.query_name = name
        seg.query_sequence = seq
        seg.flag = flag
        seg.reference_id = 0
        seg.reference_start = pos
        seg.mapping_quality = 60
        seg.cigarstring = cigar
        seg.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        return seg

    # Two distinct primary alignments that happen to share the exact same
    # locus -- this is what exercises the per-locus cache.
    primary_a = _segment("readA", 1000, "50M", "A" * 50, flag=0)
    primary_b = _segment("readB", 1000, "50M", "C" * 50, flag=0)

    # Supplementary alignments for each, elsewhere in the genome, each
    # pointing back (via SA tag) at the same 1000bp locus above.
    supp_a = _segment("readA", 2000, "20S30M", "N" * 50, flag=0x800)
    supp_a.set_tag("SA", "chr1,1001,+,50M,60,0;")
    supp_b = _segment("readB", 2050, "20S30M", "N" * 50, flag=0x800)
    supp_b.set_tag("SA", "chr1,1001,+,50M,60,0;")

    # A plain soft-clipped (non-supplementary) read, to confirm the merged
    # single-pass loop still emits a soft-clip row alongside the
    # supplementary rows above.
    plain_clip = _segment("readC", 2100, "40M10S", "G" * 40 + "T" * 10, flag=0)

    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for seg in (primary_a, primary_b, supp_a, supp_b, plain_clip):
            bam.write(seg)
    pysam.index(str(path))


def test_supplementary_rows_use_their_own_primary_sequence_not_a_cached_swap(tmp_path):
    bam_path = tmp_path / "reads.bam"
    _make_bam(bam_path)

    candidates = pd.DataFrame(
        [{"window": "w1", "chrom": "chr1", "chromStart": 1900, "chromEnd": 2200}]
    )
    candidate_path = tmp_path / "candidates.tsv"
    candidates.to_csv(candidate_path, sep="\t", index=False)

    result = find_fusion_reads(str(candidate_path), str(bam_path))

    supp_rows = result[result["chr_primary_align"] == "chr1"]
    assert set(supp_rows["read_name"]) == {"readA", "readB"}

    row_a = supp_rows[supp_rows["read_name"] == "readA"].iloc[0]
    row_b = supp_rows[supp_rows["read_name"] == "readB"].iloc[0]

    # readA's supplementary alignment must recover primary_a's sequence
    # ("A"*50), readB's must recover primary_b's ("C"*50) -- not swapped.
    assert row_a["sequence"] == "A" * 50
    assert row_b["sequence"] == "C" * 50
    assert row_a["coord_primary_align"] == 1000
    assert row_a["strand_primary_align"] == "+"

    # The plain soft-clipped read still produces its own row via the same
    # single-pass loop.
    plain_rows = result[result["read_name"] == "readC"]
    assert len(plain_rows) == 1
    assert plain_rows.iloc[0]["clipped_sequence"] == "T" * 10
    assert plain_rows.iloc[0]["chr_primary_align"] == ""

    # Supplementary reads with their own soft clip also still emit a plain
    # soft-clip row in addition to the supplementary row (unchanged from
    # before the merge).
    own_clip_rows = result[
        (result["read_name"].isin(["readA", "readB"]))
        & (result["chr_primary_align"] == "")
    ]
    assert len(own_clip_rows) == 2


def test_a_supplementary_read_does_not_truncate_later_reads_in_the_same_region(
    tmp_path,
):
    """pysam does not support two concurrent fetch() iterators on the same
    AlignmentFile handle. get_primary_sequence()'s lookup for a
    supplementary read used to run on the SAME handle driving the outer
    per-region scan -- a nested bam.fetch() call while that outer scan was
    still mid-iteration silently corrupted/reset it, so every soft-clipped
    read positioned *after* the first supplementary read in a region was
    dropped. Regression-tested with a supplementary read early in the region
    followed by many plain soft-clipped reads.
    """
    bam_path = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 200000}]}

    def _segment(name, pos, cigar, seq, flag):
        seg = pysam.AlignedSegment()
        seg.query_name = name
        seg.query_sequence = seq
        seg.flag = flag
        seg.reference_id = 0
        seg.reference_start = pos
        seg.mapping_quality = 60
        seg.cigarstring = cigar
        seg.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        return seg

    # Supplementary read early in the region, pointing (via SA tag) at a
    # primary locus far away -- the seek distance is what actually disturbs
    # the outer iterator; a nearby primary locus doesn't reliably trigger
    # the corruption in a small synthetic BAM.
    supp = _segment("suppRead", 2000, "20S30M", "N" * 50, flag=0x800)
    supp.set_tag("SA", "chr1,100001,+,50M,60,0;")
    primary = _segment("primaryRead", 100000, "50M", "A" * 50, flag=0)

    # Plain soft-clipped reads AFTER the supplementary read in position
    # order -- these are exactly what got silently dropped.
    later_clips = [
        _segment(f"laterRead{i}", 2010 + i * 10, "40M10S", "G" * 40 + "T" * 10, 0)
        for i in range(5)
    ]

    # Reads must be written in position order for BAM indexing.
    reads_by_position = sorted(
        [supp, primary, *later_clips], key=lambda r: r.reference_start
    )
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for seg in reads_by_position:
            bam.write(seg)
    pysam.index(str(bam_path))

    candidates = pd.DataFrame(
        [{"window": "w1", "chrom": "chr1", "chromStart": 1900, "chromEnd": 2100}]
    )
    candidate_path = tmp_path / "candidates.tsv"
    candidates.to_csv(candidate_path, sep="\t", index=False)

    result = find_fusion_reads(str(candidate_path), str(bam_path))

    later_rows = result[result["read_name"].str.startswith("laterRead")]
    assert len(later_rows) == 5, (
        "all 5 soft-clipped reads after the supplementary read must still "
        "be found, not truncated by a corrupted outer fetch() iterator"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
