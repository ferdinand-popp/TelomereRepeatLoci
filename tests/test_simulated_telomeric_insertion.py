from __future__ import annotations

import csv
import os
import stat
import subprocess
import sys
from pathlib import Path
import shutil

import pysam


def _write_alignment_bam(path: Path) -> None:
    header = {
        "HD": {"VN": "1.0", "SO": "coordinate"},
        "SQ": [{"SN": "1", "LN": 100000}],
    }

    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        # Add non-telomeric reads
        for read_name, pos in [
            ("read1/1", 1201),
            ("read1/2", 1251),
            ("read2/1", 1301),
            ("read2/2", 1351),
        ]:
            read = pysam.AlignedSegment()
            read.query_name = read_name.split("/")[0]
            read.flag = 99 if "/1" in read_name else 147  # Properly paired flags
            read.reference_id = 0
            read.reference_start = pos - 1
            read.mapping_quality = 60
            read.cigarstring = "50M"
            read.query_sequence = "A" * 50
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            read.next_reference_id = 0
            read.next_reference_start = pos + 50 if "/1" in read_name else pos - 50
            read.template_length = 100 if "/1" in read_name else -100
            bam.write(read)

        # Add telomeric reads
        telomeric_reads = [
            ("tel1/1", 1400, "40M10S", "TTAGGGTTAG"),
            ("tel1/2", 1450, "50M", ""),
            ("tel2/1", 1500, "35M15S", "TTAGGGTTAGGGTTA"),
            ("tel2/2", 1550, "50M", ""),
            ("tel3/1", 1600, "30M20S", "TTAGGGTTAGGGTTAGGGTT"),
            ("tel3/2", 1650, "50M", ""),
        ]
        for read_name, start0, cigar, tel_clip in telomeric_reads:
            match_len = int(cigar.split("M", 1)[0])
            read = pysam.AlignedSegment()
            read.query_name = read_name.split("/")[0]
            read.flag = 99 if "/1" in read_name else 147  # Properly paired flags
            read.reference_id = 0
            read.reference_start = start0
            read.mapping_quality = 60
            read.cigarstring = cigar
            read.query_sequence = "A" * match_len + tel_clip
            read.query_qualities = pysam.qualitystring_to_array(
                "I" * (match_len + len(tel_clip))
            )
            read.next_reference_id = 0
            read.next_reference_start = (
                start0 + 50 if "/1" in read_name else start0 - 50
            )
            read.template_length = 100 if "/1" in read_name else -100
            bam.write(read)

        # Add discordant reads (mapped mate and unmapped mate with telomeric motif)
        discordant_reads = [
            (
                "disc1",
                1700,
                "50M",
                "TTAGGGTTAGGGTTAGGG",
            ),  # Mapped mate and unmapped mate
            ("disc2", 1800, "50M", "TTAGGGTTAGGG"),  # Mapped mate and unmapped mate
        ]
        for read_name, start0, cigar, tel_clip in discordant_reads:
            # Mapped mate
            mapped_read = pysam.AlignedSegment()
            mapped_read.query_name = read_name
            mapped_read.flag = 99  # Properly paired, first in pair
            mapped_read.reference_id = 0
            mapped_read.reference_start = start0
            mapped_read.mapping_quality = 60
            mapped_read.cigarstring = cigar
            mapped_read.query_sequence = "A" * 50
            mapped_read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            mapped_read.next_reference_id = -1  # Mate is unmapped
            mapped_read.next_reference_start = 0
            mapped_read.template_length = 0
            bam.write(mapped_read)

            # Unmapped mate with telomeric motif
            unmapped_read = pysam.AlignedSegment()
            unmapped_read.query_name = read_name
            unmapped_read.flag = 141  # Unmapped, second in pair
            unmapped_read.reference_id = -1
            unmapped_read.reference_start = -1
            unmapped_read.mapping_quality = 0
            unmapped_read.cigarstring = None
            unmapped_read.query_sequence = tel_clip
            unmapped_read.query_qualities = pysam.qualitystring_to_array(
                "I" * len(tel_clip)
            )
            unmapped_read.next_reference_id = 0  # Mate is mapped
            unmapped_read.next_reference_start = start0
            unmapped_read.template_length = 0
            bam.write(unmapped_read)

        # Add unmapped reads (must be written last)
        discordant_telomeric_reads = [
            ("disc1", 1, "TTAGGGTTAGGGTTAGGG"),
            ("disc2", 1, "TTAGGGTTAGGG"),
        ]
        for read_name, reference, tel_clip in discordant_telomeric_reads:
            read = pysam.AlignedSegment()
            read.query_name = read_name
            read.flag = 141  # Unmapped, second in pair
            read.reference_id = reference
            read.reference_start = reference
            read.mapping_quality = 0
            read.cigarstring = None
            read.query_sequence = tel_clip
            read.query_qualities = pysam.qualitystring_to_array("I" * len(tel_clip))
            read.next_reference_id = 0  # Reference ID of the mapped mate
            read.next_reference_start = 1700 if "disc1" in read_name else 1800
            read.template_length = 0
            bam.write(read)

    pysam.index(str(path))


def _write_filtered_bam(path: Path) -> None:
    header = {
        "HD": {"VN": "1.0", "SO": "coordinate"},
        "SQ": [{"SN": "1", "LN": 100000}],
    }

    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        # Include only unmapped telomeric reads from discordant pairs
        discordant_telomeric_reads = [
            ("disc1", 1, "TTAGGGTTAGGGTTAGGG"),
            ("disc2", 1, "TTAGGGTTAGGG"),
        ]
        for read_name, reference, tel_clip in discordant_telomeric_reads:
            read = pysam.AlignedSegment()
            read.query_name = read_name
            read.flag = 141  # Unmapped, second in pair
            read.reference_id = reference
            read.reference_start = reference
            read.mapping_quality = 0
            read.cigarstring = None
            read.query_sequence = tel_clip
            read.query_qualities = pysam.qualitystring_to_array("I" * len(tel_clip))
            read.next_reference_id = 0  # Reference ID of the mapped mate
            read.next_reference_start = 1700 if "disc1" in read_name else 1800
            read.template_length = 0
            bam.write(read)

    pysam.index(str(path))


def _write_samtools_stub(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python3
import sys
import pysam


def main():
    args = sys.argv[1:]
    if not args or args[0] != 'view':
        raise SystemExit(2)

    exclude = 0
    i = 1
    while i < len(args) and args[i].startswith('-'):
        if args[i] == '-F':
            exclude = int(args[i + 1])
            i += 2
        else:
            i += 1

    bam_path = args[i]
    region = args[i + 1] if i + 1 < len(args) else None

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        if region:
            chrom, span = region.split(':', 1)
            start_s, end_s = span.split('-', 1)
            start = max(0, int(start_s) - 1)
            end = int(end_s)
            reads = bam.fetch(chrom, start, end)
        else:
            reads = bam.fetch(until_eof=True)

        for read in reads:
            if exclude and (read.flag & exclude):
                continue
            fields = [
                read.query_name,
                str(read.flag),
                read.reference_name or '*',
                str(read.reference_start + 1),
                str(read.mapping_quality),
                read.cigarstring or '*',
                '*',
                '0',
                '0',
                read.query_sequence or '*',
                '*',
            ]
            print('\t'.join(fields))


if __name__ == '__main__':
    main()
"""
    )
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def test_pipeline_with_simulated_discordant_reads() -> None:
    tmp_path = Path("/Users/ferdinandpopp/PycharmProjects/TelomereRepeatLoci/debug")
    repo_root = Path(__file__).resolve().parents[1]
    pid = "PID001"

    telomerehunter_dir = tmp_path / "telomerehunter_output"
    telomerehunter_dir.mkdir(parents=True, exist_ok=True)

    # Create the TelomereHunter sample directory with the exact token expected by main.py
    th_sample_dir = telomerehunter_dir / f"tumor_TelomerCnt_{pid}"
    th_sample_dir.mkdir(parents=True, exist_ok=True)

    alignment_bam = tmp_path / "tumor_input.bam"
    # place the filtered BAM inside the TelomereHunter sample dir and name it "<PID>_filtered.bam"
    filtered_bam = th_sample_dir / f"{pid}_filtered.bam"

    _write_alignment_bam(alignment_bam)
    _write_filtered_bam(filtered_bam)

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    samtools_stub = bin_dir / "samtools"
    _write_samtools_stub(samtools_stub)

    env = dict(os.environ)
    env["PATH"] = f"{bin_dir}:{env.get('PATH', '')}"

    try:
        subprocess.run(
            [
                sys.executable,
                str(repo_root / "src" / "main.py"),
                "--tumor-bam",
                str(alignment_bam),
                "--telomerehunter-dir",
                str(telomerehunter_dir),
            ],
            check=True,
            cwd=repo_root,
            env=env,
        )

    finally:
        debug_root = repo_root / "debug" / "test_simulated_telomeric_insertion"
        debug_root.mkdir(parents=True, exist_ok=True)

        output_dir = (
            telomerehunter_dir.parent / f"{telomerehunter_dir.name}_TelomereRepeatLoci"
        )

        # windows_file = (
        #     output_dir / "tables" / f"{pid}_discordant_reads_1_kb_windows.tsv"
        # )
        # candidates_file = (
        #     output_dir
        #     / "candidate_region_tables"
        #     / f"{pid}_telomere_insertions_candidate_regions.tsv"
        # )
        # extended_file = (
        #     output_dir
        #     / "candidate_region_tables"
        #     / f"{pid}_telomere_insertions_candidate_regions_extended.tsv"
        # )
        # consensus_file = (
        #     output_dir
        #     / "candidate_region_tables"
        #     / f"{pid}_telomere_insertions_candidate_regions_extended_with_consensus.tsv"
        # )
        # bed_zoomed_out = (
        #     output_dir
        #     / "plots"
        #     / "bedfiles"
        #     / "zoomed_out"
        #     / f"{pid}_telomere_insertions.bed"
        # )

        # Persist the entire output directory into debug_root for inspection
        if output_dir.exists():
            dest_dir = debug_root / output_dir.name
            if dest_dir.exists():
                shutil.rmtree(dest_dir)
            try:
                # Python 3.8+ supports dirs_exist_ok
                shutil.copytree(output_dir, dest_dir, dirs_exist_ok=True)
            except TypeError:
                # Fallback for older Python versions
                shutil.copytree(output_dir, dest_dir)
        else:
            # Ensure debug root exists even if no output produced
            debug_root.mkdir(parents=True, exist_ok=True)

        # Use the persisted copy for assertions/reads
        copied_base = debug_root / output_dir.name
        copied_windows = (
            copied_base / "tables" / f"{pid}_discordant_reads_1_kb_windows.tsv"
        )
        copied_candidates = (
            copied_base
            / "candidate_region_tables"
            / f"{pid}_telomere_insertions_candidate_regions.tsv"
        )
        copied_extended = (
            copied_base
            / "candidate_region_tables"
            / f"{pid}_telomere_insertions_candidate_regions_extended.tsv"
        )
        copied_consensus = (
            copied_base
            / "candidate_region_tables"
            / f"{pid}_telomere_insertions_candidate_regions_extended_with_consensus.tsv"
        )
        copied_bed_zoomed_out = (
            copied_base
            / "plots"
            / "bedfiles"
            / "zoomed_out"
            / f"{pid}_telomere_insertions.bed"
        )

        for path in [
            copied_windows,
            copied_candidates,
            copied_extended,
            copied_consensus,
            copied_bed_zoomed_out,
        ]:
            assert path.exists(), f"Expected persisted file {path} to exist"

        with copied_windows.open(newline="") as handle:
            windows_rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(windows_rows) == 1
        assert windows_rows[0]["window"] == "1_1000_+"
        assert windows_rows[0]["tumor_discordant_read_count"] == "3"

        with copied_extended.open(newline="") as handle:
            extended_rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(extended_rows) == 1
        assert extended_rows[0]["insertion_site"] == "1500"
        assert extended_rows[0]["reads_supporting_insertion_pos"] == "3"

        with copied_consensus.open(newline="") as handle:
            consensus_rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(consensus_rows) == 1
        assert consensus_rows[0]["consensus"]
