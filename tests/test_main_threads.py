"""Regression tests for concurrent tumor/control branch execution in main.py.

process_sample() runs the tumor and (when --control-bam is given) control
branches of the discordant-read and clipped-read steps concurrently, since
neither branch needs the other's output until a later joining step. These
tests monkeypatch subprocess.run with a fake that records call timing instead
of actually invoking the pipeline scripts, so they can assert on ordering and
overlap without needing real BAM files or scripts.
"""

from __future__ import annotations

import argparse
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import pytest

from telomererepeatloci import main


@dataclass
class Call:
    command: list
    start: float
    end: float


def _make_fake_run(calls, delay=0.03, fail_matcher=None):
    def fake_run(command, check=False, capture_output=False, text=False):
        start = time.perf_counter()
        time.sleep(delay)
        end = time.perf_counter()
        calls.append(Call(command=list(command), start=start, end=end))
        if fail_matcher is not None and fail_matcher(command):
            raise subprocess.CalledProcessError(
                returncode=1, cmd=command, stderr="boom"
            )
        return subprocess.CompletedProcess(
            args=command, returncode=0, stdout="", stderr=""
        )

    return fake_run


def _script_name(call):
    return Path(call.command[1]).name


def _is_call(command, script_name):
    return Path(command[1]).name == script_name


def _call_for(calls, script_name, contains):
    matches = [
        c for c in calls if _script_name(c) == script_name and contains in c.command
    ]
    assert len(matches) == 1, (
        f"expected exactly one {script_name} call containing {contains!r}, "
        f"got {len(matches)}"
    )
    return matches[0]


def _make_bam(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()
    return path


def _base_args(tmp_path, use_control=False):
    tel_tumor_bam = _make_bam(
        tmp_path / "tumorSample_TelomerCnt_PID1" / "tumor_filtered_intratelomeric.bam"
    )
    tumor_bam = _make_bam(tmp_path / "tumor.bam")

    args = argparse.Namespace(
        tumor_bam=str(tumor_bam),
        tel_tumor_bam=str(tel_tumor_bam),
        control_bam="",
        tel_control_bam="",
        output_dir=str(tmp_path / "out"),
        tumor_sample_name="tumor",
        control_sample_name="control",
        blacklist="no_file",
        tumor_discordant_read_lower_limit=3.0,
        control_discordant_read_upper_limit=0.0,
        consider_blacklist=False,
        reference_fasta="",
        skip_visualization=True,
        plot_min_support=2.0,
        site_window=100,
        max_tumor_noise_ratio=0.8,
        control_max_seq_distance=2,
        control_max_telo_clipped_at_site=2,
        min_control_reads_at_site=3,
        max_reads_at_site=1600,
        samtoolsbin="samtools",
    )

    if use_control:
        tel_control_bam = _make_bam(
            tmp_path
            / "controlSample_TelomerCnt_PID1"
            / "control_filtered_intratelomeric.bam"
        )
        control_bam = _make_bam(tmp_path / "control.bam")
        args.control_bam = str(control_bam)
        args.tel_control_bam = str(tel_control_bam)

    return args


def test_tumor_only_runs_sequentially_in_expected_order(tmp_path, monkeypatch):
    args = _base_args(tmp_path, use_control=False)
    calls = []
    monkeypatch.setattr(main.subprocess, "run", _make_fake_run(calls))

    main.process_sample(args, Path("/fake/scripts"))

    assert [_script_name(c) for c in calls] == [
        "find_discordant_reads.py",
        "add_mate_mapq.py",
        "count_discordant_reads.py",
        "get_candidate_regions.py",
        "find_fusion_reads.py",
        "predict_insertion_sites.py",
        "get_consensus.py",
        "assess_site_confidence.py",
        "filter_by_site_confidence.py",
        "make_bed_for_visualization.py",
    ]
    for earlier, later in zip(calls, calls[1:]):
        assert earlier.end <= later.start


def test_control_branches_run_concurrently(tmp_path, monkeypatch):
    args = _base_args(tmp_path, use_control=True)
    calls = []
    monkeypatch.setattr(main.subprocess, "run", _make_fake_run(calls, delay=0.05))

    main.process_sample(args, Path("/fake/scripts"))

    tumor_discordant = _call_for(calls, "find_discordant_reads.py", args.tel_tumor_bam)
    control_discordant = _call_for(
        calls, "find_discordant_reads.py", args.tel_control_bam
    )
    assert tumor_discordant.start < control_discordant.end
    assert control_discordant.start < tumor_discordant.end

    tumor_fusion = _call_for(calls, "find_fusion_reads.py", args.tumor_bam)
    control_fusion = _call_for(calls, "find_fusion_reads.py", args.control_bam)
    assert tumor_fusion.start < control_fusion.end
    assert control_fusion.start < tumor_fusion.end


def test_error_in_one_branch_waits_for_the_other_then_raises(tmp_path, monkeypatch):
    args = _base_args(tmp_path, use_control=True)
    calls = []

    def fail_matcher(command):
        return (
            _is_call(command, "find_discordant_reads.py")
            and args.tel_control_bam in command
        )

    monkeypatch.setattr(
        main.subprocess,
        "run",
        _make_fake_run(calls, delay=0.03, fail_matcher=fail_matcher),
    )

    with pytest.raises(RuntimeError, match="control branch failed"):
        main.process_sample(args, Path("/fake/scripts"))

    tumor_mapq_calls = [
        c
        for c in calls
        if _script_name(c) == "add_mate_mapq.py" and args.tumor_bam in c.command
    ]
    assert tumor_mapq_calls, (
        "tumor branch should have completed its discordant chain despite the "
        "control branch failing"
    )
