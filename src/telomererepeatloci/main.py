#!/usr/bin/env python3

import argparse
import concurrent.futures
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path


def _timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def parse_args():
    parser = argparse.ArgumentParser(description="Run the TelomereRepeatLoci workflow.")
    parser.add_argument("--tumor-bam", required=True, help="Required Tumor BAM file.")
    parser.add_argument(
        "--control-bam",
        default="",
        help="Optional control BAM file. If not provided, workflow runs in tumor-only mode.",
    )
    parser.add_argument(
        "--tel-tumor-bam",
        required=True,
        default="",
        help="Tumor BAM file to use for discordant-read screening.",
    )
    parser.add_argument(
        "--tel-control-bam",
        default="",
        help="Control BAM file to use for discordant-read screening.",
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help=(
            "Optional output directory. If not provided, a sibling directory named "
            "<telomerehunter-dir>_TelomereRepeatLoci is created next to the tumor "
            "TelomereHunter folder."
        ),
    )
    parser.add_argument("--tumor-sample-name", default="tumor")
    parser.add_argument("--control-sample-name", default="control")
    parser.add_argument("--blacklist", default="no_file")
    parser.add_argument("--tumor-discordant-read-lower-limit", type=float, default=3.0)
    parser.add_argument(
        "--control-discordant-read-upper-limit", type=float, default=0.0
    )
    parser.add_argument("--consider-blacklist", action="store_true")
    parser.add_argument("--reference-fasta", default="")
    parser.add_argument(
        "--skip-visualization",
        action="store_true",
        help="Skip generation of zoomed-in plots.",
    )
    parser.add_argument(
        "--plot-min-support",
        type=float,
        default=2.0,
        help=(
            "Minimum reads_supporting_insertion_pos required to include a region in "
            "plot BEDs. Default: 2."
        ),
    )
    parser.add_argument(
        "--site-window",
        type=int,
        default=100,
        help=(
            "Flank (bp) around insertion_site used to compute read depth and "
            "collect candidate clipped reads for the confidence-scoring step. "
            "Default: 100."
        ),
    )
    parser.add_argument(
        "--max-tumor-noise-ratio",
        type=float,
        default=0.8,
        help=(
            "Maximum allowed tumor_noise_ratio before a region is dropped by "
            "the site-confidence filtering step. Default: 0.8."
        ),
    )
    parser.add_argument(
        "--control-max-seq-distance",
        type=int,
        default=2,
        help=(
            "Maximum Hamming distance (over 12 bp at the breakpoint) between "
            "a control and tumor clipped sequence at the insertion site for "
            "it to still count as a germline match, used by the "
            "site-confidence filtering step. Default: 2."
        ),
    )
    parser.add_argument(
        "--control-max-telo-clipped-at-site",
        type=int,
        default=3,
        help=(
            "Maximum number of telomeric clipped reads allowed in control at "
            "the insertion site (regardless of sequence match) before a "
            "region is dropped by the site-confidence filtering step. "
            "Default: 3."
        ),
    )
    parser.add_argument(
        "--min-control-reads-at-site",
        type=int,
        default=3,
        help=(
            "Minimum control_all_reads_at_site required to trust a 'no "
            "control telomeric clips at the site' verdict as evidence of "
            "tumor specificity, used by the site-confidence filtering step. "
            "Below this, thin control coverage is dropped rather than kept "
            "on the strength of an uninformative clean result. Default: 3."
        ),
    )
    parser.add_argument(
        "--max-reads-at-site",
        type=int,
        default=1600,
        help=(
            "Maximum all_reads_at_site or control_all_reads_at_site allowed "
            "before a region is dropped as a collapsed-repeat/mapping-"
            "artifact pileup, used by the site-confidence filtering step. "
            "visualize_telomere_insertions.py silently skips plotting once "
            "either side's read count in the region reaches 3000, so "
            "candidates this deep would never render anyway. Default: 1600."
        ),
    )
    parser.add_argument("--samtoolsbin", default="samtools")

    args = parser.parse_args()

    # Build a dict of only arguments actually provided on the command line
    provided = {}
    for action in parser._actions:
        if action.dest == "help":
            continue
        value = getattr(args, action.dest)
        if value is not None:
            # Only include if it differs from default
            if value != action.default:
                provided[action.dest] = value

    print("Provided arguments:")
    for k, v in provided.items():
        print(f"{k} = {v}")

    return args


def run_command(command):
    start = time.monotonic()
    print(f"[{_timestamp()}] Running:", " ".join(command))
    subprocess.run(command, check=True)
    elapsed = time.monotonic() - start
    print(f"[{_timestamp()}] ---Done subprocess--- ({elapsed:.1f}s)")


def run_command_captured(command, label):
    """Like run_command, but buffers child stdout/stderr and prints them as one
    labeled block after the process finishes, instead of streaming them live.

    Used only for branches that run concurrently with another branch, where
    two processes writing to the same terminal at once would interleave their
    output line-by-line (or mid-line).
    """
    start = time.monotonic()
    print(f"[{_timestamp()}] [{label}] Running:", " ".join(command))
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    if result.stdout:
        print(f"[{label}] stdout:\n{result.stdout}", end="")
    if result.stderr:
        print(f"[{label}] stderr:\n{result.stderr}", end="")
    elapsed = time.monotonic() - start
    print(f"[{_timestamp()}] [{label}] ---Done subprocess--- ({elapsed:.1f}s)")


def run_concurrent_branches(branches):
    """Run independent (label, zero-arg callable) branches concurrently.

    Waits for all branches to finish before raising, and combines every
    failure into one error, since there's no cheap way to kill an in-flight
    subprocess from another thread and ThreadPoolExecutor already blocks
    until all workers finish on exit regardless.
    """
    errors = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(branches)) as executor:
        futures = {executor.submit(fn): label for label, fn in branches}
        for future in concurrent.futures.as_completed(futures):
            label = futures[future]
            try:
                future.result()
            except Exception as exc:
                errors.append((label, exc))
    if errors:
        raise RuntimeError(
            "; ".join(f"{label} branch failed: {exc}" for label, exc in errors)
        )


def get_filtered_bam(th_sample_dir):
    sample_dir = Path(th_sample_dir)
    preferred = sorted(sample_dir.glob("*_filtered_intratelomeric.bam"))
    if len(preferred) == 1:
        return preferred[0]
    if len(preferred) > 1:
        match_names = ", ".join(str(match) for match in preferred)
        raise ValueError(
            "Multiple *_filtered_intratelomeric.bam files found in "
            f"{sample_dir}: {match_names}"
        )

    matches = sorted(sample_dir.glob("*_filtered.bam"))
    if not matches:
        raise FileNotFoundError(
            f"No *_filtered_intratelomeric.bam or *_filtered.bam found in {sample_dir}"
        )
    if len(matches) > 1:
        match_names = ", ".join(str(match) for match in matches)
        raise ValueError(
            f"Multiple *_filtered.bam files found in {sample_dir}: {match_names}"
        )
    return matches[0]


def extract_pid_from_folder(folder_path):
    folder_name = Path(folder_path).name
    token = "_TelomerCnt_"
    if token not in folder_name:
        raise ValueError(
            f"Could not extract PID from folder name '{folder_name}'. "
            f"Expected pattern like '<sample>{token}<PID>'."
        )
    pid = folder_name.split(token, 1)[1]
    if not pid:
        raise ValueError(
            f"Could not extract PID from folder name '{folder_name}': empty PID."
        )
    return pid


def get_output_dir(args, tumor_th_dir):
    if args.output_dir:
        return Path(args.output_dir)
    telomerehunter_dir = tumor_th_dir.parent
    return (
        telomerehunter_dir.parent
        / f"{telomerehunter_dir.name.replace('_TelomerCnt', '')}_TelomereRepeatLoci"
    )


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def process_sample(args, scripts_dir):
    tumor_bam = Path(args.tumor_bam)
    if not tumor_bam.exists():
        raise FileNotFoundError(f"Missing tumor BAM: {tumor_bam}")

    tumor_filtered_bam = Path(args.tel_tumor_bam)
    if not tumor_filtered_bam.exists():
        raise FileNotFoundError(f"Missing telomeric tumor BAM: {tumor_filtered_bam}")

    tumor_th_dir = tumor_filtered_bam.parent
    pid = extract_pid_from_folder(tumor_th_dir)
    print(f"Detected telomere BAM dir         : {tumor_th_dir}")
    print(f"Detected PID                      : {pid}")

    use_control = bool(args.control_bam)
    control_bam = None
    control_filtered_bam = None

    if use_control:
        control_bam = Path(args.control_bam)
        if not control_bam.exists():
            raise FileNotFoundError(f"Missing control BAM: {control_bam}")

        control_filtered_bam = Path(args.tel_control_bam)
        if not control_filtered_bam.exists():
            raise FileNotFoundError(
                f"Missing telomeric control BAM: {control_filtered_bam}"
            )

        control_th_dir = control_filtered_bam.parent
        control_pid = extract_pid_from_folder(control_th_dir)
        print(f"Detected control telomere BAM dir : {control_th_dir}")

        if control_pid != pid:
            raise ValueError(
                "Tumor/control PID mismatch from folder names: "
                f"tumor PID='{pid}', control PID='{control_pid}'."
            )

    output_dir = get_output_dir(args, tumor_th_dir)

    tables_dir = output_dir / "tables"
    clipped_dir = output_dir / "clipped_reads"
    candidate_dir = output_dir / "candidate_region_tables"
    bed_zoomed_out_dir = output_dir / "plots" / "bedfiles" / "zoomed_out"
    bed_zoomed_in_dir = output_dir / "plots" / "bedfiles" / "zoomed_in"
    plot_zoomed_in_dir = output_dir / "plots" / "zoomed_in"

    for path in [
        tables_dir,
        clipped_dir,
        candidate_dir,
        bed_zoomed_out_dir,
        bed_zoomed_in_dir,
        plot_zoomed_in_dir,
    ]:
        ensure_dir(path)

    # Discordant reads: tumor and (if present) control are independent until
    # count_discordant_reads.py needs both, so run them concurrently when
    # there's a control branch to overlap with.
    def discordant_chain(filtered_bam, bam, discordant_path, mapq_path, run):
        run(
            [
                sys.executable,
                str(scripts_dir / "find_discordant_reads.py"),
                "-i",
                str(filtered_bam),
                "-o",
                str(discordant_path),
            ]
        )
        run(
            [
                sys.executable,
                str(scripts_dir / "add_mate_mapq.py"),
                "-i",
                str(discordant_path),
                "-b",
                str(bam),
                "-o",
                str(mapq_path),
            ]
        )

    tumor_discordant = (
        tables_dir / f"{pid}_{args.tumor_sample_name}_discordant_reads.tsv"
    )
    tumor_discordant_with_mapq = (
        tables_dir
        / f"{pid}_{args.tumor_sample_name}_discordant_reads_filtered_with_mapq.tsv"
    )
    control_discordant = (
        tables_dir / f"{pid}_{args.control_sample_name}_discordant_reads.tsv"
    )
    control_discordant_with_mapq = Path("NULL")

    if use_control:
        control_discordant_with_mapq = (
            tables_dir
            / f"{pid}_{args.control_sample_name}_discordant_reads_filtered_with_mapq.tsv"
        )
        run_concurrent_branches(
            [
                (
                    "tumor",
                    lambda: discordant_chain(
                        tumor_filtered_bam,
                        tumor_bam,
                        tumor_discordant,
                        tumor_discordant_with_mapq,
                        run=lambda cmd: run_command_captured(cmd, "tumor"),
                    ),
                ),
                (
                    "control",
                    lambda: discordant_chain(
                        control_filtered_bam,
                        control_bam,
                        control_discordant,
                        control_discordant_with_mapq,
                        run=lambda cmd: run_command_captured(cmd, "control"),
                    ),
                ),
            ]
        )
    else:
        discordant_chain(
            tumor_filtered_bam,
            tumor_bam,
            tumor_discordant,
            tumor_discordant_with_mapq,
            run=run_command,
        )

    windows = tables_dir / f"{pid}_discordant_reads_1_kb_windows.tsv"
    run_command(
        [
            sys.executable,
            str(scripts_dir / "count_discordant_reads.py"),
            "-t",
            str(tumor_discordant_with_mapq),
            "-c",
            str(control_discordant_with_mapq),
            "-b",
            args.blacklist,
            "-o",
            str(windows),
        ]
    )

    candidates = candidate_dir / f"{pid}_telomere_insertions_candidate_regions.tsv"
    run_command(
        [
            sys.executable,
            str(scripts_dir / "get_candidate_regions.py"),
            str(windows),
            str(candidates),
            str(args.tumor_discordant_read_lower_limit),
            str(args.control_discordant_read_upper_limit),
            str(args.consider_blacklist),
        ]
    )

    # Tumor-centric downstream steps. Control's find_fusion_reads.py (also
    # checked against the same candidate windows, so later steps can tell
    # whether control shows the same clipped-telomere signal) isn't needed
    # again until assess_site_confidence.py, so it can run concurrently with
    # tumor's whole find_fusion -> predict -> consensus chain.
    clipped = clipped_dir / f"{pid}_{args.tumor_sample_name}_clipped_reads.tsv"
    control_clipped = Path("NULL")
    extended = (
        candidate_dir / f"{pid}_telomere_insertions_candidate_regions_extended.tsv"
    )
    extended_with_consensus = (
        candidate_dir
        / f"{pid}_telomere_insertions_candidate_regions_extended_with_consensus.tsv"
    )

    def step_find_fusion(bam, out_path, run):
        run(
            [
                sys.executable,
                str(scripts_dir / "find_fusion_reads.py"),
                str(candidates),
                str(bam),
                str(out_path),
            ]
        )

    def step_predict(run):
        run(
            [
                sys.executable,
                str(scripts_dir / "predict_insertion_sites.py"),
                str(candidates),
                str(clipped),
                str(tumor_discordant_with_mapq),
                str(extended),
            ]
        )

    def step_consensus(run):
        cmd = [
            sys.executable,
            str(scripts_dir / "get_consensus.py"),
            str(extended),
            str(clipped),
            str(extended_with_consensus),
        ]
        if args.reference_fasta:
            cmd.extend(["--reference", args.reference_fasta])
        run(cmd)

    if use_control:
        control_clipped = (
            clipped_dir / f"{pid}_{args.control_sample_name}_clipped_reads.tsv"
        )

        def tumor_locus_chain():
            run_tumor = lambda cmd: run_command_captured(cmd, "tumor")  # noqa: E731
            step_find_fusion(tumor_bam, clipped, run=run_tumor)
            step_predict(run=run_tumor)
            step_consensus(run=run_tumor)

        def control_locus_chain():
            step_find_fusion(
                control_bam,
                control_clipped,
                run=lambda cmd: run_command_captured(cmd, "control"),
            )

        run_concurrent_branches(
            [("tumor", tumor_locus_chain), ("control", control_locus_chain)]
        )
    else:
        step_find_fusion(tumor_bam, clipped, run=run_command)
        step_predict(run=run_command)
        step_consensus(run=run_command)

    extended_with_confidence = (
        candidate_dir
        / f"{pid}_telomere_insertions_candidate_regions_extended_with_confidence.tsv"
    )
    run_command(
        [
            sys.executable,
            str(scripts_dir / "assess_site_confidence.py"),
            str(extended_with_consensus),
            str(clipped),
            str(control_clipped),
            str(tumor_bam),
            str(control_bam) if use_control else "NULL",
            str(extended_with_confidence),
            "--site-window",
            str(args.site_window),
        ]
    )

    extended_with_confidence_filtered = (
        candidate_dir
        / f"{pid}_telomere_insertions_candidate_regions_extended_with_confidence_filtered.tsv"
    )
    run_command(
        [
            sys.executable,
            str(scripts_dir / "filter_by_site_confidence.py"),
            str(extended_with_confidence),
            str(extended_with_confidence_filtered),
            "--max-tumor-noise-ratio",
            str(args.max_tumor_noise_ratio),
            "--control-max-seq-distance",
            str(args.control_max_seq_distance),
            "--control-max-telo-clipped-at-site",
            str(args.control_max_telo_clipped_at_site),
            "--min-insertion-support",
            str(args.plot_min_support),
            "--min-control-reads-at-site",
            str(args.min_control_reads_at_site),
            "--max-reads-at-site",
            str(args.max_reads_at_site),
        ]
    )

    bed_zoomed_out = bed_zoomed_out_dir / f"{pid}_telomere_insertions.bed"
    bed_zoomed_in = bed_zoomed_in_dir / f"{pid}_telomere_insertions.bed"
    run_command(
        [
            sys.executable,
            str(scripts_dir / "make_bed_for_visualization.py"),
            str(extended_with_confidence_filtered),
            str(bed_zoomed_out),
            str(bed_zoomed_in),
            pid,
            "--min-support",
            str(args.plot_min_support),
        ]
    )

    if not args.skip_visualization:
        visualize_cmd = [
            sys.executable,
            str(scripts_dir / "visualize_telomere_insertions.py"),
            "--tumor",
            str(tumor_bam),
            "--ref",
            args.reference_fasta,
            "--bed",
            str(bed_zoomed_in),
            "--samtoolsbin",
            args.samtoolsbin,
            "--colored_reads_tumor",
            str(tumor_discordant_with_mapq),
            "--clipped_reads_tumor",
            str(clipped),
            "--prefix",
            f"{plot_zoomed_in_dir}/",
            "--outfile",
            str(plot_zoomed_in_dir / f"{pid}_done.txt"),
        ]
        if use_control:
            visualize_cmd.extend(
                [
                    "--control",
                    str(control_bam),
                    "--colored_reads_control",
                    str(control_discordant_with_mapq),
                ]
            )
        run_command(visualize_cmd)


def get_version_from_package():
    try:
        from telomererepeatloci.version import __version__

        return __version__
    except ImportError:
        return "unknown - check pyproject.toml file"


def main():
    print(f"TelomereRepeatLoci - version {get_version_from_package()}")
    args = parse_args()
    scripts_dir = Path(__file__).resolve().parent
    start = time.monotonic()
    print(f"[{_timestamp()}] --- Processing sample ---")
    process_sample(args, scripts_dir)
    elapsed = time.monotonic() - start
    print(f"[{_timestamp()}] --- Done processing sample --- ({elapsed:.1f}s)")


if __name__ == "__main__":
    main()
