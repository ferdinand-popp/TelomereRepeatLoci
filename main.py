#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the TelomereRepeatLoci workflow without Snakemake/YAML."
    )
    parser.add_argument("--tumor-bam", required=True)
    parser.add_argument("--control-bam", required=True)
    parser.add_argument(
        "--telomerehunter-dir",
        required=True,
        help=(
            "Parent directory containing TelomereHunter output folders, e.g. "
            "<telomerehunter-dir>/tumor_TelomerCnt_<PID> and "
            "<telomerehunter-dir>/control_TelomerCnt_<PID>."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help=(
            "Optional output directory. If not provided, a sibling directory named "
            "<telomerehunter-dir>_TelomereRepeatLoci is created."
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
        "--run-visualization",
        action="store_true",
        help="Generate zoomed-in plots.",
    )
    parser.add_argument("--samtoolsbin", default="samtools")
    return parser.parse_args()


def run_command(command):
    print("Running:", " ".join(command))
    subprocess.run(command, check=True)
    print("---Done subprocess---")


def get_filtered_bam(th_sample_dir):
    matches = sorted(Path(th_sample_dir).glob("*_filtered.bam"))
    if not matches:
        raise FileNotFoundError(f"No *_filtered.bam found in {th_sample_dir}")
    if len(matches) > 1:
        match_names = ", ".join(str(match) for match in matches)
        raise ValueError(
            f"Multiple *_filtered.bam files found in {th_sample_dir}: {match_names}"
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


def detect_th_sample_dirs(telomerehunter_dir):
    base = Path(telomerehunter_dir)
    if not base.exists():
        raise FileNotFoundError(f"Missing telomerehunter-dir: {base}")
    if not base.is_dir():
        raise NotADirectoryError(f"--telomerehunter-dir is not a directory: {base}")

    tumor_dirs = []
    control_dirs = []

    for child in base.iterdir():
        if not child.is_dir():
            continue
        name_lower = child.name.lower()
        if "tumor_telomercnt_" in name_lower:
            tumor_dirs.append(child)
        elif "control_telomercnt_" in name_lower:
            control_dirs.append(child)

    if len(tumor_dirs) != 1:
        names = ", ".join(str(p) for p in sorted(tumor_dirs)) or "none"
        raise ValueError(
            f"Expected exactly 1 tumor TelomereHunter folder under {base}, "
            f"found {len(tumor_dirs)}: {names}"
        )
    if len(control_dirs) != 1:
        names = ", ".join(str(p) for p in sorted(control_dirs)) or "none"
        raise ValueError(
            f"Expected exactly 1 control TelomereHunter folder under {base}, "
            f"found {len(control_dirs)}: {names}"
        )

    tumor_dir = tumor_dirs[0]
    control_dir = control_dirs[0]

    tumor_pid = extract_pid_from_folder(tumor_dir)
    control_pid = extract_pid_from_folder(control_dir)

    if tumor_pid != control_pid:
        raise ValueError(
            "Tumor/control PID mismatch from TelomereHunter folder names: "
            f"tumor PID='{tumor_pid}', control PID='{control_pid}'."
        )

    return tumor_dir, control_dir, tumor_pid


def get_output_dir(args):
    if args.output_dir:
        return Path(args.output_dir)
    telomerehunter_dir = Path(args.telomerehunter_dir)
    return telomerehunter_dir.parent / f"{telomerehunter_dir.name}_TelomereRepeatLoci"


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def process_sample(args, scripts_dir):
    tumor_bam = Path(args.tumor_bam)
    control_bam = Path(args.control_bam)
    if not tumor_bam.exists():
        raise FileNotFoundError(f"Missing tumor BAM: {tumor_bam}")
    if not control_bam.exists():
        raise FileNotFoundError(f"Missing control BAM: {control_bam}")

    tumor_th_dir, control_th_dir, pid = detect_th_sample_dirs(args.telomerehunter_dir)
    tumor_filtered_bam = get_filtered_bam(tumor_th_dir)
    control_filtered_bam = get_filtered_bam(control_th_dir)

    output_dir = get_output_dir(args)

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

    # Tumor discordant reads
    tumor_discordant = (
        tables_dir / f"{pid}_{args.tumor_sample_name}_discordant_reads.tsv"
    )
    run_command(
        [
            sys.executable,
            str(scripts_dir / "find_discordant_reads.py"),
            "-i",
            str(tumor_filtered_bam),
            "-o",
            str(tumor_discordant),
        ]
    )

    tumor_discordant_with_mapq = (
        tables_dir
        / f"{pid}_{args.tumor_sample_name}_discordant_reads_filtered_with_mapq.tsv"
    )
    run_command(
        [
            sys.executable,
            str(scripts_dir / "add_mate_mapq.py"),
            "-i",
            str(tumor_discordant),
            "-b",
            str(tumor_bam),
            "-o",
            str(tumor_discordant_with_mapq),
        ]
    )

    # Control discordant reads
    control_discordant = (
        tables_dir / f"{pid}_{args.control_sample_name}_discordant_reads.tsv"
    )
    run_command(
        [
            sys.executable,
            str(scripts_dir / "find_discordant_reads.py"),
            "-i",
            str(control_filtered_bam),
            "-o",
            str(control_discordant),
        ]
    )

    control_discordant_with_mapq = (
        tables_dir
        / f"{pid}_{args.control_sample_name}_discordant_reads_filtered_with_mapq.tsv"
    )
    run_command(
        [
            sys.executable,
            str(scripts_dir / "add_mate_mapq.py"),
            "-i",
            str(control_discordant),
            "-b",
            str(control_bam),
            "-o",
            str(control_discordant_with_mapq),
        ]
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

    # Tumor-centric downstream steps
    clipped = clipped_dir / f"{pid}_{args.tumor_sample_name}_clipped_reads.tsv"
    run_command(
        [
            sys.executable,
            str(scripts_dir / "find_fusion_reads.py"),
            str(candidates),
            str(tumor_bam),
            str(clipped),
        ]
    )

    extended = (
        candidate_dir / f"{pid}_telomere_insertions_candidate_regions_extended.tsv"
    )
    run_command(
        [
            sys.executable,
            str(scripts_dir / "predict_insertion_sites.py"),
            str(candidates),
            str(clipped),
            str(tumor_discordant_with_mapq),
            str(extended),
        ]
    )

    extended_with_consensus = (
        candidate_dir
        / f"{pid}_telomere_insertions_candidate_regions_extended_with_consensus.tsv"
    )
    cmd = [
        sys.executable,
        str(scripts_dir / "get_consensus.py"),
        str(extended),
        str(clipped),
        str(extended_with_consensus),
    ]
    if args.reference_fasta:
        cmd.extend(["--reference", args.reference_fasta])
    run_command(cmd)

    bed_zoomed_out = bed_zoomed_out_dir / f"{pid}_telomere_insertions.bed"
    bed_zoomed_in = bed_zoomed_in_dir / f"{pid}_telomere_insertions.bed"
    run_command(
        [
            sys.executable,
            str(scripts_dir / "make_bed_for_visualization.py"),
            str(extended),
            str(bed_zoomed_out),
            str(bed_zoomed_in),
            pid,
        ]
    )

    if args.run_visualization:
        visualize_cmd = [
            sys.executable,
            str(scripts_dir / "visualize_telomere_insertions.py"),
            "--tumor",
            str(tumor_bam),
            "--control",
            str(control_bam),
            "--ref",
            args.reference_fasta,
            "--bed",
            str(bed_zoomed_in),
            "--samtoolsbin",
            args.samtoolsbin,
            "--colored_reads_tumor",
            str(tumor_discordant_with_mapq),
            "--colored_reads_control",
            str(control_discordant_with_mapq),
            "--clipped_reads_tumor",
            str(clipped),
            "--prefix",
            f"{plot_zoomed_in_dir}/",
            "--outfile",
            str(plot_zoomed_in_dir / f"{pid}_done.txt"),
        ]
        run_command(visualize_cmd)


def get_version_from_package():
    try:
        from importlib.metadata import PackageNotFoundError, version

        try:
            return version("telomererepeatloci")
        except PackageNotFoundError:
            return "unknown - check pyproject.toml file"
    except ImportError:
        return "unknown - check pyproject.toml file"


def main():
    print(f"TelomereRepeatLoci - version {get_version_from_package()}")
    args = parse_args()
    scripts_dir = Path(__file__).resolve().parent / "src"
    print("--- Processing sample ---")
    process_sample(args, scripts_dir)


if __name__ == "__main__":
    main()
