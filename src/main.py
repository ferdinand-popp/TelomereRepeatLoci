#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Run the TelomereRepeatLoci workflow.")
    parser.add_argument("--tumor-bam", required=True, help="Required Tumor BAM file.")
    parser.add_argument(
        "--control-bam",
        default="",
        help="Optional control BAM file. If not provided, workflow runs in tumor-only mode.",
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

    # Derive the tumor TelomereHunter directory and PID directly from the BAM path.
    # Expected layout: <telomerehunter-dir>/tumor_TelomerCnt_<PID>/<bam-file>
    tumor_th_dir = tumor_bam.parent
    pid = extract_pid_from_folder(tumor_th_dir)
    print(f"Detected tumor TelomereHunter dir : {tumor_th_dir}")
    print(f"Detected PID                      : {pid}")

    tumor_filtered_bam = get_filtered_bam(tumor_th_dir)

    use_control = bool(args.control_bam)
    control_bam = None
    control_filtered_bam = None

    if use_control:
        control_bam = Path(args.control_bam)
        if not control_bam.exists():
            raise FileNotFoundError(f"Missing control BAM: {control_bam}")

        # Derive the control TelomereHunter directory from the control BAM path.
        # Expected layout: <telomerehunter-dir>/control_TelomerCnt_<PID>/<bam-file>
        control_th_dir = control_bam.parent
        control_pid = extract_pid_from_folder(control_th_dir)
        print(f"Detected control TelomereHunter dir: {control_th_dir}")

        if control_pid != pid:
            raise ValueError(
                "Tumor/control PID mismatch from TelomereHunter folder names: "
                f"tumor PID='{pid}', control PID='{control_pid}'."
            )

        control_filtered_bam = get_filtered_bam(control_th_dir)

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

    # Optional control discordant reads
    control_discordant_with_mapq = Path("NULL")
    if use_control:
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
    scripts_dir = Path(__file__).resolve().parent
    print("--- Processing sample ---")
    process_sample(args, scripts_dir)


if __name__ == "__main__":
    main()
