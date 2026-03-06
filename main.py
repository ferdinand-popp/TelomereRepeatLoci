#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the TelomereRepeatLoci workflow without Snakemake/YAML."
    )
    parser.add_argument("--results-per-pid-dir", required=True)
    parser.add_argument("--telomerehunter-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--pids",
        default="all",
        help='Space-separated PIDs (e.g. "PID1 PID2") or "all".',
    )
    parser.add_argument("--tumor-sample-name", default="tumor")
    parser.add_argument("--control-sample-name", default="control")
    parser.add_argument(
        "--with-control",
        action="store_true",
        help="Include control sample processing.",
    )
    parser.add_argument("--bam-suffix", default="_merged.mdup.bam")
    parser.add_argument("--blacklist", default="no_file")
    parser.add_argument("--tumor-discordant-read-lower-limit", type=float, default=3.0)
    parser.add_argument(
        "--control-discordant-read-upper-limit", type=float, default=0.0
    )
    parser.add_argument("--consider-blacklist", action="store_true")
    parser.add_argument("--reference-fasta", default="")
    parser.add_argument(
        "--run-telomerehunter",
        action="store_true",
        help="Run telomerehunter if intratelomeric BAM outputs are missing.",
    )
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


def pid_list(args):
    results_dir = Path(args.results_per_pid_dir)
    if args.pids == "all":
        return sorted(
            d.name
            for d in results_dir.iterdir()
            if d.is_dir() and not d.name.startswith(".")
        )
    return [x for x in args.pids.split() if x]


def alignment_bam_path(results_dir, pid, sample, suffix):
    return Path(results_dir) / pid / "alignment" / f"{sample}_{pid}{suffix}"


def intratel_bam_path(telomerehunter_dir, pid, sample):
    return (
        Path(telomerehunter_dir)
        / pid
        / f"{sample}_TelomerCnt_{pid}"
        / f"{pid}_filtered_intratelomeric.bam"
    )


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def process_pid(args, scripts_dir, pid):
    samples = [args.tumor_sample_name]
    if args.with_control:
        samples.append(args.control_sample_name)

    input_bams = {}
    for sample in samples:
        bam = alignment_bam_path(args.results_per_pid_dir, pid, sample, args.bam_suffix)
        if not bam.exists():
            raise FileNotFoundError(
                f"Missing alignment BAM for PID '{pid}', sample '{sample}': {bam}"
            )
        input_bams[sample] = bam

    intratel_bams = {
        sample: intratel_bam_path(args.telomerehunter_dir, pid, sample)
        for sample in samples
    }

    if args.run_telomerehunter and not all(
        path.exists() for path in intratel_bams.values()
    ):
        cmd = [
            "telomerehunter",
            "-p",
            pid,
            "-o",
            args.telomerehunter_dir,
            "-ibt",
            str(input_bams[args.tumor_sample_name]),
        ]
        if args.with_control:
            cmd.extend(["-ibc", str(input_bams[args.control_sample_name]), "-pl"])
        cmd.extend(["-pff", "all"])
        run_command(cmd)

    for sample, intratel_bam in intratel_bams.items():
        if not intratel_bam.exists():
            raise FileNotFoundError(
                f"Missing intratelomeric BAM for PID '{pid}', sample '{sample}': {intratel_bam}. "
                "Either provide it or use --run-telomerehunter."
            )

    tables_dir = Path(args.output_dir) / "tables"
    clipped_dir = Path(args.output_dir) / "clipped_reads"
    candidate_dir = Path(args.output_dir) / "candidate_region_tables"
    bed_zoomed_out_dir = Path(args.output_dir) / "plots" / "bedfiles" / "zoomed_out"
    bed_zoomed_in_dir = Path(args.output_dir) / "plots" / "bedfiles" / "zoomed_in"
    plot_zoomed_in_dir = Path(args.output_dir) / "plots" / "zoomed_in"

    for path in [
        tables_dir,
        clipped_dir,
        candidate_dir,
        bed_zoomed_out_dir,
        bed_zoomed_in_dir,
        plot_zoomed_in_dir,
    ]:
        ensure_dir(path)

    for sample in samples:
        discordant = tables_dir / f"{pid}_{sample}_discordant_reads.tsv"
        run_command(
            [
                sys.executable,
                str(scripts_dir / "find_discordant_reads.py"),
                "-i",
                str(intratel_bams[sample]),
                "-o",
                str(discordant),
            ]
        )

        discordant_with_mapq = (
            tables_dir / f"{pid}_{sample}_discordant_reads_filtered_with_mapq.tsv"
        )
        run_command(
            [
                sys.executable,
                str(scripts_dir / "add_mate_mapq.py"),
                "-i",
                str(discordant),
                "-b",
                str(input_bams[sample]),
                "-o",
                str(discordant_with_mapq),
            ]
        )

    discordant_tumor = (
        tables_dir
        / f"{pid}_{args.tumor_sample_name}_discordant_reads_filtered_with_mapq.tsv"
    )
    discordant_control = (
        tables_dir
        / f"{pid}_{args.control_sample_name}_discordant_reads_filtered_with_mapq.tsv"
        if args.with_control
        else Path("NULL")
    )
    windows = tables_dir / f"{pid}_discordant_reads_1_kb_windows.tsv"
    run_command(
        [
            sys.executable,
            str(scripts_dir / "count_discordant_reads.py"),
            "-t",
            str(discordant_tumor),
            "-c",
            str(discordant_control),
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

    clipped_read_samples = [args.tumor_sample_name]
    if args.run_visualization and args.with_control:
        clipped_read_samples.append(args.control_sample_name)

    for sample in clipped_read_samples:
        clipped = clipped_dir / f"{pid}_{sample}_clipped_reads.tsv"
        run_command(
            [
                sys.executable,
                str(scripts_dir / "find_fusion_reads.py"),
                str(candidates),
                str(input_bams[sample]),
                str(clipped),
            ]
        )

    tumor_clipped = clipped_dir / f"{pid}_{args.tumor_sample_name}_clipped_reads.tsv"
    extended = (
        candidate_dir / f"{pid}_telomere_insertions_candidate_regions_extended.tsv"
    )
    run_command(
        [
            sys.executable,
            str(scripts_dir / "predict_insertion_sites.py"),
            str(candidates),
            str(tumor_clipped),
            str(discordant_tumor),
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
        str(tumor_clipped),
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
            str(input_bams[args.tumor_sample_name]),
            "--ref",
            args.reference_fasta,
            "--bed",
            str(bed_zoomed_in),
            "--samtoolsbin",
            args.samtoolsbin,
            "--colored_reads_tumor",
            str(discordant_tumor),
            "--clipped_reads_tumor",
            str(tumor_clipped),
            "--prefix",
            f"{plot_zoomed_in_dir}/",
            "--outfile",
            str(plot_zoomed_in_dir / f"{pid}_done.txt"),
        ]
        if args.with_control:
            visualize_cmd.extend(
                [
                    "--control",
                    str(input_bams[args.control_sample_name]),
                    "--colored_reads_control",
                    str(
                        tables_dir
                        / f"{pid}_{args.control_sample_name}_discordant_reads_filtered_with_mapq.tsv"
                    ),
                    "--clipped_reads_control",
                    str(
                        clipped_dir
                        / f"{pid}_{args.control_sample_name}_clipped_reads.tsv"
                    ),
                ]
            )
        run_command(visualize_cmd)


def main():
    args = parse_args()
    scripts_dir = Path(__file__).resolve().parent / "src"
    pids = pid_list(args)
    if not pids:
        raise ValueError("No PIDs found to process.")

    for pid in pids:
        print(f"--- Processing PID {pid} ---")
        process_pid(args, scripts_dir, pid)


if __name__ == "__main__":
    main()
