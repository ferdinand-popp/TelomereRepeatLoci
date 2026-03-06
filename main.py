#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the TelomereRepeatLoci workflow without Snakemake/YAML."
    )
    parser.add_argument("--input-bam", required=True)
    parser.add_argument("--telomerehunter-dir", required=True)
    parser.add_argument(
        "--output-dir",
        default="",
        help=(
            "Optional output directory. If not provided, a sibling directory named "
            "<telomerehunter-dir>_TelomereRepeatLoci is created."
        ),
    )
    parser.add_argument("--tumor-sample-name", default="tumor")
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


def get_intratelomeric_bam(telomerehunter_dir):
    matches = sorted(Path(telomerehunter_dir).glob("*_filtered_intratelomeric.bam"))
    if not matches:
        raise FileNotFoundError(
            f"No *_filtered_intratelomeric.bam found in {telomerehunter_dir}"
        )
    if len(matches) > 1:
        raise ValueError(
            f"Multiple intratelomeric BAM files found in {telomerehunter_dir}: {matches}"
        )
    return matches[0]


def get_output_dir(args):
    if args.output_dir:
        return Path(args.output_dir)
    telomerehunter_dir = Path(args.telomerehunter_dir)
    return telomerehunter_dir.parent / f"{telomerehunter_dir.name}_TelomereRepeatLoci"


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def process_sample(args, scripts_dir):
    input_bam = Path(args.input_bam)
    if not input_bam.exists():
        raise FileNotFoundError(f"Missing input BAM: {input_bam}")

    intratel_bam = get_intratelomeric_bam(args.telomerehunter_dir)
    pid = intratel_bam.stem.removesuffix("_filtered_intratelomeric")

    tables_dir = get_output_dir(args) / "tables"
    clipped_dir = get_output_dir(args) / "clipped_reads"
    candidate_dir = get_output_dir(args) / "candidate_region_tables"
    bed_zoomed_out_dir = get_output_dir(args) / "plots" / "bedfiles" / "zoomed_out"
    bed_zoomed_in_dir = get_output_dir(args) / "plots" / "bedfiles" / "zoomed_in"
    plot_zoomed_in_dir = get_output_dir(args) / "plots" / "zoomed_in"

    for path in [
        tables_dir,
        clipped_dir,
        candidate_dir,
        bed_zoomed_out_dir,
        bed_zoomed_in_dir,
        plot_zoomed_in_dir,
    ]:
        ensure_dir(path)

    sample = args.tumor_sample_name
    discordant = tables_dir / f"{pid}_{sample}_discordant_reads.tsv"
    run_command(
        [
            sys.executable,
            str(scripts_dir / "find_discordant_reads.py"),
            "-i",
            str(intratel_bam),
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
            str(input_bam),
            "-o",
            str(discordant_with_mapq),
        ]
    )

    discordant_tumor = (
        tables_dir
        / f"{pid}_{args.tumor_sample_name}_discordant_reads_filtered_with_mapq.tsv"
    )
    discordant_control = Path("NULL")
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

    clipped = clipped_dir / f"{pid}_{args.tumor_sample_name}_clipped_reads.tsv"
    run_command(
        [
            sys.executable,
            str(scripts_dir / "find_fusion_reads.py"),
            str(candidates),
            str(input_bam),
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
            str(input_bam),
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
        run_command(visualize_cmd)


def main():
    args = parse_args()
    scripts_dir = Path(__file__).resolve().parent / "src"
    print("--- Processing sample ---")
    process_sample(args, scripts_dir)


if __name__ == "__main__":
    main()
