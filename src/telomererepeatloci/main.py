#!/usr/bin/env python3

import argparse
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from pipeline.tables import WINDOWS_COLUMNS, read_tsv, write_tsv
from . import (
    add_mate_mapq,
    count_discordant_reads,
    find_discordant_reads,
    find_fusion_reads,
    get_candidate_regions,
    get_consensus,
    make_bed_for_visualization,
    predict_insertion_sites,
    visualize_telomere_insertions,
)


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
    parser.add_argument("--samtoolsbin", default="samtools")

    args = parser.parse_args()

    provided = {}
    for action in parser._actions:
        if action.dest == "help":
            continue
        value = getattr(args, action.dest)
        if value is not None and value != action.default:
            provided[action.dest] = value

    print("Provided arguments:")
    for k, v in provided.items():
        print(f"{k} = {v}")

    return args


def _log(message):
    print(f"[main] {message}")


def _run_timed_step(step_name, fn, *args, **kwargs):
    started = time.perf_counter()
    _log(f"Starting step: {step_name}")
    result = fn(*args, **kwargs)
    elapsed = time.perf_counter() - started
    _log(f"Finished step: {step_name} in {elapsed:.2f}s")
    return result


def _run_parallel_named_steps(named_steps):
    started = time.perf_counter()
    names = ", ".join(name for name, _, _, _ in named_steps)
    _log(f"Starting parallel steps: {names}")
    results = {}
    with ThreadPoolExecutor(max_workers=len(named_steps)) as executor:
        future_to_name = {
            executor.submit(fn, *args, **kwargs): name
            for name, fn, args, kwargs in named_steps
        }
        for future in as_completed(future_to_name):
            name = future_to_name[future]
            results[name] = future.result()
            _log(f"Finished step: {name}")
    elapsed = time.perf_counter() - started
    _log(f"Finished parallel block in {elapsed:.2f}s")
    return results


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


def process_sample(args):
    workflow_started = time.perf_counter()
    tumor_bam = Path(args.tumor_bam)
    if not tumor_bam.exists():
        raise FileNotFoundError(f"Missing tumor BAM: {tumor_bam}")

    tumor_filtered_bam = Path(args.tel_tumor_bam)
    if not tumor_filtered_bam.exists():
        raise FileNotFoundError(f"Missing telomeric tumor BAM: {tumor_filtered_bam}")

    tumor_th_dir = tumor_filtered_bam.parent
    pid = extract_pid_from_folder(tumor_th_dir)
    _log(f"Detected telomere BAM dir: {tumor_th_dir}")
    _log(f"Detected PID: {pid}")

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
        _log(f"Detected control telomere BAM dir: {control_th_dir}")

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

    tumor_discordant = (
        tables_dir / f"{pid}_{args.tumor_sample_name}_discordant_reads.tsv"
    )
    tumor_discordant_with_mapq = (
        tables_dir
        / f"{pid}_{args.tumor_sample_name}_discordant_reads_filtered_with_mapq.tsv"
    )

    control_discordant_with_mapq = Path("NULL")
    if use_control:
        control_discordant = (
            tables_dir / f"{pid}_{args.control_sample_name}_discordant_reads.tsv"
        )
        control_discordant_with_mapq = (
            tables_dir
            / f"{pid}_{args.control_sample_name}_discordant_reads_filtered_with_mapq.tsv"
        )
        _run_parallel_named_steps(
            [
                (
                    f"find_discordant_reads:{args.tumor_sample_name}",
                    find_discordant_reads.run,
                    (str(tumor_filtered_bam), str(tumor_discordant)),
                    {},
                ),
                (
                    f"find_discordant_reads:{args.control_sample_name}",
                    find_discordant_reads.run,
                    (str(control_filtered_bam), str(control_discordant)),
                    {},
                ),
            ]
        )
        _run_parallel_named_steps(
            [
                (
                    f"add_mate_mapq:{args.tumor_sample_name}",
                    add_mate_mapq.add_mate_mapq_file,
                    (
                        str(tumor_discordant),
                        str(tumor_bam),
                        str(tumor_discordant_with_mapq),
                    ),
                    {},
                ),
                (
                    f"add_mate_mapq:{args.control_sample_name}",
                    add_mate_mapq.add_mate_mapq_file,
                    (
                        str(control_discordant),
                        str(control_bam),
                        str(control_discordant_with_mapq),
                    ),
                    {},
                ),
            ]
        )
    else:
        _run_timed_step(
            f"find_discordant_reads:{args.tumor_sample_name}",
            find_discordant_reads.run,
            str(tumor_filtered_bam),
            str(tumor_discordant),
        )
        _run_timed_step(
            f"add_mate_mapq:{args.tumor_sample_name}",
            add_mate_mapq.add_mate_mapq_file,
            str(tumor_discordant),
            str(tumor_bam),
            str(tumor_discordant_with_mapq),
        )

    windows = tables_dir / f"{pid}_discordant_reads_1_kb_windows.tsv"
    windows_df = _run_timed_step(
        "count_discordant_reads",
        count_discordant_reads.compute_windows,
        str(tumor_discordant_with_mapq),
        str(control_discordant_with_mapq),
        args.blacklist,
        str(windows),
    )
    write_tsv(windows_df, windows, WINDOWS_COLUMNS)
    _log(f"Wrote windows table: {windows}")

    candidates = candidate_dir / f"{pid}_telomere_insertions_candidate_regions.tsv"
    candidates_df = _run_timed_step(
        "get_candidate_regions",
        get_candidate_regions.filter_candidates,
        read_tsv(windows),
        args.tumor_discordant_read_lower_limit,
        args.control_discordant_read_upper_limit,
        str(args.consider_blacklist),
    )
    write_tsv(candidates_df, candidates, list(candidates_df.columns))
    _log(f"Wrote candidate regions: {candidates}")

    clipped = clipped_dir / f"{pid}_{args.tumor_sample_name}_clipped_reads.tsv"
    clipped_df = _run_timed_step(
        f"find_fusion_reads:{args.tumor_sample_name}",
        find_fusion_reads.find_fusion_reads,
        str(candidates),
        str(tumor_bam),
    )
    write_tsv(clipped_df, clipped, find_fusion_reads.FUSION_READS_COLUMNS)
    _log(f"Wrote clipped reads: {clipped}")

    extended = (
        candidate_dir / f"{pid}_telomere_insertions_candidate_regions_extended.tsv"
    )
    extended_df, extended_fields = _run_timed_step(
        "predict_insertion_sites",
        predict_insertion_sites.predict_insertions,
        str(candidates),
        str(clipped),
        str(tumor_discordant_with_mapq),
        str(tumor_bam),
        str(control_bam) if use_control else "",
    )
    write_tsv(extended_df, extended, extended_fields)
    _log(f"Wrote extended candidate regions: {extended}")

    extended_with_consensus = (
        candidate_dir
        / f"{pid}_telomere_insertions_candidate_regions_extended_with_consensus.tsv"
    )
    consensus_df, consensus_fields = _run_timed_step(
        "get_consensus",
        get_consensus.build_consensus,
        str(extended),
        str(clipped),
        args.reference_fasta,
    )
    write_tsv(consensus_df, extended_with_consensus, consensus_fields)
    _log(f"Wrote consensus table: {extended_with_consensus}")

    bed_zoomed_out = bed_zoomed_out_dir / f"{pid}_telomere_insertions.bed"
    bed_zoomed_in = bed_zoomed_in_dir / f"{pid}_telomere_insertions.bed"
    _run_timed_step(
        "make_bed_for_visualization",
        make_bed_for_visualization.build_beds,
        str(extended),
        pid,
        str(bed_zoomed_out),
        str(bed_zoomed_in),
        args.plot_min_support,
    )

    if not args.skip_visualization:
        visualize_argv = [
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
            visualize_argv.extend(
                [
                    "--control",
                    str(control_bam),
                    "--colored_reads_control",
                    str(control_discordant_with_mapq),
                ]
            )
        visualize_args = visualize_telomere_insertions.parse_args(visualize_argv)
        _run_timed_step(
            "visualize_telomere_insertions",
            visualize_telomere_insertions.run,
            visualize_args,
        )

    elapsed = time.perf_counter() - workflow_started
    _log(f"Workflow finished successfully in {elapsed:.2f}s")


def get_version_from_package():
    try:
        from telomererepeatloci.version import __version__

        return __version__
    except ImportError:
        return "unknown - check pyproject.toml file"


def main():
    print(f"TelomereRepeatLoci - version {get_version_from_package()}")
    args = parse_args()
    _log("Processing sample")
    process_sample(args)


if __name__ == "__main__":
    main()
