from dataclasses import dataclass
from pathlib import Path

from telomererepeatloci import (
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
from pipeline.tables import read_tsv, write_tsv


@dataclass
class DebugArgs:
    tumor_bam: str
    tel_tumor_bam: str
    control_bam: str = ""
    tel_control_bam: str = ""
    output_dir: str = ""
    tumor_sample_name: str = "tumor"
    control_sample_name: str = "control"
    blacklist: str = "no_file"
    tumor_discordant_read_lower_limit: float = 3.0
    control_discordant_read_upper_limit: float = 0.0
    consider_blacklist: bool = False
    reference_fasta: str = ""
    run_visualization: bool = False
    samtoolsbin: str = "samtools"


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


def get_output_dir(args: DebugArgs, tumor_th_dir: Path) -> Path:
    if args.output_dir:
        return Path(args.output_dir)
    telomerehunter_dir = tumor_th_dir.parent
    return (
        telomerehunter_dir.parent
        / f"{telomerehunter_dir.name.replace('_TelomerCnt', '')}_TelomereRepeatLoci"
    )


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def _validate_inputs(args: DebugArgs):
    tumor_bam = Path(args.tumor_bam)
    if not tumor_bam.exists():
        raise FileNotFoundError(f"Missing tumor BAM: {tumor_bam}")

    tumor_filtered_bam = Path(args.tel_tumor_bam)
    if not tumor_filtered_bam.exists():
        raise FileNotFoundError(f"Missing telomeric tumor BAM: {tumor_filtered_bam}")

    control_bam = None
    control_filtered_bam = None
    use_control = bool(args.control_bam)

    if use_control:
        control_bam = Path(args.control_bam)
        if not control_bam.exists():
            raise FileNotFoundError(f"Missing control BAM: {control_bam}")

        control_filtered_bam = Path(args.tel_control_bam)
        if not control_filtered_bam.exists():
            raise FileNotFoundError(
                f"Missing telomeric control BAM: {control_filtered_bam}"
            )

    return tumor_bam, tumor_filtered_bam, control_bam, control_filtered_bam


def run_direct_pipeline(args: DebugArgs):
    tumor_bam, tumor_filtered_bam, control_bam, control_filtered_bam = _validate_inputs(
        args
    )

    tumor_th_dir = tumor_filtered_bam.parent
    pid = extract_pid_from_folder(tumor_th_dir)

    use_control = bool(args.control_bam)
    if use_control:
        control_th_dir = control_filtered_bam.parent
        control_pid = extract_pid_from_folder(control_th_dir)
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
    find_discordant_reads.run(str(tumor_filtered_bam), str(tumor_discordant))

    tumor_discordant_with_mapq = (
        tables_dir
        / f"{pid}_{args.tumor_sample_name}_discordant_reads_filtered_with_mapq.tsv"
    )
    add_mate_mapq.add_mate_mapq_file(
        str(tumor_discordant), str(tumor_bam), str(tumor_discordant_with_mapq)
    )

    control_discordant_with_mapq = Path("NULL")
    if use_control:
        control_discordant = (
            tables_dir / f"{pid}_{args.control_sample_name}_discordant_reads.tsv"
        )
        find_discordant_reads.run(str(control_filtered_bam), str(control_discordant))

        control_discordant_with_mapq = (
            tables_dir
            / f"{pid}_{args.control_sample_name}_discordant_reads_filtered_with_mapq.tsv"
        )
        add_mate_mapq.add_mate_mapq_file(
            str(control_discordant),
            str(control_bam),
            str(control_discordant_with_mapq),
        )

    windows = tables_dir / f"{pid}_discordant_reads_1_kb_windows.tsv"
    windows_df = count_discordant_reads.compute_windows(
        str(tumor_discordant_with_mapq),
        str(control_discordant_with_mapq),
        args.blacklist,
        str(windows),
    )
    write_tsv(windows_df, windows, count_discordant_reads.WINDOWS_COLUMNS)

    candidates = candidate_dir / f"{pid}_telomere_insertions_candidate_regions.tsv"
    windows_df_read = read_tsv(windows)
    candidates_df = get_candidate_regions.filter_candidates(
        windows_df_read,
        args.tumor_discordant_read_lower_limit,
        args.control_discordant_read_upper_limit,
        str(args.consider_blacklist),
    )
    write_tsv(candidates_df, candidates, list(windows_df_read.columns))

    clipped = clipped_dir / f"{pid}_{args.tumor_sample_name}_clipped_reads.tsv"
    clipped_df = find_fusion_reads.find_fusion_reads(str(candidates), str(tumor_bam))
    write_tsv(clipped_df, clipped, find_fusion_reads.FUSION_READS_COLUMNS)

    extended = (
        candidate_dir / f"{pid}_telomere_insertions_candidate_regions_extended.tsv"
    )
    extended_df, extended_fields = predict_insertion_sites.predict_insertions(
        str(candidates), str(clipped), str(tumor_discordant_with_mapq)
    )
    write_tsv(extended_df, extended, extended_fields)

    extended_with_consensus = (
        candidate_dir
        / f"{pid}_telomere_insertions_candidate_regions_extended_with_consensus.tsv"
    )
    consensus_df, consensus_fields = get_consensus.build_consensus(
        str(extended), str(clipped), args.reference_fasta
    )
    write_tsv(consensus_df, extended_with_consensus, consensus_fields)

    bed_zoomed_out = bed_zoomed_out_dir / f"{pid}_telomere_insertions.bed"
    bed_zoomed_in = bed_zoomed_in_dir / f"{pid}_telomere_insertions.bed"
    make_bed_for_visualization.build_beds(
        str(extended), pid, str(bed_zoomed_out), str(bed_zoomed_in)
    )

    if args.run_visualization:
        visualize_args = visualize_telomere_insertions.parse_args(
            [
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
        )
        if use_control:
            visualize_args = visualize_telomere_insertions.parse_args(
                [
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
                    "--control",
                    str(control_bam),
                    "--colored_reads_control",
                    str(control_discordant_with_mapq),
                ]
            )
        visualize_telomere_insertions.run(visualize_args)


def make_debug_args(
    tumor_bam: str,
    tel_tumor_bam: str,
    control_bam: str = "",
    tel_control_bam: str = "",
    output_dir: str = "",
    tumor_sample_name: str = "tumor",
    control_sample_name: str = "control",
    blacklist: str = "no_file",
    tumor_discordant_read_lower_limit: float = 3.0,
    control_discordant_read_upper_limit: float = 0.0,
    consider_blacklist: bool = False,
    reference_fasta: str = "",
    run_visualization: bool = True,
    samtoolsbin: str = "samtools",
):
    return DebugArgs(
        tumor_bam=tumor_bam,
        tel_tumor_bam=tel_tumor_bam,
        control_bam=control_bam,
        tel_control_bam=tel_control_bam,
        output_dir=output_dir,
        tumor_sample_name=tumor_sample_name,
        control_sample_name=control_sample_name,
        blacklist=blacklist,
        tumor_discordant_read_lower_limit=tumor_discordant_read_lower_limit,
        control_discordant_read_upper_limit=control_discordant_read_upper_limit,
        consider_blacklist=consider_blacklist,
        reference_fasta=reference_fasta,
        run_visualization=run_visualization,
        samtoolsbin=samtoolsbin,
    )
