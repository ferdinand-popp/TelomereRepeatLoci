"""
Author: Lina Sieverling
Affiliation: DKFZ Heidelberg
Aim: A Snakemake workflow to find telomere insertions
Date: Thu Aug 18 17:46:12 CEST 2016
Run: snakemake -s <Snakefile> --configfile <config.yaml> --cores <N>

Note: current Snakemake requires --cores to be set explicitly for local
execution, e.g. --cores 4 or --cores all.
"""

#---------------------------------------------------------------------------------------
# setup: imports, logging, config validation
#---------------------------------------------------------------------------------------
import logging
import os

import pandas as pd

logger = logging.getLogger("telomere_insertions")
logger.setLevel(logging.INFO)
if not logger.handlers:
    _handler = logging.StreamHandler()
    _handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    logger.addHandler(_handler)
# Don't rely on / interfere with the root logger's handlers (Snakemake configures
# its own before the Snakefile is parsed, which can silently swallow messages
# from a plain logging.basicConfig() call since basicConfig() is a no-op once
# the root logger already has handlers).
logger.propagate = False

REQUIRED_CONFIG_KEYS = [
    "samples",
    "pids",
    "telomerehunter_dir",
    "telomereinsertion_dir",
    "src_dir",
    "R_function_file",
    "blacklist",
    "sleep_sec_limit",
    "tumor_discordant_read_lower_limit",
    "reference_fasta",
]

missing_config_keys = [key for key in REQUIRED_CONFIG_KEYS if key not in config]
if missing_config_keys:
    raise ValueError(
        "Missing required config key(s): " + ", ".join(missing_config_keys)
        + ". Please check your --configfile."
    )

if len(config["samples"]) not in (1, 2):
    raise ValueError('config["samples"] must contain either 1 (tumor only) or 2 (tumor + control) entries.')


def _is_enabled_config_path(path_value):
    return path_value not in [None, "", "no_file"]


def _parse_bool_config(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    value_str = str(value).strip().lower()
    if value_str in ["1", "true", "yes", "y", "on"]:
        return True
    if value_str in ["0", "false", "no", "n", "off"]:
        return False
    return default


def _hms_to_minutes(hms):
    """Convert an 'HH:MM:SS' walltime string into an integer number of minutes
    (rounded up), for use with the modern `resources: runtime=...` directive."""
    hours, minutes, seconds = (int(x) for x in hms.split(":"))
    total_minutes = hours * 60 + minutes + (1 if seconds else 0)
    return max(total_minutes, 1)


def _mem_to_mb(mem_str):
    """Convert an old-style '150m' / '1g' memory string into an integer number
    of megabytes, for use with the modern `resources: mem_mb=...` directive."""
    mem_str = str(mem_str).strip().lower()
    if mem_str.endswith("g"):
        return int(float(mem_str[:-1]) * 1024)
    if mem_str.endswith("m"):
        return int(float(mem_str[:-1]))
    return int(mem_str)


#---------------------------------------------------------------------------------------
# load bam paths (either from an explicit TSV, or derived from results_per_pid_dir)
#---------------------------------------------------------------------------------------

def _load_bam_paths_from_tsv(tsv_file, sample_names):
    """Read pid -> {sample_name: bam_path} from a TSV using pandas.

    The pid column is taken positionally (always the first column),
    regardless of what its header is actually called (e.g. "pid",
    "patient_id", "sample_id", ...).

    The alignment bam path columns are required and expected to be named
    "path_to_<sample>_bam", e.g. for sample_names = ["tumor", "control"]:
    "path_to_tumor_bam" and "path_to_control_bam".

    The TelomereHunter intratelomeric bam path columns are optional and
    expected to be named "path_to_<sample>_intratelomeric_bam". If present,
    these are used (when skip_telomerehunter=true) instead of the default
    TelomereHunter output path pattern.

    Returns a tuple (bam_paths, intratelomeric_bam_paths), each a
    pid -> {sample_name: path} dict. intratelomeric_bam_paths only contains
    entries for the sample columns that were actually present in the TSV.
    """
    # Auto-detect the delimiter (tab or comma) rather than assuming tab, since
    # files named/passed as "tsv" sometimes turn out to actually be comma-separated.
    df = pd.read_csv(tsv_file, sep=None, engine="python", dtype=str, keep_default_na=False)
    if df.columns.empty:
        raise ValueError(f"bam_files_tsv has no header or is empty: {tsv_file}")
    df.columns = [c.strip() for c in df.columns]

    if len(df.columns) == 1:
        raise ValueError(
            f"bam_files_tsv only parsed a single column ({df.columns[0]!r}) from {tsv_file}. "
            "This usually means the file's delimiter wasn't detected correctly — "
            "double check it's actually tab- or comma-separated."
        )

    # First column holds the pid, whatever its header is named.
    pid_column = df.columns[0]

    bam_columns = {sample_name: f"path_to_{sample_name}_bam" for sample_name in sample_names}
    missing_columns = [col for col in bam_columns.values() if col not in df.columns]
    if missing_columns:
        raise ValueError(f"bam_files_tsv is missing required columns: {', '.join(missing_columns)}")

    # Optional: only include the intratelomeric-bam columns that actually exist in the TSV.
    th_bam_columns_all = {sample_name: f"path_to_{sample_name}_intratelomeric_bam" for sample_name in sample_names}
    th_bam_columns = {s: c for s, c in th_bam_columns_all.items() if c in df.columns}

    for col in [pid_column, *bam_columns.values(), *th_bam_columns.values()]:
        df[col] = df[col].str.strip()

    empty_pid_mask = df[pid_column] == ""
    for row_number in df.index[empty_pid_mask]:
        # +2: 1 to account for the header row, 1 to switch to 1-based line numbers
        logger.warning(f"Skipping row {row_number + 2} with empty pid in bam_files_tsv: {tsv_file}")
    df = df[~empty_pid_mask]

    duplicate_pids = sorted(df[pid_column][df[pid_column].duplicated()].unique())
    if duplicate_pids:
        raise ValueError(f"bam_files_tsv contains duplicate pid(s): {', '.join(duplicate_pids)} in {tsv_file}")

    df = df.set_index(pid_column)

    bam_paths = {
        pid: {sample_name: row[bam_col] for sample_name, bam_col in bam_columns.items()}
        for pid, row in df[list(bam_columns.values())].iterrows()
    }

    th_bam_paths = {}
    if th_bam_columns:
        th_bam_paths = {
            pid: {sample_name: row[bam_col] for sample_name, bam_col in th_bam_columns.items()}
            for pid, row in df[list(th_bam_columns.values())].iterrows()
        }
        logger.info(
            "bam_files_tsv provides TelomereHunter intratelomeric bam paths for: "
            + ", ".join(sorted(th_bam_columns.values()))
        )

    return bam_paths, th_bam_paths


explicit_bam_files_tsv = config.get("bam_files_tsv", "no_file")
use_explicit_bam_paths = _is_enabled_config_path(explicit_bam_files_tsv)

# Config switch:
#   skip_telomerehunter: true/false
# If true, workflow will assume TelomereHunter outputs already exist and will not run run_telomerehunter.
skip_telomerehunter = _parse_bool_config(config.get("skip_telomerehunter", False), default=False)

if use_explicit_bam_paths:
    bam_files_by_pid, th_bam_files_by_pid = _load_bam_paths_from_tsv(explicit_bam_files_tsv, config["samples"])
else:
    bam_files_by_pid = {}
    th_bam_files_by_pid = {}

TELOMEREHUNTER_DIR = config["telomerehunter_dir"]
TELOMEREINSERTION_DIR = config["telomereinsertion_dir"]
SRC_DIR = config["src_dir"]
R_FUNCTION_FILE = config["R_function_file"]
SAMPLES = config["samples"]
REFERENCE_FASTA = config["reference_fasta"]

# Modern Snakemake is stricter about wildcard resolution; constraining
# {pid} and {sample} avoids ambiguous-wildcard errors in the DAG.
wildcard_constraints:
    pid=r"[^/]+",
    sample="|".join(SAMPLES)


def get_alignment_bam(pid_name, sample_name):
    if use_explicit_bam_paths:
        return bam_files_by_pid[pid_name][sample_name]
    return f'{config["results_per_pid_dir"]}/{pid_name}/alignment/{sample_name}_{pid_name}{config["bam_suffix"]}'


def get_telomerehunter_intratelomeric_bam(pid_name, sample_name):
    # Only honor a TSV-supplied custom path when we're actually trusting
    # pre-existing TelomereHunter outputs (skip_telomerehunter=true).
    # Otherwise the pipeline itself will always write to the standard
    # TelomereHunter output location below, regardless of what the TSV says.
    if skip_telomerehunter and use_explicit_bam_paths:
        custom_path = th_bam_files_by_pid.get(pid_name, {}).get(sample_name, "")
        if custom_path:
            return custom_path
    return f"{TELOMEREHUNTER_DIR}/{pid_name}/{sample_name}_TelomerCnt_{pid_name}/{pid_name}_filtered_intratelomeric.bam"


if config["pids"] == "all":
    if use_explicit_bam_paths:
        pids_candidates = sorted(bam_files_by_pid.keys())
    else:
        pids_candidates = sorted(p for p in os.listdir(config["results_per_pid_dir"]) if not p.startswith('.'))
else:
    pids_candidates = config["pids"].split(' ')


#---------------------------------------------------------------------------------------
# remove PIDs without (existing) bam files
#---------------------------------------------------------------------------------------

def _bam_status_table(pids_list, sample_names):
    """Tidy pid/sample table with the resolved bam path and whether it exists on disk."""
    records = []
    for pid_name in pids_list:
        for sample_name in sample_names:
            bam_file = None
            if use_explicit_bam_paths:
                if pid_name not in bam_files_by_pid:
                    records.append((pid_name, sample_name, None))
                    continue
                bam_file = bam_files_by_pid[pid_name].get(sample_name, "")
                if bam_file == "":
                    records.append((pid_name, sample_name, None))
                    continue
            else:
                bam_file = get_alignment_bam(pid_name, sample_name)
            records.append((pid_name, sample_name, bam_file))

    status = pd.DataFrame(records, columns=["pid", "sample", "bam_file"])
    status["exists"] = status["bam_file"].apply(lambda p: p is not None and os.path.exists(p))
    return status


bam_status = _bam_status_table(pids_candidates, SAMPLES)
for _, row in bam_status[~bam_status["exists"]].iterrows():
    if row["bam_file"] is None:
        logger.warning(f"{row['pid']}: no BAM entry found for {row['sample']} in bam_files_tsv, skipping this pid!")
    else:
        logger.warning(
            f"{row['pid']}: alignment bam file for {row['sample']} sample is missing "
            f"({row['bam_file']}), skipping this pid!"
        )

pids_with_missing_bam = set(bam_status.loc[~bam_status["exists"], "pid"])
pids = [p for p in pids_candidates if p not in pids_with_missing_bam]

# Extra validation when skipping TelomereHunter: ensure expected outputs are present
if skip_telomerehunter:
    th_records = [
        (pid_name, sample_name, get_telomerehunter_intratelomeric_bam(pid_name, sample_name))
        for pid_name in pids
        for sample_name in SAMPLES
    ]
    th_status = pd.DataFrame(th_records, columns=["pid", "sample", "bam_file"])
    th_status["exists"] = th_status["bam_file"].apply(os.path.exists)
    for _, row in th_status[~th_status["exists"]].iterrows():
        logger.warning(
            f"{row['pid']}: TelomereHunter output missing for {row['sample']} "
            f"({row['bam_file']}), skipping this pid!"
        )
    pids_with_missing_th = set(th_status.loc[~th_status["exists"], "pid"])
    pids = [p for p in pids if p not in pids_with_missing_th]

if not pids:
    raise ValueError(
        "No PIDs remain after validation — see the WARNING messages above for why each "
        "PID was dropped (missing alignment bam, missing TelomereHunter output, etc.), "
        "and check your config / bam_files_tsv paths."
    )


#------------------------------------------------------------------
# rule all
#------------------------------------------------------------------

if len(SAMPLES) == 2:
    pids_control = pids
else:
    pids_control = []

localrules: all

rule all:
    input:
        [get_telomerehunter_intratelomeric_bam(pid, "tumor") for pid in pids],
        [get_telomerehunter_intratelomeric_bam(pid, "control") for pid in pids_control],
        expand(TELOMEREINSERTION_DIR + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended_with_consensus.tsv', pid=pids),
        expand(TELOMEREINSERTION_DIR + '/plots/zoomed_in/{pid}_done.txt', pid=pids)


#------------------------------------------------------------------
# run telomerehunter
#------------------------------------------------------------------

if len(SAMPLES) == 2:
    input_list = [
        lambda wildcards: get_alignment_bam(wildcards.pid, SAMPLES[0]),
        lambda wildcards: get_alignment_bam(wildcards.pid, SAMPLES[1])
    ]
    output_list = [
        TELOMEREHUNTER_DIR + '/{pid}/' + SAMPLES[0] + '_TelomerCnt_{pid}/{pid}_filtered_intratelomeric.bam',
        TELOMEREHUNTER_DIR + '/{pid}/' + SAMPLES[1] + '_TelomerCnt_{pid}/{pid}_filtered_intratelomeric.bam'
    ]
    telomerehunter_shell_extra = "-ibc {input[1]}"
    telomerehunter_threads = 2
elif len(SAMPLES) == 1:
    input_list = [
        lambda wildcards: get_alignment_bam(wildcards.pid, SAMPLES[0]),
        lambda wildcards: get_alignment_bam(wildcards.pid, SAMPLES[0])
    ]
    output_list = [
        TELOMEREHUNTER_DIR + '/{pid}/' + SAMPLES[0] + '_TelomerCnt_{pid}/{pid}_filtered_intratelomeric.bam'
    ]
    telomerehunter_shell_extra = ""
    telomerehunter_threads = 1

if not skip_telomerehunter:
    rule run_telomerehunter:
        input:
            input_list
        output:
            output_list
        threads: telomerehunter_threads
        resources:
            mem_mb=_mem_to_mb("150m"),
            runtime=_hms_to_minutes("24:00:00")
        params:
            jobname="{pid}_telomerehunter",
            sleep_sec_limit=config["sleep_sec_limit"],
            telomerehunter_dir=TELOMEREHUNTER_DIR,
            extra=telomerehunter_shell_extra
        shell:
            "sleep $((1 + RANDOM % {params.sleep_sec_limit}))s; "
            "set +u; module load Micromamba/2.0.2-0; module load R/3.4.2; set -u; "
            "time micromamba run -n telomereEnv telomerehunter2 -p {wildcards.pid} -o {params.telomerehunter_dir} -ibt {input[0]} {params.extra}-pff all"
else:
    logger.info("skip_telomerehunter=true -> run_telomerehunter rule disabled; assuming existing TelomereHunter2 outputs.")


#------------------------------------------------------------------
# find discordant reads
#------------------------------------------------------------------

rule find_discordant_reads:
    input:
        lambda wildcards: get_telomerehunter_intratelomeric_bam(wildcards.pid, wildcards.sample)
    output:
        TELOMEREINSERTION_DIR + '/tables/{pid}_{sample}_discordant_reads.tsv'
    resources:
        mem_mb=_mem_to_mb("150m"),
        runtime=_hms_to_minutes("0:59:00")
    params:
        jobname="{pid}_find_discordant_reads_{sample}",
        sleep_sec_limit=config["sleep_sec_limit"],
        src_dir=SRC_DIR
    shell:
        "sleep $((1 + RANDOM % {params.sleep_sec_limit}))s; "
        "set +u; module load Micromamba/2.0.2-0; set -u; "
        "micromamba run -n telomereEnv python {params.src_dir}/find_discordant_reads.py -i {input} -o {output}"


#------------------------------------------------------------------
# add mate mapping quality
#------------------------------------------------------------------

rule add_mate_mapq:
    input:
        discordant_reads=TELOMEREINSERTION_DIR + '/tables/{pid}_{sample}_discordant_reads.tsv',
        bam=lambda wildcards: get_alignment_bam(wildcards.pid, wildcards.sample)
    output:
        TELOMEREINSERTION_DIR + '/tables/{pid}_{sample}_discordant_reads_filtered_with_mapq.tsv'
    resources:
        mem_mb=_mem_to_mb("100m"),
        runtime=_hms_to_minutes("50:00:00")
    params:
        jobname="{pid}_add_mate_mapq_{sample}",
        sleep_sec_limit=config["sleep_sec_limit"],
        src_dir=SRC_DIR
    shell:
        "sleep $((1 + RANDOM % {params.sleep_sec_limit}))s; "
        "set +u; module load Micromamba/2.0.2-0; set -u; "
        "micromamba run -n telomereEnv python {params.src_dir}/add_mate_mapq.py -i {input.discordant_reads} -b {input.bam} -o {output}"


#------------------------------------------------------------------
# count discordant reads
#------------------------------------------------------------------

paired_t_c_flag = False

if len(SAMPLES) == 2:
    input_list = [
        TELOMEREINSERTION_DIR + '/tables/{pid}_' + SAMPLES[0] + '_discordant_reads_filtered_with_mapq.tsv',
        TELOMEREINSERTION_DIR + '/tables/{pid}_' + SAMPLES[1] + '_discordant_reads_filtered_with_mapq.tsv'
    ]
    tumor_input = input_list[0]
    control_input = input_list[1]
    paired_t_c_flag = True
elif len(SAMPLES) == 1:
    input_list = [TELOMEREINSERTION_DIR + '/tables/{pid}_' + SAMPLES[0] + '_discordant_reads_filtered_with_mapq.tsv']
    tumor_input = input_list[0]
    control_input = "NULL"

if not os.path.exists(config["blacklist"]) and not paired_t_c_flag:
    logger.warning("Please provide paired tumor-control samples or a blacklist, otherwise no proper filtering for false positives is possible!")

rule count_discordant_reads:
    input:
        input_list
    output:
        windowTable=TELOMEREINSERTION_DIR + '/tables/{pid}_discordant_reads_1_kb_windows.tsv'
    resources:
        mem_mb=_mem_to_mb("100m"),
        runtime=_hms_to_minutes("0:59:00")
    params:
        blacklist=config["blacklist"],
        jobname="{pid}_count_discordant_reads",
        r_function_file=R_FUNCTION_FILE,
        src_dir=SRC_DIR,
        tumor=tumor_input,
        control=control_input
    shell:
        "R --no-save --slave --args -t {params.tumor} -c {params.control} -b {params.blacklist} "
        "-o {output.windowTable} -f {params.r_function_file} < {params.src_dir}/count_discordant_reads.R"


#------------------------------------------------------------------
# get candidate regions
#------------------------------------------------------------------

rule get_candidate_regions:
    input:
        windowTable=TELOMEREINSERTION_DIR + '/tables/{pid}_discordant_reads_1_kb_windows.tsv'
    output:
        candidateRegions=TELOMEREINSERTION_DIR + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions.tsv'
    resources:
        mem_mb=_mem_to_mb("100m"),
        runtime=_hms_to_minutes("0:15:00")
    params:
        jobname="{pid}_get_candidate_regions",
        r_function_file=R_FUNCTION_FILE,
        src_dir=SRC_DIR,
        tumor_discordant_read_lower_limit=config["tumor_discordant_read_lower_limit"],
        control_discordant_read_upper_limit=config["control_discordant_read_upper_limit"],
        consider_blacklist="True" if _is_enabled_config_path(config.get("blacklist")) else "False"
    shell:
        "R --no-save --slave --args {input.windowTable} {output.candidateRegions} "
        "{params.tumor_discordant_read_lower_limit} {params.control_discordant_read_upper_limit} "
        "{params.consider_blacklist} {params.r_function_file} < {params.src_dir}/get_candidate_regions.R"


#------------------------------------------------------------------
# predict insertion sites
#------------------------------------------------------------------

if len(SAMPLES) == 2:
    bam = lambda wildcards: get_alignment_bam(wildcards.pid, wildcards.sample)
elif len(SAMPLES) == 1:
    bam = lambda wildcards: get_alignment_bam(wildcards.pid, SAMPLES[0])

rule find_fusion_reads:
    input:
        candidateRegions=TELOMEREINSERTION_DIR + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions.tsv',
        bam=bam
    output:
        TELOMEREINSERTION_DIR + '/clipped_reads/{pid}_{sample}_clipped_reads.tsv'
    resources:
        mem_mb=_mem_to_mb("1g"),
        runtime=_hms_to_minutes("100:00:00")
    params:
        jobname="{pid}_find_fusion_reads",
        r_function_file=R_FUNCTION_FILE,
        src_dir=SRC_DIR
    shell:
        "R --no-save --slave --args {input.candidateRegions} {input.bam} {output} "
        "{params.r_function_file} < {params.src_dir}/find_fusion_reads.R"

if len(SAMPLES) == 2:

    rule predict_insertion_sites:
        input:
            candidateRegions=TELOMEREINSERTION_DIR + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions.tsv',
            clippedReads=TELOMEREINSERTION_DIR + '/clipped_reads/{pid}_' + SAMPLES[0] + '_clipped_reads.tsv',
            discordantReads=TELOMEREINSERTION_DIR + '/tables/{pid}_' + SAMPLES[0] + '_discordant_reads_filtered_with_mapq.tsv',
            clippedReadsControl=TELOMEREINSERTION_DIR + '/clipped_reads/{pid}_' + SAMPLES[1] + '_clipped_reads.tsv',
            tumorBam=lambda wildcards: get_alignment_bam(wildcards.pid, SAMPLES[0]),
            controlBam=lambda wildcards: get_alignment_bam(wildcards.pid, SAMPLES[1])
        output:
            TELOMEREINSERTION_DIR + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended.tsv'
        resources:
            mem_mb=_mem_to_mb("100m"),
            runtime=_hms_to_minutes("0:10:00")
        params:
            jobname="{pid}_predict_insertion_sites",
            r_function_file=R_FUNCTION_FILE,
            src_dir=SRC_DIR
        message: "--- {wildcards.pid}: predict insertion sites ---"
        shell:
            "R --no-save --slave --args "
            "--candidate_region_file {input.candidateRegions} "
            "--clipped_reads_file {input.clippedReads} "
            "--discordant_read_file {input.discordantReads} "
            "--outfile {output} "
            "--function_file {params.r_function_file} "
            "--bamfile_tumor {input.tumorBam} "
            "--bamfile_control {input.controlBam} "
            "--clipped_reads_control_file {input.clippedReadsControl} "
            "< {params.src_dir}/predict_insertion_sites.R"

elif len(SAMPLES) == 1:

    rule predict_insertion_sites:
        input:
            candidateRegions=TELOMEREINSERTION_DIR + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions.tsv',
            clippedReads=TELOMEREINSERTION_DIR + '/clipped_reads/{pid}_' + SAMPLES[0] + '_clipped_reads.tsv',
            discordantReads=TELOMEREINSERTION_DIR + '/tables/{pid}_' + SAMPLES[0] + '_discordant_reads_filtered_with_mapq.tsv',
            tumorBam=lambda wildcards: get_alignment_bam(wildcards.pid, SAMPLES[0])
        output:
            TELOMEREINSERTION_DIR + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended.tsv'
        resources:
            mem_mb=_mem_to_mb("100m"),
            runtime=_hms_to_minutes("0:10:00")
        params:
            jobname="{pid}_predict_insertion_sites",
            r_function_file=R_FUNCTION_FILE,
            src_dir=SRC_DIR
        message: "--- {wildcards.pid}: predict insertion sites ---"
        shell:
            "R --no-save --slave --args "
            "--candidate_region_file {input.candidateRegions} "
            "--clipped_reads_file {input.clippedReads} "
            "--discordant_read_file {input.discordantReads} "
            "--outfile {output} "
            "--function_file {params.r_function_file} "
            "--bamfile_tumor {input.tumorBam} "
            "< {params.src_dir}/predict_insertion_sites.R"

#------------------------------------------------------------------
# get consensus sequence of insertion (and bases in sequence microhomology)
#------------------------------------------------------------------

rule get_consensus:
    input:
        candidateRegions=TELOMEREINSERTION_DIR + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended.tsv',
        clippedReads=TELOMEREINSERTION_DIR + '/clipped_reads/{pid}_' + SAMPLES[0] + '_clipped_reads.tsv'
    output:
        TELOMEREINSERTION_DIR + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended_with_consensus.tsv'
    resources:
        mem_mb=_mem_to_mb("500m"),
        runtime=_hms_to_minutes("0:10:00")
    params:
        jobname="{pid}_get_consensus",
        r_function_file=R_FUNCTION_FILE,
        src_dir=SRC_DIR,
        reference_fasta=REFERENCE_FASTA
    message: "--- {wildcards.pid}: get consensus ---"
    shell:
        "R --no-save --slave --args {input.candidateRegions} {input.clippedReads} {output} "
        "{params.reference_fasta} {params.r_function_file} < {params.src_dir}/get_consensus.R"


#------------------------------------------------------------------
# make plots
#------------------------------------------------------------------

rule make_bed_for_visualization:
    input:
        candidateRegions=TELOMEREINSERTION_DIR + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended.tsv'
    output:
        outfile1=TELOMEREINSERTION_DIR + '/plots/bedfiles/zoomed_out/{pid}_telomere_insertions.bed',
        outfile2=TELOMEREINSERTION_DIR + '/plots/bedfiles/zoomed_in/{pid}_telomere_insertions.bed',
        outfile3=TELOMEREINSERTION_DIR + '/plots/bedfiles/flagged/{pid}_telomere_insertions_review_flagged.bed'
    resources:
        mem_mb=_mem_to_mb("100m"),
        runtime=_hms_to_minutes("0:10:00")
    params:
        jobname="{pid}_make_bed",
        src_dir=SRC_DIR
    shell:
        "R --no-save --slave --args {input.candidateRegions} {output.outfile1} {output.outfile2} {output.outfile3} "
        "{wildcards.pid} < {params.src_dir}/make_bed_for_visualization.R"


if len(SAMPLES) == 2:

    rule visualize_zoomed_in:
        input:
            bed=TELOMEREINSERTION_DIR + '/plots/bedfiles/zoomed_in/{pid}_telomere_insertions.bed',
            review_bed=TELOMEREINSERTION_DIR + '/plots/bedfiles/flagged/{pid}_telomere_insertions_review_flagged.bed',
            tumor_bam=lambda wildcards: get_alignment_bam(wildcards.pid, SAMPLES[0]),
            control_bam=lambda wildcards: get_alignment_bam(wildcards.pid, SAMPLES[1]),
            discordant_reads_tumor=TELOMEREINSERTION_DIR + '/tables/{pid}_' + SAMPLES[0] + '_discordant_reads_filtered_with_mapq.tsv',
            discordant_reads_control=TELOMEREINSERTION_DIR + '/tables/{pid}_' + SAMPLES[1] + '_discordant_reads_filtered_with_mapq.tsv',
            clipped_reads_tumor=TELOMEREINSERTION_DIR + '/clipped_reads/{pid}_' + SAMPLES[0] + '_clipped_reads.tsv',
            clipped_reads_control=TELOMEREINSERTION_DIR + '/clipped_reads/{pid}_' + SAMPLES[1] + '_clipped_reads.tsv'
        output:
            TELOMEREINSERTION_DIR + '/plots/zoomed_in/{pid}_done.txt'
        resources:
            mem_mb=_mem_to_mb("3g"),
            runtime=_hms_to_minutes("10:00:00")
        params:
            jobname="{pid}_visualize_zoomed_in",
            sleep_sec_limit=config["sleep_sec_limit"],
            src_dir=SRC_DIR,
            reference_fasta=REFERENCE_FASTA,
            prefix=TELOMEREINSERTION_DIR + "/plots/zoomed_in/",
            review_prefix=TELOMEREINSERTION_DIR + "/plots/flagged/"
        shell:
            """
            sleep $((1 + RANDOM % {params.sleep_sec_limit}))s
            set +u; module load Micromamba/2.0.2-0; set -u
            micromamba run -n telomereEnv python {params.src_dir}/visualize_telomere_insertions.py \
                --control {input.control_bam} \
                --tumor {input.tumor_bam} \
                --ref {params.reference_fasta} \
                --bed {input.bed} \
                --samtoolsbin samtools \
                --colored_reads_tumor {input.discordant_reads_tumor} \
                --colored_reads_control {input.discordant_reads_control} \
                --clipped_reads_tumor {input.clipped_reads_tumor} \
                --clipped_reads_control {input.clipped_reads_control} \
                --prefix {params.prefix} \
                --review_bed {input.review_bed} \
                --review_prefix {params.review_prefix} \
                --outfile {output}
            """

elif len(SAMPLES) == 1:

    rule visualize_zoomed_in:
        input:
            bed=TELOMEREINSERTION_DIR + '/plots/bedfiles/zoomed_in/{pid}_telomere_insertions.bed',
            review_bed=TELOMEREINSERTION_DIR + '/plots/bedfiles/flagged/{pid}_telomere_insertions_review_flagged.bed',
            tumor_bam=lambda wildcards: get_alignment_bam(wildcards.pid, SAMPLES[0]),
            discordant_reads_tumor=TELOMEREINSERTION_DIR + '/tables/{pid}_' + SAMPLES[0] + '_discordant_reads_filtered_with_mapq.tsv',
            clipped_reads_tumor=TELOMEREINSERTION_DIR + '/clipped_reads/{pid}_' + SAMPLES[0] + '_clipped_reads.tsv'
        output:
            TELOMEREINSERTION_DIR + '/plots/zoomed_in/{pid}_done.txt'
        resources:
            mem_mb=_mem_to_mb("3g"),
            runtime=_hms_to_minutes("10:00:00")
        params:
            jobname="{pid}_visualize_zoomed_in",
            sleep_sec_limit=config["sleep_sec_limit"],
            src_dir=SRC_DIR,
            reference_fasta=REFERENCE_FASTA,
            prefix=TELOMEREINSERTION_DIR + "/plots/zoomed_in/",
            review_prefix=TELOMEREINSERTION_DIR + "/plots/flagged/"
        shell:
            """
            sleep $((1 + RANDOM % {params.sleep_sec_limit}))s
            set +u; module load Micromamba/2.0.2-0; set -u
            micromamba run -n telomereEnv python {params.src_dir}/visualize_telomere_insertions.py \
                --tumor {input.tumor_bam} \
                --ref {params.reference_fasta} \
                --bed {input.bed} \
                --samtoolsbin samtools \
                --colored_reads_tumor {input.discordant_reads_tumor} \
                --clipped_reads_tumor {input.clipped_reads_tumor} \
                --prefix {params.prefix} \
                --review_bed {input.review_bed} \
                --review_prefix {params.review_prefix} \
                --outfile {output}
            """
