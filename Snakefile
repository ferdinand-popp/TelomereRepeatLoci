"""
Author: Lina Sieverling
Affiliation: DKFZ Heidelberg
Aim: A Snakemake workflow to find telomere insertions
Date: Thu Aug 18 17:46:12 CEST 2016
Run: snakemake -s <Snakefile> --configfile <config.yaml> 
Run Example: 

source activate telomereEnv
snakemake -s /home/sieverli/Code/telomere_insertion_analysis/snakemake_telomere_insertions/Snakefile --configfile /abi/data/sieverling/projects/NB_Telomeres/src/config_snakemake_telomere_insertions.ya[...]

"""
#---------------------------------------------------------------------------------------
# get PIDs
#---------------------------------------------------------------------------------------
from os import listdir
import csv
import os


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


def _load_bam_paths_from_tsv(tsv_file, sample_names):
    """Read pid -> {sample_name: bam_path} from a TSV.

    The pid column is taken positionally (always the first column),
    regardless of what its header is actually called (e.g. "pid",
    "patient_id", "sample_id", ...).

    The bam path columns are expected to be named "path_to_<sample>_bam",
    e.g. for sample_names = ["tumor", "control"] the required columns are
    "path_to_tumor_bam" and "path_to_control_bam".
    """
    bam_paths = {}

    with open(tsv_file, "r") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("bam_files_tsv has no header or is empty: " + tsv_file)

        # First column holds the pid, whatever its header is named.
        pid_column = reader.fieldnames[0]

        bam_columns = {sample_name: f"path_to_{sample_name}_bam" for sample_name in sample_names}
        missing_columns = [col for col in bam_columns.values() if col not in reader.fieldnames]
        if missing_columns:
            raise ValueError("bam_files_tsv is missing required columns: " + ", ".join(missing_columns))

        # start=2 so the first data row (after header) is reported as line 2
        for row_number, row in enumerate(reader, start=2):
            pid = row[pid_column].strip()
            if pid == "":
                print(f"Skipping row {row_number} with empty pid in bam_files_tsv: {tsv_file}")
                continue
            bam_paths[pid] = {}
            for sample_name, bam_column in bam_columns.items():
                bam_paths[pid][sample_name] = row[bam_column].strip()

    return bam_paths


explicit_bam_files_tsv = config.get("bam_files_tsv", "no_file")
use_explicit_bam_paths = _is_enabled_config_path(explicit_bam_files_tsv)

# Config switch:
#   skip_telomerehunter: true/false
# If true, workflow will assume TelomereHunter outputs already exist and will not run run_telomerehunter.
skip_telomerehunter = _parse_bool_config(config.get("skip_telomerehunter", False), default=False)

if use_explicit_bam_paths:
    bam_files_by_pid = _load_bam_paths_from_tsv(explicit_bam_files_tsv, config["samples"])
else:
    bam_files_by_pid = {}

TELOMEREHUNTER_DIR = config["telomerehunter_dir"]
TELOMEREINSERTION_DIR = config["telomereinsertion_dir"]
SRC_DIR = config["src_dir"]
R_FUNCTION_FILE = config["R_function_file"]
SAMPLES = config["samples"]

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
    return f"{TELOMEREHUNTER_DIR}/{pid_name}/{sample_name}_TelomerCnt_{pid_name}/{pid_name}_filtered_intratelomeric.bam"


if config["pids"] == "all":
    if use_explicit_bam_paths:
        pids = sorted(bam_files_by_pid.keys())
    else:
        pids = sorted([i for i in listdir(config["results_per_pid_dir"]) if not i.startswith('.')])
else:
    pids = config["pids"].split(' ')


#---------------------------------------------------------------------------------------
# remove PIDs without bam files
#---------------------------------------------------------------------------------------

pids_remove = []

for pid_name in pids:
    for sample_name in SAMPLES:
        bam_file = None
        if use_explicit_bam_paths:
            if pid_name not in bam_files_by_pid:
                print(f"{pid_name}: no BAM entry found in bam_files_tsv, skipping this pid!")
                pids_remove.append(pid_name)
                break
            bam_file = bam_files_by_pid[pid_name].get(sample_name, "")
            if bam_file == "":
                print(f"{pid_name}: BAM path for {sample_name} is missing in bam_files_tsv, skipping this pid!")
                pids_remove.append(pid_name)
                break
        else:
            bam_file = get_alignment_bam(pid_name, sample_name)

        if not os.path.exists(bam_file):
            print(f"{pid_name}: alignment bam file for {sample_name} sample is missing, skipping this pid!")
            pids_remove.append(pid_name)
            break

# Extra validation when skipping TelomereHunter:
# ensure expected TelomereHunter outputs are present
if skip_telomerehunter:
    for pid_name in pids:
        for sample_name in SAMPLES:
            th_bam = get_telomerehunter_intratelomeric_bam(pid_name, sample_name)
            if not os.path.exists(th_bam):
                print(f"{pid_name}: TelomereHunter output missing for {sample_name} ({th_bam}), skipping this pid!")
                pids_remove.append(pid_name)
                break

pids = [x for x in pids if x not in pids_remove]


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
        expand(TELOMEREHUNTER_DIR + '/{pid}/tumor_TelomerCnt_{pid}/{pid}_filtered_intratelomeric.bam', pid=pids),
        expand(TELOMEREHUNTER_DIR + '/{pid}/control_TelomerCnt_{pid}/{pid}_filtered_intratelomeric.bam', pid=pids_control),
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
    telomerehunter_shell_extra = "-ibc {input[1]} -pl "
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
            "set +u; module load Micromamba/2.0.2-0; micromamba activate telomereEnv; set -u; "
            "module load R/3.4.2; "
            "time telomerehunter -p {wildcards.pid} -o {params.telomerehunter_dir} -ibt {input[0]} {params.extra}-pff all"
else:
    print("skip_telomerehunter=true -> run_telomerehunter rule disabled; assuming existing TelomereHunter outputs.")


#------------------------------------------------------------------
# find discordant reads
#------------------------------------------------------------------

rule find_discordant_reads:
    input:
        TELOMEREHUNTER_DIR + '/{pid}/{sample}_TelomerCnt_{pid}/{pid}_filtered_intratelomeric.bam',
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
        "set +u; module load Micromamba/2.0.2-0; micromamba activate telomereEnv; set -u; "
        "python {params.src_dir}/find_discordant_reads.py -i {input} -o {output}; "
        "set +u; micromamba deactivate; set -u"


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
        "set +u; module load Micromamba/2.0.2-0; micromamba activate telomereEnv; set -u; "
        "python {params.src_dir}/add_mate_mapq.py -i {input.discordant_reads} -b {input.bam} -o {output}; "
        "set +u; micromamba deactivate; set -u"


#------------------------------------------------------------------
# count discordant reads
#------------------------------------------------------------------

paired_t_c_flag = False

if len(SAMPLES) == 2:
    input_list = [
        TELOMEREINSERTION_DIR + '/tables/{pid}_' + SAMPLES[0] + '_discordant_reads_filtered_with_mapq.tsv',
        TELOMEREINSERTION_DIR + '/tables/{pid}_' + SAMPLES[1] + '_discordant_reads_filtered_with_mapq.tsv'
    ]
    tumor_input = "{input[0]}"
    control_input = "{input[1]}"
    paired_t_c_flag = True
elif len(SAMPLES) == 1:
    input_list = [TELOMEREINSERTION_DIR + '/tables/{pid}_' + SAMPLES[0] + '_discordant_reads_filtered_with_mapq.tsv']
    tumor_input = "{input[0]}"
    control_input = "NULL"

if not os.path.exists(config["blacklist"]) and not paired_t_c_flag:
    print("Please provide paired tumor-control samples or a blacklist, otherwise no proper filtering for false positives is possible!")

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
        "R-3.2.2 --no-save --slave --args -t {params.tumor} -c {params.control} -b {params.blacklist} "
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
        control_discordant_read_upper_limit=config["control_discordant_read_upper_limit"]
    shell:
        "R-3.2.2 --no-save --slave --args {input.windowTable} {output.candidateRegions} "
        "{params.tumor_discordant_read_lower_limit} {params.control_discordant_read_upper_limit} "
        "{params.r_function_file} < {params.src_dir}/get_candidate_regions.R"


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
        "R-3.2.2 --no-save --slave --args {input.candidateRegions} {input.bam} {output} "
        "{params.r_function_file} < {params.src_dir}/find_fusion_reads.R"


rule predict_insertion_sites:
    input:
        candidateRegions=TELOMEREINSERTION_DIR + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions.tsv',
        clippedReads=TELOMEREINSERTION_DIR + '/clipped_reads/{pid}_' + SAMPLES[0] + '_clipped_reads.tsv',
        discordantReads=TELOMEREINSERTION_DIR + '/tables/{pid}_' + SAMPLES[0] + '_discordant_reads_filtered_with_mapq.tsv'
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
        "R-3.2.2 --no-save --slave --args {input.candidateRegions} {input.clippedReads} {input.discordantReads} "
        "{output} {params.r_function_file} < {params.src_dir}/predict_insertion_sites.R"


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
        src_dir=SRC_DIR
    message: "--- {wildcards.pid}: get consensus ---"
    shell:
        "R-3.2.2 --no-save --slave --args {input.candidateRegions} {input.clippedReads} {output} "
        "{params.r_function_file} < {params.src_dir}/get_consensus.R"


#------------------------------------------------------------------
# make plots
#------------------------------------------------------------------

rule make_bed_for_visualization:
    input:
        candidateRegions=TELOMEREINSERTION_DIR + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended.tsv'
    output:
        outfile1=TELOMEREINSERTION_DIR + '/plots/bedfiles/zoomed_out/{pid}_telomere_insertions.bed',
        outfile2=TELOMEREINSERTION_DIR + '/plots/bedfiles/zoomed_in/{pid}_telomere_insertions.bed'
    resources:
        mem_mb=_mem_to_mb("100m"),
        runtime=_hms_to_minutes("0:10:00")
    params:
        jobname="{pid}_make_bed",
        src_dir=SRC_DIR
    shell:
        "R-3.2.2 --no-save --slave --args {input.candidateRegions} {output.outfile1} {output.outfile2} "
        "{wildcards.pid} < {params.src_dir}/make_bed_for_visualization.R"


if len(SAMPLES) == 2:

    rule visualize_zoomed_in:
        input:
            bed=TELOMEREINSERTION_DIR + '/plots/bedfiles/zoomed_in/{pid}_telomere_insertions.bed',
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
            prefix=TELOMEREINSERTION_DIR + "/plots/zoomed_in/"
        shell:
            """
            sleep $((1 + RANDOM % {params.sleep_sec_limit}))s
            set +u; module load Micromamba/2.0.2-0; micromamba activate telomereEnv; set -u
            python {params.src_dir}/visualize_telomere_insertions.py \
                --control {input.control_bam} \
                --tumor {input.tumor_bam} \
                --ref /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/1KGRef/hs37d5.fa \
                --bed {input.bed} \
                --samtoolsbin samtools-1.3.1 \
                --colored_reads_tumor {input.discordant_reads_tumor} \
                --colored_reads_control {input.discordant_reads_control} \
                --clipped_reads_tumor {input.clipped_reads_tumor} \
                --clipped_reads_control {input.clipped_reads_control} \
                --prefix {params.prefix} \
                --outfile {output}
            """

elif len(SAMPLES) == 1:

    rule visualize_zoomed_in:
        input:
            bed=TELOMEREINSERTION_DIR + '/plots/bedfiles/zoomed_in/{pid}_telomere_insertions.bed',
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
            prefix=TELOMEREINSERTION_DIR + "/plots/zoomed_in/"
        shell:
            """
            sleep $((1 + RANDOM % {params.sleep_sec_limit}))s
            set +u; module load Micromamba/2.0.2-0; micromamba activate telomereEnv; set -u
            python {params.src_dir}/visualize_telomere_insertions.py \
                --tumor {input.tumor_bam} \
                --ref /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/1KGRef/hs37d5.fa \
                --bed {input.bed} \
                --samtoolsbin samtools-1.3.1 \
                --colored_reads_tumor {input.discordant_reads_tumor} \
                --clipped_reads_tumor {input.clipped_reads_tumor} \
                --prefix {params.prefix} \
                --outfile {output}
            """
