# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Snakemake workflow (`Snakefile`) that detects telomere insertion loci (telomeric sequence
inserted into non-telomeric genomic locations) from WGS BAM files. It chains Python and R
scripts under `src/` into a DAG of per-PID (patient/sample ID) rules, driven by TelomereHunter
output plus discordant/clipped read analysis.

## Commands

Setup the conda/micromamba environment (`telomereEnv`, includes python 3.11, R 4.4, pysam,
samtools, snakemake, and required R packages):

```bash
./setup_micromamba_telomereEnv.sh
```

Run the workflow (copy `config_snakemake_TelomereRepeatLoci_example.yaml` and fill in the
placeholder paths first):

```bash
snakemake -s $REPO_DIR/Snakefile --configfile $REPO_DIR/config_snakemake_TelomereRepeatLoci_example.yaml --cores 1
```

Individual rules can be targeted with `--until <rule_name>` or by requesting a specific output
path, which is useful when iterating on a single stage of the pipeline. Use `-n` (dry run) to
check the DAG before actually launching jobs, especially after editing wildcard logic in the
Snakefile.

There is no test runner or lint config currently checked into the repo (`tests/`, `src/pipeline/`,
and `src/telomererepeatloci/` are present but empty except for `__pycache__` — treat them as
placeholders, not as the real test/package layout).

## Architecture

### Two input modes (Snakefile)

- **TSV mode (preferred)**: `bam_files_tsv` points at a tab/comma-separated file with a `pid`
  column plus `path_to_<sample>_bam` columns (and optionally `path_to_<sample>_intratelomeric_bam`
  for skipping TelomereHunter). Parsed by `_load_bam_paths_from_tsv`.
- **Legacy directory mode**: `bam_files_tsv: no_file`. BAM paths are derived from
  `results_per_pid_dir/{pid}/alignment/{sample}_{pid}{bam_suffix}`.

`get_alignment_bam()` and `get_telomerehunter_intratelomeric_bam()` are the two path-resolution
functions that abstract over these modes — almost every rule's `input:` goes through one of them.

PIDs are validated up front: `_bam_status_table` checks that alignment BAMs exist, and (if
`skip_telomerehunter: true`) that the expected TelomereHunter intratelomeric BAM outputs exist
too. PIDs that fail validation are dropped with a warning rather than failing the whole run.

### Tumor-only vs tumor+control

`config["samples"]` has either 1 entry (`[tumor]`) or 2 (`[tumor, control]`). This flag fans out
through the Snakefile — `run_telomerehunter`, `count_discordant_reads`, `predict_insertion_sites`,
and `visualize_zoomed_in` each have separate rule definitions for the 1-sample and 2-sample case,
selected with `if len(SAMPLES) == 2: ... elif len(SAMPLES) == 1: ...` blocks. When modifying one
of these rules, check whether the parallel single-sample version also needs the change.

### Rule pipeline (in dependency order)

1. `run_telomerehunter` — runs TelomereHunter on the raw BAM(s) to produce
   `*_filtered_intratelomeric.bam`. Skippable via `skip_telomerehunter: true` if these outputs
   already exist (optionally at custom paths supplied via the TSV).
2. `find_discordant_reads` (`src/find_discordant_reads.py`) — pulls read pairs where one mate is
   intratelomeric and the other isn't.
3. `add_mate_mapq` (`src/add_mate_mapq.py`) — annotates those reads with mate mapping quality
   from the original alignment BAM.
4. `count_discordant_reads` (`src/count_discordant_reads.R`) — bins discordant reads into 1kb
   genomic windows, tumor vs. control.
5. `get_candidate_regions` (`src/get_candidate_regions.R`) — filters windows into candidate
   insertion regions using `tumor_discordant_read_lower_limit` /
   `control_discordant_read_upper_limit` and an optional blacklist of excluded 1kb regions.
6. `find_fusion_reads` (`src/find_fusion_reads.R`) — finds soft-clipped/fusion reads at candidate
   regions, per sample.
7. `predict_insertion_sites` (`src/predict_insertion_sites.R`) — combines candidate regions,
   clipped reads, and discordant reads to predict exact breakpoint positions; has separate
   tumor-only and tumor+control code paths (control-sample read counting/ratios at the predicted
   site only apply in 2-sample mode).
8. `get_consensus` (`src/get_consensus.R`) — builds a consensus sequence for each predicted
   insertion and reports sequence microhomology.
9. `make_bed_for_visualization` (`src/make_bed_for_visualization.R`) — writes zoomed-in/zoomed-out
   BED files per PID.
10. `visualize_zoomed_in` (`src/visualize_telomere_insertions.py`) — renders read-level plots
    around each predicted insertion (needs a reference FASTA and `samtools`).

`src/functions.R` holds small shared R helpers (e.g. `compareNA`) sourced by the R scripts above.

### Shared config keys

All rules read from the same `config` dict validated at the top of the Snakefile
(`REQUIRED_CONFIG_KEYS`) — `telomerehunter_dir`, `telomereinsertion_dir`, `src_dir`,
`R_function_file`, `blacklist`, `sleep_sec_limit`, `tumor_discordant_read_lower_limit`, etc. The
`sleep $((1 + RANDOM % {sleep_sec_limit}))s` at the start of several shell blocks staggers
micromamba environment activation across concurrently submitted cluster jobs — this is
intentional, not incidental jitter.

### Outputs

Final per-PID outputs (also the `rule all` targets):
- `{telomereinsertion_dir}/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended_with_consensus.tsv`
- `{telomereinsertion_dir}/plots/zoomed_in/{pid}_done.txt`

## Notes when editing the Snakefile

- `wildcard_constraints` pins `{sample}` to the configured sample names and `{pid}` to
  non-slash strings — needed for unambiguous DAG resolution in modern Snakemake.
- Resource directives use `mem_mb=`/`runtime=` (minutes) via the `_mem_to_mb`/`_hms_to_minutes`
  helpers, which convert from the old-style `"150m"`/`"1g"` and `"HH:MM:SS"` strings still used
  in `params`/comments elsewhere in the file — keep using those helpers rather than hardcoding
  raw minute/MB values.
