# TelomereRepeatLoci

Snakemake workflow to detect telomere repeat loci from WGS BAM files.

## What changed

This workflow now supports a **file-path input mode** so you can run it with BAM paths you provide directly (instead of scanning PID directories).

## Setup on NPC

1. Clone this repository on NPC.
2. Make sure these tools are available in your environment/module setup:
   - Snakemake
   - Python with `pysam`
   - R (with required packages used by scripts in `src/`)
   - samtools
   - TelomereHunter
3. Copy and edit the example config:
   - `config_snakemake_TelomereRepeatLoci_example.yaml`

## Input file list (recommended mode)

Create a TSV file (for example `bam_files.tsv`) with explicit BAM paths.
Values must be **tab-separated**.

### Tumor-control mode

Header:

```text
pid	tumor_bam	control_bam
```

Example:

```text
PID001	/npc/path/PID001_tumor.bam	/npc/path/PID001_control.bam
```

### Tumor-only mode

If your config uses only one sample (`samples: [tumor]`), provide:

Header:

```text
pid	tumor_bam
```

Example:

```text
PID001	/npc/path/PID001_tumor.bam
```

## Config notes

In `config_snakemake_TelomereRepeatLoci_example.yaml`:

- Set `bam_files_tsv` to your TSV path (recommended).
- Keep `pids: all` to run all PIDs from the TSV, or set `pids` to a space-separated subset.
- `results_per_pid_dir` and `bam_suffix` are only needed for legacy directory-based mode (`bam_files_tsv: no_file`).

## Run

```bash
snakemake -s $REPO_DIR/Snakefile \
  --configfile $REPO_DIR/config_snakemake_TelomereRepeatLoci_example.yaml \
  --cores 1
```

Use your usual NPC cluster options/profile as needed.

## Outputs

Main outputs are written under `telomereinsertion_dir` from your config, including:

- `candidate_region_tables/*_extended_with_consensus.tsv`
- `plots/zoomed_in/*_done.txt`

## Legacy mode (still supported)

If you set `bam_files_tsv: no_file`, the workflow keeps the original directory-based lookup:

`results_per_pid_dir/{PID}/alignment/{sample}_{PID}{bam_suffix}`
