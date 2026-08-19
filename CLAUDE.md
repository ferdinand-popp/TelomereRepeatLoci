# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Python command-line tool (`telomere-repeat-loci`, entry point `src/telomererepeatloci/main.py`)
that detects telomere repeat loci (telomeric sequence inserted into non-telomeric genomic
locations) from WGS BAM/CRAM files. It chains a fixed sequence of standalone Python scripts under
`src/telomererepeatloci/` as subprocesses to go from TelomereHunter output through discordant/
clipped-read analysis to a final annotated candidate-region table and IGV-like plots.

There used to be a Snakemake + R implementation of this workflow; it has been fully replaced by
this Python CLI. There is no `Snakefile`, no R scripts, and no YAML pipeline config on `main`
anymore — do not reintroduce assumptions from that era. The old Snakemake/R version is preserved
on the `copilot/create-readme-and-update-setup` branch for reference, not for active development.

## Commands

Install dependencies (uv is the expected package manager):

```bash
uv sync
```

Run the workflow (tumor+control):

```bash
uv run telomere-repeat-loci \
  --tumor-bam /path/to/tumor.bam \
  --control-bam /path/to/control.bam \
  --tel-tumor-bam /path/to/tumor_filtered_intratelomeric.bam \
  --tel-control-bam /path/to/control_filtered_intratelomeric.bam \
  --reference-fasta /path/to/reference.fa
```

`--tumor-bam`/`--tel-tumor-bam` are required; everything else (control, blacklist, thresholds,
output dir, visualization) is optional — see `src/telomererepeatloci/main.py:parse_args` for the
full flag list and defaults, or the README's "Command-line options" table.

TelomereHunter (or TelomereHunter2) is **not invoked by this tool** — `--tel-tumor-bam` /
`--tel-control-bam` must already exist before running `telomere-repeat-loci`; it only screens
whatever filtered intratelomeric BAM you point it at.

Lint/format/tests:

```bash
uv run ruff check --fix .
uv run ruff format .
uv run pytest -v   # WIP coverage — see tests/
```

## Architecture

### Orchestration (`main.py`)

`parse_args()` defines the CLI surface, then `process_sample()` runs each pipeline step as a
`subprocess` call (via `run_command()`), passing the previous step's output file as the next
step's input. There is no DAG/scheduler — it's a straight-line sequence, with an `if use_control:`
branch that duplicates the tumor-only discordant-read steps for the control sample.

Key helpers in `main.py`:
- `get_filtered_bam()` — locates a TelomereHunter-style `*_filtered_intratelomeric.bam` (falling
  back to `*_filtered.bam`) in a given TelomereHunter sample directory. Not used by
  `process_sample` directly (callers pass `--tel-tumor-bam`/`--tel-control-bam` explicitly), but
  used to validate/discover TelomereHunter output layouts.
- `extract_pid_from_folder()` — derives the PID from a `<sample>_TelomerCnt_<PID>` directory name;
  raises if the tumor and control PIDs (from their respective TelomereHunter dirs) don't match.
- `get_output_dir()` — if `--output-dir` isn't given, derives
  `<telomerehunter-dir>_TelomereRepeatLoci` as a sibling of the tumor TelomereHunter directory.

### Tumor-only vs tumor+control

Controlled by whether `--control-bam` is supplied. `use_control = bool(args.control_bam)` gates:
control discordant-read extraction, the control-sample steps passed into
`count_discordant_reads.py`, and the `--control`/`--colored_reads_control` flags forwarded to
`visualize_telomere_insertions.py`. When `use_control` is false, `control_discordant_with_mapq` is
the literal string `"NULL"`, which downstream R-derived-but-now-Python scripts treat as "no control
data".

### Pipeline steps (in order, each a standalone script under `src/telomererepeatloci/`)

1. **`find_discordant_reads.py`** — scans an intratelomeric BAM with pysam for reads where the
   mate is mapped elsewhere (not intratelomeric); extracts read name + mate chrom/position. Run
   once per sample (tumor, and control if present).
2. **`add_mate_mapq.py`** — looks up each mate's strand and mapping quality from the *original*
   alignment BAM (not the intratelomeric one) and appends them; drops reads mapping to decoy
   sequences.
3. **`count_discordant_reads.py`** — bins discordant reads (mapq > 30) into strand-specific,
   overlapping 1 kb windows (500 bp step) across the genome, tumor vs. control counts side by
   side; marks windows against `--blacklist` if one is supplied. The windows TSV also carries
   per-window read-name sets (`_tumor_read_names`, `_control_read_names`).
4. **`get_candidate_regions.py`** — filters windows into candidate regions using
   `--tumor-discordant-read-lower-limit` / `--control-discordant-read-upper-limit`, optionally
   dropping blacklisted windows (`--consider-blacklist`); uses the read-name sets to fuse
   overlapping windows on the same chrom/strand when their supporting reads overlap, then drops
   the read-name columns from the output.
5. **`find_fusion_reads.py`** — for each candidate region (tumor sample only), finds soft-clipped
   reads directly and hard-clipped/supplementary-alignment reads via `samtools view`, merges them,
   and counts TTAGGG/CCCTAA repeats in the clipped sequence.
6. **`predict_insertion_sites.py`** — combines candidate regions, clipped reads, and discordant
   reads to predict the exact breakpoint (`insertion_site`), requiring clipping orientation and
   position to be consistent with the discordant-read strand/median position; counts unique-cigar
   clipped reads per position to filter mapping artifacts. See the README's "Planned
   improvements" section for known follow-up work on this script (breakpoint clustering,
   confidence scoring, ambiguous-call reporting — not yet implemented).
7. **`get_consensus.py`** — builds a per-position consensus telomere sequence from clipped reads
   (majority base, else "N") and computes microhomology against the reference genome
   (`--reference-fasta`) at each locus.
8. **`make_bed_for_visualization.py`** — writes zoomed-out and zoomed-in BED files per PID,
   filtering by `--plot-min-support` (minimum `reads_supporting_insertion_pos`).
9. **`visualize_telomere_insertions.py`** — renders IGV-like read/coverage plots around each
   predicted insertion using pysam directly (the `--samtoolsbin` flag is accepted for
   compatibility but isn't load-bearing); highlights discordant reads and telomeric vs.
   non-telomeric clipped bases. Skippable with `--skip-visualization`.

`src/pipeline/tables.py` holds shared TSV column-name constants (`WINDOWS_COLUMNS`,
`FUSION_READS_COLUMNS`, etc.) and `read_tsv`/`write_tsv`/`sanitize_tsv_values` helpers for
consistent, null-byte-stripped TSV I/O — used by the pipeline scripts.

### Outputs

All outputs live under `<output-dir>/` (see `get_output_dir()`):
- `tables/` — discordant-read tables and 1kb-window counts
- `clipped_reads/` — per-sample clipped/fusion-read tables
- `candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended_with_consensus.tsv`
  — the final annotated result table
- `plots/bedfiles/{zoomed_out,zoomed_in}/` — BED files driving visualization
- `plots/zoomed_in/{pid}_done.txt` plus per-locus plot images — visualization output (unless
  `--skip-visualization`)

### Blacklists

`blacklists/PCAWG_tel_ins_blacklist.tsv` and `blacklists/NB_tel_ins_blacklist.tsv` are
cohort-derived false-positive window lists (see `blacklists/README`). An external hg19
problematic-region blacklist (ENCODE/Boyle-Lab) can also be used — see the README's "Running the
workflow" section for the download command.

## Notes when editing

- Coordinate columns written anywhere in the Python workflow are 0-based, half-open
  (pysam/BED-style) — keep new code consistent with this.
- `src/pipeline/` and `src/telomererepeatloci/` are both real, in-use packages (registered in
  `pyproject.toml`'s `[tool.setuptools] packages`) — not placeholders.
- `tests/` has real content (`simple_debug.py`, `debug_wrapper.py`,
  `test_predict_insertion_sites.py`) but coverage is a work in progress; don't assume the test
  suite is exhaustive.
- When changing a pipeline script's CLI (positional args or flags), update the corresponding
  `run_command([...])` call in `main.py:process_sample` — argument order there must match the
  script's `argparse` definition exactly since most steps take positional args, not flags.
