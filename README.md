# TelomereRepeatLoci
*Python command-line workflow for detection of telomere repeat loci from WGS data*

This Python command-line workflow detects telomere repeat loci within cancer genomes from WGS data. The input are BAM or CRAM files from a tumor and a control sample (if available). In the first step, telomeric reads are extracted using the tool [TelomereHunter](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2851-0). From the extracted telomeric reads, discordant reads are retrieved, where one mate is intratelomeric and the other mate is mapped to the chromosome. In regions with discordant reads, it then searches for clipped reads to find the precise position of the inserted telomere sequence.

<p align="center">
  <img src="resources/images/telomere_repeat_locus_schematic.png" alt="Detection of telomere repeat loci" width="700" />
</p>

<p align="center">
  <img src="resources/images/telomere_repeat_locus_example.png" alt="Example of telomere repeat locus" width="500" />
</p>

If you are using the workflow, please cite:


> **TelomereHunter – in silico estimation of telomere content and composition from cancer genomes** <br>
Lars Feuerbach, Lina Sieverling, Katharina I. Deeg, Philip Ginsbach, Barbara Hutter, Ivo Buchhalter, Paul A. Northcott, Sadaf S. Mughal, Priya Chudasama, Hanno Glimm, Claudia Scholl, Peter Lichter, Stefan Fröhling, Stefan M. Pfister, David T. W. Jones, Karsten Rippe & Benedikt Brors <br>
*BMC Bioinformaticsvolume 20, Article number: 272 (2019)*


> **Alternative lengthening of telomeres in childhood neuroblastoma from genome to proteome** <br>
Sabine A. Hartlieb, Lina Sieverling, Michal Nadler-Holly, Matthias Ziehm, Umut H. Toprak, Carl Herrmann, Naveed Ishaque, Konstantin Okonechnikov, Moritz Gartlgruber, Young-Gyu Park, Elisa Maria Wecht, Kai-Oliver Henrich, Larissa Savelyeva, Carolina Rosswog, Matthias Fischer, Barbara Hero, David T.W. Jones, Elke Pfaff, Olaf Witt, Stefan M. Pfister, Jan Koster, Richard Volckmann, Katharina Kiesel, Karsten Rippe, Sabine Taschner-Mandl, Peter Ambros, Benedikt Brors, Matthias Selbach, Lars Feuerbach, Frank Westermann <br>
*under revision*

The workflow was also used in the following publication (where telomere repeat loci were termed "telomere insertions"):

> **Genomic footprints of activated telomere maintenance mechanisms in cancer** <br>
Lina Sieverling, Chen Hong, Sandra D. Koser, Philip Ginsbach, Kortine Kleinheinz, Barbara Hutter, Delia M. Braun, Isidro Cortés-Ciriano, Ruibin Xi, Rolf Kabbe, Peter J. Park, Roland Eils, Matthias Schlesner, PCAWG-Structural Variation Working Group, Benedikt Brors, Karsten Rippe, David T. W. Jones, Lars Feuerbach & PCAWG Consortium <br>
*Nature Communications volume 11, Article number: 733 (2020)*



---


### Detailed description of individual steps in the workflow

<p align="center">
  <img src="resources/images/TelomereRepeatLoci_workflow.png" alt="TelomereRepeatLoci workflow" width="700" />
</p>

#### 1. Run TelomereHunter

  Information on TelomereHunter can be found in the [publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2851-0).
  If intratelomeric BAM outputs are already present (e.g. precomputed with TelomereHunter2), the workflow reuses them and skips rerunning TelomereHunter.
  

#### 2. Find candidate regions with discordant reads

Note: The windows TSV retains per-window discordant read-name sets (`_tumor_read_names`, `_control_read_names`) so that read counts stay exact across the window-merge step below. `get_candidate_regions.py` itself no longer does any window fusion — it only thresholds the already-merged windows and drops the read-name-set columns from the candidate output.

The python script `find_discordant_reads.py` goes through the intratelomeric read BAM file produced by TelomereHunter with the module [pysam](https://pysam.readthedocs.io/en/latest/), which is a wrapper for [SAMtools](http://www.htslib.org). Reads that fulfill the following criteria are considered discordant intratelomeric reads: 1) mate is mapped and its reference ID is known, 2) mate is not an intratelomeric read. For each discordant intratelomeric read, the read name as well as the chromosome and position of the mate is extracted from the QNAME, RNEXT and PNEXT fields of the SAM format, respectively. The results are saved in a table. In the script `add_mate_mapq.py`, the strand and the mapping quality of the primary alignment of each chromosomal mate from the discordant read table are retrieved from the alignment BAM file and added to the table — reported as the mate's own aligned position (not the position of the mate's own mate, which loops back to the original discordant read). To keep this lookup fast, the mate BAM-fetch calls are sorted by `(chrom, position)` before running, and the reference length per chromosome is cached, so access moves forward through the coordinate-sorted BAM instead of jumping randomly. Reads mapping to decoy sequences are removed. Until this point, the scripts are run individually for the tumor and the control sample. `count_discordant_reads.py` summarizes the number of discordant reads in the tumor and control sample. For this, the genome is split into strand-specific, non-overlapping 1 kb windows (floor-based on read position). For each window, the number of discordant reads with a mapping quality of over 30 is counted. If a blacklist of false positive regions is provided, windows contained in the blacklist are marked. Adjacent same-chrom/strand windows are then merged into a single region whenever *both* have nonzero tumor discordant-read support, **before** the tumor/control thresholds below are applied — mirroring the original R implementation — so a real locus whose support straddles a window boundary isn't undercounted below the threshold. The script `get_candidate_regions.py` filters the (already-merged) list of windows to get candidate regions of somatic telomere repeat loci. Candidate regions must contain a minimum number of discordant reads in the tumor sample (set to 3 and 4 for the PCAWG and neuroblastoma analysis, respectively) and a maximum number of discordant reads in the control sample (usually 0). If specified by the user, windows contained in the blacklist are removed. This step is especially important to rule out false positives if no control sample is available.

#### 3. Find precise locus with clipped reads

For each candidate region obtained in the previous step, clipped reads that span the telomere repeat locus junction site are searched for with `find_fusion_reads.py`. First, the script searches for soft-clipped sequences. For this, all reads in the candidate region +/- 300 bp are extracted, including the read name, sequence, position, cigar and flag. The reads are then filtered and only those containing an "S" in the cigar string are kept. Moreover, the end position of the clipped sequence is extracted. Next, hard-clipped reads are obtained by searching for supplementary alignments in the candidate region +/- 300 bp. If the candidate region is on the (+) strand, supplementary alignments are extracted with `samtools view -f 2048 -F 16`, i.e. reads that are supplementary alignments and not on the reverse strand. For those on the (-) strand, the command "samtools view -f 2064" was used, i.e. reads that are supplementary alignments and on the reverse strand. In contrast to soft-clipped reads, the SAM format does not contain the clipped sequence of supplementary alignments in the SEQ field. Therefore, the full sequence must be retrieved from the primary alignment of the read. For this, the position and strand of the primary alignment is obtained from the SA tag of the supplementary alignment; this lookup runs against a *separate* pysam `AlignmentFile` handle from the one driving the region scan, since htslib does not support two concurrent `fetch()` iterators sharing one handle. Distinct primary loci are resolved once each (stopping as soon as every needed read at that locus is found), rather than caching every read seen at a locus, so a repeat-collapsed/high-depth primary locus doesn't blow up memory or silently fall back to truncated evidence. If supplementary and primary alignments are on opposite strands, the sequence is reverse complemented. All information on soft- and hard-clipped sequences is then merged into one table (written incrementally, one row at a time, to bound peak memory on high-depth regions). By taking the length of the clipped sequences into account — for a hard-clipped read, correctly consuming the "H" cigar op once the primary read's full sequence has been substituted in, not just "S" — the clipped parts of the read sequences are obtained. For each read, the number of TTAGGG and CCCTAA repeats in the clipped sequence are counted. The position of the clipped sequence, i.e. whether sequences were clipped in the upstream or downstream end of the read alignment, is inferred from the cigars.
The exact position of the telomere repeat locus is obtained from the position of the clipped reads by `predict_insertion_sites.py`. For this, only reads that contain at least one telomeric repeat in the clipped sequence are taken into account. If the discordant reads map to the (+) strand, the clipped parts of the reads need to be at the end of the aligned read. If the discordant reads map to the (-) strand, clipping needs to occur at the start of the reads. Moreover, the clipping position needs to be downstream or upstream of the median discordant read positions, respectively. Finally, a frequency table of the number of clipped reads ending or starting at different positions, respectively, is calculated. Here, only clipped reads with unique cigars at each position are counted. This filter was included because mapping artifacts were observed where all clipped reads mapped to exactly the same position. For each candidate region, the total number of clipped reads supporting the telomere repeat locus, the orientation of the telomere sequence (TTAGGG or CCCTAA on the forward strand) and the total number of TTAGGG and CCCTAA counts in the fusion reads is reported.

##### Planned improvements for insertion-site prediction

These updates target `predict_insertion_sites.py` to make breakpoint selection more robust and to surface additional QC signals in the output table.

- Cluster soft-clip positions within a small tolerance (e.g., +/-5 bp) and use the cluster median as `insertion_site`.
- Report `insertion_site_spread_bp` (max-min in cluster) as a simple uncertainty metric.
- Replace the median-mate-position cutoff with a discordant support window (e.g., q10-q90 with padding), and require the clip cluster to fall near that interval.
- Add table-only filtering and weighting: `min_mapq=30` when available and `min_clipped_len=15` as a soft filter based on clipped-sequence length.
- Keep off-orientation clips, but downweight them and report `support_expected_orientation` and `support_unexpected_orientation`.
- Improve tie handling: report `ambiguous_insertion_site` plus `insertion_site_candidates` instead of dropping calls.
- Add a simple `insertion_confidence` score from cluster support, unique cigars, telomere motif counts, spread, and second-best support.

Minimal new output fields:

- insertion_site
- insertion_site_spread_bp
- reads_supporting_insertion_pos
- unique_cigars_supporting
- ambiguous_insertion_site
- insertion_site_candidates
- insertion_confidence

#### 4. Construct telomeric sequences at the telomere repeat loci
From the clipped sequences at the telomere repeat loci, the telomere sequences flanking each locus are reconstructed with `get_consensus.py`. For each position in the clipped sequences, the frequency of each base is calculated. If a base has a frequency of at least 0.65, this base is used for the consensus sequence. Otherwise, it is set to "N".
Assuming that the telomere sequence at the repeat locus consists exclusively of t-type repeats, microhomology between the reference genome and the telomere sequence can be determined. For this, the reference genome sequence 20 bp upstream of the telomere repeat locus is extracted. The t-type telomere repeat of the inserted telomere sequence that is closest to the locus is extended and each base pair is compared to that of the reference genome. Every match is counted as a base pair of sequence homology between the reference genome and the telomere sequence. As soon as a base pair does not match, the microhomology is disrupted and further homology is not considered. If the bases upstream of the first t-type repeat in the inserted telomere sequence do not match an incomplete t-type repeat, the microhomology cannot be determined and is set to "?". The information on the telomere consensus sequence and the base pairs of microhomology is added to the telomere repeat locus table.

<p align="center">
  <img src="resources/images/microhomology_examples.jpeg" alt="Microhomology examples" width="400" />
</p>


#### 5. Make IGV-like plot
To rule out remaining false positives, each telomere repeat locus should be checked manually. To facilitate this process, Integrative Genomics Viewer (IGV)-like plots of each telomere repeat locus are made. The script `make_bed_for_visualization.py` makes tables in BED format that contain the reference genome start and end positions used for the plots, which are 100 bp up- and downstream of the telomere repeat locus. This table is then used as input for the script visualize_telomere_insertions.py. The script was adapted from [here](https://github.com/DKFZ-ODCF/IndelCallingWorkflow/blob/master/resources/analysisTools/indelCallingWorkflow/visualize.py). Given the alignment BAM files of the tumor and control sample, the script generates a PDF file for each genomic region in the input BED file, in which the reads surrounding the telomere repeat loci are displayed in the tumor and in the control sample. Moreover, panels with the coverage in the region are plotted. Several new features were added to the original script: the discordant reads obtained in previous steps of the TelomereRepeatLoci pipeline are highlighted, hard- clipped bases are obtained from the primary alignments and displayed, non-telomeric clipped bases are transparent, while telomeric clipped bases remain opaque. With the resulting images, tumor and control sample can easily be compared and artifact-prone regions, e.g. with a lot of clipped reads, can be identified.

#### 6. Site-level confidence diagnostics

For each candidate region with a predicted `insertion_site`, `assess_site_confidence.py` adds a set of
diagnostic columns to help spot two common sources of false positives: generically noisy/repetitive loci,
and germline/artifactual signal that is also present in the control. These columns are added on top of the
consensus table — nothing is dropped or filtered by this step.

`find_fusion_reads.py` is now also run against the control BAM (mirroring the existing tumor step), producing
a control clipped-reads table, so control's clipped reads at each site can be examined directly rather than
only contributing to the discordant-read-count filter from step 2.

For the tumor sample, reads covering a `--site-window` bp flank (default 100 bp, matching the zoomed-in plot
window) around `insertion_site` are fetched from the original `--tumor-bam` (duplicates excluded), giving
`all_reads_at_site`. Among the clipped reads already found for that candidate region, those whose own clip
position (the "end" of a downstream/`+` clip, or "start" of an upstream/`-` clip) lands exactly at
`insertion_site` are counted as `clipped_reads_at_site`, with the telomeric subset as `telo_clipped_reads_at_site`.
`tumor_noise_ratio = (clipped_reads_at_site - telo_clipped_reads_at_site) / all_reads_at_site` — a high value means
most of the clipping activity at that exact base is unrelated to a telomere junction, i.e. a generically messy
locus (including chromosome-end fade-outs, where `all_reads_at_site` naturally drops).

For the control sample (if `--control-bam` is given), the same depth count gives `control_all_reads_at_site`,
and clipped reads at the exact `insertion_site`, matching the same strand-implied orientation as the tumor call
and telomeric, are counted as `control_telo_clipped_at_insertion_site`. If both tumor and control have such
reads, their clipped sequences are compared with a Hamming distance over the 12 bases starting at the breakpoint
and reading away from it — oriented consistently for both strands (for a downstream/`+` clip that's the clip's
first 12 bases as-is; for an upstream/`-` clip it's the *last* 12 bases of the clip, reversed, since the
breakpoint-adjacent base there is the last character of the raw clipped substring). The minimum distance over
all tumor/control read pairs is reported as `control_min_seq_distance_to_tumor` — a low value alongside a
non-zero `control_telo_clipped_at_insertion_site` is a strong germline/artifact signal, since it means control
shows essentially the same clip at the same base.

Because telomere repeats are just `TTAGGG`/`CCCTAA` repeated, this distance alone is not very discriminating
(most telomeric clips look somewhat alike within a few bases) — read it together with the raw count, not alone.

Result file:
```
<output-dir>/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended_with_confidence.tsv
```
This is the consensus table (step 4) with the columns above appended.

#### 7. Drop low-confidence regions by threshold

`filter_by_site_confidence.py` always runs on the table from step 6 and writes a filtered table alongside
it (the unfiltered table from step 6 is also always kept, so nothing is lost). A region is dropped if:

- it has no predicted `insertion_site` at all — there is no locus to review or plot, so it can't be
  carried into the "table to use", or
- `reads_supporting_insertion_pos` is below `--plot-min-support` (default `2`, the same threshold
  `make_bed_for_visualization.py` uses to decide whether a region gets a BED entry/plot at all) — this
  keeps the filtered table from carrying rows that can never be reviewed because they'll never be
  plotted, or
- `tumor_noise_ratio` is present and exceeds `--max-tumor-noise-ratio` (default `0.8`), or
- control shows telomeric clipped reads at the insertion site (`control_telo_clipped_at_insertion_site > 0`)
  **and either** their sequence is close enough to tumor's to count as a match
  (`control_min_seq_distance_to_tumor <= --control-max-seq-distance`, default `2`), **or** there are simply
  too many of them regardless of sequence similarity (`control_telo_clipped_at_insertion_site >
  --control-max-telo-clipped-at-site`, default `3`).

For regions that do have an `insertion_site`, blank/missing diagnostics (e.g. no control BAM) never cause
a drop on their own — only a computed value that actually exceeds a threshold does. These defaults are a
starting point, not validated thresholds — tune them via
`--max-tumor-noise-ratio`/`--control-max-seq-distance`/`--control-max-telo-clipped-at-site` against your
own data by comparing the kept/dropped regions against manual review (step 8) and the plots.

Result file:
```
<output-dir>/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended_with_confidence_filtered.tsv
```

#### 8. Manual review of predicted insertions

The pipeline does not classify candidate loci as true or false positives on its own — a human still has
to look at each one. This step is manual and happens outside the CLI, but the outputs below are what
you review and where the final calls end up.

**What to look at.** `make_bed_for_visualization.py` reads from the confidence-filtered table (step 7's
output), so only regions that survived the confidence filter get a BED entry and a plot. For each PID, it
writes two BED files under `plots/bedfiles/`: `zoomed_in/{pid}_telomere_insertions.bed` (±100 bp around
`insertion_site`, used to render the plots) and `zoomed_out/{pid}_telomere_insertions.bed` (±500 bp,
generated for a wider manual look in IGV/another BAM viewer if the 100 bp plot isn't enough context — it
is not itself rendered by `visualize_telomere_insertions.py`). Both BEDs only include regions with at
least `--plot-min-support`
supporting reads (default 2).

For every region in the zoomed-in BED, `visualize_telomere_insertions.py` renders one plot at
`plots/zoomed_in/{pid}_{chrom}_{insertion_site}.pdf`. Once all plots for a PID finish, a
`plots/zoomed_in/{pid}_done.txt` marker file is written; its presence is how you (or a wrapper script)
can tell that visualization for a sample completed rather than crashed partway through. Absence of
`_done.txt` after a run means visualization did not finish and some loci may be missing plots.

Review each plot the same way the pipeline's design intends: compare the tumor and control read/coverage
tracks, check whether the highlighted discordant reads and telomeric (opaque) vs. non-telomeric
(transparent) clipped bases look like a real fusion junction, and flag loci with the artifact patterns
called out above (e.g. many clipped reads collapsing onto the exact same position, or clipping/coverage
signal that is also present in the control).

**Where results end up.** `review_insertion_plots.ipynb` (repo root) drives this review: point its
`RESULTS_ROOT` at the parent directory holding one `<sample>_TelomereRepeatLoci` output dir per patient,
and it discovers every patient's confidence-filtered table and zoomed-in plots, flattens them into one
Pass/Fail/note review queue across all patients, and — on its final cell — left-joins your annotations
back onto each patient's own table by `(chrom, insertion_site)`, writing a sibling
`*_confidence_filtered_annotated.tsv` per patient (never overwriting the pipeline's own output). Progress
is resumable: re-running the notebook skips plots already annotated in its `insertion_site_manual_annotations.csv`.
A small `insertion_site_summary_per_pid.csv` with pass counts per PID is written alongside it. The result table to carry forward — pre-annotation, or with your `manual_decision`/`manual_note`
columns if you've run the notebook — is:

```
<output-dir>/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended_with_confidence_filtered.tsv
```

with low-confidence regions already dropped per step 7's thresholds. The unfiltered
`..._extended_with_confidence.tsv` (step 6) is also always kept, in case a region you expected to see
was dropped and you want to check why before adjusting the thresholds.

This is the final, fully annotated table: one row per predicted telomere repeat locus, combining the
candidate-region/discordant-read counts, the predicted `insertion_site` and `reads_supporting_insertion_pos`,
`sum_TTAGGG_count`/`sum_CCCTAA_count` and `repeat_forward` (insertion orientation), the
`consensus`/`flanking_seq`/`bp_microhomology` columns from the consensus step, and the
`tumor_noise_ratio`/`control_telo_clipped_at_insertion_site`/`control_min_seq_distance_to_tumor` diagnostic
columns from step 6 (see above — already used by step 7's thresholds to drop rows, but still worth
eyeballing alongside the plots to judge whether a threshold needs tuning). Match a plot back to its row in this table via `chrom` and
`insertion_site` (the same values used in the plot's filename and the zoomed-in BED's `pos` column).

---

## Running the workflow

Install from GitHub (no PyPI release yet):
```bash
pip install git+https://github.com/ferdinand-popp/TelomereRepeatLoci.git
```

Or with uv:
```bash
uv pip install git+https://github.com/ferdinand-popp/TelomereRepeatLoci.git
```

Using uv (python package manager) is recommended to run the workflow. After cloning the repository, run `uv sync` to install the required dependencies.
```bash
git clone https://github.com/ferdinand-popp/TelomereRepeatLoci.git
uv sync
```

The workflow is now started directly via Python (no Snakemake/YAML config required):

```bash
uv run telomere-repeat-loci \
  --tumor-bam /path/to/tumor_input.bam \
  --control-bam /path/to/control_input.bam \
  --tel-tumor-bam /path/to/tumor_intratelomeric.bam \
  --tel-control-bam /path/to/control_intratelomeric.bam \
  --blacklist /path/to/blacklist.tsv \
  --tumor-discordant-read-lower-limit 3 \
  --control-discordant-read-upper-limit 0 \
  --consider-blacklist \
  --reference-fasta /path/to/reference.fa
```

The `--blacklist` file lists windows to exclude as likely false positives. Cohort-derived blacklists
(`PCAWG_tel_ins_blacklist.tsv`, `NB_tel_ins_blacklist.tsv`) are provided in `blacklists/`; see
`blacklists/README` for details. For an external hg19 blacklist of problematic genomic regions,
download the ENCODE/Boyle-Lab blacklist:
```bash
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz
```

Minimal single-sample run (reference FASTA is still required for microhomology analysis and visualization):

```bash
uv run telomere-repeat-loci \
  --tumor-bam /path/to/tumor_input.bam \
  --tel-tumor-bam /path/to/tumor_intratelomeric.bam \
  --reference-fasta /path/to/reference.fa
```

By default, output files are written to a new sibling directory outside the provided
TelomereHunter output directory:
`<telomerehunter-dir>_TelomereRepeatLoci`.
You can still override this with `--output-dir`.

## Manual annotations of insertion candidates

Use the Jupyter Notebook to screen and annotate all sites from a run and get the results documented
```bash
uv sync --group notebook

uv run jupyter lab
```
Use the exposed endpoint with jupyter lab and adapt the paths to the plot files and the _confidence_filtered.tsv files where results should be documented. Then use the PASS or FAIL buttons or fail with comment.


### Command-line options

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--tumor-bam` | Yes | — | Tumor BAM/CRAM file |
| `--tel-tumor-bam` | Yes | — | TelomereHunter-filtered intratelomeric tumor BAM, used for discordant-read screening |
| `--control-bam` | No | `""` | Control BAM/CRAM; omit to run in tumor-only mode |
| `--tel-control-bam` | No | `""` | Filtered intratelomeric control BAM |
| `--output-dir` | No | derived | Output directory; defaults to `<telomerehunter-dir>_TelomereRepeatLoci` next to the tumor TelomereHunter folder |
| `--tumor-sample-name` | No | `tumor` | Label used in output file names for the tumor sample |
| `--control-sample-name` | No | `control` | Label used in output file names for the control sample |
| `--blacklist` | No | `no_file` | TSV of 1 kb windows to exclude as likely false positives |
| `--consider-blacklist` | No | off | Apply the blacklist when filtering candidate regions |
| `--tumor-discordant-read-lower-limit` | No | `3.0` | Minimum discordant reads (mapq > 30) required in the tumor sample per candidate region |
| `--control-discordant-read-upper-limit` | No | `0.0` | Maximum discordant reads allowed in the control sample per candidate region |
| `--reference-fasta` | No | `""` | Reference FASTA; required for microhomology analysis and visualization (unless `--skip-visualization`) |
| `--skip-visualization` | No | off | Skip generation of zoomed-in IGV-like plots |
| `--plot-min-support` | No | `2.0` | Minimum `reads_supporting_insertion_pos` required to include a region in plot BEDs |
| `--site-window` | No | `100` | Flank (bp) around `insertion_site` used for the site-level confidence diagnostics (depth and clipped-read collection) |
| `--max-tumor-noise-ratio` | No | `0.8` | Max allowed `tumor_noise_ratio` before a region is dropped by the confidence-filtering step |
| `--control-max-seq-distance` | No | `2` | Max Hamming distance (12 bp at the breakpoint) between control/tumor clipped sequences to count as a germline match, used by the confidence-filtering step |
| `--control-max-telo-clipped-at-site` | No | `3` | Max telomeric clipped reads allowed in control at the insertion site regardless of sequence match, used by the confidence-filtering step |
| `--samtoolsbin` | No | `samtools` | Path/name of the samtools binary (kept for compatibility; visualization uses pysam directly) |

## Notes

- `main.py`'s `parse_args()` defines the CLI surface, then `process_sample()` runs each pipeline
  step above as a subprocess in a fixed, straight-line sequence (no DAG/scheduler), passing each
  step's output file as the next step's input. An `if use_control:` branch duplicates the
  tumor-only discordant-read/clipped-read steps for the control sample when `--control-bam` is
  given; when it isn't, several downstream scripts receive the literal string `"NULL"` in place of
  a control file path, which they treat as "no control data" rather than a missing-file error.
- `src/pipeline/tables.py` holds shared TSV column-name constants (`WINDOWS_COLUMNS`,
  `FUSION_READS_COLUMNS`, etc.) and `read_tsv`/`write_tsv`/`sanitize_tsv_values` helpers for
  consistent, null-byte-stripped TSV I/O, used by all the pipeline scripts above.
- If you change a pipeline script's CLI (positional args or flags), update the corresponding
  `run_command([...])` call in `main.py:process_sample` to match — argument order there must line
  up exactly with the script's own `argparse` definition, since most steps take positional args,
  not flags.
- There used to be a Snakemake + R implementation of this workflow; it has been fully replaced by
  this Python CLI (no `Snakefile`, no R scripts, no YAML pipeline config on `main`). The old
  Snakemake/R version is preserved on the `legacy/r-snakemake-workflow` branch for reference, not
  for active development.
- `--tel-tumor-bam` and `--tel-control-bam` can be any BAM you want to screen (not limited to TelomereHunter outputs).
- Discordant read screening uses non-overlapping, floor-based 1 kb windows, with adjacent
  same-chrom/strand windows merged before thresholding when both have tumor support (see step 2
  above).
- When `--control-bam` is given, the discordant-read and clipped-read steps run the tumor and
  control branches concurrently (up to 2 threads) instead of back-to-back, since neither branch
  needs the other's output until a later joining step. Each branch's subprocess output is
  captured and printed as one labeled block (`[tumor]`/`[control]`) after it finishes, instead of
  interleaved live, to keep logs readable. There's no flag for this — cap available cores at
  job-submission time (e.g. a scheduler's `--cpus-per-task`) if you want to force single-core
  execution.
- All coordinate columns written by the Python workflow are 0-based, half-open (pysam/BED-style).
- Visualization uses pysam directly; the `--samtoolsbin` flag is kept for compatibility.
- run tests with `uv run pytest -v` -> WIP
- `uv run ruff check --fix .`
- `uv run ruff format .`

## Testcase
Download the sample CRAM (tumor BAM input), its index, and the reference FASTA (needed for visualization and microhomology) from 1000genomes:

```bash
mkdir -p data
curl -L -C - \
  -o data/HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram \
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00152/alignment/HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram

curl -L -C - \
  -o data/HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram.crai \
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00152/alignment/HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram.crai

curl -L -C - \
  -o data/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

Then run TelomereHunter2 on the CRAM to generate the intratelomeric BAM files and use those as inputs to `telomere-repeat-loci` (see `tests/simple_debug.py` for the expected paths).

Example [TelomereHunter2](https://github.com/ualbertalab/TelomereHunter2) command (adjust paths as needed):

```bash
uv run telomerehunter2 \
  -ibt data/HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram \
  -p HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage \
  -o results/ \
  -b hg38
```

Example run (tumor-only, lowered lower read limit for regions on small testing file --plot-min-support 2):

```bash
uv run telomere-repeat-loci \
  --tumor-bam data/HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram \
  --tel-tumor-bam results/HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage/tumor_TelomerCnt_HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage/HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage_filtered_intratelomeric.bam \
  --reference-fasta data/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --plot-min-support 2
```

Result file `/results/.../candidate_region_tables/..._telomere_insertions_candidate_regions_extended_with_confidence_filtered.tsv` should have regions and plots should be generated for ChrX and Chr22
