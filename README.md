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

Note: The windows TSV now retains per-window discordant read-name sets (`_tumor_read_names`, `_control_read_names`). `get_candidate_regions.py` uses these sets to fuse overlapping windows on the same chrom/strand when the supporting reads overlap, then drops the read-name columns from the candidate output.

The python script `find_discordant_reads.py` goes through the intratelomeric read BAM file produced by TelomereHunter with the module [pysam](https://pysam.readthedocs.io/en/latest/), which is a wrapper for [SAMtools](http://www.htslib.org). Reads that fulfill the following criteria are considered discordant intratelomeric reads: 1) mate is mapped and its reference ID is known, 2) mate is not an intratelomeric read. For each discordant intratelomeric read, the read name as well as the chromosome and position of the mate is extracted from the QNAME, RNEXT and PNEXT fields of the SAM format, respectively. The results are saved in a table. In the script `add_mate_mapq.py`, the strand and the mapping quality of the primary alignment of each chromosomal mate from the discordant read table are retrieved from the alignment BAM file and added to the table. Reads mapping to decoy sequences are removed. Until this point, the scripts are run individually for the tumor and the control sample. `count_discordant_reads.py` summarizes the number of discordant reads in the tumor and control sample. For this, the genome is split into strand-specific, 1 kb windows with a 500 bp step. For each window, the number of discordant reads with a mapping quality of over 30 is counted. If a blacklist of false positive regions is provided, windows contained in the blacklist are marked. The overlapping 1 kb windows already account for discordant reads that fall near window boundaries, so neighboring windows are not merged. The script `get_candidate_regions.py` filters the list of windows to get candidate regions of somatic telomere repeat loci. Candidate regions must contain a minimum number of discordant reads in the tumor sample (set to 3 and 4 for the PCAWG and neuroblastoma analysis, respectively) and a maximum number of discordant reads in the control sample (usually 0). If specified by the user, windows contained in the blacklist are removed. This step is especially important to rule out false positives if no control sample is available.

#### 3. Find precise locus with clipped reads

For each candidate region obtained in the previous step, clipped reads that span the telomere repeat locus junction site are searched for with `find_fusion_reads.py`. First, the script searches for soft-clipped sequences. For this, all reads in the candidate region +/- 300 bp are extracted, including the read name, sequence, position, cigar and flag. The reads are then filtered and only those containing an "S" in the cigar string are kept. Moreover, the end position of the clipped sequence is extracted. Next, hard-clipped reads are obtained by searching for supplementary alignments in the candidate region +/- 300 bp. If the candidate region is on the (+) strand, supplementary alignments are extracted with `samtools view -f 2048 -F 16`, i.e. reads that are supplementary alignments and not on the reverse strand. For those on the (-) strand, the command "samtools view -f 2064" was used, i.e. reads that are supplementary alignments and on the reverse strand. In contrast to soft-clipped reads, the SAM format does not contain the clipped sequence of supplementary alignments in the SEQ field. Therefore, the full sequence must be retrieved from the primary alignment of the read. For this, the position and strand of the primary alignment is obtained from the SA tag of the supplementary alignment. `samtools view` is used on the alignment BAM file to extract reads in the region of the primary alignment, which are further filtered by read name and strand to obtain the read sequence of the primary alignment. If supplementary and primary alignments are on opposite strands, the sequence is reverse complemented. All information on soft- and hard-clipped sequences is then merged into one table. By taking the length of the clipped sequences into account, the clipped parts of the read sequences are obtained. For each read, the number of TTAGGG and CCCTAA repeats in the clipped sequence are counted. The position of the clipped sequence, i.e. whether sequences were clipped in the upstream or downstream end of the read alignment, is inferred from the cigars.
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

## Notes

- The scripts in `src/` are orchestrated by `main.py`.
- Legacy R helper scripts were removed; the workflow now uses Python scripts only.
- `--tel-tumor-bam` and `--tel-control-bam` can be any BAM you want to screen (not limited to TelomereHunter outputs).
- Discordant read screening uses overlapping 1 kb windows with a 500 bp step.
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

Result file `/results/.../candidate_region_tables/..._telomere_insertions_candidate_regions_extended_with_consensus.tsv` should have regions and 
