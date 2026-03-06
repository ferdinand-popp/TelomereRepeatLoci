# Author: Lina Sieverling

#!/usr/bin/python

# Usage: python /home/sieverli/Code/telomere_insertion_analysis/snakemake_telomere_insertions/src/add_mate_mapq.py \
#                -i <*_tumor_discordant_reads.tsv> \
#                -b <*_merged.mdup.bam> \
#                -o <*_discordant_reads_filtered_with_mapq.tsv>


# Description: - parses tables with telomere insertion reads
#              - skips all reads where mates are mapped to decoy sequences
#              - retrieves mapping quality of mates from original BAM file and adds it to output table
#              - if mate is not found, mapping quality is empty


import sys
import getopt
import numpy
import pysam

# ----------------------------------------------------------------
# read command line args
# ----------------------------------------------------------------
myopts, args = getopt.getopt(sys.argv[1:], "i:b:o:")

for opt, arg in myopts:
    if opt == "-i":
        telomere_insertion_table_file = arg
    elif opt == "-b":
        alignment_bam_file = arg
    elif opt == "-o":
        outfile_path = arg
    else:
        print("Usage: %s -i input_table -b bam_file -o output_file" % sys.argv[0])


# list of chromosomes accepted for output
chromosome_list = [str(i) for i in range(1, 22 + 1)] + ["X", "Y"]

#########################################################################################################################################

# ----------------------------------------------------------------
# read in table containing read names from telomere insertions
# ----------------------------------------------------------------

telomere_insertion_table = numpy.genfromtxt(
    telomere_insertion_table_file,
    skip_header=1,
    delimiter="\t",
    dtype=str,  # robust across numpy versions
    encoding="utf-8",
    comments=None,
)

# If the input has only a single data line, genfromtxt returns 1D.
# Normalize to 2D so the loop below works consistently.
if telomere_insertion_table.ndim == 1:
    telomere_insertion_table = numpy.atleast_2d(telomere_insertion_table)

# ----------------------------------------------------------------
# get mapping quality and strand of mate from original BAM file
# ----------------------------------------------------------------

output = "\t".join(
    ["read_name", "mate_chr", "mate_position", "mate_mapq", "mate_strand"]
)

# Open BAM once, query many times (requires BAM index .bai/.csi)
bam = pysam.AlignmentFile(alignment_bam_file, "rb")

for read in telomere_insertion_table:
    # Expect at least 3 columns: read_name, chromosome, position
    if len(read) < 3:
        continue

    read_name = str(read[0]).strip()
    chromosome = str(read[1]).strip()
    position = str(read[2]).strip()

    # skip mates mapped to decoy sequences
    if chromosome not in chromosome_list:
        continue

    # pysam uses 0-based, half-open intervals for fetch.
    # Input position here is expected to be 1-based (samtools region syntax),
    # so convert to 0-based coordinates.
    try:
        pos1 = int(position)
    except ValueError:
        mapq = ""
        strand = ""
        read_list = [read_name, chromosome, position, mapq, strand]
        output += "\n" + "\t".join(read_list)
        continue

    start0 = pos1 - 1
    end0 = pos1  # fetch one base

    mapq = ""
    strand = ""

    # Equivalent to: samtools view -F 2304 <bam> chr:pos-pos | grep "<read_name>"
    # -F 2304 filters out: 0x100 (secondary) + 0x800 (supplementary)
    try:
        for aln in bam.fetch(chromosome, start0, end0):
            # skip secondary and supplementary alignments (0x100 and 0x800)
            if aln.is_secondary or aln.is_supplementary:
                continue

            # mimic grep for read_name (exact qname match)
            if aln.query_name != read_name:
                continue

            # Found the mate alignment at that locus
            mapq = str(aln.mapping_quality)

            # Equivalent to checking flag & 0x10 for reverse strand
            strand = "-" if aln.is_reverse else "+"
            break
    except (ValueError, OSError):
        # ValueError can happen if contig not present in BAM header
        # OSError can happen for missing index, etc.
        mapq = ""
        strand = ""

    read_list = [read_name, chromosome, position, mapq, strand]
    output += "\n" + "\t".join(read_list)

bam.close()

# ----------------------------------------------------------------
# write output
# ----------------------------------------------------------------

outfile = open(outfile_path, "w")
outfile.write(output)
outfile.close()
