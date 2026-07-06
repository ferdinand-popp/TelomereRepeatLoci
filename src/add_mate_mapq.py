#!/usr/bin/env python
# Author: Lina Sieverling (original), optimized version
#
# Usage: python add_mate_mapq_optimized.py \
#                -i <*_tumor_discordant_reads.tsv> \
#                -b <*_merged.mdup.bam> \
#                -o <*_discordant_reads_filtered_with_mapq.tsv>
#
# Description: - parses tables with telomere insertion reads
#              - skips all reads where mates are mapped to decoy sequences
#              - retrieves mapping quality of mates from original BAM file and adds it to output table
#              - if mate is not found, mapping quality is empty
#
# Speedup notes vs. original:
#   - Uses pysam for indexed random-access BAM lookups instead of spawning
#     `samtools view | grep` once per row (this was the dominant cost -
#     one subprocess fork/exec per input row).
#   - Requires the BAM to be indexed (.bai present) - pysam will error out
#     helpfully if not.
#   - Reads the input table with pandas rather than numpy.genfromtxt.
#   - Collects output rows in a list and writes once, instead of repeated
#     string concatenation.

import argparse
import pandas as pd
import pysam

# chromosomes accepted for output
CHROMOSOME_SET = {str(i) for i in range(1, 23)} | {"X", "Y"}


def get_mate_info(bam, chromosome, position):
    """
    Look up mapq/strand for a read at a given chromosome:position,
    skipping secondary/supplementary alignments (equivalent to `-F 2304`).
    Returns (mapq, strand) as strings, empty if not found.
    """
    # pysam fetch is 0-based, half-open; samtools view CLI region is 1-based inclusive.
    # to match "chr:pos-pos" from the original command, use pos-1 .. pos
    for aln in bam.fetch(chromosome, position - 1, position):
        if aln.flag & 2304:  # secondary or supplementary
            continue
        mapq = str(aln.mapping_quality)
        strand = "-" if aln.is_reverse else "+"
        return mapq, strand
    return "", ""


def main():
    parser = argparse.ArgumentParser(description="Add mate mapq/strand to telomere insertion table.")
    parser.add_argument("-i", dest="table_file", required=True, help="input discordant reads tsv")
    parser.add_argument("-b", dest="bam_file", required=True, help="indexed bam file")
    parser.add_argument("-o", dest="outfile_path", required=True, help="output tsv")
    args = parser.parse_args()

    # read input table (no header assumptions beyond skipping the header row)
    table = pd.read_csv(args.table_file, sep="\t", header=0, dtype=str)

    out_rows = [["read_name", "mate_chr", "mate_position", "mate_mapq", "mate_strand"]]

    with pysam.AlignmentFile(args.bam_file, "rb") as bam:
        for row in table.itertuples(index=False):
            read_name = row[0]
            chromosome = row[1]
            position = int(row[2])

            # skip mates mapped to decoy sequences
            if chromosome not in CHROMOSOME_SET:
                continue

            mapq, strand = "", ""
            for aln in bam.fetch(chromosome, position - 1, position):
                if aln.flag & 2304:  # skip secondary/supplementary
                    continue
                if aln.query_name != read_name:
                    continue
                mapq = str(aln.mapping_quality)
                strand = "-" if aln.is_reverse else "+"
                break

            out_rows.append([read_name, chromosome, str(position), mapq, strand])

    with open(args.outfile_path, "w") as outfile:
        outfile.write("\n".join("\t".join(r) for r in out_rows))


if __name__ == "__main__":
    main()