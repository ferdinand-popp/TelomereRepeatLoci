#!/usr/bin/env python3
# Author: Lina Sieverling

import sys
import getopt
import os
import numpy
import pysam

telomere_insertion_table_file = None
alignment_bam_file = None
outfile_path = None

myopts, args = getopt.getopt(sys.argv[1:], "i:b:o:")
for opt, arg in myopts:
    if opt == "-i":
        telomere_insertion_table_file = arg
    elif opt == "-b":
        alignment_bam_file = arg
    elif opt == "-o":
        outfile_path = arg

if not (telomere_insertion_table_file and alignment_bam_file and outfile_path):
    print(
        f"Usage: {sys.argv[0]} -i input_table -b bam_file -o output_file",
        file=sys.stderr,
    )
    sys.exit(2)

chromosome_list = [str(i) for i in range(1, 22 + 1)] + ["X", "Y"]
WINDOW_BP = 5000
FALLBACK_CONTIG_SCAN = True

telomere_insertion_table = numpy.genfromtxt(
    telomere_insertion_table_file,
    skip_header=1,
    delimiter="\t",
    dtype=str,
    encoding="utf-8",
    comments=None,
)

if telomere_insertion_table.size == 0:
    rows = numpy.empty((0, 0), dtype=str)
elif telomere_insertion_table.ndim == 1:
    rows = numpy.atleast_2d(telomere_insertion_table)
else:
    rows = telomere_insertion_table

# --- ensure BAM can be opened ---
try:
    with pysam.AlignmentFile(alignment_bam_file, "rb") as bam_check:
        if bam_check.closed:
            raise OSError(f"Failed to open BAM file: {alignment_bam_file}")
        has_idx = bam_check.has_index()
except OSError as e:
    print(f"[ERROR] Could not open BAM or check index: {e}", file=sys.stderr)
    sys.exit(1)

# If no index: sort first, then index, then use sorted BAM
bam_for_processing = alignment_bam_file
if not has_idx:
    base, ext = os.path.splitext(alignment_bam_file)
    sorted_bam = f"{base}.sorted{ext if ext else '.bam'}"

    print(f"[INFO] BAM index missing for {alignment_bam_file}", file=sys.stderr)
    print(f"[INFO] Sorting BAM -> {sorted_bam}", file=sys.stderr)
    try:
        pysam.sort("-o", sorted_bam, alignment_bam_file)
    except Exception as e:
        print(f"[ERROR] Failed to sort BAM: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Creating index for sorted BAM -> {sorted_bam}", file=sys.stderr)
    try:
        pysam.index(sorted_bam)
    except Exception as e:
        print(
            f"[ERROR] Failed to create BAM index for sorted BAM: {e}", file=sys.stderr
        )
        sys.exit(1)

    bam_for_processing = sorted_bam

# main processing
try:
    with pysam.AlignmentFile(bam_for_processing, "rb") as bam:
        if bam.closed:
            raise OSError(f"Failed to open BAM file: {bam_for_processing}")

        if not bam.has_index():
            raise OSError(
                f"BAM index still unavailable after sort/index attempt: {bam_for_processing}"
            )

        contigs = set(bam.references)

        output_lines = [
            "\t".join(
                ["read_name", "mate_chr", "mate_position", "mate_mapq", "mate_strand"]
            )
        ]

        for read in rows:
            if len(read) < 3:
                continue

            read_name = str(read[0]).strip()
            chromosome = str(read[1]).strip()
            position = str(read[2]).strip()

            if chromosome not in chromosome_list:
                continue

            mapq = ""
            strand = ""

            if chromosome not in contigs:
                output_lines.append(
                    "\t".join([read_name, chromosome, position, mapq, strand])
                )
                continue

            try:
                pos1 = int(position)  # expected 1-based
            except ValueError:
                output_lines.append(
                    "\t".join([read_name, chromosome, position, mapq, strand])
                )
                continue

            contig_len = bam.get_reference_length(chromosome)
            start0 = max(0, pos1 - 1 - WINDOW_BP)
            end0 = min(contig_len, pos1 + WINDOW_BP)

            found = False
            try:
                for aln in bam.fetch(chromosome, start0, end0):
                    if aln.is_secondary or aln.is_supplementary:
                        continue
                    if aln.query_name != read_name:
                        continue
                    mapq = str(aln.mapping_quality)
                    strand = "-" if aln.is_reverse else "+"
                    found = True
                    break

                if (not found) and FALLBACK_CONTIG_SCAN:
                    for aln in bam.fetch(chromosome):
                        if aln.is_secondary or aln.is_supplementary:
                            continue
                        if aln.query_name != read_name:
                            continue
                        mapq = str(aln.mapping_quality)
                        strand = "-" if aln.is_reverse else "+"
                        break

            except (ValueError, OSError):
                mapq = ""
                strand = ""

            output_lines.append(
                "\t".join([read_name, chromosome, position, mapq, strand])
            )

except OSError as e:
    print(f"[ERROR] Could not open/read BAM during processing: {e}", file=sys.stderr)
    sys.exit(1)

with open(outfile_path, "w", encoding="utf-8") as outfile:
    outfile.write("\n".join(output_lines))
