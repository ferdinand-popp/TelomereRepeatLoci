# Author: originally written by Philip Ginsbach (supervized by Ivo Buchhalter), adjusted by Lina Sieverling

#!/usr/bin/python

# This script generates a png and pdf file for each entry in a bed file. The file
# displays the reads in a bam file around the bed entry (+-WINDOW_SIZE).

# only plots non-duplicate reads

# Usage: set +u; source activate telomereEnv; set -u;
#        python ~/Code/telomere_insertion_analysis/snakemake_telomere_insertions/src/visualize_telomere_insertions.py \
#        --control <alignment bam file of control sample> \
#        --tumor  <alignment bam file of tumor sample>  \
#        --ref <fasta reference genome>  \
#        --bed <BED file of regions to plot>  \
#        --samtoolsbin <samtools version to use> \
#        --colored_reads_tumor <*/tables/*_tumor_discordant_reads_filtered_with_mapq.tsv> \
#        --colored_reads_control <*/tables/*_control_discordant_reads_filtered_with_mapq.tsv> \
#        --clipped_reads_tumor <*/clipped_reads/*_tumor_clipped_reads.tsv> \
#        --clipped_reads_control <*/clipped_reads/*_control_clipped_reads.tsv> \
#        --prefix <output directory> \
#        --outfile <output log file>


import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plot
from matplotlib import gridspec
import pandas as pd
import pysam


basepair_colors_axis = {
    "A": "#009600",
    "C": "#3030fe",
    "G": "#d17105",
    "T": "#ff0000",
    "N": "#00ffff",
}
basepair_colors = {
    "A": "#4DE34D",
    "C": "#7D7DFF",
    "G": "#FFBE52",
    "T": "#FF4D4D",
    "N": "#00ffff",
}


def parse_args(argv=None):
    argument_parser = argparse.ArgumentParser(
        description="This script generates a pdf file for each entry in a bed file."
    )
    argument_parser.add_argument(
        "--control",
        metavar="FILE",
        type=str,
        default=None,
        help="input bam file of the control",
    )
    argument_parser.add_argument(
        "--tumor",
        metavar="FILE",
        type=str,
        required=True,
        help="input bam file of the tumor",
    )
    argument_parser.add_argument(
        "--ref",
        metavar="FILE",
        type=str,
        required=True,
        help="input reference genome file (fastq format)",
    )
    argument_parser.add_argument(
        "--bed", metavar="FILE", type=str, default=None, help="input bed file"
    )
    argument_parser.add_argument(
        "--annotations",
        metavar="FILE",
        type=str,
        default=None,
        help="annotation track indexed with tabix",
    )
    argument_parser.add_argument(
        "--prefix",
        metavar="PREFIX",
        type=str,
        default="./",
        help="target directory and file name prefix for generated output files",
    )
    argument_parser.add_argument(
        "--samtoolsbin",
        metavar="N",
        type=str,
        default="samtools",
        help="(unused) kept for compatibility",
    )
    argument_parser.add_argument(
        "--tabixbin",
        metavar="N",
        type=str,
        default="tabix",
        help="(unused) kept for compatibility",
    )
    argument_parser.add_argument(
        "--colored_reads_tumor",
        type=str,
        default=None,
        help="a table with reads to color",
    )
    argument_parser.add_argument(
        "--colored_reads_control",
        type=str,
        default=None,
        help="a table with reads to color",
    )
    argument_parser.add_argument(
        "--clipped_reads_tumor",
        type=str,
        default=None,
        help="a table with all clipped read sequences in tumor sample",
    )
    argument_parser.add_argument(
        "--clipped_reads_control",
        type=str,
        default=None,
        help="a table with all clipped read sequences in control sample",
    )
    argument_parser.add_argument(
        "--outfile",
        type=str,
        default=None,
        help="path to dummy output file (needed for snakemake)",
    )
    return argument_parser.parse_args(argv)


class ReferenceBuffer(object):
    def __init__(self, filename, chromosome):
        self.fasta = pysam.FastaFile(filename)
        self.chromosome = chromosome
        self.offset = 0
        self.sequence = ""

    def __getitem__(self, pos):
        if len(self.sequence) > pos - self.offset >= 0:
            return self.sequence[pos - self.offset]

        region_start = max(0, pos - 1000)
        region_end = pos + 1001
        seq = self.fasta.fetch(self.chromosome, region_start, region_end).upper()
        self.offset = region_start
        self.sequence = seq
        return self.sequence[pos - self.offset]


_TABIX_HANDLE = None


def get_tabix_handle(path):
    global _TABIX_HANDLE
    if _TABIX_HANDLE is None:
        _TABIX_HANDLE = pysam.TabixFile(path)
    return _TABIX_HANDLE


def get_annotations(region, annotations_path):
    if not annotations_path:
        return None
    tabix = get_tabix_handle(annotations_path)
    chrom, span = region.split(":", 1)
    start_s, end_s = span.split("-", 1)
    try:
        start = int(start_s)
        end = int(end_s)
    except ValueError:
        return None
    candidates = [chrom, chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"]
    for chrom_name in candidates:
        try:
            rows = list(tabix.fetch(chrom_name, max(0, start - 1), end))
        except ValueError:
            rows = []
        if rows:
            annotations = []
            for line in rows:
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                name = parts[3] if len(parts) > 3 else ""
                try:
                    ann_start = int(parts[1])
                    ann_end = int(parts[2])
                except ValueError:
                    continue
                annotations.append([name, ann_start, ann_end])
            return annotations
    return None


def parse_cigar(cigar, pos):
    cigar_struct = [["", None]]

    for char in cigar:
        if "0" <= char <= "9":
            if not isinstance(cigar_struct[-1][0], str):
                cigar_struct.append(["", None])

            cigar_struct[-1][0] = cigar_struct[-1][0] + char

        else:
            cigar_struct[-1][0] = int(cigar_struct[-1][0])
            cigar_struct[-1][1] = char

    if cigar_struct[-1][1] is None:
        del cigar_struct[-1]

    abs_pos = pos
    rel_pos = 0

    if cigar_struct[0][1] == "S":
        abs_pos -= cigar_struct[0][0]

    cigar_struct2 = []

    cigar_start = True

    for n, t in cigar_struct:
        if t in ["M", "S"] or (t == "H" and not cigar_start):
            for m in range(n):
                cigar_struct2.append((t, abs_pos + m, rel_pos + m))

            abs_pos += n
            rel_pos += n

        elif t == "D":
            for m in range(n):
                cigar_struct2.append((t, abs_pos + m, rel_pos))

            abs_pos += n

        elif t == "I":
            cigar_struct2.append((t, abs_pos, rel_pos + m))

            for m in range(1, n):
                cigar_struct2.append(("i", abs_pos, rel_pos + m))

            rel_pos += n

        elif t == "H" and cigar_start:
            for m in range(n):
                cigar_struct2.append((t, abs_pos - n + m, rel_pos + m))

        cigar_start = False

    return cigar_struct2


def plot_histogram(cigars, ax):
    points_deleted = []
    points_clipped = []
    points_unclipped = []

    for cigar in cigars:
        for c_type, abs_pos, rel_pos in cigar:
            if c_type in ["M"]:
                points_deleted.append(abs_pos)
                points_clipped.append(abs_pos)
                points_unclipped.append(abs_pos)
            elif c_type in "D":
                points_clipped.append(abs_pos)
                points_unclipped.append(abs_pos)
            elif c_type in ["S"]:
                points_unclipped.append(abs_pos)
            elif c_type in ["H"]:
                points_unclipped.append(abs_pos)

    points_deleted = sorted(points_deleted)
    points_clipped = sorted(points_clipped)
    points_unclipped = sorted(points_unclipped)

    min_unclipped = min(points_unclipped)
    max_unclipped = max(points_unclipped)

    histogram_deleted = [0] * (1 + max_unclipped - min_unclipped)
    histogram_clipped = [0] * (1 + max_unclipped - min_unclipped)
    histogram_unclipped = [0] * (1 + max_unclipped - min_unclipped)

    for p in points_deleted:
        histogram_deleted[p - min_unclipped] += 1

    for p in points_clipped:
        histogram_clipped[p - min_unclipped] += 1

    for p in points_unclipped:
        histogram_unclipped[p - min_unclipped] += 1

    original_x = range(min_unclipped, max_unclipped + 1)
    subdivided_x = [a + b for a in original_x for b in [-0.33, 0.33]][1:-2]

    smooth_unclipped = [x for x in histogram_unclipped for y in [0, 1]][1:-2]
    smooth_clipped = [x for x in histogram_clipped for y in [0, 1]][1:-2]
    smooth_deleted = [x for x in histogram_deleted for y in [0, 1]][1:-2]

    ax.fill_between(subdivided_x, 0, smooth_unclipped, color="lightblue")
    ax.fill_between(subdivided_x, 0, smooth_clipped, color="darkblue")
    ax.fill_between(subdivided_x, 0, smooth_deleted, color="blue")
    ax.set_ylim(ymin=0)


def plot_cigars(
    cigars,
    sequences,
    reverses,
    read_names,
    read_flags,
    colored_reads,
    ax,
    reference_function,
):
    def _base_color(sequence, rel_pos):
        if 0 <= rel_pos < len(sequence):
            base = sequence[rel_pos]
        else:
            base = "N"
        return basepair_colors.get(base, basepair_colors["N"])

    patches = []

    read_colors = []

    right_limits = [0] * len(cigars)

    read_strands = ["-" if flag & 0x10 else "+" for flag in read_flags]

    for cigar, sequence, reverse, read_name, read_strand in zip(
        cigars, sequences, reverses, read_names, read_strands
    ):
        for j in range(len(cigars)):
            if right_limits[j] <= cigar[0][1]:
                line = j
                break

        clipped_cigar = [(t, a, r) for (t, a, r) in cigar if t != "S"]
        right_limits[line] = cigar[-1][1] + 10

        # get color of read
        if read_name + "_" + read_strand in colored_reads:
            read_color = "#999999"  # "#99ccff"
        elif any(
            read_name in r for r in colored_reads
        ):  # read_name partial match colored_read:
            read_color = "#bfbfbf"  # "#B3E6FF"
        else:
            read_color = "#f2f2f2"  # "#e6e6e6"

        read_colors.append(read_color)

        # check if sequence contains telomere repeat
        c_types = [x[0] for x in cigar]
        if "S" in c_types or "H" in c_types:
            if "TTAGGG" in sequence or "CCCTAA" in sequence:
                alpha_clipped = 1
            else:
                alpha_clipped = 0.3

        ax.barh(
            -line + 0.5,
            clipped_cigar[-1][1] - clipped_cigar[0][1] + 1,
            height=1,
            left=clipped_cigar[0][1] - 0.5,
            color=read_color,
            linewidth=0,
        )

        for (
            c_type,
            abs_pos,
            rel_pos,
        ) in cigar:  # this loop is slow for clipped reads => can create problems when there are too many clipped reads in window
            if c_type == "D":
                ax.barh(
                    -line + 0.5,
                    1,
                    height=1,
                    left=abs_pos - 0.5,
                    color="#505050",
                    linewidth=0,
                )
            # elif c_type == 'M' and sequence[rel_pos] != reference_function[abs_pos]:
            #  ax.barh( -line, 1, height=1, left=abs_pos-0.5, color=basepair_colors[sequence[rel_pos]], linewidth=0, alpha=0.5 )      # color of point mutations
            elif c_type == "S":
                ax.barh(
                    -line + 0.5,
                    1,
                    height=1,
                    left=abs_pos - 0.5,
                    color=_base_color(sequence, rel_pos),
                    linewidth=0,
                    alpha=alpha_clipped,
                )  # color of softclipped bases
            elif c_type == "H":
                ax.barh(
                    -line + 0.5,
                    1,
                    height=1,
                    left=abs_pos - 0.5,
                    color=_base_color(sequence, rel_pos),
                    linewidth=0,
                    alpha=alpha_clipped,
                )  # color of hardclipped bases

        for c_type, abs_pos, rel_pos in cigar:
            if c_type == "I":
                ax.barh(
                    -line + 0.5,
                    0.6,
                    height=1,
                    left=abs_pos - 0.8,
                    fill=False,
                    linewidth=1,
                )

        # ax.barh( -line, cigar[-1][1]-cigar[0][1]+1, height=1, left=cigar[0][1]-0.5, fill=False, linewidth=0 ) # linewidth changed to 0

        # add arrow for read direction
        if reverse:
            patches.append(
                matplotlib.patches.Polygon(
                    [
                        (cigar[0][1] - 0.5, -line),
                        (cigar[0][1] - 0.5, 1 - line),
                        (cigar[0][1] - 2, 0.5 - line),
                    ]
                )
            )
            # patches.append( matplotlib.patches.Polygon( [(cigar[0][1]-0.5,-0.5-line),(cigar[0][1]-0.5,0.5-line),(cigar[0][1]-2,-line)] ) )
            # patches.append( matplotlib.patches.Polygon( [(cigar[-1][1]+0.5,1-line),(cigar[-1][1]+1,1-line),(cigar[-1][1]+0.5,0.5-line)] ) )
            # patches.append( matplotlib.patches.Polygon( [(cigar[-1][1]+0.5,-line),(cigar[-1][1]+1,-line),(cigar[-1][1]+0.5,0.5-line)] ) )
        else:
            patches.append(
                matplotlib.patches.Polygon(
                    [
                        (cigar[-1][1] + 0.5, -line),
                        (cigar[-1][1] + 0.5, 1 - line),
                        (cigar[-1][1] + 2, 0.5 - line),
                    ]
                )
            )
            # patches.append( matplotlib.patches.Polygon( [(cigar[-1][1]+0.5,-0.5-line),(cigar[-1][1]+0.5,0.5-line),(cigar[-1][1]+2,0-line)] ) )
            # patches.append( matplotlib.patches.Polygon( [(cigar[0][1]-0.5,1-line),(cigar[0][1]-1,1-line),(cigar[0][1]-0.5,0.5-line)] ) )
            # patches.append( matplotlib.patches.Polygon( [(cigar[0][1]-0.5,-line),(cigar[0][1]-1,-line),(cigar[0][1]-0.5,0.5-line)] ) )

    # collection = matplotlib.collections.PatchCollection( patches, linewidths=0.5, edgecolors="black", facecolors="yellow" )
    collection = matplotlib.collections.PatchCollection(
        patches, linewidths=0, edgecolors=None, facecolors=read_colors
    )
    ax.add_collection(collection)

    ax.set_ylim(bottom=(1 - len([limit for limit in right_limits if limit > 0])), top=1)
    ax.set_yticks([])


def plot_region(
    parsed_arguments,
    clipped_reads_tumor_dict,
    clipped_reads_control_dict,
    region_chrom,
    region_center,
    region_left,
    region_right,
    plot_title,
):
    region_string = "%s:%i-%i" % (region_chrom, region_left, region_right)

    samtools_reads1 = []
    if parsed_arguments.control is not None:
        samtools_reads1 = read_bam_region(parsed_arguments.control, region_string)

    samtools_reads2 = read_bam_region(parsed_arguments.tumor, region_string)

    annotations = get_annotations(region_string, parsed_arguments.annotations)

    if len(samtools_reads1) < 3000 and len(samtools_reads2) < 3000:
        # print( "region %s annotations: " % region_string, annotations )

        reference_buffer = ReferenceBuffer(parsed_arguments.ref, region_chrom)

        # set figure size
        # fig = plot.figure(figsize=(19.2, 10.8))
        plot.figure(figsize=(35, 20))

        if annotations:
            grid = gridspec.GridSpec(
                9,
                3,
                height_ratios=[1, 5, 15, 2, 5, 15, 2, 1, 1],
                hspace=0,
                width_ratios=[1, 42, 1],
                wspace=0,
                left=0,
                right=1,
                bottom=0,
                top=1,
            )
            ax = [plot.subplot(grid[i, 1]) for i in [1, 2, 4, 5, 7]]

        else:
            grid = gridspec.GridSpec(
                7,
                3,
                height_ratios=[1, 3, 15, 2, 3, 30, 1],
                hspace=0,  #
                width_ratios=[1, 42, 1],
                wspace=0,
                left=0,
                right=1,
                bottom=0,
                top=1,
            )
            ax = [plot.subplot(grid[i, 1]) for i in [1, 2, 4, 5]]

        # plot control
        if parsed_arguments.control:
            if parsed_arguments.colored_reads_control:
                colored_reads_control = getColoredReads(
                    parsed_arguments.colored_reads_control, region_chrom
                )
            else:
                colored_reads_control = []

            plot_histogram(
                [
                    parse_cigar(read[5], int(read[3]))
                    for read in samtools_reads1
                    if read[5] != "*"
                ],
                ax[0],
            )

            plot_cigars(
                [
                    parse_cigar(read[5], int(read[3]))
                    for read in samtools_reads1
                    if read[5] != "*"
                ],
                [
                    get_sequence(
                        read, parsed_arguments.control, clipped_reads_control_dict
                    )
                    for read in samtools_reads1
                    if read[5] != "*"
                ],
                [
                    bool(int(read[1]) & 0x10)
                    for read in samtools_reads1
                    if read[5] != "*"
                ],
                [read[0] for read in samtools_reads1 if read[5] != "*"],
                [int(read[1]) for read in samtools_reads1 if read[5] != "*"],
                colored_reads_control,
                ax[1],
                reference_buffer,
            )

        # plot tumor
        if parsed_arguments.colored_reads_tumor:
            colored_reads_tumor = getColoredReads(
                parsed_arguments.colored_reads_tumor, region_chrom
            )
        else:
            colored_reads_tumor = []

        plot_histogram(
            [
                parse_cigar(read[5], int(read[3]))
                for read in samtools_reads2
                if read[5] != "*"
            ],
            ax[2],
        )

        plot_cigars(
            [
                parse_cigar(read[5], int(read[3]))
                for read in samtools_reads2
                if read[5] != "*"
            ],
            [
                get_sequence(read, parsed_arguments.tumor, clipped_reads_tumor_dict)
                for read in samtools_reads2
                if read[5] != "*"
            ],
            [bool(int(read[1]) & 0x10) for read in samtools_reads2 if read[5] != "*"],
            [read[0] for read in samtools_reads2 if read[5] != "*"],
            [int(read[1]) for read in samtools_reads2 if read[5] != "*"],
            colored_reads_tumor,
            ax[3],
            reference_buffer,
        )

        for axis in ax:
            axis.set_xlim(xmin=region_left - 0.5, xmax=region_right + 0.5)

        ax[0].set_xticks([])
        ax[2].set_xticks([])

        ax[0].yaxis.set_tick_params(labelleft=True, labelright=True)
        ax[2].yaxis.set_tick_params(labelleft=True, labelright=True)

        visible_basepairs = [
            reference_buffer[i] for i in range(region_left, region_right + 1)
        ]

        ax[1].xaxis.set_tick_params(width=0)
        ax[1].set_xticks([i for i in range(region_left, region_right + 1)])
        ax[1].xaxis.set_ticklabels(visible_basepairs)

        for tick in ax[1].get_xticklabels():
            tick.set_color(basepair_colors_axis[tick._text])

        ax[3].xaxis.set_tick_params(width=0)
        ax[3].set_xticks([i for i in range(region_left, region_right + 1)])
        ax[3].xaxis.set_ticklabels(visible_basepairs)

        for tick in ax[3].get_xticklabels():
            tick.set_color(basepair_colors_axis[tick._text])

        for axis in ax[:4]:
            axis.axvline(region_center - 0.5, color="black", linewidth=0.5)
            axis.axvline(region_center + 0.5, color="black", linewidth=0.5)
            # for x in range( 10, parsed_arguments.window, 10 ):
            #  axis.axvline( region_center - x - 0.5, color="black", linewidth=0.25 )
            #  axis.axvline( region_center + x + 0.5, color="black", linewidth=0.25 )

        ax[0].set_title("%s - control" % plot_title)
        ax[2].set_title("%s - tumor" % plot_title)

        ax[0].ticklabel_format(style="plain", axis="x", useOffset=False)
        ax[2].ticklabel_format(style="plain", axis="x", useOffset=False)

        if annotations:
            ax[4].set_title(
                "annotations from " + parsed_arguments.annotations.split("/")[-1]
            )
            ax[4].set_xticks([])
            ax[4].set_yticks([])
            ax[4].axis("off")

            for ann in annotations:
                ann[1] = max(ann[1], region_left)
                ann[2] = min(ann[2], region_right)

                ax[4].barh(
                    0, ann[2] - ann[1] + 1, height=1, left=ann[1] - 0.5, color="#c8c8c8"
                )
                ax[4].text(
                    float(ann[1] + ann[2]) / 2.0, 0.5, ann[0], ha="center", va="center"
                )
    else:
        print("too many reads for plotting")
        print("reads tumor: " + str(len(samtools_reads2)))
        print("reads control: " + str(len(samtools_reads1)))


def parse_sa_tag(tag_value):
    if not tag_value:
        return None
    entry = tag_value.split(";")[0]
    parts = entry.split(",")
    if len(parts) < 4:
        return None
    chrom = parts[0]
    try:
        pos = int(parts[1])
    except ValueError:
        return None
    strand = parts[2]
    return chrom, pos, strand


def get_sequence(read, bamfile, clipped_read_dict=None):
    if "H" not in read[5]:
        return read[9]

    if clipped_read_dict:
        flag = int(read[1])
        if flag & 0x40:
            read_1_2 = "READ1"
        elif flag & 0x80:
            read_1_2 = "READ2"
        else:
            read_1_2 = ""
        return clipped_read_dict.get(read[0] + "_" + read_1_2, "")

    read_name = read[0]
    flag = int(read[1])
    strand_supp = "-" if flag & 0x10 else "+"
    is_read1 = bool(flag & 0x40)
    is_read2 = bool(flag & 0x80)

    tags = read[11:]
    sa_values = [t.split("SA:Z:", 1)[1] for t in tags if t.startswith("SA:Z:")]
    sa_tag = parse_sa_tag(sa_values[0]) if sa_values else None
    if not sa_tag:
        return read[9]

    chrom, pos, strand_primary = sa_tag
    sequence = read_primary_alignment(
        bamfile, read_name, chrom, pos, is_read1, is_read2
    )
    if strand_supp != strand_primary:
        sequence = getReverseComplement(sequence)
    return sequence


# get the reverse complement of a DNA Sequence
def getReverseComplement(sequence):
    sequence_temp = sequence.replace("A", "1")
    sequence_temp = sequence_temp.replace("C", "2")
    sequence_temp = sequence_temp.replace("G", "3")
    sequence_temp = sequence_temp.replace("T", "4")

    sequence_temp2 = sequence_temp.replace("1", "T")
    sequence_temp2 = sequence_temp.replace("2", "G")
    sequence_temp2 = sequence_temp.replace("3", "C")
    sequence_temp2 = sequence_temp.replace("4", "A")

    sequence_reverse_complement = sequence_temp2[::-1]  # reverses a string

    return sequence_reverse_complement


def getColoredReads(colored_read_file, chrom):
    df = pd.read_csv(
        colored_read_file,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        na_filter=False,
    )
    for col in ["read_name", "mate_chr", "mate_mapq", "mate_strand"]:
        if col not in df.columns:
            df[col] = ""
    df["mate_mapq"] = pd.to_numeric(df["mate_mapq"], errors="coerce").fillna(0)
    df = df[(df["mate_mapq"] > 30) & (df["mate_chr"] == chrom)]
    return (df["read_name"] + "_" + df["mate_strand"]).tolist()


def getClippedSequences(clipped_read_file):
    # makes a dictionary of clipped reads (read id (read name and read1/2) and sequence of entire read)
    try:
        df = pd.read_csv(
            clipped_read_file,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            na_filter=False,
        )
    except (FileNotFoundError, PermissionError):
        return None

    for col in ["read_name", "read_1_2", "sequence"]:
        if col not in df.columns:
            df[col] = ""
    read_ids = df["read_name"] + "_" + df["read_1_2"]
    return dict(zip(read_ids.tolist(), df["sequence"].tolist()))


def read_bam_region(bam_path, region_string):
    chrom, span = region_string.split(":", 1)
    start_s, end_s = span.split("-", 1)
    start = max(0, int(start_s) - 1)
    end = int(end_s)
    reads = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            if read.flag & 1024:
                continue
            tags = [
                f"{tag}:{value_type}:{value}"
                for tag, value, value_type in read.get_tags(with_value_type=True)
            ]
            reads.append(
                [
                    read.query_name,
                    str(read.flag),
                    read.reference_name or "*",
                    str(read.reference_start + 1),
                    str(read.mapping_quality),
                    read.cigarstring or "*",
                    "*",
                    "0",
                    "0",
                    read.query_sequence or "*",
                    "*",
                ]
                + tags
            )
    return reads


def read_primary_alignment(bamfile, read_name, chrom, pos, is_read1, is_read2):
    start = max(0, pos - 1)
    end = pos
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        candidates = [
            read
            for read in bam.fetch(chrom, start, end)
            if read.query_name == read_name
            and not read.is_supplementary
            and not read.is_secondary
        ]
    if is_read1:
        candidates = [read for read in candidates if read.is_read1]
    if is_read2:
        candidates = [read for read in candidates if read.is_read2]
    if not candidates:
        return ""
    if len(candidates) > 1:
        exact = [read for read in candidates if read.reference_start == start]
        if exact:
            candidates = exact
    return candidates[0].query_sequence or ""


def run(parsed_arguments):
    if not parsed_arguments.ref:
        raise FileNotFoundError(
            "Missing reference FASTA: provide --ref with a valid path."
        )
    if not os.path.exists(parsed_arguments.ref):
        raise FileNotFoundError(f"Reference FASTA not found: {parsed_arguments.ref}")
    if parsed_arguments.clipped_reads_tumor:
        clipped_reads_tumor_dict = getClippedSequences(
            parsed_arguments.clipped_reads_tumor
        )
    else:
        clipped_reads_tumor_dict = None

    if parsed_arguments.clipped_reads_control:
        clipped_reads_control_dict = getClippedSequences(
            parsed_arguments.clipped_reads_control
        )
    else:
        clipped_reads_control_dict = None

    if parsed_arguments.bed:
        for line in open(parsed_arguments.bed, "r"):
            if line[:1] != "#":
                region_chrom = line.split("\t")[0]
                region_left = int(line.split("\t")[1])
                region_right = int(line.split("\t")[2])
                region_center = int(line.split("\t")[3])
                pid = line.rstrip().split("\t")[4]

                if os.path.isfile(
                    "%s%s_%s_%i.png"
                    % (parsed_arguments.prefix, pid, region_chrom, region_center)
                ):
                    continue

                plot_title = "%s %s:%s" % (pid, region_chrom, region_center)

                plot_region(
                    parsed_arguments,
                    clipped_reads_tumor_dict,
                    clipped_reads_control_dict,
                    region_chrom,
                    region_center,
                    region_left,
                    region_right,
                    plot_title,
                )

                # plot.savefig( "%s%s_%s_%i.png" % ( parsed_arguments.prefix, pid, region_chrom, region_center ) )
                plot.savefig(
                    "%s%s_%s_%i.pdf"
                    % (parsed_arguments.prefix, pid, region_chrom, region_center)
                )
                plot.clf()
                plot.cla()
                plot.close()

    if parsed_arguments.outfile:
        with open(parsed_arguments.outfile, "w", encoding="utf-8") as outfile:
            outfile.write("complete")


def main(argv=None):
    parsed_arguments = parse_args(argv)
    run(parsed_arguments)


if __name__ == "__main__":
    main()
