#!/usr/bin/env python3

import argparse
from collections import defaultdict

import pandas as pd
import pysam

from pipeline.tables import read_tsv, write_tsv


EMPTY_VALUES = {"", "NA", "NaN", "nan", "None", None}
SITE_WINDOW = 100
DISTANCE_SEQ_LENGTH = 12
DUPLICATE_FLAG = 1024


def parse_int(value):
    if value in EMPTY_VALUES:
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def is_true(value):
    return str(value).strip().lower() in {"true", "t", "1"}


def reverse_seq(seq):
    return seq[::-1]


def split_clipped(clipped_sequence):
    if not clipped_sequence:
        return []
    return [part.strip() for part in str(clipped_sequence).split(",") if part.strip()]


def clip_column_for_direction(expected_pos_fusion):
    if expected_pos_fusion == "downstream":
        return "end"
    if expected_pos_fusion == "upstream":
        return "start"
    return None


def junction_oriented_clip(row):
    """First base = the base immediately at the breakpoint, for either direction."""
    direction = row.get("expected_pos_fusion")
    if direction not in {"downstream", "upstream"}:
        return ""
    parts = split_clipped(row.get("clipped_sequence"))
    if not parts:
        return ""
    if direction == "downstream":
        return parts[-1]
    return reverse_seq(parts[0])


def reads_clipped_at_site(
    rows, insertion_site, require_direction=None, telomeric_only=False
):
    """Reads whose own directional clip coordinate lands exactly at insertion_site."""
    matches = {}
    for row in rows:
        direction = row.get("expected_pos_fusion")
        clip_col = clip_column_for_direction(direction)
        if clip_col is None:
            continue
        if require_direction is not None and direction != require_direction:
            continue
        if parse_int(row.get(clip_col)) != insertion_site:
            continue
        if telomeric_only and not is_true(row.get("part_telomere")):
            continue
        matches[row.get("read_name", "")] = row
    return matches


def hamming_distance(seq_a, seq_b, length):
    return sum(1 for i in range(length) if seq_a[i] != seq_b[i])


def min_hamming_distance(rows_a, rows_b):
    seqs_a = [
        s[:DISTANCE_SEQ_LENGTH]
        for s in (junction_oriented_clip(r) for r in rows_a.values())
        if len(s) >= DISTANCE_SEQ_LENGTH
    ]
    seqs_b = [
        s[:DISTANCE_SEQ_LENGTH]
        for s in (junction_oriented_clip(r) for r in rows_b.values())
        if len(s) >= DISTANCE_SEQ_LENGTH
    ]
    if not seqs_a or not seqs_b:
        return None
    best = None
    for seq_a in seqs_a:
        for seq_b in seqs_b:
            dist = hamming_distance(seq_a, seq_b, DISTANCE_SEQ_LENGTH)
            if best is None or dist < best:
                best = dist
    return best


def chr_aliases(chrom):
    value = str(chrom).strip()
    if value.startswith("chr"):
        return [value, value[3:]]
    return [value, f"chr{value}"]


def resolve_contig(chrom, contigs):
    for alias in chr_aliases(chrom):
        if alias in contigs:
            return alias
    return None


def count_reads_covering_site(bam, chrom, site, window):
    if site is None:
        return None
    contig = resolve_contig(chrom, set(bam.references))
    if contig is None:
        return None
    contig_len = bam.get_reference_length(contig)
    start0 = max(0, site - window)
    end0 = min(contig_len, site + window)
    count = 0
    for read in bam.fetch(contig, start0, end0):
        if read.flag & DUPLICATE_FLAG:
            continue
        if read.reference_start is None or read.reference_end is None:
            continue
        if read.reference_start <= site < read.reference_end:
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_regions_file")
    parser.add_argument("tumor_clipped_reads_file")
    parser.add_argument("control_clipped_reads_file")
    parser.add_argument("tumor_bam")
    parser.add_argument("control_bam")
    parser.add_argument("outfile")
    parser.add_argument("--site-window", type=int, default=SITE_WINDOW)
    args = parser.parse_args()

    output_df, output_fields = assess_confidence(
        args.candidate_regions_file,
        args.tumor_clipped_reads_file,
        args.control_clipped_reads_file,
        args.tumor_bam,
        args.control_bam,
        args.site_window,
    )
    write_tsv(output_df, args.outfile, output_fields)


def assess_confidence(
    candidate_regions_file: str,
    tumor_clipped_reads_file: str,
    control_clipped_reads_file: str,
    tumor_bam_path: str,
    control_bam_path: str,
    site_window: int = SITE_WINDOW,
):
    candidate_df = read_tsv(candidate_regions_file)
    candidate_rows = candidate_df.to_dict("records")
    candidate_fields = list(candidate_df.columns)

    tumor_clipped_by_window = defaultdict(list)
    for row in read_tsv(tumor_clipped_reads_file).to_dict("records"):
        tumor_clipped_by_window[row.get("window", "")].append(row)

    have_control_clips = (
        bool(control_clipped_reads_file) and control_clipped_reads_file != "NULL"
    )
    control_clipped_by_window = defaultdict(list)
    if have_control_clips:
        for row in read_tsv(control_clipped_reads_file).to_dict("records"):
            control_clipped_by_window[row.get("window", "")].append(row)

    have_control_bam = bool(control_bam_path) and control_bam_path != "NULL"

    new_fields = [
        "all_reads_at_site",
        "clipped_reads_at_site",
        "telo_clipped_reads_at_site",
        "tumor_noise_ratio",
        "control_all_reads_at_site",
        "control_telo_clipped_at_insertion_site",
        "control_min_seq_distance_to_tumor",
    ]
    output_fields = candidate_fields + [
        f for f in new_fields if f not in candidate_fields
    ]

    tumor_bam = pysam.AlignmentFile(tumor_bam_path, "rb")
    control_bam = (
        pysam.AlignmentFile(control_bam_path, "rb") if have_control_bam else None
    )

    try:
        for region in candidate_rows:
            for field in new_fields:
                region[field] = ""

            insertion_site = parse_int(region.get("insertion_site"))
            strand = region.get("strand", "")
            chrom = region.get("chrom", "")
            window = region.get("window", "")
            if insertion_site is None or strand not in {"+", "-"}:
                continue

            expected_dir = "downstream" if strand == "+" else "upstream"

            tumor_rows = tumor_clipped_by_window.get(window, [])
            tumor_any_dir_at_site = reads_clipped_at_site(tumor_rows, insertion_site)
            tumor_telo_any_dir_at_site = reads_clipped_at_site(
                tumor_rows, insertion_site, telomeric_only=True
            )
            region["clipped_reads_at_site"] = str(len(tumor_any_dir_at_site))
            region["telo_clipped_reads_at_site"] = str(len(tumor_telo_any_dir_at_site))

            all_reads_at_site = count_reads_covering_site(
                tumor_bam, chrom, insertion_site, site_window
            )
            region["all_reads_at_site"] = (
                "" if all_reads_at_site is None else str(all_reads_at_site)
            )
            if all_reads_at_site:
                noise_ratio = (
                    len(tumor_any_dir_at_site) - len(tumor_telo_any_dir_at_site)
                ) / all_reads_at_site
                region["tumor_noise_ratio"] = f"{noise_ratio:.4f}"

            if have_control_bam:
                control_all_reads_at_site = count_reads_covering_site(
                    control_bam, chrom, insertion_site, site_window
                )
                region["control_all_reads_at_site"] = (
                    ""
                    if control_all_reads_at_site is None
                    else str(control_all_reads_at_site)
                )

            if have_control_clips:
                control_rows = control_clipped_by_window.get(window, [])
                control_matching = reads_clipped_at_site(
                    control_rows,
                    insertion_site,
                    require_direction=expected_dir,
                    telomeric_only=True,
                )
                region["control_telo_clipped_at_insertion_site"] = str(
                    len(control_matching)
                )

                tumor_matching = reads_clipped_at_site(
                    tumor_rows,
                    insertion_site,
                    require_direction=expected_dir,
                    telomeric_only=True,
                )
                distance = min_hamming_distance(tumor_matching, control_matching)
                region["control_min_seq_distance_to_tumor"] = (
                    "" if distance is None else str(distance)
                )
    finally:
        tumor_bam.close()
        if control_bam is not None:
            control_bam.close()

    output_df = pd.DataFrame(candidate_rows, columns=output_fields)
    return output_df, output_fields


if __name__ == "__main__":
    main()
