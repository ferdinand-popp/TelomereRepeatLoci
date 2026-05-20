#!/usr/bin/env python3
# Author: Lina Sieverling (extended debug/fix version)

import argparse
import os
import sys
from typing import Iterable, List

import pandas as pd
import pysam

from telomererepeatloci.pipeline.tables import (
    DISCORDANT_READS_WITH_MAPQ_COLUMNS,
    read_tsv,
    write_tsv,
)


WINDOW_BP = 5000
FALLBACK_CONTIG_SCAN = True
CHROMOSOME_LIST = [str(i) for i in range(1, 22 + 1)] + ["X", "Y"]


def _parse_args(argv: List[str]):
    parser = argparse.ArgumentParser(
        prog=os.path.basename(argv[0]),
        description=(
            "Add mate MAPQ/strand info for reads listed in a telomere insertion table."
        ),
    )
    parser.add_argument(
        "-i", "--input-table", required=True, dest="telomere_insertion_table_file"
    )
    parser.add_argument("-b", "--bam", required=True, dest="alignment_bam_file")
    parser.add_argument("-o", "--output", required=True, dest="outfile_path")
    return parser.parse_args(argv[1:])


def norm_read_name(name: str) -> str:
    value = str(name).strip()
    if value.endswith("/1") or value.endswith("/2"):
        value = value[:-2]
    return value


def chr_aliases(chrom: str) -> List[str]:
    value = str(chrom).strip()
    if value.startswith("chr"):
        base = value[3:]
        return [value, base]
    return [value, f"chr{value}"]


def resolve_contig(chrom: str, contigs_set: Iterable[str]):
    for alias in chr_aliases(chrom):
        if alias in contigs_set:
            return alias
    return None


def load_table(path: str) -> pd.DataFrame:
    df = read_tsv(path)
    required = {"read_name", "mate_chr", "mate_position"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            "Input table missing required columns: " + ", ".join(sorted(missing))
        )
    df = df.copy()
    df["read_name"] = df["read_name"].astype(str).str.strip()
    df["mate_chr"] = df["mate_chr"].astype(str).str.strip()
    df["mate_position"] = df["mate_position"].astype(str).str.strip()
    df["_read_name_norm"] = df["read_name"].map(norm_read_name)
    df["_chrom_base"] = df["mate_chr"].str.replace(r"^chr", "", regex=True)
    df["_pos0"] = pd.to_numeric(df["mate_position"], errors="coerce").astype("Int64")
    return df


def ensure_bam_index(path: str) -> str:
    try:
        with pysam.AlignmentFile(path, "rb") as bam_check:
            if bam_check.closed:
                raise OSError(f"Failed to open BAM file: {path}")
            has_idx = bam_check.has_index()
    except OSError as exc:
        raise OSError(f"Could not open BAM or check index: {exc}") from exc

    if has_idx:
        return path

    base, ext = os.path.splitext(path)
    sorted_bam = f"{base}.sorted{ext if ext else '.bam'}"
    pysam.sort("-o", sorted_bam, path)
    pysam.index(sorted_bam)
    return sorted_bam


def add_mate_mapq_records(
    table_df: pd.DataFrame,
    bam_path: str,
    chromosome_list: Iterable[str] = CHROMOSOME_LIST,
    window_bp: int = WINDOW_BP,
    fallback_contig_scan: bool = FALLBACK_CONTIG_SCAN,
) -> pd.DataFrame:
    if table_df.empty:
        return pd.DataFrame(columns=DISCORDANT_READS_WITH_MAPQ_COLUMNS)

    bam_for_processing = ensure_bam_index(bam_path)
    output_rows = []

    with pysam.AlignmentFile(bam_for_processing, "rb") as bam:
        contigs = set(bam.references)

        for record in table_df.to_dict("records"):
            read_name_raw = record.get("read_name", "")
            read_name_norm = record.get("_read_name_norm", "")
            chromosome_in = record.get("mate_chr", "")
            chrom_base = record.get("_chrom_base", "")
            pos0 = record.get("_pos0")

            mate_chr = ""
            mate_pos = ""
            mate_mapq = ""
            mate_strand = ""
            status = "read_not_found"

            if chrom_base not in chromosome_list:
                status = "chr_not_allowed"
                output_rows.append(
                    {
                        "read_name": read_name_raw,
                        "mate_chr": mate_chr,
                        "mate_position": mate_pos,
                        "mate_mapq": mate_mapq,
                        "mate_strand": mate_strand,
                        "status": status,
                    }
                )
                continue

            chrom_resolved = resolve_contig(chromosome_in, contigs)
            if chrom_resolved is None:
                status = "chr_not_in_bam"
                output_rows.append(
                    {
                        "read_name": read_name_raw,
                        "mate_chr": mate_chr,
                        "mate_position": mate_pos,
                        "mate_mapq": mate_mapq,
                        "mate_strand": mate_strand,
                        "status": status,
                    }
                )
                continue

            if pd.isna(pos0):
                status = "bad_pos"
                output_rows.append(
                    {
                        "read_name": read_name_raw,
                        "mate_chr": mate_chr,
                        "mate_position": mate_pos,
                        "mate_mapq": mate_mapq,
                        "mate_strand": mate_strand,
                        "status": status,
                    }
                )
                continue

            contig_len = bam.get_reference_length(chrom_resolved)
            start0 = max(0, int(pos0) - window_bp)
            end0 = min(contig_len, int(pos0) + window_bp)

            found = False
            try:
                for aln in bam.fetch(chrom_resolved, start0, end0):
                    if aln.is_secondary or aln.is_supplementary:
                        continue
                    if norm_read_name(aln.query_name) != read_name_norm:
                        continue

                    mate_mapq = str(aln.mapping_quality)
                    mate_strand = "-" if aln.is_reverse else "+"

                    if aln.mate_is_unmapped:
                        mate_chr = ""
                        mate_pos = ""
                        status = "ok_window_mate_unmapped"
                    else:
                        try:
                            mate_chr = bam.get_reference_name(aln.next_reference_id)
                        except Exception:
                            mate_chr = ""
                        mate_pos = (
                            str(aln.next_reference_start)
                            if aln.next_reference_start is not None
                            and aln.next_reference_start >= 0
                            else ""
                        )
                        status = "ok_window"
                    found = True
                    break

                if (not found) and fallback_contig_scan:
                    for aln in bam.fetch(chrom_resolved):
                        if aln.is_secondary or aln.is_supplementary:
                            continue
                        if norm_read_name(aln.query_name) != read_name_norm:
                            continue

                        mate_mapq = str(aln.mapping_quality)
                        mate_strand = "-" if aln.is_reverse else "+"

                        if aln.mate_is_unmapped:
                            mate_chr = ""
                            mate_pos = ""
                            status = "ok_fallback_contig_mate_unmapped"
                        else:
                            try:
                                mate_chr = bam.get_reference_name(aln.next_reference_id)
                            except Exception:
                                mate_chr = ""
                            mate_pos = (
                                str(aln.next_reference_start)
                                if aln.next_reference_start is not None
                                and aln.next_reference_start >= 0
                                else ""
                            )
                            status = "ok_fallback_contig"
                        found = True
                        break

                if not found:
                    status = "read_not_found"

            except (ValueError, OSError):
                mate_chr = ""
                mate_pos = ""
                mate_mapq = ""
                mate_strand = ""
                status = "fetch_error"

            output_rows.append(
                {
                    "read_name": read_name_raw,
                    "mate_chr": mate_chr,
                    "mate_position": mate_pos,
                    "mate_mapq": mate_mapq,
                    "mate_strand": mate_strand,
                    "status": status,
                }
            )

    return pd.DataFrame(output_rows)


def add_mate_mapq_file(
    input_table: str, bam_path: str, output_path: str
) -> pd.DataFrame:
    df = load_table(input_table)
    output_df = add_mate_mapq_records(df, bam_path)
    write_tsv(output_df, output_path, DISCORDANT_READS_WITH_MAPQ_COLUMNS)
    return output_df


def main():
    try:
        args = _parse_args(sys.argv)
        add_mate_mapq_file(
            args.telomere_insertion_table_file,
            args.alignment_bam_file,
            args.outfile_path,
        )
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise


if __name__ == "__main__":
    main()
