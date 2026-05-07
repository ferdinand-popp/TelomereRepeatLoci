#!/usr/bin/env python3
# Author: Lina Sieverling (extended debug/fix version)

import sys
import os
import argparse
import pandas as pd
import pysam

telomere_insertion_table_file = None
alignment_bam_file = None
outfile_path = None


def _parse_args(argv):
    p = argparse.ArgumentParser(
        prog=os.path.basename(argv[0]),
        description="Add mate MAPQ/strand info for reads listed in a telomere insertion table.",
    )
    p.add_argument(
        "-i", "--input-table", required=True, dest="telomere_insertion_table_file"
    )
    p.add_argument("-b", "--bam", required=True, dest="alignment_bam_file")
    p.add_argument("-o", "--output", required=True, dest="outfile_path")
    return p.parse_args(argv[1:])


args = _parse_args(sys.argv)
telomere_insertion_table_file = args.telomere_insertion_table_file
alignment_bam_file = args.alignment_bam_file
outfile_path = args.outfile_path

chromosome_list = [str(i) for i in range(1, 22 + 1)] + ["X", "Y"]
WINDOW_BP = 5000
FALLBACK_CONTIG_SCAN = True


def norm_read_name(x: str) -> str:
    """Normalize read names so table and BAM names can match robustly."""
    x = str(x).strip()
    if x.endswith("/1") or x.endswith("/2"):
        x = x[:-2]
    return x


def chr_aliases(chrom: str):
    """Return possible chromosome aliases, e.g. 1 <-> chr1."""
    c = str(chrom).strip()
    if c.startswith("chr"):
        base = c[3:]
        return [c, base]
    return [c, f"chr{c}"]


def resolve_contig(chrom: str, contigs_set):
    """Resolve table chromosome name to a contig present in BAM."""
    for a in chr_aliases(chrom):
        if a in contigs_set:
            return a
    return None


# Replace numpy.genfromtxt() loading with pandas
try:
    df = pd.read_csv(
        telomere_insertion_table_file,
        sep="\t",
        header=0,
        dtype=str,
        keep_default_na=False,
        na_filter=False,
    )
except Exception as e:
    print(f"[ERROR] Could not read input table: {e}", file=sys.stderr)
    sys.exit(1)

# ---- pandas best-practice: validate + normalize once, then iterate efficiently ----
required = {"read_name", "mate_chr", "mate_position"}
missing = required - set(df.columns)
if missing:
    print(
        f"[ERROR] Input table missing required columns: {', '.join(sorted(missing))}",
        file=sys.stderr,
    )
    sys.exit(1)

# normalize/derive columns used by the BAM lookup
df = df.copy()
df["read_name"] = df["read_name"].astype(str).str.strip()
df["mate_chr"] = df["mate_chr"].astype(str).str.strip()
df["mate_position"] = df["mate_position"].astype(str).str.strip()

df["_read_name_norm"] = df["read_name"].map(norm_read_name)
df["_chrom_base"] = df["mate_chr"].str.replace(r"^chr", "", regex=True)

# parse position once; keep NaN for bad values
df["_pos0"] = pd.to_numeric(df["mate_position"], errors="coerce").astype("Int64")

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
        contigs = set(bam.references)

        def emit(
            outfile, read_name_raw, mate_chr, mate_pos, mate_mapq, mate_strand, status
        ):
            outfile.write(
                "\t".join(
                    [
                        str(read_name_raw),
                        str(mate_chr or ""),
                        str(mate_pos or ""),
                        str(mate_mapq or ""),
                        str(mate_strand or ""),
                        str(status),
                    ]
                )
                + "\n"
            )

        with open(outfile_path, "w", encoding="utf-8") as outfile:
            outfile.write(
                "\t".join(
                    [
                        "read_name",
                        "mate_chr",
                        "mate_position",
                        "mate_mapq",
                        "mate_strand",
                        "status",
                    ]
                )
                + "\n"
            )

            if df.empty:
                pass
            else:
                for r in df.to_dict("records"):
                    read_name_raw = r.get("read_name", "")
                    read_name_norm = r.get("_read_name_norm", "")
                    chromosome_in = r.get("mate_chr", "")
                    chrom_base = r.get("_chrom_base", "")
                    pos0 = r.get("_pos0")

                    mate_chr = ""
                    mate_pos = ""
                    mate_mapq = ""
                    mate_strand = ""
                    status = "read_not_found"

                    if chrom_base not in chromosome_list:
                        emit(
                            outfile,
                            read_name_raw,
                            mate_chr,
                            mate_pos,
                            mate_mapq,
                            mate_strand,
                            "chr_not_allowed",
                        )
                        continue

                    chrom_resolved = resolve_contig(chromosome_in, contigs)
                    if chrom_resolved is None:
                        emit(
                            outfile,
                            read_name_raw,
                            mate_chr,
                            mate_pos,
                            mate_mapq,
                            mate_strand,
                            "chr_not_in_bam",
                        )
                        continue

                    if pd.isna(pos0):
                        emit(
                            outfile,
                            read_name_raw,
                            mate_chr,
                            mate_pos,
                            mate_mapq,
                            mate_strand,
                            "bad_pos",
                        )
                        continue

                    contig_len = bam.get_reference_length(chrom_resolved)
                    start0 = max(0, int(pos0) - WINDOW_BP)
                    end0 = min(contig_len, int(pos0) + WINDOW_BP)

                    found = False
                    try:
                        # 1) windowed fetch
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
                                    mate_chr = bam.get_reference_name(
                                        aln.next_reference_id
                                    )
                                except Exception:
                                    mate_chr = ""
                                # mate_position now emitted as 0-based to match input convention
                                mate_pos = (
                                    str(aln.next_reference_start)
                                    if aln.next_reference_start is not None
                                    and aln.next_reference_start >= 0
                                    else ""
                                )
                                status = "ok_window"

                            found = True
                            break

                        # 2) optional whole-contig fallback
                        if (not found) and FALLBACK_CONTIG_SCAN:
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
                                        mate_chr = bam.get_reference_name(
                                            aln.next_reference_id
                                        )
                                    except Exception:
                                        mate_chr = ""
                                    # mate_position now emitted as 0-based to match input convention
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

                    emit(
                        outfile,
                        read_name_raw,
                        mate_chr,
                        mate_pos,
                        mate_mapq,
                        mate_strand,
                        status,
                    )

except OSError as e:
    print(f"[ERROR] Could not open/read BAM during processing: {e}", file=sys.stderr)
    sys.exit(1)
