from pathlib import Path

import pandas as pd

DISCORDANT_READS_COLUMNS = ["read_name", "mate_chr", "mate_position"]
DISCORDANT_READS_WITH_MAPQ_COLUMNS = [
    "read_name",
    "mate_chr",
    "mate_position",
    "mate_mapq",
    "mate_strand",
    "status",
]
WINDOWS_COLUMNS = [
    "PID",
    "window",
    "chrom",
    "chromStart",
    "chromEnd",
    "strand",
    "tumor_discordant_read_count",
    "control_discordant_read_count",
    "blacklisted",
    "_tumor_read_names",
    "_control_read_names",
]
FUSION_READS_COLUMNS = [
    "window",
    "read_name",
    "read_1_2",
    "start",
    "end",
    "cigar",
    "chr_primary_align",
    "coord_primary_align",
    "strand_primary_align",
    "sequence",
    "clipped_sequence",
    "part_telomere",
    "TTAGGG_count",
    "CCCTAA_count",
    "expected_pos_fusion",
]
BED_COLUMNS = ["#chrom", "chromStart", "chromEnd", "pos", "pid"]


def read_tsv(path):
    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        na_filter=False,
    )


def write_tsv(df, path, columns):
    path_obj = Path(path)
    if path_obj.exists():
        path_obj.unlink()
    out = df.copy()
    for col in columns:
        if col not in out.columns:
            out[col] = ""
    out = out[columns]
    out = sanitize_tsv_values(out)
    out.to_csv(path, sep="\t", index=False, encoding="utf-8")


def sanitize_tsv_values(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].map(_sanitize_value)
    return df


def _sanitize_value(value: object) -> str:
    if value is None:
        return ""
    text = str(value)
    if "\x00" in text:
        text = text.replace("\x00", "")
    return text
