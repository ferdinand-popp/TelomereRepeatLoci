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
    out = df.copy()
    for col in columns:
        if col not in out.columns:
            out[col] = ""
    out = out[columns]
    out.to_csv(path, sep="\t", index=False)
