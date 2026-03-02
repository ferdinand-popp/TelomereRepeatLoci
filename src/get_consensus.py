#!/usr/bin/env python3

import argparse
import csv
import re

import pysam


EMPTY_VALUES = {"", "NA", "NaN", "nan", "None", None}


def parse_int(value):
    if value in EMPTY_VALUES:
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def reverse_seq(seq):
    return seq[::-1]


def split_clipped(clipped_sequence):
    if not clipped_sequence:
        return []
    return [part.strip() for part in str(clipped_sequence).split(",") if part.strip()]


def consensus_from_sequences(sequences):
    if not sequences:
        return ""
    max_len = max(len(s) for s in sequences)
    consensus = []
    for i in range(max_len):
        bases = [s[i] for s in sequences if i < len(s)]
        if not bases:
            continue
        read_count = len(bases)
        freqs = {
            "A": bases.count("A") / read_count,
            "C": bases.count("C") / read_count,
            "G": bases.count("G") / read_count,
            "T": bases.count("T") / read_count,
            "N": bases.count("N") / read_count,
        }
        base, max_freq = max(freqs.items(), key=lambda x: x[1])
        consensus.append(base if max_freq >= 0.65 else "N")
    return "".join(consensus)


def microhomology(seq, consensus, repeat_forward, strand):
    if not consensus or repeat_forward in EMPTY_VALUES:
        return "?"

    pos_match = [m.start() for m in re.finditer(r"TTAGGG|CCCTAA", consensus)]
    if not pos_match:
        return "?"

    repeat_junction_wrong = False

    if strand == "+":
        indice_match = pos_match[0] + 1
        junction_repeat_ref = repeat_forward[: 6 - indice_match + 1]
        junction_repeat_tel = repeat_forward[6 - indice_match + 1 : 6]
        if junction_repeat_tel != consensus[: indice_match - 1]:
            repeat_junction_wrong = True
        repeats_ref = (repeat_forward * 3 + junction_repeat_ref)[::-1]
        seq_homology = seq[::-1]
    elif strand == "-":
        last = pos_match[-1]
        indice_match = len(consensus) + 1 - (last + 6)
        junction_repeat_tel = repeat_forward[:indice_match]
        junction_repeat_ref = repeat_forward[indice_match:6]
        if junction_repeat_tel != consensus[last + 6 : len(consensus) + 1]:
            repeat_junction_wrong = True
        repeats_ref = junction_repeat_ref + repeat_forward * 3
        seq_homology = seq
    else:
        return "?"

    if indice_match > 6 or repeat_junction_wrong:
        return "?"

    cntr = 0
    for i in range(min(20, len(seq_homology), len(repeats_ref))):
        if seq_homology[i] == repeats_ref[i]:
            cntr += 1
        else:
            break
    return str(cntr)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_regions_file")
    parser.add_argument("clipped_read_file")
    parser.add_argument("outfile")
    parser.add_argument("--reference", default="")
    args = parser.parse_args()

    with open(args.candidate_regions_file, newline="") as handle:
        candidate_reader = csv.DictReader(handle, delimiter="\t")
        candidate_rows = list(candidate_reader)
        candidate_fields = candidate_reader.fieldnames or []

    with open(args.clipped_read_file, newline="") as handle:
        clipped_rows = list(csv.DictReader(handle, delimiter="\t"))

    clipped_by_window = {}
    for row in clipped_rows:
        clipped_by_window.setdefault(row.get("window", ""), []).append(row)

    fasta = pysam.FastaFile(args.reference) if args.reference else None

    extra_fields = ["consensus", "flanking_seq", "bp_microhomology"]
    out_fields = candidate_fields + [f for f in extra_fields if f not in candidate_fields]

    for region in candidate_rows:
        window = region.get("window", "")
        insertion_site = parse_int(region.get("insertion_site"))
        strand = region.get("strand", "")

        region["consensus"] = ""
        region["flanking_seq"] = ""
        region["bp_microhomology"] = ""

        if insertion_site is None:
            continue

        start_end = "end" if strand == "+" else "start"
        candidates = [
            r for r in clipped_by_window.get(window, []) if parse_int(r.get(start_end)) == insertion_site
        ]

        sequences = []
        for r in candidates:
            parts = split_clipped(r.get("clipped_sequence"))
            if not parts:
                continue
            if start_end == "end":
                sequences.append(parts[-1])
            else:
                sequences.append(reverse_seq(parts[0]))

        consensus = consensus_from_sequences(sequences)
        if strand == "-":
            consensus = reverse_seq(consensus)
        region["consensus"] = consensus

        flanking = ""
        if fasta is not None:
            chrom = region.get("chrom", "")
            chrom_chr = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
            try:
                if strand == "+":
                    flanking = fasta.fetch(chrom_chr, max(0, insertion_site - 21), insertion_site - 1)
                elif strand == "-":
                    flanking = fasta.fetch(chrom_chr, insertion_site - 1, insertion_site + 19)
            except Exception:
                flanking = ""
        region["flanking_seq"] = flanking

        region["bp_microhomology"] = microhomology(
            flanking,
            consensus,
            region.get("repeat_forward"),
            strand,
        )

    if fasta is not None:
        fasta.close()

    with open(args.outfile, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(candidate_rows)


if __name__ == "__main__":
    main()
