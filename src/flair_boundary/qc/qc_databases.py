#!/usr/bin/env python3
"""QC report for the four boundary databases.

Reports counts, peak widths, chromosome coverage, strand balance, and a
small sequence-composition sample for FANTOM, refTSS, PolyA_DB BED, and
the PolyA_DB PAS annotation.

Usage:
    python -m flair_boundary.qc.qc_databases \\
        --fantom TSS_db/FANTOM_TSS_human.bed \\
        --reftss TSS_db/refTSS_v4.1_human_coordinate.hg38.bed.txt \\
        --polyadb TTS_db/polyadb.hg38.weighted.bed \\
        --polyadb-pas TTS_db/polyAdb_human.PAS.txt \\
        --genome data/ref/GRCh38.primary_assembly.genome.fa \\
        --output qc/database_qc_report.md
"""
from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd


CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}


def _peak_width_summary(widths: np.ndarray) -> dict:
    return {
        "n": int(len(widths)),
        "min": int(widths.min()) if len(widths) else 0,
        "p10": int(np.percentile(widths, 10)) if len(widths) else 0,
        "p25": int(np.percentile(widths, 25)) if len(widths) else 0,
        "median": int(np.percentile(widths, 50)) if len(widths) else 0,
        "p75": int(np.percentile(widths, 75)) if len(widths) else 0,
        "p90": int(np.percentile(widths, 90)) if len(widths) else 0,
        "p99": int(np.percentile(widths, 99)) if len(widths) else 0,
        "max": int(widths.max()) if len(widths) else 0,
        "mean": float(widths.mean()) if len(widths) else 0.0,
    }


def _composition_sample(
    df: pd.DataFrame,
    genome_fa: str,
    n_sample: int = 1000,
    upstream: int = 200,
    downstream: int = 200,
    seed: int = 42,
) -> dict:
    """Sample N random positives, fetch ±200bp window, summarise base composition."""
    if not genome_fa or not Path(genome_fa).exists():
        return {}
    import pysam
    fa = pysam.FastaFile(genome_fa)
    rng = np.random.default_rng(seed)
    n = min(n_sample, len(df))
    sampled = df.sample(n=n, random_state=seed)

    total = Counter()
    cpg = 0
    seqs_taken = 0
    for _, row in sampled.iterrows():
        chrom = row["chrom"]
        try:
            chrom_len = fa.get_reference_length(chrom)
        except KeyError:
            continue
        pos = int(row["pos"])
        s, e = pos - upstream, pos + downstream
        if s < 0 or e > chrom_len:
            continue
        seq = fa.fetch(chrom, s, e).upper()
        if len(seq) != upstream + downstream:
            continue
        if row.get("strand") == "-":
            seq = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        for b in seq:
            total[b] += 1
        cpg += seq.count("CG")
        seqs_taken += 1

    fa.close()
    if seqs_taken == 0:
        return {}
    bases = sum(total.values())
    return {
        "n_sequences": seqs_taken,
        "window_bp": upstream + downstream,
        "frac_A": total.get("A", 0) / bases,
        "frac_C": total.get("C", 0) / bases,
        "frac_G": total.get("G", 0) / bases,
        "frac_T": total.get("T", 0) / bases,
        "frac_N": total.get("N", 0) / bases,
        "GC_content": (total.get("G", 0) + total.get("C", 0)) / bases,
        "CpG_per_kb": 1000.0 * cpg / bases,
    }


def qc_fantom(bed_file: str, genome_fa: Optional[str]) -> dict:
    # Skip the leading 'track ...' header line, then parse as BED9
    df = pd.read_csv(
        bed_file, sep="\t", header=None, skiprows=1,
        names=["chrom", "start", "end", "name", "score", "strand",
               "thickStart", "thickEnd", "rgb"],
        dtype={"start": int, "end": int, "thickStart": int, "thickEnd": int},
    )
    # Sanity check: how often does thickStart genuinely diverge from start?
    n_thick_diff = int((df["thickStart"] != df["start"]).sum())

    df["pos"] = df["thickStart"]
    widths = (df["end"] - df["start"]).clip(lower=1).values

    chrom_counts = df["chrom"].value_counts()
    canonical = df[df["chrom"].isin(CANONICAL_CHROMS)]

    out = {
        "source": "FANTOM CAGE peaks",
        "file": bed_file,
        "n_peaks_total": int(len(df)),
        "n_peaks_canonical": int(len(canonical)),
        "n_peaks_other_contigs": int(len(df) - len(canonical)),
        "strand_balance": dict(df["strand"].value_counts()),
        "width_bp": _peak_width_summary(widths),
        "top_chroms": chrom_counts.head(5).to_dict(),
        "chrom_with_zero": sorted(CANONICAL_CHROMS - set(chrom_counts.index)),
        "summit_field": "thickStart (peak summit)",
        "thickStart_diff_from_start": n_thick_diff,
    }
    if genome_fa:
        out["sequence_composition_sampled"] = _composition_sample(
            canonical[["chrom", "pos", "strand"]], genome_fa
        )
    return out


def qc_reftss(bed_file: str, genome_fa: Optional[str]) -> dict:
    df = pd.read_csv(bed_file, sep="\t")
    df = df.rename(columns={"chromosome": "chrom"})
    df["pos"] = ((df["start"].astype(int) + df["end"].astype(int)) // 2).astype(int)
    widths = (df["end"].astype(int) - df["start"].astype(int)).clip(lower=1).values

    canonical = df[df["chrom"].isin(CANONICAL_CHROMS)]
    out = {
        "source": "refTSS v4.1",
        "file": bed_file,
        "n_peaks_total": int(len(df)),
        "n_peaks_canonical": int(len(canonical)),
        "strand_balance": dict(df["strand"].value_counts()),
        "width_bp": _peak_width_summary(widths),
        "summit_field": "interval midpoint (no summit annotated)",
    }
    if genome_fa:
        out["sequence_composition_sampled"] = _composition_sample(
            canonical[["chrom", "pos", "strand"]], genome_fa
        )
    return out


def qc_fantom_cat(bed_file: str, genome_fa: Optional[str]) -> dict:
    """QC a generic merged-CAGE-cluster BED file (e.g. FANTOM5-hg38 fair peaks).

    Each row is a merged CAGE cluster. Width = `end - start` is the
    per-cluster width target for the TSS regressor.
    """
    df = pd.read_csv(
        bed_file, sep="\t", header=None,
        usecols=[0, 1, 2, 3, 4, 5],
        names=["chrom", "start", "end", "name", "score", "strand"],
        dtype={"chrom": str, "start": int, "end": int,
               "name": str, "score": str, "strand": str},
    )
    # Cluster summit is best approximated by the midpoint (the BED has no
    # explicit summit field after liftOver).
    df["pos"] = ((df["start"].astype(int) + df["end"].astype(int)) // 2).astype(int)
    widths = (df["end"] - df["start"]).clip(lower=1).values

    canonical = df[df["chrom"].isin(CANONICAL_CHROMS)]
    out = {
        "source": "Merged CAGE cluster BED (width target source)",
        "file": bed_file,
        "n_peaks_total": int(len(df)),
        "n_peaks_canonical": int(len(canonical)),
        "n_peaks_other_contigs": int(len(df) - len(canonical)),
        "strand_balance": dict(df["strand"].value_counts()),
        "width_bp": _peak_width_summary(widths),
        "summit_field": "interval midpoint",
    }
    if genome_fa:
        out["sequence_composition_sampled"] = _composition_sample(
            canonical[["chrom", "pos", "strand"]], genome_fa
        )
    return out


def qc_polyadb(bed_file: str, pas_file: Optional[str], genome_fa: Optional[str]) -> dict:
    df = pd.read_csv(
        bed_file, sep="\t", header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
    )
    df["pos"] = np.where(df["strand"] == "+", df["end"] - 1, df["start"])
    widths = (df["end"] - df["start"]).clip(lower=1).values

    canonical = df[df["chrom"].isin(CANONICAL_CHROMS)]
    out = {
        "source": "PolyA_DB",
        "file": bed_file,
        "n_peaks_total": int(len(df)),
        "n_peaks_canonical": int(len(canonical)),
        "strand_balance": dict(df["strand"].value_counts()),
        "width_bp": _peak_width_summary(widths),
    }

    if pas_file and Path(pas_file).exists():
        pas_df = pd.read_csv(pas_file, sep="\t")
        # Distribution of PAS signal types
        sig_col = "PAS Signal" if "PAS Signal" in pas_df.columns else None
        if sig_col is not None:
            out["pas_signal_distribution"] = dict(
                pas_df[sig_col].fillna("missing").value_counts()
            )
        loc_col = "Intron/exon location"
        if loc_col in pas_df.columns:
            out["intron_exon_distribution"] = dict(
                pas_df[loc_col].fillna("missing").value_counts().head(10)
            )

    if genome_fa:
        out["sequence_composition_sampled"] = _composition_sample(
            canonical[["chrom", "pos", "strand"]], genome_fa
        )
    return out


def write_markdown_report(qcs: List[dict], output_path: str) -> None:
    out = ["# Boundary Database QC Report", ""]
    out.append(f"_Generated by `flair_boundary.qc.qc_databases`._")
    out.append("")
    for qc in qcs:
        if not qc:
            continue
        out.append(f"## {qc['source']}")
        out.append("")
        out.append(f"- **File:** `{qc['file']}`")
        out.append(f"- **Peaks (total):** {qc['n_peaks_total']:,}")
        out.append(f"- **Peaks on canonical chroms (chr1–22, X, Y, M):** "
                   f"{qc['n_peaks_canonical']:,}")
        if qc.get("n_peaks_other_contigs"):
            out.append(f"- **Peaks on other contigs:** {qc['n_peaks_other_contigs']:,}")
        if qc.get("strand_balance"):
            out.append(f"- **Strand balance:** {qc['strand_balance']}")
        if qc.get("summit_field"):
            out.append(f"- **Boundary coordinate:** {qc['summit_field']}")
        if "thickStart_diff_from_start" in qc:
            n = qc["thickStart_diff_from_start"]
            total = qc["n_peaks_total"]
            out.append(
                f"- **`thickStart` differs from `start`:** {n:,} / {total:,} rows "
                f"({100*n/max(total,1):.2f}%) — if 0, the BED9 has no separate "
                f"summit annotation and `thickStart` simply duplicates `start`."
            )

        w = qc["width_bp"]
        out.append("")
        out.append("**Peak width (bp):**")
        out.append("")
        out.append("| min | p10 | p25 | median | p75 | p90 | p99 | max | mean |")
        out.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        out.append(
            f"| {w['min']} | {w['p10']} | {w['p25']} | {w['median']} | "
            f"{w['p75']} | {w['p90']} | {w['p99']} | {w['max']} | {w['mean']:.1f} |"
        )

        if qc.get("top_chroms"):
            out.append("")
            out.append(
                "**Top 5 chroms:** "
                + ", ".join(f"{c}={n:,}" for c, n in qc["top_chroms"].items())
            )
        if qc.get("chrom_with_zero"):
            out.append(
                f"- Canonical chroms with zero peaks: {qc['chrom_with_zero']}"
            )

        if qc.get("pas_signal_distribution"):
            out.append("")
            out.append("**PAS signal distribution:**")
            out.append("")
            out.append("| Signal | Count |")
            out.append("|---|---:|")
            for k, v in qc["pas_signal_distribution"].items():
                out.append(f"| {k} | {v:,} |")

        if qc.get("intron_exon_distribution"):
            out.append("")
            out.append("**Intron/exon location:**")
            out.append("")
            out.append("| Location | Count |")
            out.append("|---|---:|")
            for k, v in qc["intron_exon_distribution"].items():
                out.append(f"| {k} | {v:,} |")

        comp = qc.get("sequence_composition_sampled") or {}
        if comp:
            out.append("")
            out.append(
                f"**Sequence composition** "
                f"(strand-aware, ±{comp['window_bp']//2}bp around summit, "
                f"n={comp['n_sequences']} sampled):"
            )
            out.append("")
            out.append("| %A | %C | %G | %T | %N | GC | CpG/kb |")
            out.append("|---:|---:|---:|---:|---:|---:|---:|")
            out.append(
                f"| {comp['frac_A']*100:.1f} | {comp['frac_C']*100:.1f} | "
                f"{comp['frac_G']*100:.1f} | {comp['frac_T']*100:.1f} | "
                f"{comp['frac_N']*100:.2f} | "
                f"{comp['GC_content']*100:.1f}% | {comp['CpG_per_kb']:.1f} |"
            )

        out.append("")
        out.append("---")
        out.append("")

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text("\n".join(out))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fantom", help="FANTOM TSS BED file (raw permissive peaks)")
    p.add_argument("--fantom-cat", help="FANTOM-CAT lv3_robust CAGE_cluster BED (hg38)")
    p.add_argument("--reftss", help="refTSS BED file")
    p.add_argument("--polyadb", help="PolyA_DB BED6 file")
    p.add_argument("--polyadb-pas", help="PolyA_DB PAS annotation file (TSV)")
    p.add_argument("--genome", help="Genome FASTA for sequence composition sampling")
    p.add_argument("--output", required=True, help="Output markdown path")
    args = p.parse_args()

    qcs: List[dict] = []
    if args.fantom:
        print(f"QC: FANTOM ({args.fantom})...")
        qcs.append(qc_fantom(args.fantom, args.genome))
    if args.fantom_cat:
        print(f"QC: FANTOM-CAT ({args.fantom_cat})...")
        qcs.append(qc_fantom_cat(args.fantom_cat, args.genome))
    if args.reftss:
        print(f"QC: refTSS ({args.reftss})...")
        qcs.append(qc_reftss(args.reftss, args.genome))
    if args.polyadb:
        print(f"QC: PolyA_DB ({args.polyadb})...")
        qcs.append(qc_polyadb(args.polyadb, args.polyadb_pas, args.genome))

    write_markdown_report(qcs, args.output)
    print(f"\nWrote: {args.output}")


if __name__ == "__main__":
    main()
