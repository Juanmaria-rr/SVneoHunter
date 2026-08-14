#!/usr/bin/env python3
"""
check_strand_bias.py — is the gene strand of our candidates distributed as chance?

WHY
---
An SV-annotation pipeline that treats minus-strand transcripts with plus-strand
logic produces a systematic error: the 5'/3' assignment inverts, and the fusion
protein is built from the wrong side. The failure is silent — genes are still
named, peptides are still emitted — so it has to be looked for.

This measures the strand composition at each level of the analysis and tests it
against a background. Where the bias first appears localises its cause:

    reference catalogue only   -> inherited from upstream, not ours
    candidate peptides too     -> our generation or annotation
    only the intersection      -> something in the cross itself

BACKGROUND CHOICE MATTERS
-------------------------
The naive background is "half of all coding genes". The better one weights by
gene length, because a structural variant is likelier to break a long gene, and
if long genes were strand-skewed the naive test would misfire. Both are
reported; on GRCh38 they are close (49.3% vs 48.5%), which is itself worth
knowing.

WHAT A POSITIVE RESULT DOES AND DOES NOT MEAN
---------------------------------------------
A significant skew is evidence that something is strand-dependent. It is NOT
proof of an annotation artefact: gene selection, the reference cohort's own
composition, or a genuinely strand-correlated biology would all produce it.
Rule out the annotation directly (see docs/OPEN_QUESTIONS.md) rather than
inferring it from this test.

USAGE
-----
    PYENSEMBL_CACHE_DIR=... python tools/check_strand_bias.py \\
        --results-dir results [--reference /path/to/catalogue.tsv] \\
        [--release 115]
"""
from __future__ import annotations

import argparse
import collections
import os
import pathlib
import sys

import pandas as pd


def gene_index(release: int) -> dict:
    """gene symbol -> (strand, genomic length), for coding genes."""
    from pyensembl import EnsemblRelease
    genome = EnsemblRelease(release)
    index = {}
    for gene in genome.genes():
        if gene.biotype == "protein_coding" and gene.gene_name not in index:
            index[gene.gene_name] = (gene.strand, gene.end - gene.start)
    return index


def backgrounds(index: dict) -> tuple[float, float]:
    """(by gene count, by genomic length) fraction of coding genes on minus."""
    counts = collections.Counter(strand for strand, _ in index.values())
    by_count = counts["-"] / sum(counts.values())
    lengths = collections.Counter()
    for strand, length in index.values():
        lengths[strand] += length
    by_length = lengths["-"] / sum(lengths.values())
    return by_count, by_length


def test_set(label: str, genes, index: dict, p_null: float) -> dict:
    """Binomial test of one gene set's strand composition against the background."""
    from scipy import stats
    observed = [index[g] for g in {str(x) for x in genes if pd.notna(x)} if g in index]
    n = len(observed)
    if not n:
        return {"set": label, "n": 0}
    minus = sum(1 for strand, _ in observed if strand == "-")
    p = stats.binomtest(minus, n, p_null).pvalue
    # `median()` of an empty series is NaN, and `NaN or 0` is NaN, not 0 —
    # NaN is truthy. The same trap has bitten this codebase three times now:
    # never fall back on truthiness when the value can be NaN.
    def median_kb(strand):
        lengths = [L for st, L in observed if st == strand]
        return round(pd.Series(lengths).median() / 1000) if lengths else None

    return {"set": label, "n": n, "minus": minus, "fraction_minus": round(minus / n, 4),
            "p_value": p, "median_len_minus_kb": median_kb("-"),
            "median_len_plus_kb": median_kb("+")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--reference", default=None,
                    help="reference catalogue TSV; defaults to the configured one")
    ap.add_argument("--release", type=int, default=115)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    repo = pathlib.Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(repo / "src"))

    index = gene_index(args.release)
    by_count, by_length = backgrounds(index)
    print(f"background, coding genes on the minus strand:")
    print(f"  by gene count      {by_count:.1%}")
    print(f"  by genomic length  {by_length:.1%}   <- the fairer null: SVs break "
          f"long genes more often\n")

    rows = []
    reference = args.reference
    if reference is None:
        from svneo import config as config_mod
        cfg = config_mod.load(str(repo / "config" / "rpe1_hmf.yaml"))
        reference, gene_column = cfg.reference.peptides, cfg.reference.gene_column
    else:
        gene_column = "gene1"
    if os.path.exists(reference):
        ref = pd.read_csv(reference, sep="\t", low_memory=False)
        if gene_column in ref.columns:
            rows.append(test_set("reference catalogue", ref[gene_column], index, by_length))

    results = pathlib.Path(args.results_dir)
    for run_dir in sorted(p for p in results.iterdir() if p.is_dir()):
        peptides = next(run_dir.glob("*.all_neopeptides.txt"), None)
        if peptides:
            table = pd.read_csv(peptides, sep="\t", low_memory=False)
            genes = pd.concat([table.get("gene1"), table.get("gene2")]).dropna()
            rows.append(test_set(f"candidates: {run_dir.name}", genes, index, by_length))

    master = results / "master_table.tsv"
    if master.exists():
        table = pd.read_csv(master, sep="\t", low_memory=False)
        rows.append(test_set("matched candidates", table.get("ref_gene"), index, by_length))
        supported = table[table.rna_tier.isin(["STRONG", "SUGGESTIVE", "WEAK"])] \
            if "rna_tier" in table.columns else table.iloc[0:0]
        if len(supported):
            rows.append(test_set("RNA-supported events", supported.ref_gene, index, by_length))

    report = pd.DataFrame(rows)
    print(report.to_string(index=False))

    if args.out:
        report.to_csv(args.out, sep="\t", index=False)
        print(f"\nwrote {args.out}")

    skewed = report[(report.n > 20) & (report.p_value < 0.01)] \
        if "p_value" in report.columns else report.iloc[0:0]
    if len(skewed):
        print(f"\n{len(skewed)} set(s) deviate from the background at p < 0.01.")
        print("Where the deviation FIRST appears localises the cause; see")
        print("docs/OPEN_QUESTIONS.md before attributing it to annotation.")


if __name__ == "__main__":
    main()
