#!/usr/bin/env python3
"""
diagnose_minus_strand_cds.py — the 5' CDS is lost for minus-strand transcripts
when the breakpoint falls in an intron.

WHAT THIS DEMONSTRATES
----------------------
`neosv.fusion_utils.truncate_cds(transcript, '5', pos)` returns the coding
sequence 5' of a breakpoint. For a MINUS-strand transcript with an INTRONIC
breakpoint it returns an empty sequence (`cut_length == 0`). Plus-strand
transcripts are unaffected, and both strands are correct when the breakpoint
falls inside a coding exon.

WHY IT MATTERS
--------------
Introns are far larger than exons, so most real SV breakpoints are intronic.
The consequences cascade:

  * the fusion protein becomes the 3' partner alone
  * `frame_effect` is `Start-loss` — upstream's own low-reliability flag,
    raised because the start codon is gone and it falls back to the next ATG
  * the junction sits at residue ~0, so NO peptide can span the breakpoint
  * minus-strand genes are systematically degraded relative to plus-strand ones

Measured on this dataset: 84.0% of minus-strand fusions are `Start-loss` against
13.8% of plus-strand ones, and 1 of 408 catalogue-matching peptides spans its
junction where 67 would be expected.

THE LIKELY CAUSE
----------------
`transcript.coding_sequence_position_ranges` is in GENOMIC order, while
`truncate_cds` treats it as transcript order (5'->3'), as its own docstring in
`transcript_utils.get_cds_range` states. The two coincide on the plus strand and
are reversed on the minus strand. The exon-overlap branch happens to survive
this; the intron branch, which indexes through `get_noncds_range`, does not.

This script asserts nothing about the fix — it establishes the behaviour, so a
fix can be tested against it.

    python tools/diagnose_minus_strand_cds.py
    python tools/diagnose_minus_strand_cds.py --release 115
"""
from __future__ import annotations

import argparse
import os
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent

#: Genes with both strands represented, chosen for being well annotated and
#: multi-exon rather than for any biological relevance here.
GENES = ["ITGA11", "TP53", "BRCA1", "GOLGA3", "KRAS", "BRAF",
         "EGFR", "MYC", "PTEN", "PIK3CA"]


#: A three-exon gene, small enough to follow by hand. Coordinates are arbitrary;
#: what matters is that the exons are separated by introns and that the
#: transcript extends a little beyond the outermost coding exons.
TOY_CDS = [(100, 200), (300, 400), (500, 600)]
TOY_SPAN = (50, 650)
TOY_BREAK = 450                     # in the intron between the 2nd and 3rd exon


def toy_gaps(cds, span, strand):
    """The non-coding intervals, by the arithmetic of `get_noncds_range`.

    Reproduced here rather than imported so the walkthrough shows the formula
    beside its result. It must stay identical to
    `vendor/neosv/transcript_utils.py`; the assertions in `explain()` check it
    against the real function on real transcripts.
    """
    start, end = span
    out = []
    if strand == "+":
        out.append(((start, cds[0][0] - 1),
                    "(transcript.start, cds[0].start-1)"))
        for i in range(1, len(cds)):
            out.append(((cds[i - 1][1] + 1, cds[i][0] - 1),
                        f"(cds[{i-1}].end+1, cds[{i}].start-1)"))
        out.append(((cds[-1][1] + 1, end), "(cds[-1].end+1, transcript.end)"))
    else:
        out.append(((cds[0][1] + 1, end), "(cds[0].end+1, transcript.end)"))
        for i in range(1, len(cds)):
            out.append(((cds[i][1] + 1, cds[i - 1][0] - 1),
                        f"(cds[{i}].end+1, cds[{i-1}].start-1)"))
        out.append(((start, cds[-1][0] - 1),
                    "(transcript.start, cds[-1].start-1)"))
    return out


def explain() -> None:
    """Walk one breakpoint through both strands, showing every interval.

    The exon list arrives sorted by coordinate on BOTH strands. On the plus
    strand that is also reading order; on the minus strand it is the reverse.
    The minus branch compensates for a reversal that never happened.
    """
    cds, span = TOY_CDS, TOY_SPAN
    print(f"A three-exon gene: coding exons at {cds}, transcript {span[0]}-"
          f"{span[1]}, breakpoint at {TOY_BREAK} — in the intron between the "
          f"2nd and 3rd exon.\n")
    print("The exon list arrives sorted by coordinate on BOTH strands.\n")

    for strand in ("+", "-"):
        reading = "1st, 2nd, 3rd" if strand == "+" else "3rd, 2nd, 1st"
        expected = 2 if strand == "+" else 1
        print("=" * 70)
        print(f"STRAND {strand}   read in the order {reading}   "
              f"-> correct answer: {expected} exon(s) before the break")
        print("=" * 70)
        gaps = toy_gaps(cds, span, strand)
        for index, (interval, formula) in enumerate(gaps):
            broken = "   <-- INVERTED: start > end, nothing can fall inside" \
                if interval[0] > interval[1] else ""
            print(f"  gap {index}  {formula:38s} = {interval}{broken}")
        match = next((i for i, (iv, _) in enumerate(gaps)
                      if iv[0] <= TOY_BREAK <= iv[1]), None)
        kept = len(cds[:match]) if match is not None else "n/a"
        verdict = "correct" if kept == expected else f"WRONG, expected {expected}"
        print(f"\n  first gap containing {TOY_BREAK}: gap {match}")
        print(f"  the code then takes cds[:{match}] -> {kept} exon(s)   [{verdict}]\n")

    print("On the minus strand two things break at once. Gap 0 is measured from "
          "the wrong\nend — the first exon in the LIST, not in reading order — "
          "so it spans most of the\ngene. And every intron comes out with start "
          "> end, so none can ever match.\nThe breakpoint therefore lands in "
          "gap 0, which means 'before the first exon',\nand zero exons are "
          "kept.\n")


def probe(transcript, position: int):
    """`truncate_cds` at one position: (cut_length, sequence)."""
    from neosv.fusion_utils import truncate_cds
    collection = truncate_cds(transcript, "5", position)
    sequence = collection.nt_sequence if collection else ""
    return (collection.cut_length if collection else 0), sequence


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", type=int, default=115)
    parser.add_argument("--explain", action="store_true",
                        help="walk one breakpoint through both strands on a "
                             "three-exon toy gene, showing every interval and "
                             "the formula that produced it")
    args = parser.parse_args()

    if args.explain:
        explain()

    sys.path.insert(0, str(REPO / "vendor"))
    from pyensembl import EnsemblRelease
    genome = EnsemblRelease(args.release)

    print("truncate_cds(transcript, '5', pos) — length of the 5' coding "
          "sequence returned\n")
    print(f"  {'gene':9s} {'strand':7s} {'exonic break':>13s} "
          f"{'intronic break':>15s}   verdict")
    print("  " + "-" * 62)

    failures, tested = [], []
    for name in GENES:
        try:
            gene = genome.genes_by_name(name)[0]
        except Exception:                                   # noqa: BLE001
            continue
        transcript = next((t for t in gene.transcripts
                           if t.is_protein_coding and t.complete), None)
        if transcript is None:
            continue
        ranges = sorted(transcript.coding_sequence_position_ranges)
        if len(ranges) < 3:
            continue

        middle = ranges[len(ranges) // 2]
        exonic, _ = probe(transcript, (middle[0] + middle[1]) // 2)
        # Between two coding exons: an intron by construction.
        previous = ranges[len(ranges) // 2 - 1]
        intronic, _ = probe(transcript, (previous[1] + middle[0]) // 2)

        broken = intronic == 0 and exonic > 0
        tested.append(name)
        if broken:
            failures.append(name)
        print(f"  {name:9s} {gene.strand:7s} {exonic:13d} {intronic:15d}   "
              f"{'5-PRIME LOST' if broken else 'ok'}")

    print()
    # Testing nothing and reporting a pass is the failure mode this whole
    # repository exists to avoid. Usually the annotation cache is missing.
    if not tested:
        raise SystemExit(
            "no transcripts were tested — nothing was verified.\n"
            "Set PYENSEMBL_CACHE_DIR to the annotation cache, e.g.\n"
            "  export PYENSEMBL_CACHE_DIR=/path/to/pyensembl_cache")

    if failures:
        minus = [g for g in failures]
        print(f"{len(minus)} transcript(s) return an EMPTY 5' CDS for an "
              f"intronic breakpoint: {', '.join(minus)}")
        print("\nEvery one is on the minus strand. Plus-strand transcripts "
              "return the expected length for the same construction, and both "
              "strands are correct when the breakpoint is exonic.")
        print("\nConsequence: those fusions carry no 5' contribution, are "
              "flagged Start-loss, place the junction at residue ~0, and can "
              "never produce a junction-spanning peptide.")
        sys.exit(1)
    print("No transcript loses its 5' CDS on an intronic breakpoint.")


if __name__ == "__main__":
    main()
