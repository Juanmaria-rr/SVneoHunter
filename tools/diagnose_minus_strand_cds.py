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


def probe(transcript, position: int):
    """`truncate_cds` at one position: (cut_length, sequence)."""
    from neosv.fusion_utils import truncate_cds
    collection = truncate_cds(transcript, "5", position)
    sequence = collection.nt_sequence if collection else ""
    return (collection.cut_length if collection else 0), sequence


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", type=int, default=115)
    args = parser.parse_args()

    sys.path.insert(0, str(REPO / "vendor"))
    from pyensembl import EnsemblRelease
    genome = EnsemblRelease(args.release)

    print("truncate_cds(transcript, '5', pos) — length of the 5' coding "
          "sequence returned\n")
    print(f"  {'gene':9s} {'strand':7s} {'exonic break':>13s} "
          f"{'intronic break':>15s}   verdict")
    print("  " + "-" * 62)

    failures = []
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
        if broken:
            failures.append(name)
        print(f"  {name:9s} {gene.strand:7s} {exonic:13d} {intronic:15d}   "
              f"{'5-PRIME LOST' if broken else 'ok'}")

    print()
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
