# Shared locus, different peptide — the derivation

Structural variants recur between patients and cell lines far more often than
the neopeptides they produce do. This document derives that gap step by step,
says what each number is counted over, and separates the causes.

It exists because the gap is the easiest way to overstate a recurrence result.
"We found matching SVs between patients and the cell lines" is true here, and
means far less than it sounds: what matched was usually a *place*, not a
*consequence*.

All figures are generated from `results/master_sv.tsv` by
`tools/analyse_locus_vs_peptide.py` and refreshed into this document by
`tools/sync_event_tables.py`, so no number here is transcribed by hand.

---

## What is being modelled, and what is not

**This is a decomposition, not a statistical model.** There is no fitting, no
inference, no null. Every figure is a count of rows satisfying a stated
predicate, and each predicate reads columns that are in the shared table. That
is deliberate: the question — *why do junctions match more often than peptides?*
— is answered by partitioning the data, and a model would add assumptions where
none are needed.

What *is* modelled are the definitions. Four choices determine every number
below, and a different choice would give different figures without either being
wrong:

| Choice | Value here | What it decides |
|---|---|---|
| What counts as the same place | distance to the nearest patient breakpoint, minimum over the two breakends | the proximity ladder |
| How near is "near" | 1,000 bp | the population of the cause breakdown |
| What counts as the same mechanism | same gene *and* same SV type | the level-2 rung |
| What counts as the same consequence | exact peptide string identity | the level-1 rung |

The distances are reported unthresholded in `patient_bp_dist_bp`, so any of
these can be recomputed against a different cutoff without re-running the
pipeline.

**One property of the ladder matters more than the thresholds**: it is ordered
by *precedence*, not by distance. A junction is labelled with the strongest
level it reaches, and a gene or peptide match outranks proximity. So a row
labelled `identical_peptide` may sit far from any patient breakpoint — one does,
at 165 kb — and reading the ladder as a distance scale would be wrong. Measured
here: 4 junctions beyond 100 kb carry a gene- or peptide-level label.

---

## Step 0 — the universe

Every admitted junction across the four cell lines, panel of normals reported but
not enforced. Junction, not breakend, and not peptide: the two breakends of one
rearrangement are one row.

<!-- analysis:locus_vs_peptide_per_line -->
| Level | RPE1-TP53 | RPE1-TP53-BRCA1 | RPE1-TP53-BRCA2 | RPE1-WT |
|---|---|---|---|---|
| `identical_peptide` | 0 | 1 | 0 | 53 |
| `same_gene_and_svtype` | 3 | 3 | 5 | 316 |
| `same_gene` | 2 | 2 | 4 | 107 |
| `breakpoint_within_1kb` | 3 | 5 | 3 | 382 |
| `breakpoint_within_10kb` | 20 | 35 | 19 | 1,167 |
| `breakpoint_within_100kb` | 32 | 41 | 42 | 2,615 |
| `not_seen_in_patients` | 28 | 26 | 27 | 2,530 |
| **total** | **88** | **113** | **100** | **7,170** |
<!-- /analysis:locus_vs_peptide_per_line -->

RPE1-WT dominates by 63-81x depending on the line compared. Its call set is the whole background genome,
while the derived lines carry only what was acquired relative to their parent, so
no count below should be compared across lines without that in view.

## Step 1 — distance to the nearest patient breakpoint

For each breakend, the distance to the closest breakpoint in the patient
catalogue on the same chromosome; the junction takes the **smaller** of its two.

Taking the minimum is a choice in favour of sensitivity: it asks "is either end
near a patient event?". Whether *both* ends are near is a separate and stronger
question, carried in `patient_bp_both_within_10kb`.

<!-- analysis:locus_vs_peptide_distance -->
| Distance to nearest patient breakpoint | Junctions | % |
|---|---|---|
| 0 (same coordinate) | 105 | 1.4 |
| 1 - 1,000 | 407 | 5.4 |
| 1,001 - 10,000 | 1,451 | 19.4 |
| 10,001 - 100,000 | 2,893 | 38.7 |
| > 100,000 | 2,615 | 35.0 |
| no patient breakpoint on that chromosome | 0 | 0.0 |
<!-- /analysis:locus_vs_peptide_distance -->

Every junction has a defined distance — the catalogue covers enough of the genome
that no admitted junction is on a chromosome without one. Proximity alone is
therefore nearly meaningless at the wide end: **65% of junctions land within
100 kb of some patient breakpoint**, which is a statement about breakpoint
density more than about recurrence.

## Steps 2–3 — the evidence ladder

Each junction is assigned the strongest level it reaches, by precedence:

| Level | Requires | What it demonstrates |
|---|---|---|
| `identical_peptide` | a candidate peptide string present verbatim in the catalogue | a shared **consequence** |
| `same_gene_and_svtype` | the same gene broken, by the same class of event | a shared **mechanism** |
| `same_gene` | the same gene broken | a shared **target** |
| `breakpoint_within_Nkb` | proximity only | a shared **locus** |
| `not_seen_in_patients` | none of the above within 100 kb | — |

<!-- analysis:locus_vs_peptide -->
| Level reached | Junctions | % |
|---|---|---|
| **`identical_peptide`** | 54 | 0.7 |
| `same_gene_and_svtype` | 327 | 4.4 |
| `same_gene` | 115 | 1.5 |
| `breakpoint_within_1kb` | 393 | 5.3 |
| `breakpoint_within_10kb` | 1,241 | 16.6 |
| `breakpoint_within_100kb` | 2,730 | 36.5 |
| `not_seen_in_patients` | 2,611 | 34.9 |
| *total* | 7,471 | |

| Junctions within 1,000 bp of a patient breakpoint | n |
|---|---|
| produce no peptide at all — cannot match by construction | 354 |
| **produce peptides, none matches** | **135** |
| produce peptides and at least one matches | 23 |
| *total* | 512 |

| At the *same coordinate* as a patient breakpoint | n |
|---|---|
| junctions | 105 |
| of which same SV type as the patient event | 92 |
| of which carry inserted bases at the junction | 77 (median 7 bp) |
| of which produce peptides | 91 |
| **of which share a peptide** | **19** |
<!-- /analysis:locus_vs_peptide -->

**The headline: 4,806 junctions reach a locus or mechanism level against 54
reaching an identical peptide — 89:1.**

## Steps 4–5 — why the gap exists

Restricting to the 512 junctions within 1 kb of a patient breakpoint, three
different things prevent a peptide match. They are not interchangeable, and
collapsing them into one "did not match" hides the only interesting case:

| | n | What it means |
|---|---|---|
| produce no peptide at all | 354 | **Could not have matched.** Intergenic, intronic without protein consequence, or producing sequence identical to the wild type. |
| produce peptides, none matches | 135 | **Same place, different sequence.** The case this document is about. |
| produce peptides, at least one matches | 23 | |

Most near-misses are not near-misses at all: 69% of junctions beside a patient
breakpoint never had a peptide to compare.

## Steps 6–8 — what makes the sequence differ

The strongest form of "same place" is landing on the *exact* coordinate. Of
those, most are also the same class of event — and still do not share a peptide.

The mechanism that separates them is the sequence inserted at the junction:

| At a shared coordinate, producing peptides | median `insert_len` |
|---|---|
| peptides do **not** match | 6 bp |
| peptides match | 2 bp |

73% of same-coordinate junctions carry an insertion at all (median 7 bp). Two
breaks at one position with different inserted bases build different fusion
sequences — and an insertion whose length is not a multiple of three shifts the
reading frame, so everything downstream translates to an unrelated sequence
rather than a similar one. Observed at distance 0 with no match: ANO9 (310
inserted bases), PLAT (301), THSD7B (95), CNTNAP2 (63).

Two further mechanisms appear in the same subset:

- **the SV class differs at one locus** — ANO9 and PLAT are duplications where
  the patient event is a deletion; same breakpoint, different consequence
- **the lesion size differs** while the breakpoint coincides — the junction is
  in the same place but joins different sequence

## Step 9 — even a match is mostly not a match

Of the peptides a *matching* junction produces, the median fraction found in the
catalogue is **0.040**; the minimum is 0.003 — PPP1R12A produces 302 candidate
peptides of which 1 is in the catalogue.

This is expected. An 8–11mer sliding window over a fusion protein yields dozens
of overlapping sequences, and the catalogue holds only those that passed its own
filters. But it is the quantitative reason the level-1 cross is exact sequence
identity rather than proximity: **locus-level overlap almost never implies
peptide-level overlap**, and the two must never be quoted interchangeably.

---

## The full derivation

Every step, its population, its operation and the columns it reads. Each row can
be reproduced from `master_sv.tsv` with the named columns alone.

<!-- analysis:locus_vs_peptide_steps -->
| Step | Population | Operation | Columns read | n |
|---|---|---|---|---|
| 0 | every admitted junction, all lines | rows of master_sv.tsv | `—` | **7,471** |
| 1 | the same junctions | distance to the nearest patient breakpoint, min over the two breakends | `patient_bp_dist_bp` | **7,471** |
| 2 | the same junctions | assign the highest evidence level reached, by precedence | `patient_evidence` | **7,471** |
| 3a | all junctions | reach a locus or mechanism level, but not a peptide | `patient_evidence` | **4,806** |
| 3b | all junctions | reach an identical peptide | `patient_evidence` | **54** |
| 4 | all junctions | restrict to within 1,000 bp of a patient breakpoint | `patient_bp_dist_bp` | **512** |
| 5a | step 4 | produce no candidate peptide — cannot match | `n_peptides_generated` | **354** |
| 5b | step 4 | produce peptides, none identical to the catalogue | `n_peptides_generated, n_peptides_matched` | **135** |
| 5c | step 4 | produce peptides and at least one matches | `n_peptides_matched` | **23** |
| 6 | all junctions | restrict to the SAME coordinate as a patient breakpoint | `patient_bp_dist_bp` | **105** |
| 7a | step 6 | also the same SV type as the patient event | `svtype, patient_svtype` | **92** |
| 7b | step 6 | carry inserted bases | `insert_len` | **77** |
| 8a | step 6 | produce peptides | `n_peptides_generated` | **91** |
| 8b | step 8a | share a peptide | `n_peptides_matched` | **19** |
<!-- /analysis:locus_vs_peptide_steps -->

---

## What this does not establish

- **Nothing about causation.** A shared locus may reflect a fragile site, a
  mappability artefact, or convergent biology; this counts co-occurrence and
  does not distinguish them. The catalogue's own `is_cfs` column flags common
  fragile sites and is carried through untouched for exactly this reason.
- **Nothing about the derived lines.** RPE1-TP53, -BRCA1 and -BRCA2 contribute
  301 junctions between them against RPE1-WT's 7,170. Their rungs are single
  digits and should not be read as rates.
- **No claim that 1 kb is the right window.** It is a reporting choice. The
  unthresholded distance is in the table.
- **No correction for breakpoint density.** Junctions and patient breakpoints
  are both non-uniform along the genome, so proximity at the 100 kb scale is
  partly a statement about where breaks happen at all. The permutation null in
  the main analysis addresses this for peptide matches; it is not applied to the
  proximity rungs.

## Reproducing

```bash
python tools/build_master_table.py --runs <run> --suffix <s>   # the input table
python tools/analyse_locus_vs_peptide.py                       # every figure here
python tools/sync_event_tables.py                              # refresh this file
```

Figures come from `patient_evidence`, `patient_bp_dist_bp`, `insert_len`,
`n_peptides_generated`, `n_peptides_matched`, `svtype` and `patient_svtype`. All
are defined in `results/master_sv_column_dictionary.tsv`.
