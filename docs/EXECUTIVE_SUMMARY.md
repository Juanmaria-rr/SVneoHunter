# Executive summary — finding and validating SV-derived neopeptides

How a catalogue of patient neopeptides is tested for recurrence in a set of
samples, what each filter measures and why its threshold sits where it does, and
what the four RPE1 cell lines actually showed.

This document is in two halves. **Part A** is the method: every stage, every
criterion, every threshold with its justification, and what each one does *not*
capture. **Part B** is the result, one section per cell line, with the cascade in
numbers and the genes named. Part A can be read without Part B; Part B cannot be
read without Part A, because a count means nothing without the criterion that
produced it.

Code references are `file.py` · `function()` · `L<line>`, and are verified
mechanically by `tools/check_docs.py`, which fails if a line reference, a
constant or a function name has drifted.

---

## The question

> A cohort of patients yields neopeptides created by structural variants. Do the
> same peptides arise in an independent sample — and if they do, is that
> recurrence real, or an artefact?

An apparent match can arise from five things, only the first of which is
interesting:

1. a genuine shared genomic event
2. a common germline polymorphism present in both, and in most people
3. mapping noise or a caller artefact at a recurrent site
4. low-complexity sequence that matches by composition rather than by descent
5. coincidental convergence of two short strings

The pipeline exists to separate (1) from (2)–(5), and to make the attrition at
every step explicit rather than to report a single surviving number.

### The two inputs

Nothing in `src/svneo` names a cohort, tissue or cell line. A run is defined by:

```
REFERENCE   a catalogue of peptides to look for   (any cohort, any caller)
SAMPLES     anything with an SV VCF               (+ RNA and DNA if available)
```

For the analysis reported here:

| | |
|---|---|
| Reference catalogue | `HMF_curated_2856` — 2,856 unique peptides across 379 genes |
| Catalogue SV types | 2,576 DEL, 241 DUP, 39 TRA |
| Catalogue peptide lengths | 8-mer 182, 9-mer 1,386, 10-mer 879, 11-mer 409 |
| Samples | RPE1-WT (parental), RPE1-TP53, RPE1-TP53-BRCA1, RPE1-TP53-BRCA2 |
| SV caller | PURPLE/GRIDSS2 |
| Annotation | GRCh38, Ensembl 115 — **the same release both sides** |
| Population reference | gnomAD-SV v4.1 |

The Ensembl release matching on both sides is not a detail. Peptide sequences
depend on the annotation; crossing catalogues built on different releases
compares strings that were never comparable, and no liftover is performed
anywhere in this pipeline.

---

# Part A — the method

## Overview of the cascade

| Stage | Module | Question |
|---|---|---|
| 1 | `vcf.py` | Which SV records are admissible, and what is the real lesion size? |
| 2 | `generators/` | What neopeptides could this sample's SVs produce? |
| 3–5 | `cross.py` | Which catalogue peptides match, how tightly, and how many *events* is that? |
| 6 | `confidence.py` | Is the SV call believable, and is it private to the sample? |
| 7–8 | `rna.py` | Is the locus transcribed, and do reads **cross** the junction? |
| 9 | `synthesis.py` | What survives, and which sample explains it? |
| — | `null_model.py` | How many matches would chance alone produce? |

**Every threshold lives in one module**, `src/svneo/criteria.py`, one constant
per decision, each with the reasoning beside it. A run writes the full manifest
into its own `summary.json`, so no number can outlive the thresholds that
produced it.

### The counting unit changes four times

This is the single most common way to misread the results, so it is stated
before any figure:

```
VCF records  →  breakends  →  junctions  →  peptide rows  →  genomic events
   17,953        14,340        1,143         20,710             13
```

A junction has two breakends, so pairing roughly halves the count. One junction
produces many peptides, because an 8–11mer sliding window over a fusion protein
yields up to ~40 overlapping sequences. And many peptides from one junction are
**one finding**, not many: in one observed case 21 apparent matches were a single
53 bp deletion, a 21-fold overstatement. `events` is the reported unit.

---

## Stage 1 — Admission · `src/svneo/vcf.py`

`read_breakends()` L75 · `admit()` L107 · `event_size()` L151 · `pair_junctions()` L170

### 1.1 Caller FILTER

| | |
|---|---|
| **Measures** | whether the SV caller flagged the record as passing its own checks |
| **How** | read from the VCF `FILTER` column |
| **Threshold** | somatic: keep `FILTER ∈ {PASS}` (`SOMATIC_KEEP_FILTERS`); germline: all records pass by construction |
| **Why** | records the caller rejected carry no usable evidence, and it made that judgement with more information than we have |
| **Does not capture** | nothing about whether the SV is *private* — that is stage 6 |

Two rejection classes matter. `PON` means the breakpoint was seen in a panel of
normal samples. `INFERRED` means the caller deduced the breakend from a
copy-number transition with **no read support at all** — there is no junction
sequence to translate, so such a record cannot yield a peptide even in principle.

### 1.2 Paired breakends

| | |
|---|---|
| **Measures** | whether both ends of the junction are known |
| **How** | the ALT of a mated breakend carries a bracketed partner coordinate, `N[chr5:123456[` |
| **Threshold** | ALT contains `[` or `]` |
| **Why** | a junction peptide is built by joining two sequences; one coordinate has nothing to join to |
| **Does not capture** | this is a **limitation, not a quality judgement** — a single breakend may be a perfectly real rearrangement |

### 1.3 Panel of normals

| | |
|---|---|
| **Measures** | how many unrelated normal samples show the same breakpoint |
| **How** | `PON_COUNT` INFO field (germline), or the caller's `FILTER=PON` verdict (somatic) |
| **Threshold** | `PON_MAX = 10`, i.e. keep `PON_COUNT < 10`. Absent means not in the panel, which passes (`PON_ABSENT_MEANS = 0`) |
| **Why** | a breakpoint in thousands of unrelated normals is a recurrent artefact or a common polymorphism — either way not this sample's own |
| **Does not capture** | a panel has **false negatives for inherited variation**: observed candidates with `PON_COUNT` 1 and 3 carried population frequencies of 0.287 and 0.263 |

**This is an absolute count, not a frequency**, and only means something against
the panel size. The largest count observed in this panel is ~12,000, so `< 10`
means *seen in under 0.08% of normals* — strict. The same number against a
50-sample panel would mean `< 20%`, which would be permissive. The pipeline
therefore reports `pon_fraction` alongside the raw count, and runs an unfiltered
branch so the filter's cost is measured rather than assumed.

**The panel reaches a candidate twice**, and this is where a sensitivity branch
has to be built carefully:

| Route | Germline call set | Somatic call set |
|---|---|---|
| Admission | `PON_COUNT` is an INFO field this pipeline thresholds → `pon_max: null` suffices | the caller already decided and wrote `FILTER=PON`; those records are dropped for not being `PASS` **before** `pon_max` is read → needs `admit_panel_filtered: true` |
| Privacy (stage 6) | `is_private` combines the panel with population frequency → `pon_in_privacy: false` removes the panel's vote | same |

A branch setting only `pon_max: null` is therefore a **no-op on somatic call
sets**: it produces output identical to the filtered branch under an unfiltered
label. Measured on these lines: 82, 58 and 96 records per derived line carry
`FILTER=PON` and `pon_max` could never have reached any of them.

### 1.4 Event size is type-aware

`event_size()` L151 computes the **lesion**, not the breakend span, because the
two differ:

| Type | Size |
|---|---|
| DEL | `span − 1` bases deleted |
| DUP | `span` — the duplicated tract |
| INS | the inserted length, which the span misrepresents entirely |
| BND | `None` — inter-chromosomal events have no defined size |

Using the span made 248 insertions look like 1 bp events in an earlier
implementation. `None` is propagated as `None`, never as 0: an undefined size is
not a small one.

---

## Stage 2 — Peptide generation · `src/svneo/generators/`

The only pluggable stage, isolated behind the contract in `base.py`.

Peptides are generated by [NeoSV](https://github.com/ysbioinfo/NeoSV) (Shi, Jing
& Xi, *Genome Biol* 24:169, 2023), vendored under its MIT licence. The junction
sequence is reconstructed, checked for frame, translated, and cut into 8–11mers
(`PEPTIDE_LENGTHS`). **Peptides already present in the wild-type protein of
either transcript are subtracted**, so what remains is what the SV created.

| | |
|---|---|
| **Measures** | whether an admitted junction yields any peptide at all |
| **Threshold** | at least one peptide survives wild-type subtraction |
| **Why** | most junctions fall in intergenic space, in introns without changing the protein, or reproduce the wild-type sequence — none can produce a neoantigen by definition |
| **Does not capture** | nothing about presentation; see below |

**The vendored copy carries exactly one patch**, required for correctness rather
than convenience. pyensembl changed its coding-sequence accessor at v2.3.13;
NeoSV's compensating `− 3` then shifts the 3′ side of a fusion one amino acid out
of frame. Unpatched with pyensembl 2.10.1 the tool emits 170 peptides where the
patched version emits 64, **with none in common** — the frame-shifted sequences
survive wild-type subtraction as spurious neopeptides. The bug was identified and
first corrected in [NeoSV-Trace](https://github.com/winterga/NeoSV-Trace) by
Greyson Wintergerst; this backend reproduces that fork's output exactly (64 of 64
peptides, identical `(peptide, sv_id)` keys), verified by
`tools/compare_generators.py`.

### MHC binding prediction is deliberately not applied

`RUN_MHC_PREDICTION = False`. Whether the sample's HLA would present a peptide
has no bearing on whether the peptide **exists** in it, and a sequence-identity
test is HLA-independent by construction. Filtering on predicted binding here
would discard true sequence matches and make the result depend on a licensed
predictor. Presentability is a separate question, asked later and only of
survivors — and it is not answered in this document.

`_neosv_extensions.py` adds what NeoSV does not provide and this analysis needs:
SV identifiers recovered from the VCF by coordinate, the junction offset within
the fusion protein, and a candidate table written **before** any MHC step.

---

## Stages 3–5 — Cross, QC, event collapse · `src/svneo/cross.py`

`level1_identical()` L88 · `level2_gene_svtype()` L143 · `level3_proximity()` L176 · `sequence_qc()` L290 · `to_events()` L312

### 3.1 Three levels of evidence, which are not interchangeable

| Level | Match on | What it demonstrates |
|---|---|---|
| 1 | identical peptide sequence | a shared **consequence** |
| 2 | gene + SV type | a shared **mechanism** |
| 3 | breakpoint within 1/10/50/100 kb (`PROXIMITY_WINDOWS_KB`) | a shared **locus** only |

Level 1 is the recurrence question itself, and the match is **exact** — no
alignment, no similarity score. A near-match is a different peptide that would be
presented differently or not at all; anything looser answers a different
question.

`master_sv.tsv` carries `patient_evidence` naming the strongest level each
junction reaches (`identical_peptide`, `same_gene_and_svtype`, `same_gene`,
`breakpoint_within_1kb|10kb|100kb`, `not_seen_in_patients`) alongside the
unthresholded distance in `patient_bp_dist_bp`.

### 3.2 Gene concordance

| | |
|---|---|
| **Measures** | whether the gene broken in this sample is the gene the catalogue attributes the peptide to |
| **How** | compared against **both** breakends — a junction joins two genes and the catalogue may annotate either side |
| **Threshold** | the catalogue's gene matches either breakend's gene |
| **Why** | short peptides can coincide between unrelated genomes by chance, and coincidence would not preferentially land in the same gene |
| **Does not capture** | a **genuine convergent peptide from a different gene is excluded by construction** — a deliberate trade of sensitivity for specificity that must be stated whenever "credible" counts are quoted |

Testing only the first breakend loses real matches: 8 of 64 in one run, all
translocations.

### 3.3 Sequence complexity

| | |
|---|---|
| **Measures** | whether the peptide is compositionally trivial |
| **How** | Shannon entropy, longest homopolymer run, distinct residues, single-residue fraction |
| **Threshold** | `LC_MIN_SHANNON_ENTROPY = 2.0`, `LC_MAX_HOMOPOLYMER_RUN = 4`, `LC_MIN_DISTINCT_RESIDUES = 3`, `LC_MAX_SINGLE_AA_FRACTION = 0.5` |
| **Why** | a low-complexity string matches by composition rather than by descent; `FFFFFFFFF` will match something somewhere |
| **Does not capture** | low-complexity regions can be genuinely immunogenic; this is a specificity choice |

### 3.4 Self test

Peptides present in the normal proteome are removed: a sequence the healthy
proteome already contains is not neo-anything. Without a proteome the test is
recorded as `NA`, never as passed.

### 3.5 Collapse to events

`to_events()` deduplicates on `EVENT_DEDUP_COLUMN = sample_sv_id`. This is the
most consequential reduction in the pipeline and the one most often skipped: many
peptides from one junction are one finding. A related trap is deduplicating on a
per-peptide value such as an odds ratio, which leaves window peptides separate.

### 3.6 A criterion that is deliberately OFF

`REQUIRE_SPANS_JUNCTION = False`. Whether a peptide crosses the breakpoint is
**annotated** (`spans_junction`) but does not filter. Applying it as a filter
discarded 63 of 64 catalogue matches in one run — 13 events collapsed to 1 —
because a sliding window legitimately produces peptides on either side of the
junction whose novelty comes from the reading frame, not from straddling the
breakpoint. The column is there to be used in interpretation; it is not a gate.

---

## Stage 6 — Confidence and privacy · `src/svneo/confidence.py`

`annotate_confidence()` L40 · `annotate_gnomad()` L129 · `annotate_privacy()` L61 · `_overlap_join()` L153

### 6.1 SV call confidence

| | |
|---|---|
| **Measures** | whether the SV call itself is believable |
| **Threshold** | `HC_MIN_SEGMAPQ = 30`, `HC_MIN_VF = 5`, `HC_MIN_QUAL = 20`, `HC_MIN_SV_SIZE = 100` |
| **Why** | a breakend in repetitive sequence, supported by few fragments, or below the size floor is where callers make their mistakes |
| **Applied?** | **Reported, not enforced.** No event is dropped for failing it |

`TRUST_COPY_NUMBER_FIELDS = False`: copy-number fields are excluded from the
judgement because they depend on a purity/ploidy fit that is unreliable for
clonal samples.

### 6.2 Privacy needs two independent filters

| Filter | Threshold | Rationale |
|---|---|---|
| Panel of normals | `PON_MAX = 10` | recurrent site in unrelated normals |
| Population frequency | `GNOMAD_MAX_AF = 0.001` on `GNOMAD_AF_FIELD = popmax` | inherited variation the panel misses |

An event is private only if it clears **both**. This is not redundancy: the panel
has false negatives for inherited variation, and the two observed cases that
motivated it (`PON_COUNT` 1 with AF 0.287; `PON_COUNT` 3 with AF 0.263) pass the
panel while being documented common polymorphisms.

Matching against gnomAD-SV uses **reciprocal overlap** at
`GNOMAD_RECIPROCAL_OVERLAP = 0.5`, the community default for SV matching.

**All nine ancestry groups are carried into every table** —
`gnomad_af_afr, ami, amr, asj, eas, fin, mid, nfe, sas` — plus
`gnomad_af_popmax` and the population it came from. Which group is the right
reference depends on the donor's ancestry, **which the pipeline does not know**.
`GNOMAD_AF_FIELD` selects only what the automatic verdict uses; every group's
frequency is in the table so any row can be re-judged without re-running
anything. A global average hides the spread that matters: one event here runs
from 0.591 (sas) to 0.905 (afr).

Without a gnomAD resource, every `gnomad_af_*` column is `NA` and `privacy_note`
records `"gnomAD not evaluated — PON only"`. **Not evaluated is not clean.**

---

## Stages 7–8 — Expression and junction evidence · `src/svneo/rna.py`

`expression()` L39 · `isofox_context()` L105 · `select_test()` L217 · `count_at_breakend()` L268 · `rna_tier()` L338

This is the stage on which two independent workflows disagreed, so its rules are
stated as rules. Three events called STRONG by the other workflow, over the same
BAM:

```
NTRK1   BND       coverage bp1=3, bp2=8,913,  crossing reads=0   ->  NONE
GOLGA3  DEL 2 bp  coverage 468/467,           crossing reads=0   ->  test invalid
ITGA11  DEL 220bp coverage 2,372/2,367,       crossing reads=31  ->  STRONG
```

Three STRONG, zero survivors: two had nothing crossing the junction, and the
third is a common germline polymorphism.

### Rule 1 — only junction-crossing reads tier an event

Coverage is context, not evidence. A breakpoint inside a highly expressed gene
has thousands of reads whether or not the junction exists.
`COVERAGE_TIERS_EVENTS = False`.

A soft-clip shows a read *ends* at a breakpoint, not that it crosses.
`SOFTCLIPS_TIER_EVENTS = False`. An inter-chromosomal candidate once scored "5
split reads" on clipped reads whose supplementary alignments landed 3–156 Mb from
the partner breakend — no read connected the two loci at all.
`SA_PARTNER_TOLERANCE = 1000` bp is now required of an SA tag.

### Rule 2 — never `max()` coverage across breakends

`COVERAGE_SUMMARY = "min"`. An event is in transcribed territory only if **both**
ends are. Taking the maximum turned NTRK1's silent chr1 locus into "8,913 reads"
— the count belonged to its chr8 partner.

### Rule 3 — no valid test means UNTESTABLE, never NONE

`NONE` asserts "tested and negative". The test is chosen from event geometry by
`select_test()`:

| Geometry | Test | Why |
|---|---|---|
| inter-chromosomal | `chimeric` | needs reads mapping to both loci |
| intra-chromosomal, size ≥ `MIN_TESTABLE_GAP_SIZE = 20` | `sizegap` | a CIGAR `N` gap matching the deletion |
| insertion-driven, `insert_len ≥ MIN_INSERT_LEN = 10` | `insertion` | the span is not the lesion |
| anything smaller | `none` → **UNTESTABLE** | aligners do not emit `N` skips below their minimum intron, so a 2 bp deletion *cannot* produce the signal |

Scoring an untestable event `NONE` fabricates a negative result. Gap matching
allows `NGAP_POSITION_TOLERANCE = 15` bp and `NGAP_SIZE_TOLERANCE = 10` bp.

### Rule 4 — count fragments, not alignment records

`COUNT_UNIQUE_FRAGMENTS = True`, `MIN_READ_MAPQ = 20`. One fragment can emit
several alignment records: support that read as 7 alignments was 3 fragments, one
contributing 4 records at MAPQ 3 — STAR's multi-mapping value. The MAPQ floor
matches GRIDSS2's own threshold for treating a read as unmapped.

### Evidence tiers

| Tier | Crossing fragments |
|---|---|
| STRONG | ≥ `STRONG_MIN_JUNCTION_READS = 5` (aligned with pVACfuse) |
| SUGGESTIVE | ≥ `SUGGESTIVE_MIN_JUNCTION_READS = 3` |
| WEAK | ≥ 1 |
| NONE | 0, with a valid test |
| UNTESTABLE | no valid test |
| NA | no RNA data for this sample |

**A negative RNA result is not absence of the lesion.** The gene may not be
expressed in that sample, and nonsense-mediated decay of an aberrant transcript
is a real possibility.

`EVALUATE_ALL_JUNCTIONS = True`: stages 6–8 run on every admitted junction, not
only on those that matched, so `master_sv.tsv` describes the whole call set
rather than leaving 99.3% of its rows empty.

---

## The null model · `src/svneo/null_model.py`

Raw match counts are not comparable between samples whose candidate universes
differ by orders of magnitude — 20,710 candidates against 54. Every sample
therefore reports a **per-1,000-candidate rate** and a permutation null.

`NULL_STRATEGY = "shuffle_residues"`, `NULL_PERMUTATIONS = 1000`: each candidate
peptide's residues are shuffled, preserving length and amino-acid composition,
and the catalogue intersection is recomputed.

**The null is measurably too lax, and by construction.** Shuffling residues
destroys the sequence but keeps composition, so it estimates how often a peptide
of this composition hits the catalogue by chance — not how often an unrelated
*real* genomic rearrangement would. A real-sequence null would be higher, so the
enrichment figures are upper bounds. They are reported because the direction is
informative, not because the multiplier is exact.

An enrichment ratio computed against a null mean at or below the resolution of
1,000 permutations is not quotable at all: dividing by approximately zero
produces an arbitrary number. Where that happens below, the per-1,000 rate is
given instead.

---

# Part B — results

## All four lines at a glance

Both branches for every line. `PON10` enforces the panel of normals; `noPON`
reports it without filtering on it, at both the admission and the privacy step.

| | WT<br>PON10 | WT<br>noPON | TP53<br>PON10 | TP53<br>noPON | BRCA1<br>PON10 | BRCA1<br>noPON | BRCA2<br>PON10 | BRCA2<br>noPON |
|---|---|---|---|---|---|---|---|---|
| SV records | 17,953 | 17,953 | 350 | 350 | 303 | 303 | 234 | 234 |
| Caller `FILTER=PON` admitted | — | — | — | 82 | — | 58 | — | 96 |
| Admitted breakends | 2,286 | 14,340 | 102 | 176 | 172 | 226 | 114 | 200 |
| Junctions | 1,143 | 7,170 | 51 | 88 | 86 | 113 | 57 | 100 |
| Junctions yielding peptides | 354 | 2,738 | 15 | 31 | 22 | 35 | 19 | 41 |
| Candidate peptides | 20,710 | 63,036 | 64 | 78 | 536 | 572 | 54 | 148 |
| **Distinct matches** | **64** | **288** | **0** | **0** | **21** | **21** | **0** | **0** |
| Gene-concordant | 64 | 376 | — | — | 21 | 21 | — | — |
| Credible | 57 | 331 | — | — | 21 | 21 | — | — |
| **Events** | **13** | **46** | **0** | **0** | **1** | **1** | **0** | **0** |
| Private | 9 | 24 | — | — | 0 | 0 | — | — |
| High-confidence | 4 | 15 | — | — | 0 | 0 | — | — |
| RNA-supported | 1 | 3 | — | — | 0 | 0 | — | — |
| **Private + HC + RNA** | **0** | **0** | — | — | **0** | **0** | — | — |

Matches per 1,000 candidates, with the permutation null:

| Run | candidates | matches | per 1,000 | null mean ± sd | p |
|---|---|---|---|---|---|
| RPE1-WT PON10 | 20,710 | 64 | 3.09 | 1.33 ± 0.83 | < 0.001 |
| RPE1-WT noPON | 63,036 | 288 | 4.57 | 2.74 ± 1.04 | < 0.001 |
| RPE1-TP53-BRCA1 PON10 | 536 | 21 | 39.18 | 0.00 | < 0.001 |
| RPE1-TP53-BRCA1 noPON | 572 | 21 | 36.71 | 0.00 | < 0.001 |
| RPE1-TP53 (both) | 64 / 78 | 0 | 0.00 | 0.00 | 1 |
| RPE1-TP53-BRCA2 (both) | 54 / 148 | 0 | 0.00 | 0.00 | 1 |

### Three things this table says

**Privacy is read against its branch, not across branches.** WT's private count
rises from 9 to 24 on `noPON` while the branch is *less* selective. That is not a
gain: on `noPON` privacy rests on population frequency alone, so a different
question is being answered, not the same one more loosely.

**The last two rows are different questions.** An event can be confidently called
and transcribed while being a germline polymorphism carried by most of the
population. `Private + HC + RNA` is the figure to quote; it is **0 everywhere**.

**Relaxing the panel adds junctions to the derived lines and no matches.** The
82, 58 and 96 caller-rejected records raise the junction counts by 73%, 31% and
75%, and change the distinct-match count by nothing at all: 0, 21, 0. Whatever
limits recurrence detection in the derived lines, it is not the panel filter.
This is a negative result, and it is the result the branch was run to obtain.

---

## RPE1-WT — the parental line

Germline call set, and by far the largest: 17,953 records against 234–350 for the
derived lines. That is expected — it is the whole background genome, not what was
acquired relative to a parent.

### The cascade, `PON10`

| Step | n | Removed |
|---|---|---|
| SV records | 17,953 | |
| after FILTER | 17,953 | 0 — germline sets are all-PASS by construction |
| paired breakends | 14,340 | −3,613 single breakends |
| `PON_COUNT < 10` | 2,286 | **−12,054** |
| → junctions | 1,143 | pairing halves the count |
| producing a peptide | 354 | −789 hit no coding transcript |
| candidate peptides | 20,710 | from 24,548 rows |
| identical to catalogue | 64 | |
| gene-concordant | 64 | −0 |
| high-complexity, non-self | 57 | −7 low-complexity |
| **→ events** | **13** | 57 peptides collapse to 13 loci |

The panel filter removes 84% of paired breakends. Its measured cost, recovered
from the unfiltered branch: **53 junctions that did match the catalogue**.

### The 13 events

| Gene | Type | Size | Peptides | PON | gnomAD popmax | Private | HC | RNA test | min cov | Crossing | Tier |
|---|---|---|---|---|---|---|---|---|---|---|---|
| EEF1A1 | DEL | 48 | 3 | — | — | yes | no | sizegap | 36 | 0 | NONE |
| AGMO | DUP | 124 | 3 | 3 | 0.356 (fin) | **no** | yes | insertion | 2 | 0 | NONE |
| PTPRN2 | DUP | 67 | 4 | — | 0.006 (amr) | **no** | no | insertion | 0 | 0 | NONE |
| CSMD1 | DEL | 1,256 | 10 | 3 | 0.015 (nfe) | **no** | yes | sizegap | 0 | 0 | NONE |
| NTRK1 | BND | — | 1 | — | — | yes | no | chimeric | 4 | 0 | NONE |
| MAP2K3 | BND | — | 1 | — | — | yes | no | chimeric | 170 | 0 | NONE |
| PPP1R12A | DUP | 32 | 1 | 3 | — | yes | no | sizegap | 620 | 0 | NONE |
| GOLGA3 | DEL | 2 | 1 | 3 | — | yes | no | insertion | 456 | **2** | **WEAK** |
| RPH3AL | DEL | 131 | 23 | — | — | yes | yes | insertion | 1 | 0 | NONE |
| DLGAP1 | BND | — | 7 | — | — | yes | no | chimeric | 0 | 0 | NONE |
| VSTM2B | INS | 36 | 1 | 7 | — | yes | no | insertion | 0 | 0 | NONE |
| CMSS1 | DEL | 849 | 1 | 1 | 0.460 (ami) | **no** | yes | sizegap | 27 | 0 | NONE |
| MUC4 | DUP | 3,215 | 1 | — | — | yes | no | sizegap | 2 | 0 | NONE |

**Nothing satisfies every criterion.** The three final criteria are overlapping
sets, not a chain — all three start from the 13 events and one event can fail
several:

- 4 fail privacy (AGMO, PTPRN2, CSMD1, CMSS1 — real events, but common)
- 9 fail the confidence bar
- 12 have no read crossing the junction

Note **CMSS1**: `PON_COUNT = 1` — as clean as the panel can report — but a
population frequency of 0.460 in `ami`. This is the case the two-filter privacy
rule exists for; the panel alone would have passed it.

Note also that a high `min_coverage` is not evidence. PPP1R12A sits under 620
reads and MAP2K3 under 170, with **zero** crossing the junction. Both are in
transcribed territory; neither junction is transcribed.

### The `noPON` branch

46 events instead of 13, 288 distinct matches instead of 64. Three carry RNA
support:

| Gene | Crossing reads | Coverage | PON | gnomAD popmax | Private | HC | Tier |
|---|---|---|---|---|---|---|---|
| **ITGA11** | 31 | 2,412 / 2,433 | 3,513 | **0.905 (afr)** | no | yes | STRONG |
| GOLGA3 | 2 | 456 | 3 | absent | yes | no | WEAK |
| ROBO1 | 1 | 73 | 86 | 0.609 (mid) | no | no | WEAK |

**ITGA11 is the clearest illustration in the dataset of why RNA support alone
cannot promote a candidate.** The junction is real and strongly transcribed — 31
crossing fragments over ~2,400× coverage, a high-confidence call. It is also
present in **every** gnomAD population sampled, from 0.591 (sas) to 0.905 (afr),
with `PON_COUNT` 3,513. A common germline deletion polymorphism, not a somatic
neoantigen.

The only strongly transcribed junction in this lineage is a polymorphism. That is
the finding, and it survives both branches.

---

## RPE1-TP53

Somatic call set relative to RPE1-WT. 350 records, of which the caller rejected
160 as `INFERRED` and 82 as `PON`.

| Step | PON10 | noPON |
|---|---|---|
| SV records | 350 | 350 |
| Caller `FILTER=PON` admitted | — | 82 |
| after FILTER | 108 | 190 |
| paired | 102 | 176 |
| junctions | 51 | 88 |
| producing a peptide | 15 | 31 |
| candidate peptides | 64 | 78 |
| **matches** | **0** | **0** |

**Nothing reaches the RNA stage, and that is not a null RNA result.** The chain
stops at the catalogue cross: no candidate peptide matched, so stages 6–8 had
nothing to test. Every RNA column for this line is `NA`, **not zero** — the
distinction matters, because zero would assert "tested and not transcribed".

The reason is the candidate universe, not the filters. 51 junctions produce 64
candidate peptides — against a catalogue of 2,856. Admitting the 82 panel-flagged
records raises this to 78 candidates and still yields nothing. With a universe
this small, the permutation null also produces zero matches, so `p = 1` here means
"indistinguishable from chance", not "significantly absent".

---

## RPE1-TP53-BRCA1

The only derived line with a catalogue match.

| Step | PON10 | noPON |
|---|---|---|
| SV records | 303 | 303 |
| Caller `FILTER=PON` admitted | — | 58 |
| after FILTER | 176 | 234 |
| paired | 172 | 226 |
| junctions | 86 | 113 |
| producing a peptide | 22 | 35 |
| candidate peptides | 536 | 572 |
| **distinct matches** | **21** | **21** |
| credible | 21 | 21 |
| **events** | **1** | **1** |

### The single event

| Gene | Type | Size | Peptides | PON | gnomAD popmax | Private | HC | RNA test | min cov | Crossing | Tier |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SLC9A9 | DEL | 52 | 21 | — | **0.222 (asj)** | **no** | no | sizegap | 3 | 0 | NONE |

**This event is not what its raw numbers suggest, in three separate ways.**

*The 21 matches are one finding.* All 21 peptides come from a single 52 bp
deletion, seen through the sliding window. Quoting 21 would overstate the
evidence 21-fold — this is the observed case that motivated event-level counting.

*It is line-specific within the experiment and still not private.* SLC9A9 is
absent from the other lines' call sets, which is what "specific to this line"
means experimentally. But its population frequency is 0.222 in `asj`, and absent
from the panel of normals entirely (`PON_COUNT` empty). Specific to the line and
common in the population are not mutually exclusive, and only the gnomAD filter
catches it. **The panel alone would have passed this event.**

*The per-1,000 rate is high and the enrichment ratio is meaningless.* 39.18
matches per 1,000 candidates is more than ten times RPE1-WT's 3.09. But the
permutation null mean is 0.000 — at or below the resolution of 1,000
permutations — so any enrichment ratio computed from it is a division by
approximately zero. The rate is the comparable figure; the multiplier is not.

The relaxed branch changes the denominators and nothing else: same 21 matches,
same single event, per-1,000 rate 36.71 instead of 39.18.

---

## RPE1-TP53-BRCA2

| Step | PON10 | noPON |
|---|---|---|
| SV records | 234 | 234 |
| Caller `FILTER=PON` admitted | — | 96 |
| after FILTER | 121 | 217 |
| paired | 114 | 200 |
| junctions | 57 | 100 |
| producing a peptide | 19 | 41 |
| candidate peptides | 54 | 148 |
| **matches** | **0** | **0** |

As with RPE1-TP53, **the RNA stage never ran**; its columns are `NA`, not zero.

This line shows the panel relaxation most clearly: 96 of 234 records carry
`FILTER=PON`, so admitting them nearly doubles the junction count (57 → 100) and
almost triples the candidate peptides (54 → 148). The catalogue match count stays
at zero. The panel filter was not hiding recurrent candidates here.

---

# What this analysis does not answer

- **Presentation.** Read evidence supports "the junction is transcribed". It says
  nothing about whether the peptide is processed, loaded and displayed. That
  requires immunopeptidome mass spectrometry. HLA typing (LILAC) has not been run
  on these lines, so even predicted presentability is not available.
- **Convergent peptides from a different gene.** Excluded by construction through
  gene concordance — a deliberate specificity trade.
- **Single breakends.** Dropped at admission because a junction peptide needs two
  coordinates. Some are real rearrangements.
- **Whether a negative RNA result means the lesion is absent.** It does not.
- **An unexplained strand skew.** The catalogue is 61.0% minus-strand
  (n = 379, p = 1.6 × 10⁻⁶ against a length-weighted background), while RPE1-WT's
  candidates are 43.7% (n = 174). The 5′/3′ assignment has been verified
  strand-aware — pyensembl returns exons 5′→3′, confirmed on four minus-strand
  genes — so the cause is elsewhere and unidentified. **No result here should be
  read as strand-controlled until this is resolved.** See
  [`OPEN_QUESTIONS.md`](OPEN_QUESTIONS.md).

# Interpreting the overall outcome

Reducing a four-figure candidate count to a single-digit believable set is the
**expected** outcome, not a failed analysis: a 28-team benchmark found roughly 6%
of top-ranked neoantigen predictions validate functionally. What matters is that
every order of magnitude lost is attributable to a stated criterion, which is
what the funnel and the per-filter reports provide.

The headline result across all four lines and both branches is that **no
candidate is simultaneously private, confidently called and transcribed**. The
two events that come closest fail for opposite reasons: ITGA11 has excellent read
support and is a common polymorphism; SLC9A9 is line-specific within the
experiment and is also common in the population. Neither is a recurrent
tumour-derived neoantigen.

That the catalogue peptides are found at 48–105× the composition-matched null in
RPE1-WT says the overlap is not random string coincidence. What the per-event
analysis then says is that the overlap is dominated by shared **germline**
variation — which is exactly what a recurrence test between unrelated genomes
should be expected to surface first, and exactly what the privacy filters exist
to remove.

---

# Reproducing this

The inputs are not public: the reference catalogue is access-controlled and the
sample VCFs and BAMs belong to the group that generated them. Cloning the
repository reproduces the **method**, not these figures.

```bash
python tests/test_integration.py     # the whole chain on a synthetic fixture
```

A pass means the installation computes what the reference one computes. The
fixture deliberately contains the cases that have broken this pipeline: an
insertion-driven event whose span misrepresents its size, an inter-chromosomal
junction with no defined size, a panel-of-normals record that must be rejected,
and a rejected breakend sharing a coordinate with an accepted one.

Every figure in Part B comes from `summary.json`, `credible_events.tsv` and
`stage8_rna_evidence.tsv` in the corresponding `results/<line>_<branch>/`
directory, each of which carries the complete threshold manifest that produced
it. Regenerate the per-run reports with:

```bash
python tools/build_report.py --all         # markdown, per run
python tools/build_html_report.py --all    # filter-by-filter validation pages
python tools/build_summary_html.py         # this document, as a shareable page
```
