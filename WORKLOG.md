# WORKLOG

Chronological log of work sessions. Append-only; newest entry first.
For the repository's current state see `README.md`, and for the user-facing
history of behavioural changes see `CHANGELOG.md`. This file is the narrative:
what was tried, why, and what it showed.

---

## 2026-08-14

### Logic followed

The five runs had produced numbers, but a filter-by-filter count is not the same
as a validated result: a reader cannot check `-12,054` without knowing what the
filter measures, why that criterion, and why that threshold. The goal of this
session was a per-cell-line HTML report where every count sits beside its own
justification, prototyped on RPE1-WT before extrapolating.

Two problems surfaced while building it, both about numbers that are correct but
misread:

1. **The counting unit changes three times** along the cascade (breakends →
   junctions → peptide rows → genomic events). Read as one sequence, `2,286`
   followed by `1,143` looks like an arithmetic error rather than a change of
   grain. Silent notes were not enough — the transitions are now explicit
   markers that name both units and explain the collapse.

2. **A count cannot be audited without names.** This was found by a direct
   question — *why is ITGA11 no longer a final candidate?* — that the report as
   built could not answer. ITGA11 is removed by the panel-of-normals filter,
   which acts BEFORE annotation, so the filtered run holds no record of it at
   all: the gene's absence was indistinguishable from a bug. The unfiltered
   branch exists precisely to make that cost measurable, so the report now reads
   it back and names what the threshold discarded.

### Changes made

- `tools/build_html_report.py`
  - `FILTERS`: eleven filters, each carrying `measures` / `how` / `why` /
    `threshold` / `threshold_why` next to the code that counts it, so a
    threshold cannot drift away from its justification.
  - `grain_change()`: explicit markers at the three unit transitions.
  - `unfiltered_losses()`: recovers, from the parallel `_noPON` run, the
    catalogue-matching junctions the panel filter removed — with PON_COUNT,
    gnomAD popmax and its population, and RNA support.
  - `named_losses()` / `losses()`: name the rows each later criterion removed,
    with the deciding value, for gene concordance, complexity, privacy,
    confidence and RNA.
  - RNA card distinguishes **"locus not transcribed"** (measured zero coverage)
    from **"expressed, but nothing crosses the junction"** — the distinction the
    NTRK1/GOLGA3 disagreement turned on.
  - `_text()`: single helper for absent values. `value or ""` does NOT work
    here — pandas yields `float('nan')`, NaN is truthy, and the literal string
    `nan` reaches the page. Fourth occurrence of this bug in this codebase.
  - Two CSS tokens (`--accent`, `--rule-faint`) were referenced but never
    defined; they would have failed silently. Now checked mechanically that
    every `var(--x)` used is defined.
- `docs/EXECUTIVE_SUMMARY.md`: ITGA11's frequency was quoted as `0.79` without
  naming a population — that is the global `gnomad_af`, while the pipeline
  decides on `gnomad_af_popmax` (`0.905`, afr). Now quotes popmax with its
  population and the across-population range.

### Results

Five reports in `results/reports/`. The RPE1-WT cascade:

| filter | in → out |
|---|---|
| Caller FILTER | 17,953 → 17,953 |
| Paired breakends | 17,953 → 14,340 |
| Panel of normals | 14,340 → 2,286 |
| *grain: breakends → junctions* | 2,286 → 1,143 |
| Produces a candidate peptide | 1,143 → 173 |
| Identical to catalogue peptide | 173 → 14 |
| *grain: junctions → peptide rows* | 14 → 64 |
| Gene concordance | 64 → 64 |
| Sequence complexity | 64 → 57 |
| Not in normal proteome | 57 → 57 |
| *grain: peptide rows → events* | 57 → 13 |
| Private / Confidence / RNA | 13 → 9 / 4 / 1 (overlapping sets, not a chain) |

All three verdict branches are exercised: RPE1-WT and RPE1-TP53-BRCA1 reach
"no candidate satisfies every criterion"; RPE1-TP53 and RPE1-TP53-BRCA2 reach
"nothing reaches the RNA stage" (0 catalogue matches, so RNA is `NA`, not zero);
RPE1-WT_noPON reaches "1 candidate satisfies every criterion".

**The panel filter's measured cost: 53 junctions that DID match the patient
catalogue.** The most notable is ITGA11 — a real, transcribed junction (31
crossing reads, coverage 2,412/2,433, tier STRONG) that is nonetheless present
in 59.1% (sas) to 90.5% (afr) of every gnomAD population, with PON_COUNT 3,513.
A common germline deletion polymorphism, not a somatic neoantigen. This is the
clearest illustration in the dataset of why RNA support alone cannot promote a
candidate, and why privacy is checked independently of expression.

### Open

- Repository still has **zero commits** (60 files staged). Nothing is versioned.
- Strand skew unresolved — see `docs/OPEN_QUESTIONS.md`.
- LILAC HLA typing not run; presentability layer therefore not evaluated.
- Patch 001 not yet submitted upstream to `ysbioinfo/NeoSV`.
