# WORKLOG

Chronological log of work sessions. Append-only; newest entry first.
For the repository's current state see `README.md`, and for the user-facing
history of behavioural changes see `CHANGELOG.md`. This file is the narrative:
what was tried, why, and what it showed.

---

## 2026-08-14 — reporting the panel instead of filtering on it

### Logic followed

The request was one folder per cell line where the panel of normals is reported
but does not remove candidates. The branch appeared to exist already
(`noPON`, `pon_max: null`), so this looked mechanical. It was not.

**The panel reaches a candidate by two different routes, and the existing branch
relaxed neither of them on the derived lines.**

1. *Admission.* On a germline call set the panel is an INFO field, `PON_COUNT`,
   thresholded by this pipeline — so `pon_max=None` is enough. On a somatic set
   PURPLE has already applied its own panel filter and written the verdict into
   `FILTER=PON`; admission drops those records for not being `PASS`, one step
   before `pon_max` is ever consulted. Measured: 82, 58 and 96 records per
   derived line that `pon_max` could never reach.
2. *Privacy.* `is_private` combines the panel count with the population
   frequency at stage 6. Relaxing admission alone still lets it remove exactly
   the candidates the branch was run to see.

So running the branch as configured would have produced three folders labelled
"noPON" holding output byte-identical to the filtered branch. A no-op is worse
than a missing result, because it looks like an answer.

Also found while checking the above: `docs/EXECUTIVE_SUMMARY.md` asserted "the
panel branch does not apply to [the derived lines]", which was the same mistake
written down as a conclusion.

### Changes made

- `criteria.SOMATIC_KEEP_FILTERS_REPORT_PON` — the somatic FILTER set for a
  panel-reporting branch, carrying the reasoning above. `INFERRED` stays
  excluded: it is not a panel judgement, and those breakends have no read
  support and hence no junction sequence to translate.
- `vcf.admit(..., admit_panel_filtered=)` — admits caller-rejected panel
  records and records `caller_pon_admitted` in the funnel, so the relaxation is
  visible in the output rather than inferred from a total.
- `confidence.annotate_privacy(..., pon_in_privacy=)` — the panel count is
  still computed and reported (`pass_pon`, `pon_count`, `pon_fraction`) but does
  not vote. `privacy_note` states which basis produced the verdict.
- `config.Branch` gains `admit_panel_filtered` and `pon_in_privacy`. Its
  docstring no longer claims branches differ only in admission, which is no
  longer true.
- `counts["events_private_hc_and_rna"]` — the all-criteria count.
  `events_hc_and_rna` does not include privacy; the HTML report was describing
  it as "satisfies every criterion", which overstated ITGA11 (`PON_COUNT`
  3,513, gnomAD popmax 0.905) as a surviving candidate. Both reports now
  quote the two figures separately and say why they must not be conflated.
- `tools/build_html_report.py` — a "Panel of normals — reported, not applied"
  block that takes the place the filter card would have occupied. Without it
  the cascade merely lacks a panel step, which reads as "the panel removed
  nothing" — the opposite conclusion from the same numbers.
- Three tests pinning the above (`tests/test_svneo.py`), including one that
  fails if `pon_max=None` ever starts admitting somatic `FILTER=PON` records on
  its own, since that would make the two branches indistinguishable.

### Results

Both branches were re-run for all four lines. `PON10` reproduced its previous
figures exactly (2,286 admitted, 1,143 junctions, 64 matches, 13 events,
9 private, 4 HC), which is the regression check on the refactor.

| | WT PON10 | WT noPON | TP53 PON10 | TP53 noPON | BRCA1 PON10 | BRCA1 noPON | BRCA2 PON10 | BRCA2 noPON |
|---|---|---|---|---|---|---|---|---|
| Caller `FILTER=PON` admitted | — | — | — | 82 | — | 58 | — | 96 |
| Admitted | 2,286 | 14,340 | 102 | 176 | 172 | 226 | 114 | 200 |
| Junctions | 1,143 | 7,170 | 51 | 88 | 86 | 113 | 57 | 100 |
| Distinct matches | 64 | 288 | 0 | 0 | 21 | 21 | 0 | 0 |
| Events | 13 | 46 | 0 | 0 | 1 | 1 | 0 | 0 |
| Private + HC + RNA | 0 | 0 | — | — | 0 | 0 | — | — |

**Relaxing the panel adds junctions to the derived lines and no matches.** The
junction count rises 73%, 31% and 75%, and the distinct-match count does not
move: 0, 21, 0 — identical to the filtered branch. Whatever limits recurrence
detection in the derived lines, it is not the panel filter. This is the finding
the branch was run to obtain, and it is a negative one.

WT's private count rises from 9 to 24 between branches, which is an artefact of
the question changing rather than a gain: on `noPON` privacy rests on population
frequency alone. The three RNA-supported events are unchanged — GOLGA3 (2
fragments, absent from gnomAD, but not high-confidence), ITGA11 (31 reads,
popmax 0.905) and ROBO1 (1 read, popmax 0.609) — so the earlier conclusion
holds: with the panel reported rather than enforced, the only strongly
transcribed junction in this lineage is a common germline polymorphism.

### Open

- Unchanged: strand skew, LILAC HLA typing, patch 001 upstream.

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
