# WORKLOG

Chronological log of work sessions. Append-only; newest entry first.
For the repository's current state see `README.md`, and for the user-facing
history of behavioural changes see `CHANGELOG.md`. This file is the narrative:
what was tried, why, and what it showed.

---

## 2026-09-08 — phase 2: the patched run, measured against the old one

### Logic followed

The whole analysis was re-run with patch 002 into `results_patch002/`, leaving
`results/` untouched. Replacing it would have answered the wrong question: the
useful quantity is not "what do we get now" but "how much was the defect
distorting", and that is the evidence an upstream report needs.

`tools/compare_runs.py` reports the difference on the quantities the defect was
expected to move. A figure that does NOT move is as informative as one that does.

### Results

**The defect's signature is gone, and only on the strand it affected.**
`Start-loss` falls from 84.0% to 10.1% on minus-strand fusions while the plus
strand moves 13.8% to 13.4%. A fix that repaired both equally would have meant
the diagnosis was wrong.

**The strand skew resolves in everything the pipeline generates.** RPE1-WT
candidates: 41.5% minus (p = 5e-4) to 46.8% (p = 0.36), against a 48.6%
background. Matched candidates: 74.4% (p = 7e-4) to 60.9% (p = 0.30). The only
set still deviating is the reference catalogue itself — the patient side, which
cannot be corrected from here.

**Roughly half the catalogue matches were artefactual, and none were being
missed.** RPE1-WT: 288 distinct matches to 130. RPE1-TP53-BRCA1: 21 to **zero**.
Of the 288, 130 survive, 158 vanish, and **0 are new**. The defect was
manufacturing matches, not hiding them.

**The headline conclusion is unchanged.** No candidate satisfies every criterion,
in any line, on either branch — 0 before, 0 after.

**SLC9A9 no longer exists.** RPE1-TP53-BRCA1's only event, featured in the
executive summary and in reporting to collaborators as a line-specific deletion
that was nonetheless a common polymorphism, was an artefact of the defect.

### What did not resolve

Matched peptides still almost never span their junction: 1 of 161 (0.6%) against
18.2% across the candidate universe, where before it was 1 of 387 (0.3%) against
16.5%. The ~30-fold depletion survives the patch and is therefore **not**
explained by this defect.

Most likely the reference catalogue carries the artefact still, being built with
the same tool family, so our surviving matches are the ones resembling its
non-spanning peptides. Settling it requires regenerating the catalogue with the
patched tool, which is outside this repository. Recorded as an open question.

### Changes made

- `results_patch002/` with `PROVENANCE.md` (revision, diagnostic state, what to
  expect) and `COMPARISON.md` / `.html`.
- `tools/compare_runs.py`.
- `.gitignore` excluded `results/` exactly, so a sibling output directory would
  have been committed — both carry access-controlled, patient-derived data.
  Now `results_*/` too, verified with `git check-ignore`.
- `docs/OPEN_QUESTIONS.md` question 1 marked RESOLVED and FIXED, with the
  measured effect.

### Consequence to decide

Everything published so far — the executive summary in three versions, the master
tables, the per-line reports, and what was sent to collaborators — describes
`results/`. Those figures are now known to include artefactual matches. Whether
to republish from `results_patch002/`, or to publish the comparison as the
result, is a decision about the record rather than a technical one.

## 2026-09-08 — patch 002: coding exons now reach consumers in reading order

### Logic followed

Phase 1 of the plan agreed after the root cause was found. Two things had to
happen before writing any code: correcting the record, and building a check that
could actually detect the defect.

**The recorded scope was wrong.** The previous entry said the exonic branch was
correct on both strands. It is not. That claim came from checking the returned
sequence began with `ATG` — a test that passes at any non-zero length, since the
sequence is always taken from position 0 of the coding sequence. Starting
correctly is not measuring correctly, and no probe had ever compared lengths.

### Changes made

**`tools/diagnose_minus_strand_cds.py` rewritten** to sweep a breakpoint through
every region of a transcript — each coding exon, each intron — and compare the
length returned against `expected_head()`, written from the biology rather than
from the implementation. Before the patch:

| | regions | correct |
|---|---|---|
| EGFR, PTEN, PIK3CA (plus) | 87 | 87 |
| ITGA11, TP53, BRCA1, GOLGA3, KRAS, BRAF (minus) | 188 | **0** |

Not a single minus-strand region was right. Two failure modes: an intronic
breakpoint returned an empty head, an exonic one the wrong length — too short
near the start of a transcript, **too long** near its end. The second direction
is the damaging one, since it adds sequence the gene does not contribute and the
sliding window turns it into peptides.

**`vendor/patches/002-minus-strand-cds-order`** — `get_cds_range` now sorts into
reading order (descending coordinate on the minus strand) and `get_noncds_range`
consumes it instead of reading the raw accessor. Both edits carry a `PATCHED`
comment naming the patch, following the convention of 001.

Deliberately not patched: `truncate_cds`. Its minus-strand branches are already
written for reading order and become correct once they receive it. Fixing them
individually would add a second exception on top of the first and leave the
documented invariant still false.

**`tests/test_svneo.py::test_minus_strand_cds_is_read_in_transcript_order`** —
asserts the source rather than the behaviour, so it runs without an annotation
cache. Verified in both directions: passes patched, fails on a reverted file.
Its accessor count strips comments, since the patch comments name the accessor
they replaced and counting those would make the assertion depend on prose.

### Results

| | regions correct |
|---|---|
| before | 87 / 275 (32%) |
| after | **275 / 275 (100%)** |

**Zero regression on the plus strand**: 90 probe positions compared between the
original and patched module, none changed. 189 of 194 minus-strand positions
changed; the five that did not are positions where the wrong answer happened to
coincide with the right one.

Documentation corrected throughout — `OPEN_QUESTIONS.md`, `README.md`,
`EXECUTIVE_SUMMARY.md` and `VENDOR.md` all described a single patch and an
intron-only defect.

### Not done

**Phase 2: re-running.** `results/` still holds figures produced with the defect
present, and `RPE1-WT_noPON` is from 2026-09-07 while the other seven runs are
from 2026-08-14 — a mixed state either way. The agreed approach is to run the
patched pipeline into a separate directory and compare, so the defect's impact
is measured rather than silently corrected.

What the comparison should answer: how the `frame_effect` distribution moves
(84.0% `Start-loss` on the minus strand today), whether junction-spanning
peptides recover (1 of 408 today, 67 expected), whether the strand skew
disappears — which would confirm the diagnosis — and whether any candidate now
survives every criterion.

Anticipate that apparent overlap with the catalogue may *fall*: it was built with
the same tool family and probably carries the same artefact, so correcting one
side removes a shared signature. That would be informative, not a regression.

## 2026-09-07 — the strand skew is a bug in the vendored generator

### Logic followed

Three items from the improvement list were meant to be separate investigations:
fill the empty `frameshift` column, test the strand skew on admitted breakends,
and round-trip a minus-strand fusion by hand. They turned out to be one finding.

The trigger was a question asked while reviewing the report — *how can junctions
match while the peptides differ?* — which led to checking whether the matched
peptides actually cross the breakpoint. They do not: 1 of 408, where the 16.5%
base rate predicts 67 (p = 1e-30).

### What was found, in order

**Step 1 — the skew is not in what breaks.** `check_strand_bias.py` gained
`--breakends`, which looks genes up from the admitted breakend coordinates using
pyensembl directly, touching no fusion logic. RPE1-WT's admitted breakends are
47.6% minus (p = 0.33) against a 48.6% length-weighted background; the candidate
peptides they produce are 41.5% (p = 5e-4). The deviation is introduced at
annotation. This was step 1 of the plan recorded in OPEN_QUESTIONS.md and it
eliminated the leading hypothesis.

**Step 2 — the frame column had never been populated.** The pipeline read
NeoSV's verdict as `frameshift`; upstream calls it `frame_effect`. `getattr`
returned the default for all 72,351 rows, so the column shipped empty and never
errored. With it populated: **84.0% of minus-strand fusions are `Start-loss`
against 13.8% of plus-strand ones.** `Start-loss` is upstream's own reliability
warning, not a biological class.

**Step 3 — the defect.** `truncate_cds(transcript, '5', pos)` returns an empty
5' coding sequence for a minus-strand transcript when the breakpoint falls in an
intron. Six of six minus-strand genes fail; three of three plus-strand genes
pass; both strands are correct for exonic breakpoints.

The cause: coding exon ranges arrive sorted by coordinate, which is reading order
only on the plus strand. `get_noncds_range` has a minus-strand branch that
reverses its arithmetic — but that branch assumes the list already arrives
reversed. Correcting data that was never reversed produces gap intervals with
`start > end`, which can never match, and a first interval that swallows the
whole gene (7,662,015-7,676,594 on TP53: 14.5 kb of a 14.8 kb transcript). An
intronic breakpoint therefore resolves to gap index 0, meaning "before the first
exon", so zero exons are kept and the 5' piece is empty.

### Why it survived this long

It fails silently — an inverted range raises nothing. It fails only on introns,
so a test built from exonic breakpoints passes; the first probe written here did
exactly that and reported "correct" before being redone. And the docstring of
`get_cds_range` asserts the very invariant that is violated, which is also why
an earlier entry in OPEN_QUESTIONS.md recorded exon ordering as *ruled out* —
that check was run against `transcript.exons`, not the accessor this code path
reads.

### Changes made

- `tools/diagnose_minus_strand_cds.py` — reproduces the defect, exits 1 while
  present.
- `tools/check_strand_bias.py` — `--breakends`; and it silently stopped emitting
  the matched-candidate row when `master_table.tsv` was renamed, leaving a stale
  figure in circulation. It now warns instead of skipping.
- `_neosv_extensions.py` — reads `frame_effect`; the retired `frameshift` name is
  documented as always-empty for anyone holding an older table.
- `docs/OPEN_QUESTIONS.md` — question 1 rewritten: root cause, plain-language
  mechanism with a worked three-exon example, the evidence, and the proposed fix.
- README and executive summary no longer describe the skew as unexplained.

### Not done, deliberately

**The fix.** Every figure in the repository was produced with the defect present,
so applying it changes published numbers — a decision about results, not a
refactor. `results/RPE1-WT_noPON/` was regenerated on 2026-09-07 to measure
`frame_effect`; the other seven runs are from 2026-08-14, so `results/` is in a
mixed state until that decision is taken.

Recommended: patch and run BOTH versions, so the impact of the defect is
measured rather than merely corrected — which is also the evidence upstream
would need.

## 2026-08-14 — quantifying locus-level vs peptide-level recurrence

### Logic followed

A direct question — *how can SV junctions be similar between patients and cell
lines while the neopeptides differ, and does that happen in our data?* — turned
out to have a measurable answer that the report was asserting rather than
showing. The three cross levels were described as "different kinds of evidence,
not degrees of one", which is the right claim, but nothing in the document said
by how much they diverge or why.

### Changes made

- `tools/analyse_locus_vs_peptide.py` — computes the gap and separates its
  causes. Written as a script rather than an ad-hoc query because its figures
  are quoted in a shared report.
- `tools/sync_event_tables.py` gains `<!-- analysis:NAME -->` regions, so a
  generated analysis table lives in the markdown and cannot drift from the run.
- `docs/EXECUTIVE_SUMMARY.md`: a new section, *Shared locus is not shared
  peptide*, with the tables generated by that script, plus a cross-reference
  from the three-levels table in Part A.

### Results

Across the four lines (panel reported, not enforced): **4,806 junctions reach a
locus- or mechanism-level match against 54 reaching an identical peptide, 89:1.**

Of the 512 junctions within 1 kb of a patient breakpoint, the causes separate
cleanly and are not interchangeable:

| | n |
|---|---|
| produce no peptide at all — could not match by construction | 354 |
| produce peptides, none matches — same place, different sequence | 135 |
| produce peptides and at least one matches | 23 |

The strongest subset is the 105 junctions at the *exact* coordinate of a patient
breakpoint, 92 of them also of the same SV type. Nineteen share a peptide.

**Inserted bases separate the two groups.** Among same-coordinate junctions that
produce peptides, median `insert_len` is 6 bp where no peptide matches and 2 bp
where one does; 73% of the same-coordinate set carries an insertion at all. Two
breaks at one position with different insertions build different fusion
sequences, and an insertion that is not a multiple of three shifts the frame, so
everything downstream is an unrelated sequence rather than a similar one.
Observed at distance 0 with no match: ANO9 (310 inserted bases), PLAT (301),
THSD7B (95), CNTNAP2 (63). Two further mechanisms appear in the same subset —
the SV type differing at one locus (ANO9, PLAT are duplications where the patient
event is a deletion) and the lesion size differing while the breakpoint
coincides.

**A matching junction matches few of its own peptides**: median fraction 0.040,
minimum 0.003 (PPP1R12A, 1 of 302). Expected from a sliding window, but it means
locus-level overlap almost never implies peptide-level overlap — which is the
justification for making level 1 exact sequence identity.

### Also fixed

- `build_summary_html.py` anchored only `h2`/`h3`, so a cross-reference to a
  top-level section did not resolve. Now `h1`-`h3`, with `h1` as a top-level
  contents entry; the title is excluded from the contents rather than listed as
  its first item.
- Emphasis inside a code span (`` `**x**` ``) renders the asterisks literally.
  Fixed in the generator, where the markdown is built.

## 2026-08-14 — the executive summary, rewritten in full

### Logic followed

`docs/EXECUTIVE_SUMMARY.md` gave one section to RPE1-WT and a single table row to
each of the other three lines. The request was the complete version: the method
end to end with every criterion and its justification, then the results developed
per cell line. Rewritten rather than added alongside — its old §3 was already
"stage by stage, with the code", so a second document would have duplicated the
thresholds and the two would have diverged. That failure mode is not theoretical
here: this same document asserted "the panel branch does not apply to the derived
lines", which was false.

Structure: **Part A** the method, **Part B** the results. Stated at the top that
Part B cannot be read without Part A, because a count means nothing without the
criterion that produced it.

### Changes made

- `docs/EXECUTIVE_SUMMARY.md` rewritten, 742 lines. Part A covers each stage with,
  per criterion: what it measures, how, the threshold, why that threshold, and
  **what it does not capture** — including the criteria that are deliberately OFF
  (`REQUIRE_SPANS_JUNCTION`, `SOFTCLIPS_TIER_EVENTS`, `COVERAGE_TIERS_EVENTS`,
  `RUN_MHC_PREDICTION`), since a withdrawn criterion that goes unmentioned reads
  as applied. Part B gives each line its cascade in both branches, its events
  with genes named, and what the line does and does not show.
- `tools/build_summary_html.py` — renders a `docs/` markdown file to a
  self-contained page. One source, two outputs: the markdown stays canonical and
  the page is generated, never edited.

### Results

`results/reports/EXECUTIVE_SUMMARY.html`, 56.9 KB, 38 anchored headings, 24
tables each in its own scroll container, no external resources, light and dark
themes.

Every figure was verified against the run outputs rather than transcribed from
memory. Three passes:

1. every integer in `summary.json` for all eight runs appears in the document — 0
   discrepancies
2. every `CONSTANT = value` quoted matches `criteria.py` — 0 real discrepancies
   (9 false positives from VCF fields such as `FILTER=PON`)
3. the 13-event RPE1-WT table checked row by row against `credible_events.tsv`
   joined to `stage8_rna_evidence.tsv`: gene order, SV type, privacy, confidence,
   RNA test, coverage, crossing reads and tier — all match

`check_docs.py` initially found ten stale line references, all of them line
numbers I had written from memory. Corrected with `--fix-line-numbers`.

The document's conclusions, unchanged by the rewrite but now stated with their
grounds: no candidate in any line or branch is simultaneously private,
confidently called and transcribed. The two that come closest fail for opposite
reasons — ITGA11 has 31 crossing fragments and a popmax of 0.905; SLC9A9 is
line-specific within the experiment and carries 0.222 in `asj` while being absent
from the panel entirely. Neither is a recurrent tumour-derived neoantigen, and
the second would have passed a panel-only privacy filter.

### Open

- Unchanged: strand skew, LILAC HLA typing, patch 001 upstream.

## 2026-08-14 — README audit for publication as SVneoHunter

### Logic followed

The repository is going public as `SVneoHunter`, so the README was audited
against the code rather than read through: every documented command executed,
every config key checked against the dataclasses, every documented output
checked against a real run, every link resolved. Prose claims were checked the
same way — a README that asserts something the code no longer does is worse than
one that omits it, because the reader has no reason to doubt it.

### Changes made

**Claims that were false.**

- *"Sensitivity branches differ **only** in stage-1 admission"* — untrue since
  the panel-reporting change earlier today, which added a stage-6 lever. Replaced
  with a section explaining that the panel reaches a candidate twice, why both
  levers are needed, and why `pon_max: null` alone is a no-op on somatic sets.
- *"the 12-row evidence funnel"* — it was 13, and is now 14.
- The citation section asked readers to cite "the peptide generator it vendors
  (NeoSV-Trace)". The repository vendors **NeoSV**; NeoSV-Trace is credited for
  identifying the frame bug and is explicitly not redistributed. The licence
  section said this correctly, so the two contradicted each other.

**Gaps a new user would hit.**

- No clone step and no repository name: the title said `svneo`, which is the
  Python package, not the repository. Both are now stated, with the distinction
  made explicit so `python -m svneo.run` inside `SVneoHunter/` is not surprising.
- No repository layout. Added, with a pointer to `docs/EXECUTIVE_SUMMARY.md`,
  which was not linked from the README at all.
- The `resources:` block was undocumented in full. `gnomad_sv` matters most:
  without it privacy silently degrades to the panel alone, every `gnomad_af_*`
  column is `NA`, and `privacy_note` records why — a real loss of specificity
  that a reader had no way to anticipate. Each key now has an "effect if absent"
  row, plus where to download gnomAD-SV.
- Config keys absent from the README: `isofox_prefix`, `dna_bam`, `notes`,
  `svtype_column`, `pool`, `pool_peptide_column`, `gnomad_max_af`, the `genome:`
  and `criteria:` blocks, and the two new branch keys. All documented; `pool` in
  particular gates cross level 3 entirely.
- Two tools were missing from the review table (`build_html_report.py`,
  `build_summary_html.py`).
- No example of what a run produces. Added a "What a run looks like" section with
  real console output and the funnel, including the warning that the last three
  funnel rows are overlapping sets rather than a chain.

**A usability bug found by running the documented commands.**

`cp config/template.yaml … && --dry-run` — the first thing anyone does — raised
`FileNotFoundError` behind a stack trace on the first placeholder path it met.
That reads as a broken tool rather than an unconfigured one, and revealed one bad
path per invocation. Now a `ConfigError` reports every unreadable path at once
with a line saying what to do, `run.py` prints the message without a traceback,
and the exit code is 1. Pinned by a test that loads the shipped template and
asserts more than one placeholder is named.

### Changes to the code, not only the docs

- `config.ConfigError`, and `load()` collecting all missing paths.
- `run.py` catching it and exiting with the message alone.
- `synthesis.FUNNEL_ROWS`: `events_hc_and_rna` was labelled "both — the
  strongest set". It is not — privacy is absent from it. Relabelled, and
  `events_private_hc_and_rna` added as the final row.

### Results

`results/` was regenerated from the current code so no output predates the
thresholds that produced it. Both branches ran in one invocation, which also
fixed the cross-sample aggregates: they are rewritten per invocation, so the
earlier split runs had left `null_model_rates.tsv` describing only the last
branch. It now carries all eight runs.

Verified after the rewrite: every tool, doc, config key and documented output
appears in the README; all eight relative links resolve; every table-of-contents
anchor matches a heading; `check_docs.py` clean; both test suites pass (31 unit
tests, 21 integration checks).

Enrichment over the permutation null, all eight runs:

| Run | candidates | matches | per 1,000 | null mean | enrichment |
|---|---|---|---|---|---|
| RPE1-WT PON10 | 20,710 | 64 | 3.09 | 1.30 | 49× |
| RPE1-WT noPON | 63,036 | 288 | 4.57 | 2.74 | 105× |
| RPE1-TP53-BRCA1 PON10 | 536 | 21 | 39.18 | 0.001 | — |
| RPE1-TP53-BRCA1 noPON | 572 | 21 | 36.71 | 0.000 | — |

The BRCA1 enrichment figure is not quotable: its null mean is at or below the
resolution of 1,000 permutations, so the ratio is an artefact of dividing by
approximately zero. The per-1,000 rate is the comparable number.

### Open

- Unchanged: strand skew, LILAC HLA typing, patch 001 upstream.
- `docs/EXECUTIVE_SUMMARY.md` is being rewritten as the full method-plus-results
  document; `tools/build_summary_html.py` renders it to a shareable page.

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
