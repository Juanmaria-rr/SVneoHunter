# SVneoHunter

Recurrence testing for structural-variant-derived neoantigens.

> The repository is `SVneoHunter`; the Python package inside it is `svneo`.
> Every command below is therefore `python -m svneo.run …`, run from the
> repository root.

Given a catalogue of neopeptides and a set of samples with SV calls, `svneo`
answers one question: **do those peptides reappear in those samples, and is the
recurrence real or an artefact?** An apparent match can arise from a genuine
shared genomic event, or from mapping noise, low-complexity sequence,
coincidental peptide convergence, or a common germline polymorphism. The pipeline
is built to tell them apart and to make the attrition at every step explicit.

## Contents

- [How it works](#how-it-works)
- [Repository layout](#repository-layout)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Configuration](#configuration)
- [Peptide generation (stage 2)](#peptide-generation-stage-2)
- [Outputs](#outputs)
- [What a run looks like](#what-a-run-looks-like)
- [Design rationale](#design-rationale)
- [Limitations](#limitations)
- [Reviewing a run](#reviewing-a-run)
- [Testing](#testing)
- [Licence](#licence)
- [Citation](#citation)

## How it works

The analysis is defined by two inputs, and nothing in `src/svneo` names a
specific cohort, tissue or cell line:

```
REFERENCE   a catalogue of peptides to look for   (any cohort, any caller)
SAMPLES     anything with an SV VCF               (+ RNA and DNA if available)
```

| Stage | Module | Question |
|---|---|---|
| 1 | `vcf.py` | Which SV records are admissible, and what is the real lesion size? |
| **2** | **`generators/`** | What neopeptides could this sample's SVs produce? *(pluggable)* |
| 3–5 | `cross.py` | Which reference peptides match, how tightly, and how many *events* is that? |
| 6 | `confidence.py` | Is the SV call believable, and is it private to the sample? |
| 7–8 | `rna.py` | Is the locus transcribed, and do reads **cross** the junction? |
| 9 | `synthesis.py` | What survives, and which sample explains it? |
| — | `null_model.py` | How many matches would chance alone produce? |

Samples form a lineage via `parent:`, which controls execution order (parents
first) and attribution (the earliest sample carrying a breakpoint explains it).
For unrelated samples, leave `parent: null` and the lineage collapses to a flat
list.

A sample missing a modality is recorded as `NA` for the affected stages and
continues. Missing data is never scored as a negative result.

## Repository layout

```
src/svneo/            the pipeline; nothing here names a cohort or sample
    criteria.py       EVERY threshold, one constant per decision, each with why
    config.py         the config schema (Sample, Branch, Reference, RunConfig)
    vcf.py            stage 1   admission, breakend pairing, event size
    generators/       stage 2   peptide generation (the one pluggable stage)
    cross.py          stages 3-5  catalogue cross, sequence QC, event collapse
    confidence.py     stage 6   SV confidence, PON and gnomAD privacy
    rna.py            stages 7-8  expression, junction-crossing read evidence
    synthesis.py      stage 9   funnel, cross-sample comparison, attribution
    null_model.py     how many matches chance alone would produce
    run.py            the orchestrator; entry point

vendor/neosv/         NeoSV (MIT), vendored; one documented patch
config/               template.yaml + a complete worked example
tools/                reporting and verification scripts (see Reviewing a run)
tests/                unit tests + an integration test on a synthetic fixture
docs/                 EXECUTIVE_SUMMARY.md, LOCUS_VS_PEPTIDE.md,
                      PROVENANCE.md, OPEN_QUESTIONS.md
notebooks/            an independent PySpark re-derivation of every stage
results/              run outputs (git-ignored; regenerated, not tracked)
```

Start with [`docs/EXECUTIVE_SUMMARY.md`](docs/EXECUTIVE_SUMMARY.md): the method
end to end, every criterion with its threshold and justification, and the results
per sample. `README.md` (this file) is how to run it; `WORKLOG.md` is the
chronological record of how it got here.

## Installation

```bash
git clone https://github.com/Juanmaria-rr/SVneoHunter.git
cd SVneoHunter
```

```bash
conda create -n svneo python=3.10 -y && conda activate svneo
pip install -r requirements.txt
```

The peptide generator ships with the repository (see
[Peptide generation](#peptide-generation-stage-2)), so no separate install is
needed. One reference download is required:

```bash
export PYENSEMBL_CACHE_DIR=$PWD/pyensembl_cache
pyensembl install --release 115 --species homo_sapiens     # ~3.2 GB, once
```

### Version constraints that are not optional

| Package | Constraint | Why |
|---|---|---|
| **pyensembl** | **> 2.3.13** | The vendored generator carries a patch that is correct only for releases *after* the coding-sequence accessor changed. On ≤ 2.3.13 the patch must be reverted, or every fusion is one amino acid out of frame. See `vendor/patches/001-*`. |
| Python | ≥ 3.10 | union types in annotations |
| PySpark | 3.5.x on **Java 11** | review notebook only; PySpark 4.x requires Java 17 |

Reference environment that produced the published results: Python 3.10.20,
pandas 2.3.3, pysam 0.24.0, **pyensembl 2.10.1**, PySpark 3.5.3, OpenJDK 11.0.26.
Every run records its own versions in `summary.json`.

> **The Ensembl release used to generate peptides must match the one the
> reference catalogue was built with.** Peptide sequences depend on the
> annotation; crossing catalogues built on different releases compares strings
> that were never comparable. Both sides here are Ensembl 115.

## Quick start

```bash
cp config/template.yaml config/my_run.yaml     # then edit the paths
export PYENSEMBL_CACHE_DIR=/path/to/pyensembl_cache

# 1. validate the config without computing anything
PYTHONPATH=src python -m svneo.run --config config/my_run.yaml --dry-run

# 2. run it
PYTHONPATH=src python -m svneo.run --config config/my_run.yaml

# 3. turn the outputs into tables and reports
python tools/build_master_table.py --results-dir results

# restrict the combined table to the runs a given report covers
python tools/build_master_table.py --runs SAMPLE_branch --suffix _branch
python tools/build_report.py --all          # markdown, one per sample and branch
python tools/build_html_report.py --all     # filter-by-filter validation pages
```

`PYTHONPATH=src` is required — the package is not installed, it is run from the
source tree. Use the interpreter of the environment created above; a system
Python without `pyyaml` fails at config load with `No module named 'yaml'`.

`--dry-run` checks the configuration and prints the execution plan without
reading a single BAM. Run it first: it reports **every** unreadable path at once,
so an unedited template produces a list to fix rather than one error per attempt.

```
config error in config/my_run.yaml:
these configured paths do not exist:
  reference.peptides: /path/to/reference_peptides.tsv
  samples[sample_A].vcf: /path/to/sample_A.sv.vcf.gz

If this is a freshly copied config/template.yaml, replace the
/path/to/... placeholders with real paths.
```

Options:

| Flag | Effect |
|---|---|
| `--config PATH` | the run configuration (required) |
| `--samples NAME [NAME…]` | restrict to named samples |
| `--branches NAME [NAME…]` | restrict to named sensitivity branches |
| `--out-dir PATH` | override the config's output directory |
| `--dry-run` | validate the config, print the plan, and exit |

A name that matches nothing fails loudly and lists what is available: selecting
zero samples and exiting successfully is indistinguishable from a run that
worked, and sample names use hyphens while the directories they came from often
use underscores.

One sample failing does not abort the run; its traceback is printed and the
remaining samples complete. Cross-sample outputs are rebuilt from whatever
succeeded.

> **Run every branch in one invocation.** The cross-sample files
> (`comparison.tsv`, `null_model_rates.tsv`, `attribution.tsv`) are rewritten on
> each run and describe only the samples and branches of *that* invocation. Two
> separate runs leave aggregates covering only the second.

## Configuration

`config/template.yaml` documents every field. The essentials:

```yaml
run_name: my_run
peptide_generator: neosv              # stage-2 backend
out_dir: results

genome:
  build: GRCh38                       # recorded, not enforced: no liftover
  ensembl_release: 115                # MUST match the catalogue's release

criteria:                             # optional per-run threshold overrides
  PON_MAX: 10                         # any constant from src/svneo/criteria.py

reference:
  name: my_catalogue                  # appears in the run manifest
  peptides: /path/to/reference_peptides.tsv
  peptide_column: neoantigen          # required: the peptide sequence
  gene_column: gene1                  # null disables gene concordance
  svtype_column: svtype               # null disables cross level 2
  proteome: /path/to/proteome.fa.gz   # null records the self test as NA
  pool: /path/to/catalogue_with_coords.tsv   # optional; enables cross level 3
  pool_peptide_column: neoantigen

samples:
  - name: sample_A
    vcf: /path/to/sample_A.sv.vcf.gz  # required; everything else is optional
    vcf_kind: somatic                 # germline | somatic
    parent: null                      # name of another sample, or null
    rna_bam: /path/to/sample_A.bam    # stages 7-8; NA without it
    isofox_dir: /path/to/isofox/A     # transcript-level expression context
    isofox_prefix: sample_A.isf       # filename prefix inside isofox_dir
    dna_bam: /path/to/sample_A.dna.bam  # carried for follow-up; no stage reads it yet
    notes: |                          # free text, echoed into summary.json
      Anything a reader needs to interpret this sample.
```

A sample missing a modality gets `NA` for the affected stages and continues, so
the minimum viable sample is a name and an SV VCF.

### External resources

```yaml
resources:
  pyensembl_cache: /path/to/pyensembl_cache        # required by stage 2
  gnomad_sv: /path/to/gnomad.v4.1.sv.sites.vcf.gz  # strongly recommended
  gnomad_helper: /path/to/annotator.py             # optional
  neosv_path: /path/to/an/external/neosv           # optional; vendored by default
  peptides_from:                                   # optional; skip stage 2
    sample_A: /path/to/previous/run/sample_A
```

| Key | Effect if absent |
|---|---|
| `pyensembl_cache` | stage 2 falls back to pyensembl's default location; set it explicitly so runs are reproducible across machines |
| `gnomad_sv` | **privacy degrades silently to the panel of normals alone.** Every `gnomad_af_*` column is written as `NA` and `privacy_note` records `"gnomAD not evaluated — PON only"`. A panel has false negatives for inherited variation, so this is a real loss of specificity, not a cosmetic one — see design rationale 4 |
| `gnomad_helper` | the reciprocal-overlap join runs directly; the helper only exists to reuse a site's existing annotator and its local-copy discovery |
| `neosv_path` | the vendored copy under `vendor/neosv/` is used, which is the intended default |
| `peptides_from` | stage 2 runs normally; point it at a previous run's output prefix to reproduce a run without regenerating peptides |

Download gnomAD-SV v4.1 sites (~1.7 GB) from
[gnomad.broadinstitute.org](https://gnomad.broadinstitute.org/downloads#v4-structural-variants).
It must be on the same genome build as the samples; no liftover is performed.

The three cross levels answer progressively weaker questions, and the config
controls which are available:

| Level | Needs | Question |
|---|---|---|
| 1 identical peptide | `peptide_column` | a shared **consequence** |
| 2 gene + SV type | `gene_column`, `svtype_column` | a shared **mechanism** |
| 3 breakpoint proximity | `pool` (a catalogue **with coordinates**) | a shared **locus** only |

`pool` is separate from `peptides` because the two differ in what they carry: the
curated peptide set may have no coordinates, while the larger pool has
breakpoints but is not peptide-filtered. Without `pool`, level 3 and the
`patient_evidence` proximity tiers are simply absent, not zero.

`vcf_kind` is a property of the sample, not of the filename: `germline` means
"the background genome shared with descendants" and receives the
population-frequency filters, `somatic` means "acquired relative to the parent"
and is assumed already panel-filtered by the caller.

### Sensitivity branches

A branch re-runs the same samples under a different threshold, so that any
difference in the result is attributable to that threshold and nothing else.
Branches differ in stage-1 admission and — for the panel of normals only — in
whether the panel count votes on the privacy verdict at stage 6.

```yaml
branches:
  - name: PON10
    pon_max: 10                 # keep records with PON_COUNT < 10
  - name: noPON
    pon_max: null               # this pipeline's own PON threshold: off
    admit_panel_filtered: true  # also admit somatic FILTER=PON records
    pon_in_privacy: false       # PON reported, is_private = gnomAD only
```

A branch may also set `gnomad_max_af:` to override the population-frequency
threshold for that branch alone; `null` disables the population filter. Any
threshold in `criteria.py` can likewise be overridden per run from the config's
`criteria:` block, and every override is echoed into the run manifest.

**The panel of normals reaches a candidate twice, and a branch meaning "report
the panel, do not filter on it" has to relax both:**

| Where | Germline call set | Somatic call set |
|---|---|---|
| Admission | `PON_COUNT` is an INFO field this pipeline thresholds → `pon_max: null` is enough | the caller already decided and wrote `FILTER=PON`; those records are dropped for not being `PASS` **before** `pon_max` is consulted → needs `admit_panel_filtered: true` |
| Privacy (stage 6) | `is_private` combines the panel count with population frequency → `pon_in_privacy: false` removes the panel's vote | same |

Setting only `pon_max: null` is therefore a **no-op on somatic call sets**: it
produces output identical to the filtered branch under an unfiltered label. The
panel count is measured and carried into every table in both branches; what
changes is whether it removes anything.

Thresholds live in `src/svneo/criteria.py`, one constant per decision, each with
its rationale. They can be overridden per run from the config's `criteria:`
block; every override is echoed into the run manifest, so a result can never be
traced to a threshold that is not recorded beside it.

`config/rpe1_hmf.example.yaml` is a complete worked example with every field
filled in and paths as placeholders. Local run configurations (`config/*.yaml`)
are git-ignored, since they carry absolute paths and the filenames of
access-controlled inputs.

## Peptide generation (stage 2)

Stage 2 is the only point where the pipeline depends on external
peptide-generation code, and it is isolated behind a documented contract:

```
src/svneo/generators/
    base.py           the contract every backend must satisfy
    precomputed.py    reuse a stored output, e.g. to reproduce a run
    __init__.py       backend registry
    neosv.py          the NeoSV backend — the only file that imports it
    _neosv_extensions.py   what NeoSV does not provide (this repo's code)
vendor/
    neosv/            vendored NeoSV
    LICENSE.NeoSV     MIT, Copyright (c) 2022 Yang Shi
    patches/          every modification to upstream, documented
    VENDOR.md         provenance, licence, patch policy
```

Selection is one config line: `peptide_generator: neosv`.

Peptides are generated by [NeoSV](https://github.com/ysbioinfo/NeoSV) (Shi, Jing
& Xi, *Genome Biol* 24:169, 2023), vendored under its MIT licence so the
repository runs without an external checkout or network access.

**The vendored copy carries exactly one patch**, and it is required for
correctness rather than convenience. pyensembl changed its coding-sequence
accessor at v2.3.13; NeoSV's compensating `- 3` then shifts the 3' side of a
fusion one amino acid out of frame. Unpatched with pyensembl 2.10.1 the tool
emits 170 peptides where the patched version emits 64, **with none in common** —
the frame-shifted sequences survive the wild-type subtraction as spurious
neopeptides. See [`vendor/patches/001-pyensembl-stop-codon-frame.md`](vendor/patches/001-pyensembl-stop-codon-frame.md).
The issue was identified and first corrected in
[NeoSV-Trace](https://github.com/winterga/NeoSV-Trace) by Greyson Wintergerst.

NeoSV writes only post-netMHCpan binders and carries no SV identifier, so three
things are added in `_neosv_extensions.py` — this repository's own code: SV
identifiers recovered from the VCF by coordinate, the junction offset within the
fusion protein, and a candidate-peptide table written before any MHC step. MHC
binding prediction is not part of a sequence-identity test, which is
HLA-independent by construction.

Peptide sequences themselves are upstream's, untouched. Equivalence is verified
rather than assumed: this backend reproduces NeoSV-Trace's output exactly — 64 of
64 peptides shared, identical `(peptide, sv_id)` keys, `spans_junction` agreeing
on every row. Reproduce with `tools/compare_generators.py`.

To add a different generator, write a sibling of `neosv.py`, satisfy the
contract in `base.py` (enforced by `validate_output`), and register it in
`REGISTRY`. No other module changes.

## Outputs

Per sample and branch, under `<out_dir>/<sample>[_<branch>]/`:

| File | Content |
|---|---|
| `stage1_junctions.tsv` | admitted junctions with type-aware `event_size` and `insert_len` |
| `<sample>.admitted.vcf` | the admitted call set handed to stage 2 |
| `stage3_matches.tsv` | peptide matches with concordance and QC flags |
| `stage3_gene_svtype.tsv`, `stage3_proximity.tsv` | cross levels 2 and 3 |
| `credible_events.tsv` | **the reported unit** — one row per genomic event |
| `stage7_expression.tsv`, `stage7_background.json` | expression with its background gradient |
| `stage8_rna_evidence.tsv` | per-breakend counts, test applied, evidence tier |
| `stage6_8_all_junctions.tsv` | the same evidence for **every** admitted junction, not only those that matched (`criteria.EVALUATE_ALL_JUNCTIONS`) |
| `stage3_junction_recurrence.tsv` | per junction: distance to the nearest patient breakpoint, and what was found there |
| `master_sv.tsv`, `master_peptides.tsv` | the auditable tables — values, not verdicts |
| `funnel.tsv` | the 14-row evidence funnel, fixed row order so samples are comparable line by line |
| `summary.json` | counts, null model, generator provenance, full threshold manifest |

Cross-sample, under `<out_dir>/`: `comparison.tsv`, `null_model_rates.tsv`,
`attribution.tsv`, the combined `master_sv.tsv` / `master_peptides.tsv` with their
column dictionaries, and `reports/<sample>.md`.

**Was a junction seen in patients?** `master_sv.tsv` carries `patient_evidence`,
naming the strongest level of recurrence each junction reaches — `identical_peptide`,
`same_gene_and_svtype`, `same_gene`, `breakpoint_within_1kb|10kb|100kb`, or
`not_seen_in_patients` — alongside the unthresholded distance in
`patient_bp_dist_bp`. The three levels are not the same kind of evidence: a shared
peptide is a shared consequence, a shared gene and SV type a shared mechanism, and
a nearby breakpoint only a shared locus.

**Population frequency** is carried for all nine gnomAD ancestry groups
(`gnomad_af_afr` … `gnomad_af_sas`) plus `gnomad_af_popmax`. Which one is the
right reference depends on the donor's ancestry, which the pipeline does not know;
`criteria.GNOMAD_AF_FIELD` selects what the automatic verdict uses, and every
group's frequency is in the table so it can be re-judged without re-running.

Read `comparison.tsv` on events and per-1,000-candidate rates, not on peptide
counts: a single locus seen through a sliding window can produce dozens of
matches.

## What a run looks like

Console output is one block per sample and branch, with the attrition printed as
it happens:

```
=== RPE1-WT_PON10  [germline] ===
  stage1: 17,953 records -> 2,286 admitted -> 1,143 junctions
  stage2 [neosv]: 24,548 peptide rows -> 20,710 unique
  stages3-5: 64 matches -> 57 credible -> 13 events
  stage6: 4 HC, 9 private (PON+population)
  stage8: tiers {'NONE': 12, 'WEAK': 1}
  null: 64 matches from 20,710 candidates (3.09/1,000); null 1.3 ± 0.8,
        p < 0.001, 49x over null
```

`funnel.tsv` is the same cascade as a table, with a fixed row order so two
samples can be read side by side:

| step | description | n |
|---|---|---|
| `vcf_records` | SV records in the raw VCF | 17,953 |
| `admitted` | after FILTER / pairing / population filters | 2,286 |
| `junctions` | paired breakends collapsed to junctions | 1,143 |
| `candidate_peptides` | unique candidate neopeptides | 20,710 |
| `matches_identical` | identical to a reference peptide (rows: peptide × SV) | 64 |
| `matches_unique_peptides` | distinct peptides among those matches | 64 |
| `matches_gene_concordant` | and broken in the same gene | 64 |
| `matches_credible` | and high-complexity, non-self | 57 |
| `events` | credible matches collapsed to genomic events | **13** |
| `events_private` | PON and population-frequency clean | 9 |
| `events_hc` | high-confidence SV calls | 4 |
| `events_rna_supported` | junction-crossing reads in RNA | 1 |
| `events_hc_and_rna` | high-confidence AND transcribed (privacy NOT applied) | 0 |
| `events_private_hc_and_rna` | private AND high-confidence AND transcribed | **0** |

Read the last three rows carefully: they are **overlapping sets, not a chain**.
All three start from the 13 events, and one event can fail more than one
criterion. `events_hc_and_rna` deliberately excludes privacy, so it is not the
headline figure — an event can be confidently called and transcribed while being
a germline polymorphism carried by most of the population. Quote
`events_private_hc_and_rna`.

Note also that the counting unit changes along the cascade: breakends →
junctions → peptide rows → genomic events. A peptide count read as a finding
count overstates the evidence by up to ~40×, which is why `events` is the
reported unit.

The numbers above come from access-controlled inputs and are shown to make the
output shape concrete; they are not reproducible from a clone. See
[Reproducing results](#reproducing-results).

## Design rationale

Five rules distinguish this implementation from a naive one. Each addresses a
failure mode that produces confident, wrong results, illustrated below with
observed cases.

**1. Only junction-crossing reads count as RNA evidence.** Coverage is context, and
a soft-clip shows that a read *ends* at a breakpoint, not that it crosses. An
inter-chromosomal candidate scored "5 split reads" on clipped reads whose
supplementary alignments landed 3–156 Mb from the partner breakend; no read
connected the two loci at all.

**2. Coverage is summarised across breakends with `min`, never `max`.** An event is in
transcribed territory only if both ends are. Taking the maximum reported a silent
locus carrying 3 reads as having 8,913 — the count belonged to its highly
expressed partner.

**3. An event with no valid test is `UNTESTABLE`, not `NONE`.** `NONE` asserts
"tested and negative". A 2 bp deletion cannot produce a CIGAR `N` gap, because
aligners do not emit skips shorter than their minimum intron; scoring it `NONE`
fabricates a result. The test is selected from event geometry: chimeric evidence
for inter-chromosomal junctions, gap matching for intra-chromosomal events above
the floor, and inserted-sequence matching for insertion-driven events, which the
breakend span misrepresents.

**4. Privacy requires two independent filters.** A panel of normals has false
negatives for inherited variation: candidates with `PON_COUNT` of 1 and 3 carried
population allele frequencies of 0.287 and 0.263. An event must clear both the
panel and a population-frequency threshold. A well-called, well-expressed,
read-supported deletion present in thousands of unrelated normals is real — and
disqualified, because the property that matters is not whether the call is
correct but whether it is the sample's own.

**5. Counting is done on events, never peptide rows.** One SV yields up to ~40
overlapping sliding-window peptides. In one case 21 apparent matches were a
single 53 bp deletion, a 21-fold overstatement. A related trap: deduplicating on
a per-peptide value such as an odds ratio leaves window peptides separate.

Two further points shape interpretation. Raw match counts are not comparable
between samples whose candidate universes differ by orders of magnitude, so every
sample reports a per-1,000-candidate rate and a permutation null alongside the
raw count. And a funnel that reduces a four-figure candidate count to a
single-digit set is the expected outcome, not a failure: a 28-team benchmark
found ~6% of top-ranked neoantigen predictions validate functionally. What
matters is where each order of magnitude is lost.

Full derivation, including which observation set each threshold, is in
[`docs/PROVENANCE.md`](docs/PROVENANCE.md).

## Limitations

- **Presentation is out of scope.** Read evidence supports "the junction is
  transcribed". Claiming presentation requires immunopeptidome mass spectrometry.
- **Gene concordance is part of the credibility definition**, so a genuine
  convergent peptide arising from a different gene is excluded by construction.
  This trades sensitivity for specificity and should be stated whenever
  "credible" counts are quoted.
- **A negative RNA result is not absence of the lesion.** The gene may not be
  expressed in that sample, and nonsense-mediated decay of an aberrant transcript
  is a real possibility.
- **Copy-number fields are excluded from credibility judgements** by default,
  because they depend on a caller's purity/ploidy fit that is unreliable for
  clonal samples.
- **No liftover is performed.** The reference catalogue and every sample VCF must
  be on the same genome build and annotation release.
- **A known defect in the vendored generator, not yet fixed.** For a
  minus-strand transcript whose breakpoint falls in an **intron** — most real
  breakpoints — the 5' coding segment is returned empty. The fusion becomes the
  3' partner alone, is flagged `Start-loss`, places the junction at residue ~0,
  and can never yield a junction-spanning peptide. Measured here: 84.0% of
  minus-strand fusions are `Start-loss` against 13.8% of plus-strand ones, and 1
  of 408 catalogue-matching peptides spans its junction where 67 would be
  expected. **Every figure in this repository was produced with this defect
  present.** Root cause, evidence and the proposed fix are in
  [`docs/OPEN_QUESTIONS.md`](docs/OPEN_QUESTIONS.md); reproduce with
  `python tools/diagnose_minus_strand_cds.py`.

## Reproducing results

**The inputs of the published analysis are not public.** The reference catalogue
is access-controlled, and the sample VCFs and BAMs belong to the group that
generated them. Cloning this repository therefore does not reproduce our figures;
it reproduces the *method*. Concretely:

| What | Availability |
|---|---|
| This code, thresholds, and the vendored generator | in the repository |
| gnomAD-SV v4.1 sites VCF (1.7 GB) | public download |
| Ensembl reference via pyensembl | public download |
| Reference peptide catalogue | access-controlled (cohort data agreement) |
| Sample SV calls, RNA and DNA alignments | belong to the originating group |

To confirm that an installation computes what ours does, run the integration
test. It executes the whole chain on a small synthetic dataset shipped in
`tests/data/` and asserts exact values:

```bash
python tests/test_integration.py
```

A pass means the environment reproduces the reference behaviour. A failure means
it does not — find out why before running real data. The fixture deliberately
contains the cases that have broken this pipeline: an insertion-driven event
whose span misrepresents its size, an inter-chromosomal junction with no defined
size, a panel-of-normals record that must be rejected, and a rejected breakend
sharing a coordinate with an accepted one.

To run the analysis on your own data, copy `config/template.yaml`, point it at a
peptide catalogue and a set of samples, and follow [Quick start](#quick-start).

## Reviewing a run

| Tool | Question it answers |
|---|---|
| `tools/build_master_table.py` | two tables of **values, not verdicts**: `master_sv.tsv` (one row per admitted junction) and `master_peptides.tsv` (one row per matched peptide), each with a column dictionary. `--runs`/`--suffix` restrict the combined table to the runs a given report covers |
| `tools/build_report.py` | a technical results report per sample, in markdown, generated from its outputs |
| `tools/build_html_report.py` | the same run as a **filter-by-filter validation page**: each filter with what it measures, how, why, its threshold and the threshold's justification, beside the count it removed and the names of what it removed. `--all` for every run |
| `tools/build_summary_html.py` | renders a `docs/` markdown document to a self-contained, shareable HTML page. The markdown stays canonical; the page is generated, never edited |
| `tools/build_focused_summary.py` | derives a single-branch version of a document — drops branch-comparison blocks, removes the other branch's table columns, and can drop whole sections by heading. Avoids a second hand-written document that would diverge from the first |
| `tools/sync_event_tables.py` | fills the events tables in `docs/` from the run outputs, so a results table in prose cannot drift from the run it describes. `--check` gates a commit |
| `tools/check_docs.py` | verifies the documentation still describes the code; exits 1 if stale |
| `tools/column_meanings.py` | one written explanation per master-table column, emitted into `<table>_column_dictionary.tsv`. A column with no entry is marked `UNDOCUMENTED` in the artefact and fails the test suite, so a new column cannot ship unexplained |
| `tools/analyse_locus_vs_peptide.py` | measures the gap between locus-level and peptide-level recurrence, and separates its causes. Answers "how can junctions match while peptides do not?" with numbers |
| `tools/check_strand_bias.py` | strand composition against a length-weighted background (see open questions) |
| `tools/inspect_insertion.py` | the individual reads behind an insertion call, with mapping quality, duplicate flag and in-read offset, against a background rate — is this support real or a duplicate stack? |
| `tools/compare_generators.py` | do two peptide generators produce interchangeable output? Run before swapping one. |
| `notebooks/pipeline_walkthrough_pyspark.ipynb` | every stage re-derived independently in PySpark and asserted against the pipeline. Change `SAMPLE` in the first cell to walk through any line. |

The notebook is generated from `notebooks/build_walkthrough.py`; edit that and
re-run it rather than editing the `.ipynb`.

## Testing

```bash
python tests/test_svneo.py        # unit tests: the rules that must not regress
python tests/test_integration.py  # the whole chain on the synthetic fixture
```

The suite pins the decisions above rather than pursuing coverage: each test
guards a rule whose violation produces a plausible but wrong result.

## Licence

**MIT** — see [`LICENSE`](LICENSE).

The vendored dependency, NeoSV, is also MIT (Copyright (c) 2022 Yang Shi); its
licence ships as [`vendor/LICENSE.NeoSV`](vendor/LICENSE.NeoSV) and applies to
`vendor/neosv/`. Modifications to it are limited to one documented patch under
[`vendor/patches/`](vendor/patches/).

The reading-frame bug that patch fixes was identified in
[NeoSV-Trace](https://github.com/winterga/NeoSV-Trace) by Greyson Wintergerst;
no code from that fork is redistributed here.

## Citation

If this pipeline contributes to published work, cite this repository, the
peptide generator it vendors — **NeoSV** (Shi, Jing & Xi, *Genome Biol* 24:169,
2023), see [`vendor/VENDOR.md`](vendor/VENDOR.md) — and the reference peptide
catalogue used. No code from NeoSV-Trace is redistributed here; that fork is
credited for identifying the reading-frame bug, not vendored.
