# svneo

Recurrence testing for structural-variant-derived neoantigens.

Given a catalogue of neopeptides and a set of samples with SV calls, `svneo`
answers one question: **do those peptides reappear in those samples, and is the
recurrence real or an artefact?** An apparent match can arise from a genuine
shared genomic event, or from mapping noise, low-complexity sequence,
coincidental peptide convergence, or a common germline polymorphism. The pipeline
is built to tell them apart and to make the attrition at every step explicit.

## Contents

- [How it works](#how-it-works)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Configuration](#configuration)
- [Peptide generation (stage 2)](#peptide-generation-stage-2)
- [Outputs](#outputs)
- [Design rationale](#design-rationale)
- [Limitations](#limitations)
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

## Installation

```bash
conda create -n svneo python=3.11 -y && conda activate svneo
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
PYTHONPATH=src python -m svneo.run --config config/my_run.yaml --dry-run
PYTHONPATH=src python -m svneo.run --config config/my_run.yaml
```

`PYTHONPATH=src` is required — the package is not installed, it is run from the
source tree. Use the interpreter of the environment created above; a system
Python without `pyyaml` fails at config load with `No module named 'yaml'`.

After a run, assemble the auditable candidate table:

```bash
python tools/build_master_table.py --results-dir results
```

Options:

| Flag | Effect |
|---|---|
| `--samples NAME [NAME…]` | restrict to named samples |
| `--branches NAME [NAME…]` | restrict to named sensitivity branches |
| `--out-dir PATH` | override the config's output directory |
| `--dry-run` | print the execution plan and exit |

One sample failing does not abort the run; its traceback is printed and the
remaining samples complete.

## Configuration

`config/template.yaml` documents every field. The essentials:

```yaml
run_name: my_run
peptide_generator: neosv              # stage-2 backend

reference:
  peptides: /path/to/reference_peptides.tsv
  peptide_column: neoantigen          # required
  gene_column: gene1                  # null disables gene concordance
  proteome: /path/to/proteome.fa.gz   # null records the self test as NA

samples:
  - name: sample_A
    vcf: /path/to/sample_A.sv.vcf.gz
    vcf_kind: somatic                 # germline | somatic
    parent: null
    rna_bam: ...                      # optional
    isofox_dir: ...                   # optional
```

`vcf_kind` is a property of the sample, not of the filename: `germline` means
"the background genome shared with descendants" and receives the
population-frequency filters, `somatic` means "acquired relative to the parent"
and is assumed already panel-filtered by the caller.

Sensitivity branches differ **only** in stage-1 admission, so any difference
between branches is attributable to admission and nothing else.

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
| `funnel.tsv` | the 12-row evidence funnel |
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
- **An unexplained strand skew.** Candidate genes deviate significantly from a
  length-weighted background at three levels of the analysis, in two directions.
  The 5′/3′ assignment has been verified strand-aware, so the cause is elsewhere
  and unidentified. No result here should be read as strand-controlled until it
  is resolved — see [`docs/OPEN_QUESTIONS.md`](docs/OPEN_QUESTIONS.md).

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
| `tools/build_master_table.py` | two tables of **values, not verdicts**: `master_sv.tsv` (one row per admitted junction) and `master_peptides.tsv` (one row per matched peptide) |
| `tools/build_report.py` | a technical results report per sample, generated from its outputs |
| `tools/check_docs.py` | verifies the documentation still describes the code; exits 1 if stale |
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

If this pipeline contributes to published work, cite this repository and the
peptide generator it vendors (NeoSV-Trace, see `vendor/VENDOR.md`) alongside the
reference catalogue used.
