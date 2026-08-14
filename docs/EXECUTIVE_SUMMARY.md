# Executive summary — finding and validating SV-derived neopeptides

A worked example on **RPE1-WT**, the parental near-diploid cell line, tested
against a curated catalogue of **2,856** patient-derived neopeptides.

Every number below is traceable: each section names the module, the function and
the threshold constant that produced it, and the output file it lands in.
Constants live in `src/svneo/criteria.py`; line numbers are given for
orientation and will drift, the names will not.

---

## 1. The question

Given a catalogue of neopeptides predicted from patient tumours and a sample
with its own structural-variant calls, does any of those peptides **reappear**
in the sample, and is the recurrence a shared genomic event or an artefact?

An apparent match can arise from mapping noise, low-complexity sequence,
coincidental convergence of short peptides, or a common germline polymorphism.
Separating those is the entire pipeline.

---

## 2. Result on RPE1-WT

Branch `PON10` (breakpoints seen in fewer than 10 of ~12,000 normal samples).

| # | Step | n | Where it is computed | Output |
|---|---|---|---|---|
| 0 | SV records in the raw VCF | 17,953 | `vcf.read_breakends()` | — |
| 1 | Admissible records | 2,286 | `vcf.admit()` | `<sample>.admitted.vcf` |
| 2 | Junctions | 1,143 | `vcf.pair_junctions()` | `stage1_junctions.tsv` |
| 3 | Candidate neopeptides (unique) | 20,721 | `generators/neosv.py` | `<sample>.all_neopeptides.txt` |
| 4 | Identical to a catalogue peptide (rows) | **64** | `cross.level1_identical()` | `stage3_matches.tsv` |
| 4b | …distinct peptides among them | 64 | `cross.level1_identical()` | `stage3_matches.tsv` |
| 5 | Gene-concordant | 64 | `cross.level1_identical()` | `stage3_matches.tsv` |
| 6 | Credible sequence | 57 | `cross.sequence_qc()` | `stage3_matches.tsv` |
| 7 | **Genomic events** | **13** | `cross.to_events()` | `credible_events.tsv` |
| 8 | Private (panel + population) | 9 | `confidence.annotate_privacy()` | `credible_events.tsv` |
| 9 | High-confidence SV calls | 4 | `confidence.annotate_confidence()` | `credible_events.tsv` |
| 10 | Junction-crossing reads in RNA | 1 | `rna.evaluate()` | `stage8_rna_evidence.tsv` |
| 11 | **High-confidence and RNA-supported** | **0** | `run.run_sample()` | `funnel.tsv` |

The funnel row order is fixed in `synthesis.FUNNEL_ROWS`, so every sample emits
the same rows and samples are comparable line by line. Rows 4 and 4b differ when
one peptide is produced by several SVs; on the unfiltered branch of this same
sample they are 387 and 288.

### All four lines

Every line is shown on both branches. `PON10` enforces the panel of normals;
`noPON` reports it without filtering on it — see `criteria.SOMATIC_KEEP_FILTERS_REPORT_PON`
for why that needs two levers rather than one threshold.

| | WT<br>PON10 | WT<br>noPON | TP53<br>PON10 | TP53<br>noPON | BRCA1<br>PON10 | BRCA1<br>noPON | BRCA2<br>PON10 | BRCA2<br>noPON |
|---|---|---|---|---|---|---|---|---|
| SV records | 17,953 | 17,953 | 350 | 350 | 303 | 303 | 234 | 234 |
| Caller `FILTER=PON` admitted | — | — | — | 82 | — | 58 | — | 96 |
| Admitted | 2,286 | 14,340 | 102 | 176 | 172 | 226 | 114 | 200 |
| Junctions | 1,143 | 7,170 | 51 | 88 | 86 | 113 | 57 | 100 |
| Candidate peptides | 20,710 | 63,036 | 64 | 78 | 536 | 572 | 54 | 148 |
| **Distinct matches** | **64** | **288** | **0** | **0** | **21** | **21** | **0** | **0** |
| **Events** | **13** | **46** | **0** | **0** | **1** | **1** | **0** | **0** |
| Private | 9 | 24 | — | — | 0 | 0 | — | — |
| High-confidence | 4 | 15 | — | — | 0 | 0 | — | — |
| RNA-supported | 1 | 3 | — | — | 0 | 0 | — | — |
| **Private + HC + RNA** | **0** | **0** | — | — | **0** | **0** | — | — |

Read `Private` against its branch: on `noPON` the panel count does not vote, so
privacy there rests on population frequency alone. That is why WT's private count
rises from 9 to 24 while the branch is *less* selective, not more — a different
question is being answered, not the same one more loosely.

**Relaxing the panel adds junctions to the derived lines and no matches.** The
82, 58 and 96 caller-rejected records admitted on `noPON` raise the junction
count by 73%, 31% and 75%, and the distinct-match count by nothing at all: 0, 21
and 0, identical to the filtered branch. Whatever limits recurrence detection in
the derived lines, it is not the panel filter.

Two observations the per-line numbers make plain:

**The only transcribed loci are common polymorphisms.** On the unfiltered branch
the three events with junction-crossing reads are ITGA11 (31 reads), ROBO1
(1 read) and GOLGA3 (2 fragments, absent from gnomAD). Frequencies are quoted as
`gnomad_af_popmax` with the population named, because a global average hides the
spread that matters: ITGA11 runs **0.591 (sas) to 0.905 (afr)** — never below
half in any of the nine groups — and ROBO1 peaks at **0.609 (mid)**. The two
events with real read support are common in every population sampled.

**RPE1-TP53-BRCA1's single event is not what it appeared to be.** *SLC9A9*, a
52 bp deletion carrying 21 sliding-window peptides, is absent from the other
lines' call sets — line-specific *within the experiment* — but carries a
population frequency of **0.181**. Specific to the line and a common polymorphism
are not mutually exclusive, and only the second filter catches it.

**Is 64 more than chance?** Yes, decisively. `null_model.permutation_test()`
shuffles each candidate peptide's residues, preserving length and composition,
and repeats the cross `criteria.NULL_PERMUTATIONS` (1,000) times:

| observed | expected by chance | enrichment | empirical p |
|---|---|---|---|
| 64 | 1.3 ± 0.8 | 49× | < 0.001 |

Output: `null_model_rates.tsv`. The strategy is named in
`criteria.NULL_STRATEGY` (L370).

Reducing a four-figure candidate count to a single-digit set is the expected
outcome — a 28-team benchmark found ~6% of top-ranked neoantigen predictions
validate functionally. What matters is knowing which gate removed what.

---

## 3. Stage-by-stage, with the code

### Stage 1 — Admission · `src/svneo/vcf.py`

| Operation | Function | Criteria |
|---|---|---|
| Parse every VCF record, extract evidence fields and inserted sequence | `read_breakends()` L75 | — |
| Apply FILTER / pairing / panel filters | `admit()` L107 | `SOMATIC_KEEP_FILTERS` L34, `PON_MAX` L69, `PON_ABSENT_MEANS` L73 |
| Type-aware lesion size | `event_size()` L151 | — |
| Collapse mated breakends to junctions | `pair_junctions()` L170 | — |
| Write the admitted set for stage 2 | `write_admitted_vcf()` L241 | — |
| Estimate panel size, so `PON_COUNT` can be read as a fraction | `panel_size_estimate()` L282 | `criteria.pon_fraction()` |

The breakend span is **not** the lesion size: a deletion removes `span − 1`
bases, an insertion always spans 1 while inserting many. `event_size()` encodes
this, and `insert_len` survives to stage 8 because an insertion-driven event can
only be validated by searching for the inserted sequence.

The generator applies no admission logic of its own, so it must receive
`<sample>.admitted.vcf`, never the raw file — handed the raw VCF it builds
peptides from panel and copy-number-inferred breakends.

### Stage 2 — Peptide generation · `src/svneo/generators/`

| Operation | Where | Criteria |
|---|---|---|
| Backend selection | `generators/__init__.py` `get_generator()` | config `peptide_generator:` |
| Contract every backend must satisfy | `generators/base.py` `validate_output()` | — |
| NeoSV run, stopping before MHC | `generators/neosv.py` `NeoSVGenerator.generate()` | `PEPTIDE_LENGTHS` L120, `RUN_MHC_PREDICTION` L126 |
| SV identifiers recovered from the VCF | `_neosv_extensions.build_sv_id_map()` L40 | — |
| Junction offset in the fusion protein | `_neosv_extensions.junction_indices()` L81 | — |
| Candidate table written before any MHC step | `_neosv_extensions.write_all_neopeptides()` L130 | `PEPTIDE_COLUMN` L130 |

Peptide sequences come from NeoSV (MIT, vendored at `vendor/neosv/`) with **one
patch**: `vendor/patches/001-pyensembl-stop-codon-frame`. Without it, with
pyensembl > 2.3.13, the 3′ side of every fusion is one amino acid out of frame
and the tool emits 170 spurious peptides where the patched version emits 64,
with none in common. Evidence and reproduction in `vendor/VENDOR.md`.

MHC binding prediction is deliberately excluded: a sequence-identity test is
HLA-independent, and an IC50 filter here would discard true matches.

`spans_junction` is emitted per peptide but **does not filter**
(`REQUIRE_SPANS_JUNCTION = False`, L97). After a frameshift every downstream
residue exists only because of the SV; filtering on it discarded 63 of 64 real
matches.

### Stages 3–5 — Cross, QC, event collapse · `src/svneo/cross.py`

| Operation | Function | Criteria |
|---|---|---|
| Level 1: identical peptide, with concordance flags | `level1_identical()` L84 | — |
| Level 2: same gene + same SV type | `level2_gene_svtype()` L151 | — |
| Level 3: breakpoint proximity gradient | `level3_proximity()` L164 | `PROXIMITY_WINDOWS_KB` |
| Low-complexity and self-proteome tests | `sequence_qc()` L290 | `LC_*` L108–111 |
| The credible rule | `criteria.is_credible()` | — |
| Collapse peptides to genomic events | `to_events()` L312 | `EVENT_DEDUP_COLUMN` L225, `REPORT_BOTH_EVENT_GRAINS` L232 |

Gene concordance is evaluated against **both** breakends: a junction joins two
genes and the catalogue peptide may be annotated to either. Testing only the
first side lost 8 of 64 matches, all translocations.

`credible = high-complexity AND non-self AND gene-concordant`. The cost is
stated wherever the count is quoted: a genuine convergent peptide from a
different gene is excluded by construction.

### Stage 6 — Confidence and privacy · `src/svneo/confidence.py`

| Operation | Function | Criteria |
|---|---|---|
| High-confidence verdict per event | `annotate_confidence()` L40 → `criteria.event_is_hc()` | `HC_MIN_SEGMAPQ` L240, `HC_MIN_VF` L241, `HC_MIN_QUAL` L242, `HC_MIN_SV_SIZE` L243 |
| Panel + population verdict | `annotate_privacy()` L61 → `criteria.is_private()` | `PON_MAX` L69, `GNOMAD_MAX_AF` L79 |
| Population frequency by reciprocal overlap | `annotate_gnomad()` L129, `_overlap_join()` L153 | `GNOMAD_RECIPROCAL_OVERLAP` L83 |

Two orthogonal questions kept apart. **Confidence** asks whether the call is
real; **privacy** asks whether it is the sample's own. An event can be perfect
on every evidence metric and still be disqualified for being a common
polymorphism.

Copy-number fields are excluded from every credibility judgement
(`TRUST_COPY_NUMBER_FIELDS = False`, L196): they depend on a purity/ploidy fit
that is unreliable for clonal samples.

### Stages 7–8 — Expression and junction evidence · `src/svneo/rna.py`

| Operation | Function | Criteria |
|---|---|---|
| Gene and disrupted-isoform TPM, with background gradient | `expression()` L39 | `TPM_GRADIENT` L329, `TPM_EXPRESSED` L334, `EXPRESSION_USES_BOTH_BREAKENDS` L338 |
| Choose the only valid test for the geometry | `select_test()` L217 | `MIN_TESTABLE_GAP_SIZE` L360, `MIN_INSERT_LEN` L367 |
| Count evidence at one breakend | `count_at_breakend()` L268 | `COVERAGE_WINDOW` L365, `NGAP_*` L307–308, `MIN_SOFTCLIP_LEN` L366, `MIN_READ_MAPQ` L392, `COUNT_UNIQUE_FRAGMENTS` L397 |
| Verify a supplementary alignment hits the partner | `sa_hits_partner()` L242 | `SA_PARTNER_TOLERANCE` L364 |
| Assign the evidence tier | `rna_tier()` L338 | `STRONG_MIN_JUNCTION_READS` L351, `SUGGESTIVE_MIN_JUNCTION_READS` L355 |

Three rules are enforced here and declared in the manifest so a reader cannot
mistake them:

- `SOFTCLIPS_TIER_EVENTS = False` (L348) — a soft-clip shows a read *ends* at a
  breakpoint, never that it crosses.
- `COVERAGE_TIERS_EVENTS = False` (L349) — a breakpoint inside an expressed gene
  has thousands of reads whether or not the junction exists.
- `COVERAGE_SUMMARY = "min"` (L354) — an event is in transcribed territory only
  if *both* ends are.

Test selection follows event geometry: `chimeric` for translocations, `sizegap`
for deletions above the aligner's minimum intron, `insertion` for
insertion-driven events. An event with no applicable test is `UNTESTABLE`, not
`NONE` — the latter asserts "tested and negative".

`MIN_READ_MAPQ = 20` matches GRIDSS2's unmapped threshold, the lineage that
produced these SV calls. In RNA it is not tunable: STAR emits only 0, 1, 3 and
255, with zero reads between 10 and 20 across 73,232 alignments measured.

### Stage 9 — Synthesis · `src/svneo/synthesis.py`

| Operation | Function | Criteria |
|---|---|---|
| The fixed-row funnel | `build_funnel()` L46, `FUNNEL_ROWS` | — |
| Earliest sample carrying each breakpoint | `attribute()` L53 | `ATTRIBUTION_MATCH_TOLERANCE` L435 |
| Cross-sample table with per-1,000 rates | `compare_samples()` L106 | — |
| Plain-language cautions attached to the numbers | `interpret()` L130 | — |

---

## 4. What validation looks like on one candidate

**GOLGA3** — the only event in RPE1-WT with any RNA support, and a worked
example of why a count is not evidence.

The call is an insertion misreported by breakend geometry as a 2 bp deletion:
breakend span 3, deleted length 2, **34 bases inserted**. That distinction
decides the test (`rna.select_test()` returns `insertion`, not `sizegap`); a gap
test could never confirm it.

The inserted sequence is `CGCCGGGAAGCAGGAGGGCTGGGGAGGCGGGGGG` — 82% GC, Shannon
entropy 1.40 bits, six-base homopolymer. That profile is what a chance match
looks like, so specificity was measured rather than assumed:

| Check | Result |
|---|---|
| Occurrences in 328,868 reference transcripts | **0** |
| Occurrences in 22,269 reads from 200 random loci | **0** |

| Modality | Evidence | Verdict |
|---|---|---|
| **DNA** | 28 alignments, 23 distinct in-read offsets, all MAPQ 60, no duplicates | the lesion is **real** |
| **RNA** | 7 alignments → **3 fragments**; one contributed 4 records at MAPQ 3; 2 well-mapped fragments remain, all heavily soft-clipped | transcription **unproven** (`WEAK`) |

A genuine germline insertion in the parental genome, not demonstrably
transcribed across the junction. Counting alignments rather than fragments, and
ignoring mapping quality, would have graded this `STRONG`.

Reproduce with `tools/inspect_insertion.py`, which prints every supporting read
with its MAPQ, duplicate flag, CIGAR and in-read offset, plus a background rate
over random loci.

---

## 5. What this design cannot answer

- **Presentation.** Read evidence supports "the junction is transcribed".
  Claiming "presented" needs immunopeptidome mass spectrometry.
- **Convergent peptides from a different gene**, excluded by the credibility
  definition (`criteria.is_credible()`).
- **A negative RNA result is not absence of the lesion**: the gene may not be
  expressed at all, and nonsense-mediated decay of an aberrant transcript is a
  real possibility.
- **No liftover.** Reference and samples must share genome build and annotation
  release; there is no coordinate conversion anywhere.

---

## 6. Reproducing this run

```bash
PYTHONPATH=src python -m svneo.run --config config/rpe1_hmf.yaml \
    --samples RPE1-WT --branches PON10
```

Every output directory carries `summary.json`, which contains the counts, the
null-model result, the generator version (including the vendored commit) and the
**full threshold manifest** produced by `criteria.manifest()`. A result can
therefore always be traced to the exact criteria that produced it — and cannot
outlive them.

Then assemble the auditable candidate table — one row per candidate, one column
per criterion, plus `first_failed_gate` naming the earliest gate each candidate
failed:

```bash
python tools/build_master_table.py --results-dir results
```

To review the logic rather than the numbers, `notebooks/pipeline_walkthrough_pyspark.ipynb`
re-derives every stage independently in PySpark and asserts agreement with the
pipeline at each step (16 checks). Change `SAMPLE` in its first cell to walk
through any line.

Threshold derivations, and the observation that set each one, are in
[`PROVENANCE.md`](PROVENANCE.md). The regression tests that pin them are in
`tests/test_svneo.py`.
