# Open questions

Findings that are real, reproducible, and **not yet explained**. Recorded here
rather than in the limitations section because a limitation is understood and
these are not.

---

## 1. Gene strand is skewed at three levels, in two directions

**Status:** open · **Found:** 2026-08-14 · **Reproduce:** `tools/check_strand_bias.py`

### The observation

Coding genes are split almost evenly between strands — 49.2% on minus by gene
count, 48.6% weighted by genomic length. The length-weighted figure is the fairer
null, since a structural variant is likelier to break a long gene. Against it:

| Gene set | n | on minus | p |
|---|---|---|---|
| Reference catalogue (patient side) | 379 | **60.9%** | 1.6 × 10⁻⁶ |
| Our candidate peptides (unfiltered branch) | 614 | **41.5%** | 5.1 × 10⁻⁴ |
| Candidates that match the catalogue | 43 | **74.4%** | 6.7 × 10⁻⁴ |
| Events with RNA support | 3 | 100% | 0.11 (n too small) |

Three significant deviations. The reference catalogue leans **minus**, our
candidate pool leans **plus**, and their intersection leans minus harder than
either input. That last part is what makes a simple inherited-bias explanation
insufficient: an intersection cannot ordinarily be more extreme than both of the
sets it comes from unless something in the matching is strand-dependent.

All three events with junction-crossing RNA reads are on the minus strand.

### What has been ruled out

**The 5′/3′ assignment is strand-aware, and verified.** The suspicion that
prompted this — minus-strand transcripts annotated with plus-strand coordinate
logic — does not hold for the annotation step:

- `vendor/neosv/annotation_utils.py` branches on all four strand combinations
  (`+/+`, `+/−`, `−/+`, `−/−`) crossed with the four breakend patterns, and
  truncates the terminal exon at the correct edge per strand: `end` for plus,
  `start` for minus.
- That is only correct if exons arrive in transcript order. **They do**:
  pyensembl returns `transcript.exons` 5′→3′, not in genomic order. Verified on
  GOLGA3, ITGA11, TP53 and BRCA1 — all minus-strand, all with the first exon at a
  *higher* coordinate than the last.

So the inversion artefact is not in the exon-retention logic. Something else is
strand-dependent.

### What has not been ruled out

- **Gene selection.** Which genes are broken in this sample may itself be
  strand-correlated for reasons unrelated to annotation (fragile sites,
  replication timing, repeat content).
- **The reference cohort's own composition.** The catalogue was generated
  upstream, by a different run of the same tool family. Its 60.9% is not
  something this pipeline produced.
- **Frame or sequence assembly downstream of exon retention.** Ruling out the
  5′/3′ assignment does not rule out the junction sequence being built from the
  wrong end for one strand.
- **Small numbers.** n = 43 for the matched set. The catalogue's n = 379 does not
  have that excuse.

### How to settle it

In order of decisiveness:

1. **Test the admitted breakends, before peptide generation.** If the skew is
   already present in which genes are broken, the cause is upstream of annotation
   entirely. This separates "what breaks" from "how it is annotated" and is the
   single most informative next step. `check_strand_bias.py` currently measures
   the peptide tables; extend it to `stage1_junctions.tsv` joined to gene
   annotation.
2. **Round-trip a known minus-strand fusion.** Take one annotated event on the
   minus strand, reconstruct the expected junction protein by hand from the two
   transcripts, and compare against the generated `aa_sequence`. A one-amino-acid
   or wrong-side error shows immediately.
3. **Compare strand composition of peptides that match versus those that do
   not**, within the same candidate pool. If matching itself is strand-dependent,
   the two distributions differ.
4. **Ask whether the catalogue's 60.9% is expected.** It was built from patient
   tumours with the same tool family; if that lean is present there too, it is a
   property of the method rather than of either dataset.

### Why it matters

If some strand-dependent step is systematically mis-building fusion sequences,
the peptides on that strand are wrong — and since matching is by exact peptide
string, wrong peptides would mostly fail to match, biasing *against* the affected
strand rather than producing false positives. The matched set leaning **towards**
minus is therefore the opposite of what a naive "minus-strand peptides are
corrupted" story predicts, which is itself a clue.

No result in this repository should be read as strand-controlled until this is
resolved.
