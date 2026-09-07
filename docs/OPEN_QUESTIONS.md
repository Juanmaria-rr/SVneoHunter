# Open questions

## 1. Minus-strand transcripts lose their 5' CDS on an intronic breakpoint

**RESOLVED — root cause identified 2026-09-07. The fix is not yet applied, and
every result in this repository was produced with the bug present.**

Reproduce with `python tools/diagnose_minus_strand_cds.py` (exits 1 while the
behaviour is present).

### The defect

`neosv.fusion_utils.truncate_cds(transcript, '5', pos)` returns the coding
sequence 5' of a breakpoint. For a **minus-strand transcript with an intronic
breakpoint** it returns an empty sequence, `cut_length == 0`:

| Breakpoint falls in | Plus strand | Minus strand |
|---|---|---|
| a coding exon | correct | correct |
| an intron | correct (1,006 / 253 / 1,746 nt) | **empty, 6 of 6 genes tested** |

Tested on ITGA11, TP53, BRCA1, GOLGA3, KRAS, BRAF (minus) and EGFR, PTEN,
PIK3CA (plus). Introns are far larger than exons, so most real SV breakpoints
land in them.

### Why it happens

`transcript.coding_sequence_position_ranges` is in **genomic** order, while
`truncate_cds` treats it as transcript order — `transcript_utils.get_cds_range`
states "from 5' to 3'" in its own docstring. The two coincide on the plus strand
and are reversed on the minus strand. The exon-overlap branch survives this; the
branch that handles a breakpoint outside a coding exon, which indexes through
`get_noncds_range`, does not.

### What it explains

One root cause accounts for every symptom previously listed here as separate:

- **84.0% of minus-strand fusions are flagged `Start-loss`** against 13.8% of
  plus-strand ones. `Start-loss` is upstream's own low-reliability warning: the
  start codon is gone, so it falls back to the next ATG and its comment says
  "such prediction is of low reliability, we should annotate it and remove these
  fusions when necessary."
- **The fusion protein becomes the 3' partner alone**, so the junction sits at
  residue ~0 — median `junction_aa` is 3 among catalogue-matching peptides, and
  0 for ITGA11.
- **No peptide can span the breakpoint.** 1 of 408 matching peptide rows spans
  its junction where 67 would be expected from the 16.5% base rate, p = 1e-30.
- **The strand skew appears at annotation, not at the breakends.** Measured:
  admitted breakends of RPE1-WT are 47.6% minus (p = 0.33, background 48.6%),
  while the candidate peptides they produce are 41.5% minus (p = 5e-4). Step 1
  of the old plan, now executed, and it rules out "what breaks" as the cause.
- **The catalogue leans the other way** (60.9% minus) and matched candidates
  lean harder still (74.4%). Both sides were built with the same tool family, so
  both carry the artefact; an intersection can be more extreme than either input
  when the matching itself is driven by a shared signature.

### What is NOT established

- That the matches are artefactual. A `Start-loss` fusion can still produce a
  genuinely novel protein — `is_self` is False for 403 of 408 matched peptides,
  so they are not in the normal proteome. What is established is that they are
  not junction-spanning and that their frame rests on a fallback start codon
  upstream itself calls unreliable.
- The exact index arithmetic at fault. The behaviour is reproducible and
  strand-separated; the precise off-by-one has not been isolated.
- Whether fixing it changes the headline result. It cannot be known without
  re-running, and the headline is a negative finding that a more sensitive
  minus-strand path could only add to.

### Before this was understood

The three-way skew below was recorded as an unexplained observation, with
"pyensembl returns exons 5'->3'" listed as ruled out. That was verified on
`transcript.exons` — but `coding_sequence_position_ranges`, which is what
`truncate_cds` actually reads, is genomic-ordered. The exclusion was drawn from
the wrong accessor.

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
