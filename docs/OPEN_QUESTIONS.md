# Open questions

## 1. Minus-strand transcripts get the wrong 5' CDS, in every region

**RESOLVED and FIXED.** Root cause identified 2026-09-07; patch 002 applied
2026-09-08 and its effect measured by re-running the whole analysis into
`results_patch002/` rather than overwriting `results/`.

**Measured effect of the fix** (`results_patch002/COMPARISON.md`):

| | before | after |
|---|---|---|
| `Start-loss`, minus strand | 84.0% | **10.1%** |
| `Start-loss`, plus strand | 13.8% | 13.4% |
| RPE1-WT candidates on minus | 41.5% (p = 5x10⁻⁴) | **46.8% (p = 0.36)** |
| Matched candidates on minus | 74.4% (p = 7x10⁻⁴) | **60.9% (p = 0.30)** |
| RPE1-WT distinct matches | 288 | **130** |
| RPE1-TP53-BRCA1 events | 1 (SLC9A9) | **0** |
| Candidates satisfying every criterion | 0 | 0 |

The skew resolves in everything this pipeline generates, and only on the strand
the defect touched. **No new match appeared anywhere** — 130 of 288 survive, 158
vanish, 0 are new — so the defect was manufacturing matches rather than hiding
them.

The only set still deviating is the **reference catalogue** (60.9%, p = 2x10⁻⁶),
which is the patient side and cannot be corrected from here. Related: matched
peptides still almost never span their junction (1 of 161, against 18.2% across
the candidate universe), a ~30-fold depletion the patch does **not** explain.
Both point the same way — the catalogue was built with the same tool family and
probably carries the artefact still. See question 2.

Reproduce with `python tools/diagnose_minus_strand_cds.py` (exits 1 while the
behaviour is present).

### The defect

`truncate_cds(transcript, '5', pos)` returns what a transcript contributes as the
**head** of a fusion protein: its coding sequence from the start codon to the
breakpoint. On minus-strand transcripts it returns the wrong length **in every
region tested** — not only in introns, as first recorded here.

Swept across every 5'UTR, coding exon, intron and 3'UTR of nine transcripts:

| | regions | correct |
|---|---|---|
| EGFR, PTEN, PIK3CA (plus) | 87 | **87** |
| ITGA11, TP53, BRCA1, GOLGA3, KRAS, BRAF (minus) | 188 | **0** |

Two distinct failure modes on the minus strand:

| Breakpoint in | Result |
|---|---|
| an **intron** | the head comes back **empty** — no start codon, flagged `Start-loss` |
| a **coding exon** | the length is counted from the wrong end: too short near the start of the transcript, **too long** near its end |

**The second mode is the dangerous one.** Sequence the gene does not contribute
is added to the fusion, and the sliding window turns it into peptides. Measured
on a worked example, a breakpoint inside the last coding exon returned 278 nt
where 26 were correct — about 84 residues of protein the gene never contributed.
The defect can **fabricate** candidates, not only lose them.

Reproduce with `python tools/diagnose_minus_strand_cds.py`, which sweeps every
region and exits 1 while the behaviour is present.

**Correction to an earlier version of this entry.** It stated that the exonic
branch was correct on both strands. That was verified by checking the returned
sequence began with `ATG` — a check that passes at any non-zero length, because
the sequence is always taken from position 0 of the coding sequence. Starting
correctly is not measuring correctly, and the exonic branch is wrong on the minus
strand too.

### The mechanism, in plain terms

Two facts set it up.

**DNA coordinates always count left to right.** A plus-strand gene is *read* in
that direction too. A minus-strand gene is read the other way, right to left, so
**the beginning of a minus-strand gene sits at the highest coordinate**.

**The annotation library hands over coding exons sorted by coordinate**, always.
For a plus-strand gene that order is also the reading order. For a minus-strand
gene it is the reverse of the reading order.

Take a gene with three coding exons:

```
coordinate:  100───200   300───400   500───600
exon:          [A]         [B]         [C]
                    gap         gap
```

Read on the plus strand: A, B, C — it starts at A.
Read on the minus strand: C, B, A — it starts at **C**.

Either way the library returns `[A, B, C]`.

**What the code needs.** When a breakpoint lands in an intron, the tool must know
how many exons precede it in *reading* order, so it can keep that leading piece
of the protein. To find out, it first builds the list of gaps: the stretch before
the first exon, then each intron.

**What it does instead.** `get_noncds_range` knows minus-strand genes are read
backwards, so it has a separate branch that reverses the subtractions. But that
branch assumes the exon list already arrives in reading order, `[C, B, A]`. It
arrives as `[A, B, C]`. The correction is applied to data that was never
reversed, and the gaps come out like this:

| Gap | Should be | Comes out as |
|---|---|---|
| before the first exon | above 600 | **201 → 600** |
| intron 1 | 401 → 499 | **401 → 99** |
| intron 2 | 201 → 299 | **601 → 299** |

Two things break at once. **Every intron is inverted** — "401 to 99" is not a
range, so no position can ever fall inside it. And **the first gap swallows most
of the gene**: it is measured from the end of the first exon *in the list* (A)
rather than the first in reading order (C), so instead of sitting above 600 it
collapses leftwards to cover exons B and C and both introns.

```
      50 .... 100---200 .... 300---400 .... 500---600 .... 650
              [  A  ]        [  B  ]        [  C  ]
                     |________ gap 0: 201 - 650 __________|
       outside
```

Exon A itself falls just outside, by one base, since the gap starts at its end
plus one. That changes nothing about the outcome. A breakpoint inside A takes the
exonic branch instead — which is *also* wrong on the minus strand, in the other
direction: it counts exons from the left of the list rather than from the start
of reading, so it returns the head that would belong to the mirrored breakpoint.

**The consequence.** A breakpoint at position 450 is tested against the gaps.
The introns cannot match, being inverted, so it matches the first gap — the one
that ate the gene. The tool concludes the breakpoint lies *before the first
exon*, and therefore that **zero exons precede it**. The leading piece of the
protein comes out empty.

Verified on the real TP53 transcript: the first gap is reported as
7,662,015–7,676,594, which is **14.5 kb of a 14.8 kb transcript**, and all seven
introns come back inverted.

### Why the plus strand escapes

There, coordinate order and reading order are the same, and the plus branch
reverses nothing. Every region of every plus-strand transcript tested returns the
correct length. That is also why the tool appears to work: half the genome is
processed correctly, and the other half returns output that is not obviously
absurd.

### Why it went unnoticed for so long

- **It fails silently.** An inverted range raises nothing; it simply never
  matches.
- **The obvious check passes.** The returned sequence always begins with `ATG`,
  since it is taken from position 0 of the coding sequence, so any test that
  looks at the start rather than the length reports success. Two probes written
  during this investigation did exactly that before being redone.
- **The intronic mode is silent and the exonic mode is plausible.** An empty head
  still yields a fusion — just one made entirely of the partner's tail — and a
  head of the wrong length still translates.
- **The docstring asserts the invariant that is violated.** `get_cds_range` is
  documented as returning ranges "from 5' to 3'", which is true only on the plus
  strand and is never checked.

### The fix, and its measured effect

Sort the coding ranges into reading order — descending coordinate for
minus-strand transcripts — before anything consumes them, restoring the invariant
the docstring already claims. Both `get_cds_range` and `get_noncds_range` read
the raw accessor, so both need it; `truncate_cds` is the only consumer, and its
minus-strand branches are already written for reading order and become correct
once they receive it.

Tested in memory, without modifying the vendored source:

| | regions correct |
|---|---|
| unpatched | 87 / 275 (32%) |
| patched | **275 / 275 (100%)** |

Plus-strand results are unchanged, which is the necessary condition: the fix must
repair the minus strand without disturbing the half that works.

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
