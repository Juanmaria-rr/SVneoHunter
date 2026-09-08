#!/usr/bin/env python3
"""
build_patient_catalogue.py — regenerate the patient peptide catalogue with the
minus-strand fix applied.

WHY
---
The reference catalogue was built with an unpatched NeoSV. Verified, not assumed:
its source file shows 79.2% `Start-loss` on the minus strand against 13.0% on the
plus — this repository's own signature before patch 002.

That makes the current comparison one-sided. A genuine minus-strand neoantigen is
not in the catalogue under its true sequence, so a corrected cell-line candidate
cannot match it: the recurrence is unfindable rather than absent. Every match
count reported so far is therefore a **floor**, and rebuilding the patient side
can only raise it.

WHAT THIS DOES, AND DOES NOT
----------------------------
Regenerates the **peptide layer only** — stage 1 admission and stage 2 generation
per patient, then an aggregated ranking with the same columns as the existing
`*.neoantigen_ranking.tsv`, so the two can be compared directly.

It does **not** reproduce the curation. The 2,856-peptide shortlist carries
columns from an immune-selection analysis (`or_clean`, `confirmed_ge1`/`ge5`,
`is_cfs`, `hla_pres_cov_*`) produced by the group that built it. Those cannot be
recomputed here, and this tool does not pretend to.

Nothing is overwritten: output goes to its own directory, one file per patient,
and existing files are skipped so a long run can be resumed or split.

    python tools/build_patient_catalogue.py --limit 5        # pilot, then time it
    python tools/build_patient_catalogue.py                  # the whole cohort
    python tools/build_patient_catalogue.py --aggregate-only # rebuild the ranking
"""
from __future__ import annotations

import argparse
import os
import pathlib
import sys
import time
import traceback

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

#: Columns of the existing ranking, reproduced so the two tables are comparable
#: column for column. `frameshift` keeps upstream's name: it is what the old file
#: calls the frame verdict, and renaming it here would break the comparison.
RANKING_COLUMNS = ["neoantigen", "gene1", "gene2", "svtype", "svpattern",
                   "frameshift", "n_rows", "n_patients", "n_samples"]


def patient_id(path: pathlib.Path) -> str:
    """`H00000004.purple.sv.vcf.gz` -> `H00000004`."""
    return path.name.split(".")[0]


def process(vcf_path: pathlib.Path, out_dir: pathlib.Path, cache: str,
            release: int) -> tuple[str, int, str]:
    """One patient: admit, generate peptides, write them. Returns a status."""
    from svneo import criteria, generators, vcf as vcf_mod

    name = patient_id(vcf_path)
    target = out_dir / f"{name}.peptides.tsv"
    if target.exists():
        return name, -1, "skipped (already present)"

    work = out_dir / "_work" / name
    work.mkdir(parents=True, exist_ok=True)
    try:
        breakends = vcf_mod.read_breakends(str(vcf_path))
        # Patient call sets are somatic: the caller has already applied its panel
        # filter, so admission is FILTER == PASS. Matching how the original
        # catalogue was built matters more here than any preference of ours —
        # the point is to isolate the generator change, nothing else.
        admitted, _ = vcf_mod.admit(breakends, "somatic", pon_max=None)
        if admitted.empty:
            target.write_text("\t".join(["neoantigen", "gene1", "gene2",
                                         "svtype", "frame_effect"]) + "\n")
            return name, 0, "no admitted breakends"

        admitted_vcf = work / f"{name}.admitted.vcf"
        vcf_mod.write_admitted_vcf(str(vcf_path), admitted, str(admitted_vcf))

        backend = generators.get_generator("neosv")
        generated = backend.generate(str(admitted_vcf), str(work), name,
                                     release=release, cache_dir=cache)
        peptides = generated.peptides
        keep = [c for c in ("neopeptide", "gene1", "gene2", "svtype",
                            "frame_effect") if c in peptides.columns]
        table = peptides[keep].rename(columns={"neopeptide": "neoantigen"})
        table = table.drop_duplicates()
        table.to_csv(target, sep="\t", index=False)
        return name, len(table), "ok"
    except Exception as error:                                  # noqa: BLE001
        (out_dir / "_failures").mkdir(exist_ok=True)
        (out_dir / "_failures" / f"{name}.txt").write_text(traceback.format_exc())
        return name, 0, f"FAILED: {type(error).__name__}: {error}"
    finally:
        # The intermediate VCF and generator scratch are large and worthless once
        # the peptides are written; 6,378 of them would be tens of gigabytes.
        for leftover in work.glob("*"):
            try:
                leftover.unlink()
            except OSError:
                pass


def aggregate(out_dir: pathlib.Path) -> pd.DataFrame:
    """Every patient's peptides collapsed into the ranking table."""
    files = sorted(out_dir.glob("*.peptides.tsv"))
    if not files:
        raise SystemExit(f"no per-patient files in {out_dir}")

    frames = []
    for path in files:
        try:
            table = pd.read_csv(path, sep="\t", low_memory=False)
        except pd.errors.EmptyDataError:
            continue
        if table.empty:
            continue
        table["patient"] = path.name.split(".")[0]
        frames.append(table)

    combined = pd.concat(frames, ignore_index=True)
    grouped = combined.groupby("neoantigen", dropna=True)
    ranking = grouped.agg(
        gene1=("gene1", lambda s: ";".join(sorted(set(s.dropna().astype(str))))),
        gene2=("gene2", lambda s: ";".join(sorted(set(s.dropna().astype(str))))),
        svtype=("svtype", lambda s: ";".join(sorted(set(s.dropna().astype(str))))),
        frameshift=("frame_effect",
                    lambda s: ";".join(sorted(set(s.dropna().astype(str))))),
        n_rows=("patient", "size"),
        n_patients=("patient", "nunique"),
    ).reset_index()
    # `svpattern` is not carried through this path; emitted empty so the column
    # set matches the original file rather than silently differing.
    ranking["svpattern"] = ""
    ranking["n_samples"] = ranking["n_patients"]
    return ranking[RANKING_COLUMNS].sort_values(
        ["n_patients", "n_rows"], ascending=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf-dir",
                        default="../hmf/hg38_SVs",
                        help="directory of patient *.purple.sv.vcf[.gz]")
    parser.add_argument("--out-dir", default="../results_hmf_patch002",
                        help="written here; existing per-patient files are kept")
    parser.add_argument("--limit", type=int, default=None,
                        help="process at most this many patients — use for a "
                             "pilot before committing to the full cohort")
    parser.add_argument("--release", type=int, default=115)
    parser.add_argument("--cache", default=os.environ.get("PYENSEMBL_CACHE_DIR"))
    parser.add_argument("--aggregate-only", action="store_true",
                        help="skip generation, rebuild the ranking from what is "
                             "already on disk")
    args = parser.parse_args()

    out_dir = (REPO / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not args.aggregate_only:
        vcf_dir = (REPO / args.vcf_dir).resolve()
        vcfs = sorted(list(vcf_dir.glob("*.purple.sv.vcf.gz"))
                      + list(vcf_dir.glob("*.purple.sv.vcf")))
        if not vcfs:
            raise SystemExit(f"no patient VCFs in {vcf_dir}")
        done = {p.name.split(".")[0] for p in out_dir.glob("*.peptides.tsv")}
        pending = [v for v in vcfs if patient_id(v) not in done]
        if args.limit:
            pending = pending[:args.limit]

        print(f"{len(vcfs):,} VCFs found, {len(done):,} already done, "
              f"{len(pending):,} to process now\n")

        started, peptides_total, failed = time.time(), 0, 0
        for index, path in enumerate(pending, 1):
            name, count, status = process(path, out_dir, args.cache, args.release)
            if status.startswith("FAILED"):
                failed += 1
            elif count > 0:
                peptides_total += count
            elapsed = time.time() - started
            rate = elapsed / index
            print(f"  [{index}/{len(pending)}] {name}: {status}"
                  + (f", {count:,} peptides" if count > 0 else "")
                  + f"   ({rate:.1f}s each)")
            if index == 10 or (index % 200 == 0):
                remaining = (len(vcfs) - len(done) - index) * rate
                print(f"      -- projected for the rest of the cohort: "
                      f"{remaining / 3600:.1f} h")

        print(f"\n{len(pending)} processed, {failed} failed, "
              f"{peptides_total:,} peptide rows written")
        if failed:
            print(f"tracebacks in {out_dir / '_failures'}")

    ranking = aggregate(out_dir)
    path = out_dir / "neoantigen_ranking_patch002.tsv"
    ranking.to_csv(path, sep="\t", index=False)
    print(f"\nranking: {len(ranking):,} unique peptides -> "
          f"{path.relative_to(REPO.parent)}")


if __name__ == "__main__":
    main()
