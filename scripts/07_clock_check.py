"""
07_clock_check.py
-----------------
Does this virus behave like a clock, and is the signal real or an artefact?

Dating a tree assumes that sequences sampled later are further from the root by
roughly a constant rate. That assumption has to be tested before it is used, not
after. TreeTime's root-to-tip regression does exactly that: it plots each tip's
genetic distance from the root against its sampling date and fits a line. The
slope is the substitution rate; the scatter tells you whether a clock is a fair
description of the data at all.

For DENV the expected rate is of order 5e-4 to 1e-3 substitutions per site per
year. A fitted rate far outside that range is a warning, not a result.

This script runs the regression under several scenarios, because one run cannot
distinguish a real signal from an artefact:

  wg_outgroup   whole genome, rooted on the 3I outgroup by IQ-TREE
  wg_lsq        whole genome, rooted by the clock itself (least squares)
                -> if these two agree, the root is trustworthy
  wg_no2025     whole genome with the 2025 batch removed
                -> the 2025 genomes have 80-89% coverage against 97-99% for the
                   earlier ones. Missing data lengthens branches. If the rate
                   changes a lot when they are dropped, the clock is being driven
                   by their missing data rather than by real evolution.
  wg_drop_marginal  whole genome without OZ484513.1, which passed the coverage
                   threshold by 0.003. A conclusion must not depend on one
                   borderline sequence.
  eg_lsq        E gene, clock-rooted
                -> most E-gene dates are year-only, so this is expected to be
                   worse. It is run to show that, not to be used for dating.

Nothing here produces a dated tree. This step decides whether dating is defensible.

Inputs : results/trees/{wg,eg}.treefile
         results/datasets/{wg,eg}.fasta, {wg,eg}_dates.tsv
Outputs: results/clock/<scenario>/   full TreeTime output incl. the regression plot
         results/clock/summary.tsv   rate, R^2 and outliers per scenario

Run from the project folder:
  python scripts/07_clock_check.py
"""

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import Phylo
from Bio.SeqIO.FastaIO import SimpleFastaParser

TREES = Path("results/trees")
DATA = Path("results/datasets")
OUT = Path("results/clock")

DROP_MARGINAL = ["OZ484513.1"]   # passed the whole-genome coverage cutoff by 0.003


def write_subset(stem, exclude, tag):
    """Make a pruned tree + alignment + dates file with `exclude` removed."""
    sub = OUT / "_subsets" / tag
    sub.mkdir(parents=True, exist_ok=True)
    excl = set(exclude)

    tree = Phylo.read(TREES / f"{stem}.treefile", "newick")
    for name in excl:
        try:
            tree.prune(name)
        except ValueError:
            pass
    Phylo.write(tree, sub / "tree.nwk", "newick")

    with (sub / "aln.fasta").open("w") as fh:
        for header, seq in SimpleFastaParser((DATA / f"{stem}.fasta").open()):
            acc = header.split()[0]
            if acc not in excl:
                fh.write(f">{acc}\n{seq}\n")

    d = pd.read_csv(DATA / f"{stem}_dates.tsv", sep="\t")
    d[~d["name"].isin(excl)].to_csv(sub / "dates.tsv", sep="\t", index=False)
    return sub / "tree.nwk", sub / "aln.fasta", sub / "dates.tsv"


def names_dated(dates_path, prefix):
    """Accessions in a dates file whose date starts with `prefix` (e.g. '2025')."""
    try:
        d = pd.read_csv(dates_path, sep="\t")
        return {n for n, dt in zip(d["name"], d["date"].astype(str))
                if dt.startswith(prefix)}
    except Exception:
        return set()


def run_clock(tag, tree, aln, dates, reroot=None):
    outdir = OUT / tag
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = ["treetime", "clock",
           "--tree", str(tree),
           "--aln", str(aln),
           "--dates", str(dates),
           "--outdir", str(outdir),
           "--clock-filter", "3"]        # flag tips >3 interquartile ranges off the line
    if reroot:
        cmd += ["--reroot", reroot]
    else:
        cmd += ["--keep-root"]

    print(f"\n{'=' * 62}\n  {tag}\n{'=' * 62}")
    print("  " + " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    log = (outdir / "treetime.log")
    log.write_text(res.stdout + "\n" + res.stderr)

    if res.returncode != 0:
        print(f"  FAILED (exit {res.returncode}). Last lines:")
        print("\n".join((res.stderr or res.stdout).strip().splitlines()[-12:]))
        return {"scenario": tag, "rate": None, "r2": None, "n_outliers": None,
                "status": f"failed ({res.returncode})"}

    text = res.stdout

    # TreeTime prints the numbers as '--rate:   9.231e-04' and '--r^2:  0.99',
    # then explains them in prose. The first version of this regex matched the
    # prose ("The rate is the slope...") and captured the letter 'e'. Anchor on
    # the '--' prefix and require the value to look like a number.
    NUM = r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"

    def grab(pattern):
        m = re.search(pattern, text)
        if not m:
            return None
        try:
            return float(m.group(1))
        except (TypeError, ValueError):
            return None

    rate = grab(r"--\s*rate\s*:\s*" + NUM)
    r2 = grab(r"--\s*r\^?2\s*:\s*" + NUM)

    # which tips the clock filter threw out - a flagged 2025 sequence would mean
    # the missing-data artefact is real
    flagged = []
    m = re.search(r"don't follow a loose clock.*?:\s*(.*?)(?:\n\s*\n|--|\Z)",
                  text, flags=re.S)
    if m:
        # the block lists an accession then its deviation; keep only accessions
        flagged = [x for x in re.split(r"[\s,]+", m.group(1).strip())
                   if re.match(r"^[A-Za-z]{1,3}\d{4,}(\.\d+)?$", x)]

    row = {"scenario": tag,
           "rate": rate,
           "r2": r2,
           "n_flagged": len(flagged),
           "flagged": ";".join(flagged[:20]),
           "status": "ok"}
    print(f"  rate = {rate}   R^2 = {r2}   flagged tips = {len(flagged)}")
    if flagged:
        print(f"  flagged: {', '.join(flagged[:20])}")
        # is any flagged tip from 2025? read the dates file rather than guessing
        # accession prefixes (the 2025 accessions share no common prefix)
        y2025 = names_dated(dates, "2025")
        hit = [f for f in flagged if f in y2025]
        if hit:
            print(f"  NOTE: {len(hit)} flagged tip(s) are from 2025 ({', '.join(hit)}).")
            print("        Check whether this is the missing-data artefact rather")
            print("        than real rate variation.")
    if rate is None:
        print("  (could not parse the rate - read the log directly:)")
        print("\n".join("    " + l for l in text.strip().splitlines()[-8:]))
    print(f"  full output in {outdir}/  (root_to_tip_regression plot is there)")
    return row


def main():
    if not shutil.which("treetime"):
        sys.exit("treetime not found. conda activate denv3")
    for p in [TREES / "wg.treefile", DATA / "wg.fasta", DATA / "wg_dates.tsv"]:
        if not p.exists():
            sys.exit(f"Missing {p} - run scripts/06_build_trees.sh first")

    OUT.mkdir(parents=True, exist_ok=True)
    rows = []

    # 1-2. whole genome: the two independent rootings
    rows.append(run_clock("wg_outgroup", TREES / "wg.treefile",
                          DATA / "wg.fasta", DATA / "wg_dates.tsv"))
    rows.append(run_clock("wg_lsq", TREES / "wg.treefile",
                          DATA / "wg.fasta", DATA / "wg_dates.tsv",
                          reroot="least-squares"))

    # 3. is the clock driven by the poorly covered 2025 batch?
    d = pd.read_csv(DATA / "wg_dates.tsv", sep="\t")
    y2025 = [n for n, dt in zip(d["name"], d["date"].astype(str)) if dt.startswith("2025")]
    if y2025:
        print(f"\n(dropping {len(y2025)} sequences dated 2025: {', '.join(y2025)})")
        t, a, dt = write_subset("wg", y2025, "no2025")
        rows.append(run_clock("wg_no2025", t, a, dt, reroot="least-squares"))

    # 4. does one borderline sequence matter?
    present = set(d["name"])
    marginal = [x for x in DROP_MARGINAL if x in present]
    if marginal:
        t, a, dt = write_subset("wg", marginal, "drop_marginal")
        rows.append(run_clock("wg_drop_marginal", t, a, dt, reroot="least-squares"))

    # 5. E gene, for comparison only
    if (TREES / "eg.treefile").exists():
        rows.append(run_clock("eg_lsq", TREES / "eg.treefile",
                              DATA / "eg.fasta", DATA / "eg_dates.tsv",
                              reroot="least-squares"))

    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "summary.tsv", sep="\t", index=False)

    print("\n\n================ SUMMARY ================")
    print(summary.to_string(index=False))
    print(f"\nWritten -> {OUT / 'summary.tsv'}")
    print("""
How to read this:

  * wg_outgroup vs wg_lsq        similar rate and root -> the root is trustworthy.
                                 Very different -> say so in the paper; do not pick
                                 the one you prefer.
  * wg_lsq vs wg_no2025          similar rate -> the clock is not an artefact of the
                                 2025 batch's missing data. Very different -> the
                                 2025 genomes need masking before they can be dated,
                                 and the dating must be repeated.
  * wg_lsq vs wg_drop_marginal   must be nearly identical. If not, no conclusion may
                                 rest on OZ484513.1.
  * rate                         expect roughly 5e-4 to 1e-3 subs/site/year for DENV.
                                 Far outside that -> something is wrong upstream.
  * R^2                          low R^2 with few sampling years is normal and not
                                 fatal; low R^2 with wide date spread means weak
                                 temporal signal and the dating should be presented
                                 with explicit caution or dropped.
  * flagged tips                 look at every one. A flagged Bangladeshi sequence is
                                 interesting; a flagged 2025 sequence may just be
                                 missing data.

Open results/clock/wg_lsq/*regression* and look at the plot before trusting any number.
""")


if __name__ == "__main__":
    main()
