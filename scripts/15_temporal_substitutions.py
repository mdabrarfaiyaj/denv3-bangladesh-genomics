"""
15_temporal_substitutions.py
----------------------------
Did anything actually change in the Bangladeshi lineage between 2017 and 2025?

This is the question that decides whether there is a molecular story at all. It
has to be asked carefully, because the obvious way to ask it is wrong in three
separate ways:

  1. Averaging per-year frequencies treats a year with 1 sequence the same as a
     year with 14. Counts must be POOLED (sum of carriers / sum of sequences
     that covered the position), never averaged.

  2. Coverage is wildly uneven by gene. Most Bangladeshi sequences from
     2018-2024 are E-gene fragments that do not cover NS1, NS3 or NS5 at all.
     So for the NS genes, "early vs late" is really "the 2019 whole-genome batch
     vs the 2025 whole-genome batch" - two laboratories, two pipelines, and
     9-18% missing data in the later one. That is a batch comparison wearing the
     costume of a time series.

  3. A difference of 13/14 vs 2/6 looks dramatic and is not significant once you
     test it. Every candidate gets Fisher's exact test and a Benjamini-Hochberg
     correction, because ~900 substitutions are being screened.

The direction also matters and is reported. A substitution that DISAPPEARS by
2025 means the later viruses match the (genotype III) reference where earlier
ones did not. Mass reversion to an ancestral state is not plausible biology;
it points to assembly differences or a distinct sub-clade. Only substitutions
that RISE are candidates for adaptation.

Fisher's exact test is implemented here directly so the script needs no scipy.

Input : results/mutations/bd_by_year.tsv
Output: results/mutations/temporal_test.tsv
        printed report

Run from the project folder:
  python scripts/15_temporal_substitutions.py
  python scripts/15_temporal_substitutions.py --early-max 2019 --late-min 2022
"""

import argparse
import sys
from math import comb
from pathlib import Path

import pandas as pd

IN = Path("results/mutations/bd_by_year.tsv")
OUT = Path("results/mutations/temporal_test.tsv")


def fisher_exact_two_sided(a, b, c, d):
    """
    Two-sided Fisher's exact test on [[a, b], [c, d]], by direct enumeration.
    a = early carriers, b = early non-carriers, c = late carriers, d = late non.
    Exact and dependency-free; the tables here are small.
    """
    n = a + b + c + d
    row1, row2 = a + b, c + d
    col1 = a + c

    def p_table(x):
        return (comb(row1, x) * comb(row2, col1 - x)) / comb(n, col1)

    lo = max(0, col1 - row2)
    hi = min(row1, col1)
    p_obs = p_table(a)
    # include every table at least as extreme as the observed one
    tol = 1e-12
    return min(1.0, sum(p_table(x) for x in range(lo, hi + 1)
                        if p_table(x) <= p_obs + tol))


def bh_adjust(pvals):
    """Benjamini-Hochberg FDR."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = n - rank + 1
        val = min(prev, pvals[i] * n / k)
        adj[i] = val
        prev = val
    return adj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--early-max", type=int, default=2019)
    ap.add_argument("--late-min", type=int, default=2022)
    ap.add_argument("--min-cov", type=int, default=5,
                    help="minimum sequences covering the position on EACH side")
    ap.add_argument("--alpha", type=float, default=0.05)
    args = ap.parse_args()

    if not IN.exists():
        sys.exit(f"Input not found: {IN} (run scripts/05_mutation_profile.py first)")
    d = pd.read_csv(IN, sep="\t")

    years = sorted(int(c) for c in d.columns if c.isdigit())
    early = [y for y in years if y <= args.early_max]
    late = [y for y in years if y >= args.late_min]
    print(f"early years: {early}")
    print(f"late  years: {late}")
    if not early or not late:
        sys.exit("need years on both sides of the split")

    # POOLED counts, not averaged frequencies
    d["early_n"] = d[[f"{y}_n" for y in early]].sum(axis=1)
    d["early_of"] = d[[f"{y}_of" for y in early]].sum(axis=1)
    d["late_n"] = d[[f"{y}_n" for y in late]].sum(axis=1)
    d["late_of"] = d[[f"{y}_of" for y in late]].sum(axis=1)

    print(f"\nPooled coverage across all substitutions:")
    print(f"  early: median {int(d['early_of'].median())} sequences per position "
          f"(max {int(d['early_of'].max())})")
    print(f"  late : median {int(d['late_of'].median())} sequences per position "
          f"(max {int(d['late_of'].max())})")

    ok = d[(d["early_of"] >= args.min_cov) & (d["late_of"] >= args.min_cov)].copy()
    print(f"\nSubstitutions with >= {args.min_cov} covering sequences on both "
          f"sides: {len(ok)} of {len(d)}")
    if ok.empty:
        print("\nNothing passes the coverage filter. With this data the question "
              "cannot be answered for any position - which is itself the answer.")
        return

    ok["early_freq"] = (ok["early_n"] / ok["early_of"]).round(3)
    ok["late_freq"] = (ok["late_n"] / ok["late_of"]).round(3)
    ok["direction"] = [("rose" if l > e else "fell" if l < e else "flat")
                       for e, l in zip(ok["early_freq"], ok["late_freq"])]

    ok["p"] = [fisher_exact_two_sided(int(en), int(eo - en), int(ln), int(lo - ln))
               for en, eo, ln, lo in zip(ok["early_n"], ok["early_of"],
                                         ok["late_n"], ok["late_of"])]
    ok["p_adj"] = bh_adjust(ok["p"].tolist())
    ok = ok.sort_values("p_adj")

    cols = ["substitution", "gene", "position", "early_n", "early_of",
            "early_freq", "late_n", "late_of", "late_freq", "direction",
            "p", "p_adj"]
    ok[cols].to_csv(OUT, sep="\t", index=False)

    sig = ok[ok["p_adj"] < args.alpha]
    rose = sig[sig["direction"] == "rose"]
    fell = sig[sig["direction"] == "fell"]

    print(f"\n{'=' * 70}")
    print(f"  Significant after multiple-testing correction (FDR < {args.alpha}):"
          f" {len(sig)}")
    print(f"{'=' * 70}")
    if len(sig):
        print(ok.loc[sig.index, cols].head(30).to_string(index=False))
    else:
        print("  none")

    print(f"\n  of those, ROSE by {args.late_min}+: {len(rose)}  "
          f"(candidates for adaptation)")
    print(f"           FELL by {args.late_min}+: {len(fell)}  "
          f"(later viruses match the reference where earlier ones did not)")

    print("""
How to read this:

  * Nothing significant -> the lineage shows no detectable change across the
    period. That is a real finding, and it argues the 2025-26 surge was driven
    by host population immunity rather than by viral evolution.

  * Substitutions that FELL are not evidence of adaptation. They mean the later
    sequences match the genotype III reference where the earlier ones did not.
    Across many linked positions at once, that points to assembly differences
    between the two sequencing batches, or to a distinct sub-clade - not to the
    virus reverting.

  * Substitutions that ROSE and survive correction are the only candidates
    worth putting on a structure. Before believing any of them, check that the
    position is covered in BOTH sequencing batches from 2025 (the whole genomes
    and the prM+E amplicons), using results/audit/q3_2025_cross_batch.tsv. A
    change seen in only one batch is a pipeline, not a virus.

  * Whatever the result: for the NS genes this is a comparison of the 2019
    whole-genome batch against the 2025 whole-genome batch. Two laboratories,
    two pipelines. Say so in the paper; it is the main limitation of any
    temporal claim made from this dataset.
""")


if __name__ == "__main__":
    main()
