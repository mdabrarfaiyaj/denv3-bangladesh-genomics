"""
03_lineage_geography.py
-----------------------
Before selecting sequences, find out where the Bangladeshi lineages live.

Step 2 showed that every Bangladeshi sequence from 2017-2025 is lineage 3I_A.2,
and that 3I_A.2 is mostly Bangladeshi. This script checks that properly:

  - how many 3I_A.2 sequences exist outside Bangladesh, in which countries,
    in which years, and how long they are
  - the same for the sister lineages (3I_A, 3I_A.1, 3I) that a re-introduction
    would most plausibly come from
  - which countries carry genotype I at all, so context sampling is informed
    rather than arbitrary

Nothing is filtered or written as analysis input here. This step only looks.

Input : results/metadata/metadata_clades.tsv
Output: results/metadata/lineage_geography.tsv   (country x year counts for the focal lineages)

Run from the project folder:
  python scripts/03_lineage_geography.py
"""

import sys
from pathlib import Path

import pandas as pd

IN = Path("results/metadata/metadata_clades.tsv")
OUT = Path("results/metadata/lineage_geography.tsv")

# The Bangladeshi lineage, and the lineages it sits inside. A re-introduction of
# something labelled 3I_A.2 would most likely come from one of these.
FOCAL = ["3I_A.2", "3I_A", "3I_A.1", "3I"]


def main():
    if not IN.exists():
        sys.exit(f"Input not found: {IN} (run scripts/02b_summarise_nextclade.py first)")
    df = pd.read_csv(IN, sep="\t", low_memory=False)

    # ---- 1. the focal lineage: Bangladesh vs elsewhere ----
    f = df[df["clade"] == "3I_A.2"]
    inside = (f["country"] == "Bangladesh").sum()
    outside = len(f) - inside
    print(f"== 3I_A.2 worldwide: {len(f)} sequences ==")
    print(f"Bangladesh: {inside}   elsewhere: {outside}")

    print("\nNon-Bangladeshi 3I_A.2, by country:")
    other = f[f["country"] != "Bangladesh"]
    print(other["country"].value_counts(dropna=False).to_string() or "(none)")

    print("\nNon-Bangladeshi 3I_A.2, one per line (this is the re-introduction evidence):")
    cols = [c for c in ["accession", "country", "division", "date", "date_precision",
                        "year", "length", "qc.overallStatus", "totalMissing"]
            if c in other.columns]
    if len(other):
        print(other[cols].sort_values(["year", "country"], na_position="last").to_string(index=False))

    # ---- 2. the wider genotype I group ----
    print("\n\n== Wider genotype I lineages: country x year ==")
    g = df[df["clade"].isin(FOCAL)].dropna(subset=["year"]).copy()
    g["year"] = g["year"].astype(int)
    tab = pd.crosstab([g["clade"], g["country"]], g["year"])
    # keep the printout readable: countries with at least 3 sequences
    tab = tab[tab.sum(axis=1) >= 3]
    print(tab.to_string())

    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.crosstab([g["clade"], g["country"]], g["year"]).to_csv(OUT, sep="\t")
    print(f"\nFull table written -> {OUT}")

    # ---- 3. how much usable whole-genome data exists per lineage ----
    print("\n== Whole genomes (>=10,000 bases) with a date, by lineage and country ==")
    wg = df[(df["length"] >= 10000) & df["clade"].isin(FOCAL) & df["date"].notna()]
    print(pd.crosstab(wg["clade"], wg["country"]).to_string())

    # ---- 4. quality of what we would use ----
    print("\n== QC of genotype I whole genomes ==")
    if "qc.overallStatus" in wg.columns:
        print(pd.crosstab(wg["country"], wg["qc.overallStatus"]).to_string())

    # ---- 5. the 2021 genotype III incursion, for completeness ----
    print("\n== 3III_A (the 2021 Bangladeshi incursion): where else does it occur? ==")
    a = df[df["clade"] == "3III_A"]
    print(a["country"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
