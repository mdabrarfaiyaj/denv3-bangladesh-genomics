"""
02b_summarise_nextclade.py
--------------------------
Join the Nextclade results to the metadata table and report what we have.
Nothing is filtered here. This step only answers: which lineages are present,
how good is the sequence quality, and what does Bangladesh look like by year?

Input : results/nextclade/nextclade.tsv
        results/metadata/metadata.tsv
Output: results/metadata/metadata_clades.tsv   (metadata + lineage + QC per sequence)

Run from the project folder:
  python scripts/02b_summarise_nextclade.py
"""

import sys
from pathlib import Path

import pandas as pd

NEXTCLADE = Path("results/nextclade/nextclade.tsv")
METADATA = Path("results/metadata/metadata.tsv")
OUT = Path("results/metadata/metadata_clades.tsv")


def main():
    for p in (NEXTCLADE, METADATA):
        if not p.exists():
            sys.exit(f"Input not found: {p}")

    nc = pd.read_csv(NEXTCLADE, sep="\t", low_memory=False)
    md = pd.read_csv(METADATA, sep="\t", low_memory=False)

    # Nextclade's seqName is the full FASTA header; the accession is its first field.
    nc["accession"] = nc["seqName"].astype(str).str.split().str[0]

    keep = [c for c in [
        "accession", "clade", "qc.overallScore", "qc.overallStatus",
        "alignmentStart", "alignmentEnd", "totalMissing", "totalSubstitutions",
        "totalNonACGTNs", "qc.frameShifts.status", "qc.stopCodons.status",
    ] if c in nc.columns]
    nc = nc[keep]

    df = md.merge(nc, on="accession", how="left", validate="one_to_one")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)

    # ---------------- report ----------------
    print(f"Sequences in metadata : {len(md)}")
    print(f"Rows from Nextclade   : {len(nc)}")
    unassigned = df["clade"].isna().sum() if "clade" in df else len(df)
    print(f"No lineage assigned   : {unassigned}")
    print(f"Written -> {OUT}")

    if "qc.overallStatus" in df:
        print("\n== QC status, all sequences ==")
        print(df["qc.overallStatus"].value_counts(dropna=False).to_string())

    print("\n== Top 15 lineages worldwide ==")
    print(df["clade"].value_counts().head(15).to_string())

    bd = df[df["country"] == "Bangladesh"].copy()
    print(f"\n== Bangladesh: {len(bd)} sequences ==")
    print("\nLineages found in Bangladesh (any year):")
    print(bd["clade"].value_counts(dropna=False).to_string())

    print("\nBangladesh lineage by year (rows = year, columns = lineage):")
    bdy = bd.dropna(subset=["year"])
    if not bdy.empty:
        tab = pd.crosstab(bdy["year"].astype(int), bdy["clade"].fillna("unassigned"))
        print(tab.to_string())

    print("\n== The 2025 Bangladeshi sequences, one per line ==")
    cols = [c for c in ["accession", "date", "division", "length", "clade",
                        "qc.overallStatus", "totalMissing"] if c in bd.columns]
    print(bd[bd["year"] == 2025][cols].sort_values("length", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
