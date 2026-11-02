#!/usr/bin/env python3
"""
19_position_conservation.py   (v2)
----------------------------------
How conserved is a given amino-acid position across the whole DENV-3 dataset?

Why this exists
---------------
lineage_defining.tsv says NS5:H52R is at 0.029 in the 3I_A context - 4 of 138.
That is a narrow denominator. This script asks the wider question: across all
~9,815 DENV-3 sequences, what does position 52 of NS5 actually look like?

Two very different sentences depend on the answer:
  "fixed in this lineage at a position otherwise invariant across the serotype"
  "fixed in this lineage at a position that varies across DENV-3"

The denominator trap this script is built to avoid
--------------------------------------------------
A sequence that does not COVER a position is not evidence that the position is
unchanged there. Most DENV-3 sequences here are E-gene or prM/E fragments that
never reach NS5. Nextclade reports substitutions, not per-residue coverage, so
a naive carrier count divides by thousands of sequences that could never have
shown it. Frequencies below are over COVERED sequences only; blind sequences
are counted and reported, never treated as reference.

v2 changes
----------
  * Fixed: the reference residue printed as '?' because the set of called
    wild-type residues contained numpy nan alongside the real letter, and
    `if w` lets nan through (nan is truthy). Now filtered on isinstance str.
  * Reference residue is also taken from the modal called wt rather than
    requiring unanimity, with a warning if more than one letter appears.
  * Country join is attempted on several plausible key/column names and says
    which it used, instead of silently producing 'unknown'.

Run:
  python scripts/19_position_conservation.py --positions NS5:52,NS5:139,NS5:404
  python scripts/19_position_conservation.py --positions NS5:52 --by-lineage
"""

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

# 1-based nucleotide positions on Nextclade reference NC_001475.2, from the
# dataset genome_annotation.gff3. Edit if the pinned dataset tag changes.
GENES = {
    "C": (95, 436), "prM": (437, 934), "E": (935, 2413), "NS1": (2414, 3469),
    "NS2A": (3470, 4123), "NS2B": (4124, 4513), "NS3": (4514, 6370),
    "NS4A": (6371, 6751), "2K": (6752, 6820), "NS4B": (6821, 7564),
    "NS5": (7565, 10264),
}

NEXTCLADE_TSV = "results/nextclade/nextclade.tsv"
METADATA = "results/metadata/metadata_clades.tsv"


def codon_range(gene, aa_pos):
    if gene not in GENES:
        sys.exit(f"Unknown gene {gene}. Known: {', '.join(GENES)}")
    start, end = GENES[gene]
    first = start + (aa_pos - 1) * 3
    last = first + 2
    if last > end:
        sys.exit(f"{gene}:{aa_pos} runs past the end of {gene} ({start}-{end})")
    return first, last


def parse_subs(cell):
    out = []
    if not isinstance(cell, str) or not cell.strip():
        return out
    for item in cell.split(","):
        item = item.strip()
        if ":" not in item:
            continue
        gene, _, change = item.partition(":")
        m = re.fullmatch(r"([A-Z\*])(\d+)([A-Z\*])", change)
        if m:
            out.append((gene, m.group(1), int(m.group(2)), m.group(3)))
    return out


def join_country(df, mp):
    """Try several key/column name pairs and report which worked."""
    if not mp.exists():
        print(f"  no metadata at {mp}; skipping country")
        return df
    m = pd.read_csv(mp, sep="\t", low_memory=False)
    key_candidates = ["seqName", "seqname", "accession", "name", "strain", "id"]
    ctry_candidates = ["country", "geo_country", "Country", "location",
                       "geoLocName", "geo_loc_name"]
    key = next((c for c in key_candidates if c in m.columns), None)
    ctry = next((c for c in ctry_candidates if c in m.columns), None)
    if not key or not ctry:
        print(f"  metadata columns available: {list(m.columns)[:15]}")
        print("  could not find a key and a country column; skipping country")
        return df
    sub = m[[key, ctry]].rename(columns={key: "_k", ctry: "country"})
    merged = df.merge(sub, left_on="seqName", right_on="_k", how="left")
    hit = merged["country"].notna().sum()
    if hit == 0:
        # accessions sometimes differ by version suffix
        sub["_k"] = sub["_k"].astype(str).str.split(".").str[0]
        merged = df.copy()
        merged["_base"] = merged["seqName"].astype(str).str.split().str[0] \
                                           .str.split(".").str[0]
        merged = merged.merge(sub, left_on="_base", right_on="_k", how="left")
        hit = merged["country"].notna().sum()
        print(f"  joined country on {key} (version-stripped) -> "
              f"{hit:,} matched")
    else:
        print(f"  joined country on {key}/{ctry} -> {hit:,} matched")
    return merged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--positions", required=True,
                    help="comma-separated, e.g. NS5:52,NS5:139,NS5:404")
    ap.add_argument("--nextclade", default=NEXTCLADE_TSV)
    ap.add_argument("--metadata", default=METADATA)
    ap.add_argument("--by-lineage", action="store_true")
    ap.add_argument("--min-clade-n", type=int, default=10)
    ap.add_argument("--out", default="results/structure/position_conservation.tsv")
    args = ap.parse_args()

    nc_path = Path(args.nextclade)
    if not nc_path.exists():
        sys.exit(f"Not found: {nc_path}")

    print(f"Reading {nc_path} ...")
    head = pd.read_csv(nc_path, sep="\t", nrows=0)
    cols = set(head.columns)
    need = {"seqName", "aaSubstitutions", "alignmentStart", "alignmentEnd"}
    missing = need - cols
    if missing:
        sys.exit(f"Nextclade output is missing: {sorted(missing)}\n"
                 f"Cannot determine coverage, and guessing would make every "
                 f"'absent' sequence look like evidence of conservation. "
                 f"Refusing to proceed.")
    wanted = ["seqName", "aaSubstitutions", "alignmentStart", "alignmentEnd"]
    for extra in ("clade", "qc.overallStatus"):
        if extra in cols:
            wanted.append(extra)
    df = pd.read_csv(nc_path, sep="\t", usecols=wanted, low_memory=False)
    print(f"  {len(df):,} sequences")

    df = join_country(df, Path(args.metadata))

    print("Parsing substitution calls ...")
    df["_subs"] = df["aaSubstitutions"].apply(parse_subs)

    rows = []
    for spec in [s.strip() for s in args.positions.split(",") if s.strip()]:
        gene, _, p = spec.partition(":")
        aa_pos = int(p)
        first, last = codon_range(gene, aa_pos)
        print(f"\n=== {gene}:{aa_pos}  (reference nt {first}-{last}) ===")

        covered = (df["alignmentStart"] <= first) & (df["alignmentEnd"] >= last)
        n_cov = int(covered.sum())
        n_blind = len(df) - n_cov
        print(f"  covered {n_cov:,} | blind {n_blind:,} of {len(df):,}")
        if n_cov == 0:
            print("  nothing covers this position; skipping")
            continue

        sub = df[covered].copy()

        def state(subs, gene=gene, aa_pos=aa_pos):
            for g, wt, pos, mut in subs:
                if g == gene and pos == aa_pos:
                    return mut, wt
            return None, None

        states = sub["_subs"].apply(state)
        sub["mut"] = [s[0] for s in states]
        sub["wt_called"] = [s[1] for s in states]

        # FIX (v2): nan is truthy, so `if w` let it through and wt printed '?'.
        called = [w for w in sub["wt_called"] if isinstance(w, str) and w]
        wt_counts = Counter(called)
        if not wt_counts:
            wt = "?"
            print("  no sequence carries a substitution here, so the reference "
                  "residue is not reported by Nextclade in this column")
        else:
            wt = wt_counts.most_common(1)[0][0]
            if len(wt_counts) > 1:
                print(f"  WARNING: more than one wild-type residue reported "
                      f"{dict(wt_counts)}; using the modal one")

        counts = Counter(sub["mut"].fillna("REF"))
        total = sum(counts.values())
        print(f"  reference residue (Nextclade ref NC_001475.2, genotype III): {wt}")
        print(f"  states across {total:,} covered sequences:")
        for st, n in counts.most_common():
            label = f"{wt} (reference)" if st == "REF" else st
            print(f"    {label:20s} {n:6,}  {100*n/total:6.2f}%")

        ref_n = counts.get("REF", 0)
        major = counts.most_common(1)[0]
        rows.append({
            "position": spec, "ref_residue_genotypeIII": wt,
            "n_covered": total, "n_blind": n_blind,
            "n_reference_state": ref_n,
            "pct_reference_state": round(100 * ref_n / total, 2),
            "majority_state": wt if major[0] == "REF" else major[0],
            "pct_majority_state": round(100 * major[1] / total, 2),
            "n_distinct_states": len(counts),
            "states": "; ".join(f"{(wt if k=='REF' else k)}={v}"
                                for k, v in counts.most_common()),
        })

        if args.by_lineage and "clade" in sub.columns:
            print(f"\n  by clade (>= {args.min_clade_n} covered):")
            for clade, part in sub.groupby("clade"):
                if len(part) < args.min_clade_n:
                    continue
                c = Counter(part["mut"].fillna("REF"))
                desc = ", ".join(
                    f"{(wt if k=='REF' else k)} {100*v/len(part):.0f}%"
                    for k, v in c.most_common(2))
                print(f"    {clade:16s} n={len(part):5,}  {desc}")

        if "country" in sub.columns:
            carriers = sub[sub["mut"].notna()]
            if len(carriers):
                top = Counter(carriers["country"].fillna("unknown")).most_common(8)
                print("\n  carriers by country: "
                      + ", ".join(f"{k} {v}" for k, v in top))
                if "clade" in carriers.columns:
                    print("  carriers by clade:   "
                          + ", ".join(f"{k} {v}" for k, v in
                                      Counter(carriers["clade"]).most_common(8)))

    if rows:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
        print(f"\nWrote {out}")

    print("\nNote: percentages are over COVERED sequences only.")


if __name__ == "__main__":
    main()
