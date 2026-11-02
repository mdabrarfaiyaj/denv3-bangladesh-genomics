"""
06b_inspect_trees.py
--------------------
Read the trees and answer the specific questions the project asks, instead of
squinting at a picture.

IQ-TREE was run with both --alrt and -B, so each internal node carries a label
like "97.3/99" meaning SH-aLRT / ultrafast bootstrap. A branch is solid when
SH-aLRT >= 80 AND UFBoot >= 95. Both numbers are printed for every node that
matters, because reporting one of them is how weak branches get oversold.

Questions answered:
  1. Do the 2025 Bangladeshi sequences form one cluster, and with what support?
  2. Where does PQ815901.1 (the only 2023 whole genome) sit relative to the
     2019 and 2025 sequences? This is the whole-genome persistence chain.
  3. Do the exported-looking sequences (USA 2023, Russia 2024, Qatar, etc.)
     fall inside the Bangladeshi clade?
  4. Are the 2021 3III_A sequences one cluster, and how diverse are they?
  5. Is the smallest clade containing all Bangladeshi 3I_A.2 sequences purely
     Bangladeshi, or does it contain foreign sequences too? A purely Bangladeshi
     clade spanning 2017-2025 is what persistence looks like.

Input : results/trees/{wg,eg}.treefile
        results/datasets/{wg,eg}_traits.tsv, {wg,eg}_dates.tsv
Output: printed report; results/trees/<set>_key_nodes.tsv

Run from the project folder:
  python scripts/06b_inspect_trees.py
"""

import sys
from pathlib import Path

import pandas as pd
from Bio import Phylo

TREES = Path("results/trees")
DATA = Path("results/datasets")

EXPORT_CANDIDATES = ["OR821935.1", "PP968117.1", "PX447451.1"]   # USA, Russia, Qatar


def support_of(clade):
    """Return (sh_alrt, ufboot) from an IQ-TREE node label, whatever form it took."""
    raw = clade.name if clade.name else clade.confidence
    if raw is None:
        return None, None
    s = str(raw)
    if "/" in s:
        a, _, b = s.partition("/")
        try:
            return float(a), float(b)
        except ValueError:
            return None, None
    try:
        return None, float(s)
    except ValueError:
        return None, None


def fmt_support(clade):
    a, b = support_of(clade)
    if a is None and b is None:
        return "no support value"
    if a is None:
        return f"UFBoot {b:.0f}"
    verdict = "solid" if (a >= 80 and b >= 95) else "NOT well supported"
    return f"SH-aLRT {a:.1f} / UFBoot {b:.0f}  -> {verdict}"


def tips(clade):
    return [t.name for t in clade.get_terminals()]


def mrca_report(tree, names, label, meta):
    present = [n for n in names if n in meta.index]
    missing = [n for n in names if n not in meta.index]
    found = [n for n in present if n in set(t.name for t in tree.get_terminals())]
    dropped = [n for n in present if n not in found]
    print(f"\n--- {label} ---")
    if dropped:
        print(f"  not in this tree (collapsed as identical, or filtered): {', '.join(dropped)}")
    if missing:
        print(f"  not in the dataset at all: {', '.join(missing)}")
    if len(found) < 1:
        print("  none of these sequences are in this tree")
        return None
    if len(found) == 1:
        print(f"  only one sequence present ({found[0]}); no clade to test")
        return None

    node = tree.common_ancestor(found)
    members = tips(node)
    extra = [m for m in members if m not in found]
    print(f"  sequences tested: {len(found)}")
    print(f"  smallest clade containing them holds {len(members)} tips")
    print(f"  support for that clade: {fmt_support(node)}")
    if not extra:
        print("  the clade contains ONLY these sequences -> they are monophyletic")
    else:
        print(f"  the clade also contains {len(extra)} other tips:")
        info = meta.reindex(extra)
        summary = info.groupby(["country", "clade"], dropna=False).size()
        for (country, cl), n in summary.items():
            print(f"      {n:4d}  {country} / {cl}")
        yrs = pd.to_numeric(info["year"], errors="coerce").dropna()
        if len(yrs):
            print(f"      their years: {int(yrs.min())}-{int(yrs.max())}")
    return node


def run(stem):
    tf = TREES / f"{stem}.treefile"
    if not tf.exists():
        print(f"\n({tf} not found - skipping)")
        return
    tree = Phylo.read(tf, "newick")
    traits = pd.read_csv(DATA / f"{stem}_traits.tsv", sep="\t")
    dates = pd.read_csv(DATA / f"{stem}_dates.tsv", sep="\t")
    meta = traits.merge(dates, on="name", how="left").set_index("name")
    meta["year"] = meta["date"].astype(str).str.slice(0, 4)

    leaf_names = set(t.name for t in tree.get_terminals())
    print("\n" + "=" * 70)
    print(f"  {stem.upper()} TREE: {len(leaf_names)} tips analysed "
          f"(of {len(meta)} in the dataset)")
    print("=" * 70)

    bd = meta[(meta["country"] == "Bangladesh")]
    bd_focal = bd[bd["clade"] == "3I_A.2"]
    bd_inc = bd[bd["clade"] == "3III_A"]

    # 1. the 2025 batch
    y2025 = [n for n in bd_focal.index if str(meta.at[n, "year"]) == "2025"]
    mrca_report(tree, y2025, "The 2025 Bangladeshi sequences", meta)

    # 2. the persistence chain: 2019 + 2023 + 2025
    chain = [n for n in bd_focal.index
             if str(meta.at[n, "year"]) in {"2019", "2023", "2025"}]
    node = mrca_report(tree, chain, "Persistence chain: 2019 + 2023 + 2025 Bangladeshi", meta)

    # where exactly does the 2023 genome sit?
    if "PQ815901.1" in leaf_names:
        print("\n--- PQ815901.1 (the only 2023 whole genome) ---")
        path = tree.get_path("PQ815901.1")
        if len(path) >= 2:
            parent = path[-2]
            sibs = [t for t in tips(parent) if t != "PQ815901.1"]
            print(f"  its immediate neighbours in the tree: {', '.join(sibs) if sibs else '(none)'}")
            if sibs:
                print(meta.reindex(sibs)[["country", "clade", "year"]].to_string())
            print(f"  support for that grouping: {fmt_support(parent)}")

    # 3. all Bangladeshi 3I_A.2: is the clade purely Bangladeshi?
    mrca_report(tree, list(bd_focal.index),
                "All Bangladeshi 3I_A.2 (is the clade purely Bangladeshi?)", meta)

    # 4. exports
    ex = [n for n in EXPORT_CANDIDATES if n in leaf_names]
    if ex:
        print("\n--- Possible exports from Bangladesh ---")
        for acc in ex:
            path = tree.get_path(acc)
            parent = path[-2] if len(path) >= 2 else None
            sibs = [t for t in tips(parent) if t != acc] if parent else []
            bd_sibs = [s for s in sibs if meta.at[s, "country"] == "Bangladesh"] \
                if sibs else []
            country = meta.at[acc, "country"]
            print(f"  {acc} ({country}, {meta.at[acc, 'year']}): "
                  f"{len(sibs)} immediate neighbours, {len(bd_sibs)} of them Bangladeshi")
            if parent is not None:
                print(f"      support: {fmt_support(parent)}")
            if bd_sibs:
                print(f"      -> nests with Bangladeshi sequences: likely an export")

    # 5. the 2021 incursion
    mrca_report(tree, list(bd_inc.index),
                "The 2021 Bangladeshi 3III_A incursion", meta)
    if len(bd_inc):
        inc_in_tree = [n for n in bd_inc.index if n in leaf_names]
        print(f"  {len(inc_in_tree)} of {len(bd_inc)} are present as distinct tips; "
              f"the rest were identical and collapsed.")
        print("  Identical genomes mean one introduction with no accumulated diversity -")
        print("  check metadata.tsv for a shared submitter and submission date before")
        print("  ruling out a single sequencing batch.")

    # key node table for the figures
    rows = []
    for cl in tree.get_nonterminals():
        a, b = support_of(cl)
        t = tips(cl)
        if len(t) < 2:
            continue
        info = meta.reindex(t)
        countries = info["country"].dropna().unique()
        rows.append({
            "n_tips": len(t),
            "sh_alrt": a, "ufboot": b,
            "all_bangladesh": bool(len(countries) == 1 and countries[0] == "Bangladesh"),
            "countries": ";".join(sorted(map(str, countries))[:6]),
            "clades": ";".join(sorted(info["clade"].dropna().unique())[:6]),
        })
    out = pd.DataFrame(rows).sort_values("n_tips", ascending=False)
    out.to_csv(TREES / f"{stem}_key_nodes.tsv", sep="\t", index=False)

    pure = out[(out["all_bangladesh"]) & (out["n_tips"] >= 5)]
    print(f"\n--- Purely Bangladeshi clades with 5+ tips ({len(pure)}) ---")
    if len(pure):
        print(pure[["n_tips", "sh_alrt", "ufboot", "clades"]].to_string(index=False))
    print(f"\nAll internal nodes -> {TREES / f'{stem}_key_nodes.tsv'}")


def main():
    if not (TREES / "wg.treefile").exists() and not (TREES / "eg.treefile").exists():
        sys.exit("No trees found - run scripts/06_build_trees.sh first")
    for stem in ("wg", "eg"):
        run(stem)
    print("""

Reminder on reading support:
  solid            SH-aLRT >= 80 and UFBoot >= 95
  suggestive       one of the two passes
  not supported    neither passes - the grouping is not evidence of anything

The branch that matters most is the one joining the 2025 sequences to the earlier
Bangladeshi ones. If that branch is not solid, the persistence claim has to be
stated as consistent with the data rather than demonstrated by it.
""")


if __name__ == "__main__":
    main()
