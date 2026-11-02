"""
13b_compare_runs.py
-------------------
Ask the same question of every robustness run and put the answers side by side.

The question that matters is not "is the tree identical" - it never will be, and
it does not need to be. It is:

    do the 2025 Bangladeshi sequences still fall inside a clade made up of
    Bangladeshi sequences, and is that clade well supported?

plus, as a second check,

    is the clade containing ALL Bangladeshi 3I_A.2 sequences still essentially
    Bangladeshi, or has a different context draw broken it up?

If the answer is the same across the main run, a different root, and two
different context draws, the conclusion does not depend on those choices. If any
run disagrees, that disagreement is a result and belongs in the paper.

Also reported for the outgroup run: whether IQ-TREE found a branch separating the
outgroup at all. With 3II_B it should; with 3I it cannot, because 3I is a
paraphyletic parent label.

Run from the project folder (normally called by 13_robustness.sh):
  python scripts/13b_compare_runs.py --root results/robustness
"""

import argparse
import re
from pathlib import Path

import pandas as pd
from Bio import Phylo

SOLID_ALRT, SOLID_BOOT = 80.0, 95.0


def support_of(clade):
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


def verdict(a, b):
    if a is None or b is None:
        return "no support value"
    return ("solid" if (a >= SOLID_ALRT and b >= SOLID_BOOT)
            else "NOT well supported")


def assess(tree_path, traits_path):
    """Answer both questions for one run."""
    if not Path(tree_path).exists() or not Path(traits_path).exists():
        return None
    tree = Phylo.read(tree_path, "newick")
    tr = pd.read_csv(traits_path, sep="\t")
    leaves = {t.name for t in tree.get_terminals()}

    bd = tr[(tr["country"] == "Bangladesh") & (tr["clade"] == "3I_A.2")]
    bd_names = [n for n in bd["name"] if n in leaves]
    yr_col = None  # dates live in the dates file; use accession prefix-free logic
    # 2025 sequences: read them from the dates file next to the traits file
    dates_path = Path(str(traits_path).replace("_traits.tsv", "_dates.tsv"))
    y2025 = []
    if dates_path.exists():
        d = pd.read_csv(dates_path, sep="\t")
        y2025 = [n for n, dt in zip(d["name"], d["date"].astype(str))
                 if dt.startswith("2025") and n in leaves
                 and n in set(bd["name"])]

    out = {"n_tips": len(leaves), "n_bd_focal": len(bd_names),
           "n_2025": len(y2025)}

    # Q1: the 2025 cluster
    if len(y2025) >= 2:
        node = tree.common_ancestor(y2025)
        tips = [t.name for t in node.get_terminals()]
        foreign = [t for t in tips
                   if t in set(tr["name"])
                   and tr.set_index("name").loc[t, "country"] != "Bangladesh"]
        a, b = support_of(node)
        out.update({
            "q1_clade_size": len(tips),
            "q1_foreign_in_clade": len(foreign),
            "q1_alrt": a, "q1_boot": b, "q1_verdict": verdict(a, b),
        })

    # Q2: all Bangladeshi 3I_A.2
    if len(bd_names) >= 2:
        node = tree.common_ancestor(bd_names)
        tips = [t.name for t in node.get_terminals()]
        idx = tr.set_index("name")
        foreign = [t for t in tips if t in idx.index
                   and idx.loc[t, "country"] != "Bangladesh"]
        a, b = support_of(node)
        out.update({
            "q2_clade_size": len(tips),
            "q2_foreign_in_clade": len(foreign),
            "q2_alrt": a, "q2_boot": b, "q2_verdict": verdict(a, b),
        })
    return out


def outgroup_note(log_path):
    if not Path(log_path).exists():
        return ""
    text = Path(log_path).read_text()
    if "Branch separating outgroup is not found" in text:
        return "outgroup rooting FAILED (paraphyletic outgroup)"
    if re.search(r"^\s*outgroup:", text, flags=re.M):
        return "outgroup rooting applied"
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="results/robustness")
    ap.add_argument("--main-data", default="results/datasets")
    ap.add_argument("--main-trees", default="results/trees")
    ap.add_argument("--sets", nargs="+", default=["wg"])
    args = ap.parse_args()

    root = Path(args.root)
    runs = [("main", Path(args.main_trees), Path(args.main_data), None)]
    for d in sorted(p for p in root.iterdir() if p.is_dir()) if root.exists() else []:
        runs.append((d.name, d / "trees", d / "datasets", d / "trees.log"))

    for stem in args.sets:
        rows = []
        for label, treedir, datadir, log in runs:
            res = assess(treedir / f"{stem}.treefile",
                         datadir / f"{stem}_traits.tsv")
            if res is None:
                print(f"  ({label}/{stem}: tree or traits missing - skipped)")
                continue
            res["run"] = label
            res["note"] = outgroup_note(log) if log else "clock-rooted (main run)"
            rows.append(res)
        if not rows:
            continue
        df = pd.DataFrame(rows).set_index("run")

        print(f"\n======================= {stem.upper()} =======================")
        print("\nDataset composition:")
        print(df[["n_tips", "n_bd_focal", "n_2025", "note"]].to_string())

        q1 = [c for c in df.columns if c.startswith("q1_")]
        if q1:
            print("\nQ1. Do the 2025 Bangladeshi sequences form a supported cluster?")
            print(df[q1].to_string())

        q2 = [c for c in df.columns if c.startswith("q2_")]
        if q2:
            print("\nQ2. Is the clade holding ALL Bangladeshi 3I_A.2 still "
                  "essentially Bangladeshi?")
            print(df[q2].to_string())

        # the actual verdict
        print("\n--- verdict ---")
        for key, label in (("q1_verdict", "2025 cluster"),
                           ("q2_verdict", "whole Bangladeshi clade")):
            if key not in df.columns:
                continue
            vals = df[key].dropna().unique().tolist()
            if len(vals) == 1 and vals[0] == "solid":
                print(f"  {label}: solid in ALL {len(df)} runs -> robust")
            elif all(v == "solid" for v in vals):
                print(f"  {label}: solid in every run that could be assessed")
            else:
                print(f"  {label}: DIFFERS between runs -> {dict(df[key].dropna())}")
                print(f"    Report this. A conclusion that changes with the root "
                      f"or the context draw is not a conclusion.")

    print("""
How to use this:
  * All runs 'solid' -> the result does not depend on the root or on which
    context sequences were drawn. Say so in the paper, with these numbers.
  * The outgroup run is the one that matters most: it is the only check on
    rooting, since the main analysis is clock-rooted out of necessity.
  * Small differences in clade size are expected - a different context draw puts
    different foreign sequences near the Bangladeshi clade. What must not change
    is whether the Bangladeshi sequences group together, and the support.
""")


if __name__ == "__main__":
    main()
