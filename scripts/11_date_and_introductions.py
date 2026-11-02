"""
11_date_and_introductions.py
----------------------------
Put dates on the tree, and count introductions into Bangladesh.

Two questions, in order:

  1. WHEN did this lineage arrive in Bangladesh?
     TreeTime converts branch lengths into calendar time using the clock
     established in step 7 (7.57e-4 subs/site/year, R^2 0.96). The number that
     matters is the date of the most recent common ancestor (tMRCA) of the
     Bangladeshi sequences, WITH its confidence interval. A tMRCA without an
     interval is not a result - the interval is usually wide with 20 genomes.

  2. HOW MANY separate introductions does the tree support?
     TreeTime's mugration model treats location as a character evolving along
     the tree, infers a location for every internal node, and then an
     introduction is simply a branch whose parent is inferred elsewhere and
     whose child is inferred to be in Bangladesh.

Both are run at two levels of location detail, because they fail differently:
  country  - fine-grained, but many countries have one sequence, which makes
             the model confident about things it cannot know
  region   - coarse (Bangladesh / South Asia other / Southeast Asia / East Asia
             / elsewhere), more trustworthy
Report both. If they disagree, the region-level answer is the safer one.

Honesty constraints built in:
  - the clock-filtered outliers stay as tips but TreeTime excludes them from
    rate estimation (--clock-filter 3), as in step 7
  - --confidence gives intervals on every node date; nodes are reported with them
  - ambiguous dates (YYYY-XX-XX) are handled by TreeTime as ranges, not guesses
  - introductions are counted from the inferred states, and the per-node
    probability is reported, so a weakly supported introduction can be seen as such

Inputs : results/trees/wg.treefile      (and eg.treefile)
         results/datasets/wg.fasta, wg_dates.tsv, wg_traits.tsv
Outputs: results/dating/<set>/                TreeTime's own output
         results/dating/<set>/mugration_<attr>/
         results/dating/summary.tsv           tMRCA and introduction counts

Run from the project folder:
  python scripts/11_date_and_introductions.py
  python scripts/11_date_and_introductions.py --sets wg --attrs region
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import Phylo

TREES = Path("results/trees")
DATA = Path("results/datasets")
OUT = Path("results/dating")

HOME = "Bangladesh"


def run(cmd, log_path):
    print("  " + " ".join(str(c) for c in cmd))
    res = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    log_path.write_text(res.stdout + "\n==== stderr ====\n" + res.stderr)
    if res.returncode != 0:
        print(f"  FAILED (exit {res.returncode}); last lines:")
        for line in (res.stderr or res.stdout).strip().splitlines()[-12:]:
            print("    " + line)
    return res


def strip_annotations(nexus_path, out_newick):
    """
    Pull the tree string out of a nexus file and remove [&key=value] comments,
    so Bio.Phylo can read the topology. The node states are recovered separately
    from TreeTime's confidence.csv, which is more reliable than parsing comments.
    """
    text = nexus_path.read_text()
    # TreeTime writes 'Tree tree1=(...)' - capital T, no spaces around '=', and
    # no [&R] prefix. The first version of this regex was case-sensitive and
    # silently matched nothing, which looked like a missing file.
    # The tree name must not be allowed to contain '=' or '(': the tree string
    # itself has no whitespace, so a greedy \S* swallows the whole tree and then
    # backtracks to the LAST '=' in the line (inside an annotation such as
    # date=2011.5), capturing only the tail. Restrict the name instead.
    # Leading whitespace is [ \t] only, so the match cannot start on a prior line.
    m = re.search(r"^[ \t]*tree\s+[^=\s(]*[ \t]*=[ \t]*(?:\[&[^\]]*\][ \t]*)?(.+?);[ \t]*$",
                  text, flags=re.M | re.S | re.I)
    if not m:
        # say what was actually found, instead of reporting nothing at all
        head = next((l for l in text.splitlines()
                     if l.strip().lower().startswith("tree")), None)
        print(f"    {nexus_path.name}: found the file but could not parse a tree line.")
        if head:
            print(f"    the line looks like: {head[:160]}")
        return None
    tree_str = re.sub(r"\[&[^\]]*\]", "", m.group(1))
    out_newick.write_text(tree_str.strip() + ";\n")
    return out_newick


def read_tree_any(outdir, tag):
    """TreeTime's tree output, whatever it called it, as a Bio.Phylo tree."""
    for name in (f"{tag}.nexus", f"{tag}.nwk", f"{tag}.newick"):
        p = outdir / name
        if not p.exists():
            continue
        if p.suffix == ".nexus":
            clean = strip_annotations(p, outdir / f"{tag}.clean.nwk")
            if clean:
                return Phylo.read(clean, "newick"), p
        else:
            return Phylo.read(p, "newick"), p
    return None, None


# ---------------------------------------------------------------------------
# Node identity, done properly.
#
# TreeTime carries IQ-TREE's support labels through as node NAMES, so a tree can
# contain dozens of nodes all called "100/100". Any lookup keyed on node name is
# therefore ambiguous and silently returns the wrong node. The first version of
# this script did exactly that: it reported a tMRCA of 2018.34 for a clade whose
# oldest tip is from mid-2017 (impossible), and zero introductions (also
# impossible). Nothing warned, because every lookup "succeeded".
#
# So nodes are identified by their TOPOLOGY - the set of tips below them - and
# their attributes are read from the same structure that defines them:
#   dates        -> auspice_tree.json, which nests num_date with each node
#   mugration    -> the annotated newick, parsed with its [&key=value] comments
# ---------------------------------------------------------------------------

def load_auspice(path):
    """auspice_tree.json -> nested node dicts with tips, date and interval."""
    import json
    raw = json.loads(path.read_text())
    root = raw.get("tree", raw)

    def walk(n):
        attrs = n.get("node_attrs", {})
        nd = attrs.get("num_date", {})
        ci = nd.get("confidence") or [None, None]
        node = {
            "name": n.get("name", ""),
            "date": nd.get("value"),
            "lower": ci[0], "upper": ci[1],
            "bad": (attrs.get("bad_branch", {}) or {}).get("value") == "Yes",
            "children": [walk(c) for c in n.get("children", []) or []],
        }
        node["tips"] = (frozenset([node["name"]]) if not node["children"]
                        else frozenset().union(*[c["tips"] for c in node["children"]]))
        return node

    return walk(root)


def mrca_by_topology(node, targets):
    """
    Deepest node whose tip set contains every target tip. Identified by topology,
    so duplicate node names cannot confuse it.
    """
    if not targets <= node["tips"]:
        return None
    for c in node["children"]:
        deeper = mrca_by_topology(c, targets)
        if deeper is not None:
            return deeper
    return node


def _split_annotation(s):
    """Split 'a="x,y",b=1' on commas that are not inside quotes."""
    out, buf, q = [], [], None
    for ch in s:
        if q:
            if ch == q:
                q = None
            buf.append(ch)
        elif ch in "\"'":
            q = ch
            buf.append(ch)
        elif ch == ",":
            out.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    if buf:
        out.append("".join(buf))
    d = {}
    for item in out:
        if "=" in item:
            k, _, v = item.partition("=")
            d[k.strip().lstrip("&")] = v.strip().strip("\"'")
    return d


class _NewickParser:
    """Minimal newick parser that KEEPS [&key=value] annotations per node."""

    def __init__(self, s):
        self.s, self.i = s, 0

    def _peek(self):
        return self.s[self.i] if self.i < len(self.s) else ""

    def parse(self):
        node = self._node()
        return node

    def _node(self):
        node = {"name": "", "ann": {}, "children": []}
        if self._peek() == "(":
            self.i += 1
            while True:
                node["children"].append(self._node())
                if self._peek() == ",":
                    self.i += 1
                    continue
                if self._peek() == ")":
                    self.i += 1
                break
        start = self.i
        while self._peek() and self._peek() not in "(),:;[":
            self.i += 1
        node["name"] = self.s[start:self.i].strip().strip("'\"")
        # annotations and branch lengths can appear in either order, repeatedly.
        # NOTE the explicit truthiness test: _peek() returns "" at end of input,
        # and "" in "[:" is True in Python, so without it this loop never ends.
        while self._peek() and self._peek() in "[:":
            if self._peek() == "[":
                j = self.s.find("]", self.i)
                if j < 0:
                    break
                node["ann"].update(_split_annotation(self.s[self.i + 1:j]))
                self.i = j + 1
            else:
                self.i += 1
                while self._peek() and self._peek() not in "(),;[":
                    self.i += 1
        node["tips"] = (frozenset([node["name"]]) if not node["children"]
                        else frozenset().union(*[c["tips"] for c in node["children"]]))
        return node


def parse_annotated_tree(path):
    """Read a nexus/newick tree keeping its [&...] annotations."""
    text = path.read_text()
    m = re.search(r"^[ \t]*tree\s+[^=\s(]*[ \t]*=[ \t]*(?:\[&[^\]]*\][ \t]*)?(.+?);[ \t]*$",
                  text, flags=re.M | re.S | re.I)
    s = (m.group(1) if m else text.strip().rstrip(";")).strip()
    return _NewickParser(s).parse()


def read_treetime_table(path, sep="\t"):
    """
    TreeTime writes dates.tsv with comment lines AND a header that itself starts
    with '#' (e.g. '#node<TAB>numeric date<TAB>...'). Passing comment='#' to
    pandas therefore eats the header and promotes the first data row to column
    names. Find the real header explicitly instead.
    """
    lines = path.read_text().splitlines()
    header_idx = None
    for i, line in enumerate(lines):
        low = line.lower()
        if ("numeric" in low and "date" in low) or low.lstrip("#").startswith("node"):
            header_idx = i
            break
    if header_idx is None:
        return pd.read_csv(path, sep=sep)
    header = [h.strip().lstrip("#").strip()
              for h in lines[header_idx].split(sep)]
    body = [l for l in lines[header_idx + 1:]
            if l.strip() and not l.lstrip().startswith("#")]
    rows = [l.split(sep) for l in body]
    rows = [r + [""] * (len(header) - len(r)) for r in rows]
    df = pd.DataFrame([r[:len(header)] for r in rows], columns=header)
    # convert what is numeric, leave the rest as text. pandas removed
    # errors="ignore", so do it per column instead of relying on it.
    for c in df.columns[1:]:
        conv = pd.to_numeric(df[c], errors="coerce")
        if conv.notna().any():
            df[c] = conv
    return df


def node_states_from_confidence(path):
    """
    confidence.csv from `treetime mugration`: one row per node, one column per
    state, values are probabilities. Returns {node: (best_state, probability)}.
    """
    if not path.exists():
        return {}
    df = read_treetime_table(path, sep=",")
    name_col = df.columns[0]
    states = [c for c in df.columns[1:]]
    out = {}
    for _, r in df.iterrows():
        vals = {s: pd.to_numeric(r[s], errors="coerce") for s in states}
        vals = {k: v for k, v in vals.items() if v == v}
        if not vals:
            continue
        best = max(vals, key=vals.get)
        out[str(r[name_col])] = (best, float(vals[best]))
    return out


def tmrca_of(auspice_root, names):
    """
    Date of the MRCA of `names`, with its interval, taken from the auspice tree
    where each node carries its own num_date. The node is found by topology, so
    duplicate support-label names are irrelevant.
    """
    targets = frozenset(names) & auspice_root["tips"]
    if len(targets) < 2:
        return None
    node = mrca_by_topology(auspice_root, targets)
    if node is None:
        return None
    return {
        "n_tips_in_clade": len(node["tips"]),
        "n_targets_found": len(targets),
        "date": node["date"],
        "lower": node["lower"],
        "upper": node["upper"],
        "clade_is_pure": targets == node["tips"],
    }


def count_introductions(root, attr, home=HOME):
    """
    A branch whose parent is inferred elsewhere and whose child is inferred to be
    in `home` is one introduction. States are read from the node's own inline
    annotation, so nothing depends on node names being unique.
    """
    intros = []

    def conf_of(node):
        for k, v in node["ann"].items():
            if k.lower().startswith(attr.lower()) and "conf" in k.lower():
                try:
                    return float(v)
                except ValueError:
                    return None
        return None

    def walk(node):
        p = node["ann"].get(attr)
        for child in node["children"]:
            c = child["ann"].get(attr)
            if p and c and p != home and c == home:
                intros.append({
                    "parent_state": p,
                    "child_state": c,
                    "child_name": child["name"],
                    "child_prob": conf_of(child),
                    "child_is_tip": not child["children"],
                    "tips_below": len(child["tips"]),
                    "example_tips": ";".join(sorted(child["tips"])[:4]),
                })
            walk(child)

    walk(root)
    return pd.DataFrame(intros)


def do_set(stem, attrs, args, rows):
    tf = TREES / f"{stem}.treefile"
    for p in (tf, DATA / f"{stem}.fasta", DATA / f"{stem}_dates.tsv",
              DATA / f"{stem}_traits.tsv"):
        if not p.exists():
            print(f"  missing {p} - skipping {stem}")
            return

    outdir = OUT / stem
    outdir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 70)
    print(f"  DATING: {stem}")
    print("=" * 70)
    run(["treetime",
         "--tree", tf,
         "--aln", DATA / f"{stem}.fasta",
         "--dates", DATA / f"{stem}_dates.tsv",
         "--outdir", outdir,
         "--reroot", "least-squares",
         "--clock-filter", "3",
         "--confidence",
         "--covariation"], outdir / "treetime.log")

    # dates come from the auspice json, where each node carries its own date;
    # the nexus is kept only as the input to mugration
    auspice = outdir / "auspice_tree.json"
    nexus = outdir / "timetree.nexus"
    if not auspice.exists() or not nexus.exists():
        print(f"  missing auspice_tree.json or timetree.nexus in {outdir}; "
              f"see treetime.log")
        return
    aus = load_auspice(auspice)
    print(f"  timetree read from auspice_tree.json "
          f"({len(aus['tips'])} tips, root {aus['date']:.1f})")

    traits = pd.read_csv(DATA / f"{stem}_traits.tsv", sep="\t")
    bd_focal = traits[(traits["country"] == HOME) &
                      (traits["clade"] == "3I_A.2")]["name"].tolist()
    bd_inc = traits[(traits["country"] == HOME) &
                    (traits["clade"] == "3III_A")]["name"].tolist()

    # tip dates, for the impossibility check below
    dates_lookup = {}

    def collect(node):
        if not node["children"] and node["date"] is not None:
            dates_lookup[node["name"]] = node["date"]
        for c in node["children"]:
            collect(c)

    collect(aus)

    for label, names in (("Bangladeshi 3I_A.2 (the persisting lineage)", bd_focal),
                         ("Bangladeshi 3III_A (the 2021 incursion)", bd_inc)):
        t = tmrca_of(aus, names)
        print(f"\n  tMRCA, {label}")
        if not t:
            print("    could not be determined")
            continue
        ci = (f" [{t['lower']:.2f} - {t['upper']:.2f}]"
              if t["lower"] is not None and t["upper"] is not None else "")
        print(f"    {t['date']:.2f}{ci}")
        print(f"    clade holds {t['n_tips_in_clade']} tips, of which "
              f"{t['n_targets_found']} are the target sequences"
              f"{' (clade is purely them)' if t['clade_is_pure'] else ''}")
        print("    Report the interval, never the point estimate alone.")
        # sanity check: an MRCA cannot be younger than the tips below it
        if t["date"] is not None:
            tip_dates = []
            for n in names:
                d = dates_lookup.get(n)
                if d is not None:
                    tip_dates.append(d)
            if tip_dates and t["date"] > min(tip_dates) + 0.01:
                print(f"    WARNING: this tMRCA ({t['date']:.2f}) is LATER than the "
                      f"oldest target tip ({min(tip_dates):.2f}), which is impossible. "
                      f"Do not use this number.")
        rows.append({"dataset": stem, "quantity": f"tMRCA {label}",
                     "value": t["date"], "lower": t["lower"], "upper": t["upper"],
                     "n": t["n_tips_in_clade"]})

    # ---- introductions, at each level of location detail ----
    for attr in attrs:
        if attr not in traits.columns:
            print(f"\n  traits file has no '{attr}' column - skipping")
            continue
        mdir = outdir / f"mugration_{attr}"
        mdir.mkdir(parents=True, exist_ok=True)
        states_file = mdir / "states.csv"
        traits[["name", attr]].rename(columns={"name": "name"}).to_csv(
            states_file, index=False)

        print(f"\n  mugration on '{attr}'")
        run(["treetime", "mugration",
             "--tree", nexus,
             "--states", states_file,
             "--attribute", attr,
             "--confidence",
             "--outdir", mdir], mdir / "mugration.log")

        ann_path = mdir / "annotated_tree.nexus"
        if not ann_path.exists():
            print(f"    no annotated_tree.nexus in {mdir}; see mugration.log")
            continue
        mroot = parse_annotated_tree(ann_path)

        keys = sorted({k for n in [mroot] for k in n["ann"]})
        print(f"    annotations on the root node: {keys if keys else '(none)'}")
        if attr not in keys:
            # report what the attribute is actually called rather than silently
            # counting zero introductions
            def any_keys(node, acc):
                acc.update(node["ann"].keys())
                for c in node["children"]:
                    any_keys(c, acc)
                return acc
            found = sorted(any_keys(mroot, set()))
            print(f"    '{attr}' not found as an annotation. Keys present: {found}")
            cand = [k for k in found if k.lower() == attr.lower()]
            if not cand:
                continue
            attr = cand[0]
            print(f"    using '{attr}' instead")

        intros = count_introductions(mroot, attr)

        # Some TreeTime versions write only the state inline, not its confidence.
        # Fill the probability from confidence.csv, but ONLY for nodes whose name
        # is unique in that file - internal nodes often carry IQ-TREE support
        # labels like "100/100", which are not unique and must not be matched.
        conf = mdir / "confidence.csv"
        if conf.exists() and len(intros):
            cdf = read_treetime_table(conf, sep=",")
            counts = cdf.iloc[:, 0].astype(str).value_counts()
            unique_names = set(counts[counts == 1].index)
            lookup = {}
            for _, r in cdf.iterrows():
                nm = str(r.iloc[0])
                if nm not in unique_names:
                    continue
                vals = {c: pd.to_numeric(r[c], errors="coerce")
                        for c in cdf.columns[1:]}
                vals = {k: v for k, v in vals.items() if v == v}
                if vals:
                    lookup[nm] = max(vals.values())
            intros["child_prob"] = [
                lookup.get(str(n)) if p is None or pd.isna(p) else p
                for n, p in zip(intros.get("child_name", intros["example_tips"]),
                                intros["child_prob"])]

        intros.to_csv(mdir / "introductions.tsv", sep="\t", index=False)
        probs = (pd.to_numeric(intros["child_prob"], errors="coerce")
                 if len(intros) else pd.Series(dtype=float))
        n_known = int(probs.notna().sum())
        n_sure = int((probs >= 0.9).sum())
        if n_known == 0:
            print(f"    introductions into {HOME}: {len(intros)} total "
                  f"(per-node probabilities NOT available - this TreeTime version "
                  f"writes the state inline but not its confidence, so do not read "
                  f"this as weak support)")
        else:
            print(f"    introductions into {HOME}: {len(intros)} total, "
                  f"{n_sure} of {n_known} with probability >= 0.9")
        if len(intros):
            print(intros.sort_values("tips_below", ascending=False)
                  .head(10).to_string(index=False))
        else:
            print("    Zero introductions means the model inferred Bangladesh at "
                  "every ancestor of every Bangladeshi tip, INCLUDING the root - "
                  "which usually means the root state is wrong, not that there "
                  "was no introduction. Check the root annotation above.")
        rows.append({"dataset": stem, "quantity": f"introductions ({attr})",
                     "value": len(intros), "lower": None, "upper": None,
                     "n": n_sure})

        gtr = mdir / "GTR.txt"
        if gtr.exists():
            print(f"    migration model written to {gtr}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", nargs="+", default=["wg", "eg"])
    ap.add_argument("--attrs", nargs="+", default=["region", "country"])
    args = ap.parse_args()

    if not shutil.which("treetime"):
        sys.exit("treetime not found. conda activate denv3")
    OUT.mkdir(parents=True, exist_ok=True)

    rows = []
    for stem in args.sets:
        do_set(stem, args.attrs, args, rows)

    if rows:
        s = pd.DataFrame(rows)
        s.to_csv(OUT / "summary.tsv", sep="\t", index=False)
        print("\n\n================ SUMMARY ================")
        print(s.to_string(index=False))
        print(f"\nWritten -> {OUT / 'summary.tsv'}")

    print("""
Reading this:

  * The tMRCA of the Bangladeshi 3I_A.2 clade is the arrival date of the
    persisting lineage. Expect something around 2015-2017 given that Malaysian
    and Chinese 3I_A.2 sequences date from 2015. A confidence interval spanning
    several years is normal with 20 genomes - report it, do not hide it.

  * Introductions counted at 'region' are more trustworthy than at 'country'.
    Many countries contribute one sequence, which makes the country model
    assert sources it cannot really distinguish.

  * One introduction with most of the Bangladeshi tips below it is the
    signature of a single establishment followed by local spread - which is
    what persistence means. Several introductions each with few tips below
    would mean repeated seeding instead.

  * The 3III_A tMRCA will be very recent and its interval tight, because those
    eighteen genomes are identical. That is an artefact of one sequencing batch
    (decisions log D2), not evidence of an explosive outbreak.

  * Compare the wg and eg answers. The eg tree has far more Bangladeshi tips but
    much weaker support and mostly year-only dates, so treat its tMRCA as a
    sanity check on the whole-genome estimate rather than as an independent one.
""")


if __name__ == "__main__":
    main()
