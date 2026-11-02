#!/usr/bin/env bash
#
# 10_confirm_chimera.sh
# ---------------------
# Confirm (or refute) that PQ815901.1 has two evolutionary histories.
#
# The window scan in script 09 showed PQ815901.1 matching the 3III_A consensus
# far better than its own 3I_A.2 lineage in three blocks, and 3I_A.2 everywhere
# else. If that is real, it makes a prediction that can be tested directly:
#
#     a tree built from ONLY the suspect regions should place PQ815901.1 with
#     the genotype III sequences, while a tree built from the REST of the genome
#     places it with the Bangladeshi genotype I sequences.
#
# A sequence cannot genuinely belong in two places, so if both trees put it in
# the expected clade with good support, the two histories are real. This is the
# phylogenetic-incongruence test, and it is the standard way recombination is
# demonstrated. GARD or 3SEQ would quantify breakpoints more precisely, but the
# incongruence test is what makes the case.
#
# Control: the same two trees are built with other Bangladeshi genomes present.
# They should NOT move between the trees. If everything moves, the partition is
# simply too short to resolve anything and the test says nothing.
#
# Input : results/datasets/wg.fasta, wg_traits.tsv
# Output: results/chimera/suspect.treefile, rest.treefile  + a printed comparison
#
# Run from the project folder:  bash scripts/10_confirm_chimera.sh

set -uo pipefail

FOCUS="${FOCUS:-PQ815901.1}"
OUT="results/chimera"
IN="results/datasets/wg.fasta"
TRAITS="results/datasets/wg_traits.tsv"

# The suspect blocks, read off the window scan in script 09. Adjust if the scan
# on your run puts them elsewhere - they are deliberately not hardcoded deep in
# the script.
SUSPECT_RANGES="${SUSPECT_RANGES:-1501-2250,3751-4250,4501-5250}"

[ -f "$IN" ] || { echo "Missing $IN - run scripts/04_select_sequences.py"; exit 1; }
mkdir -p "$OUT"
IQ=$(command -v iqtree2 2>/dev/null || command -v iqtree 2>/dev/null || true)
[ -n "$IQ" ] || { echo "iqtree not found - conda activate denv3"; exit 1; }

echo "Focus sequence : $FOCUS"
echo "Suspect blocks : $SUSPECT_RANGES"

python - "$IN" "$OUT" "$SUSPECT_RANGES" <<'PY'
import sys
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser

aln_path, outdir, ranges = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
blocks = []
for tok in ranges.split(","):
    a, b = tok.split("-")
    blocks.append((int(a), int(b)))

seqs = {}
with open(aln_path) as fh:
    for h, s in SimpleFastaParser(fh):
        seqs[h.split()[0]] = s.upper()
L = len(next(iter(seqs.values())))

keep = [False] * L
for a, b in blocks:
    for i in range(a - 1, min(b, L)):
        keep[i] = True

n_sus = sum(keep)
print(f"  suspect partition : {n_sus} bases")
print(f"  rest partition    : {L - n_sus} bases")

with (outdir / "suspect.fasta").open("w") as f1, (outdir / "rest.fasta").open("w") as f2:
    for name, s in seqs.items():
        sus = "".join(c for c, k in zip(s, keep) if k)
        rest = "".join(c for c, k in zip(s, keep) if not k)
        # skip sequences with almost no data in a partition; they add noise
        if sum(1 for c in sus if c in "ACGT") >= 0.5 * n_sus:
            f1.write(f">{name}\n{sus}\n")
        if sum(1 for c in rest if c in "ACGT") >= 0.5 * (L - n_sus):
            f2.write(f">{name}\n{rest}\n")
PY

for PART in suspect rest; do
  echo
  echo "=============================================================="
  echo "  tree from the $PART partition"
  echo "=============================================================="
  "$IQ" -s "$OUT/$PART.fasta" --prefix "$OUT/$PART" --seqtype DNA \
        -m GTR+F+I+G4 -B 1000 --alrt 1000 -T AUTO --redo \
        > "$OUT/$PART.iqtree.out" 2>&1 \
    && echo "  done" || { echo "  FAILED - see $OUT/$PART.iqtree.out"; tail -12 "$OUT/$PART.iqtree.out"; }
done

echo
echo "=============================================================="
echo "  Where does each sequence sit in the two trees?"
echo "=============================================================="

python - "$OUT" "$TRAITS" "$FOCUS" <<'PY'
import sys
from pathlib import Path
import pandas as pd
from Bio import Phylo

out, traits_path, focus = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
tr = pd.read_csv(traits_path, sep="\t").set_index("name")

def genotype(clade):
    s = str(clade)
    if s.startswith("3III"):
        return "III"
    if s.startswith("3II"):
        return "II"
    if s.startswith("3I"):
        return "I"
    return "?"

def neighbours(tree, name, k=6):
    """The smallest clade around `name` holding at least k other tips."""
    path = tree.get_path(name)
    for node in reversed(path[:-1]):
        tips = [t.name for t in node.get_terminals() if t.name != name]
        if len(tips) >= k:
            sup = node.name if node.name else node.confidence
            return tips, sup
    return [], None

rows = []
trees = {}
for part in ("suspect", "rest"):
    f = out / f"{part}.treefile"
    if f.exists():
        trees[part] = Phylo.read(f, "newick")

if len(trees) < 2:
    sys.exit("  need both trees; check the iqtree output files")

# the focus sequence, plus controls: other Bangladeshi genotype I genomes
controls = [n for n in tr.index
            if tr.at[n, "country"] == "Bangladesh"
            and str(tr.at[n, "clade"]).startswith("3I_")][:8]
watch = [focus] + [c for c in controls if c != focus]

for name in watch:
    rec = {"sequence": name, "own_clade": tr.at[name, "clade"] if name in tr.index else "?"}
    for part, tree in trees.items():
        leaves = {t.name for t in tree.get_terminals()}
        if name not in leaves:
            rec[f"{part}_sits_with"] = "(absent)"
            rec[f"{part}_support"] = None
            continue
        tips, sup = neighbours(tree, name)
        gts = [genotype(tr.at[t, "clade"]) for t in tips if t in tr.index]
        dom = max(set(gts), key=gts.count) if gts else "?"
        frac = gts.count(dom) / len(gts) if gts else 0
        rec[f"{part}_sits_with"] = f"genotype {dom} ({frac:.0%} of {len(gts)} neighbours)"
        rec[f"{part}_support"] = sup
    rows.append(rec)

df = pd.DataFrame(rows)
print(df.to_string(index=False))
df.to_csv(out / "partition_placement.tsv", sep="\t", index=False)

def placed(cell):
    """
    The genotype a sequence sits with, or None if it could not be measured.
    'genotype ?' and '(absent)' mean the neighbour search found nothing - that is
    a failure to measure, NOT a change of genotype. An earlier version of this
    script counted them as changes and so reported controls moving when none did.
    """
    s = str(cell).split("(")[0].strip()
    if not s or s in ("genotype ?", "(absent)", "nan"):
        return None
    return s


f = df[df["sequence"] == focus]
if len(f):
    a_raw = f.iloc[0].get("suspect_sits_with", "")
    b_raw = f.iloc[0].get("rest_sits_with", "")
    a, b = str(a_raw), str(b_raw)
    pa, pb = placed(a_raw), placed(b_raw)
    moved = pa is not None and pb is not None and pa != pb
    others = df[df["sequence"] != focus]
    comparable, others_moved = 0, 0
    for _, r in others.iterrows():
        x, y = placed(r.get("suspect_sits_with")), placed(r.get("rest_sits_with"))
        if x is None or y is None:
            continue                      # unmeasurable - excluded, not counted
        comparable += 1
        if x != y:
            others_moved += 1
    unmeasurable = len(others) - comparable
    if unmeasurable:
        print(f"  ({unmeasurable} control(s) could not be placed in one partition "
              f"and are excluded from the comparison)")
    print(f"\n  {focus}: suspect -> {a}")
    print(f"  {' ' * len(focus)}  rest    -> {b}")
    print(f"  control sequences that also changed genotype: "
          f"{others_moved} of {comparable} comparable")
    if moved and others_moved == 0:
        print("""
  VERDICT: the focus sequence changes genotype between partitions and no control
  does. Two histories in one genome - confirmed incongruence. Exclude it from the
  dated tree, and describe it as a probable recombinant or chimeric assembly
  (say which only if you can rule the other out).""")
    elif moved:
        print("""
  VERDICT: the focus sequence moves, but so do controls. The partitions may be
  too short to resolve placement reliably. Treat as inconclusive and run GARD or
  3SEQ before claiming recombination - but still exclude the sequence from dating,
  because an unexplained long branch and clock-outlier status are reason enough.""")
    else:
        print("""
  VERDICT: the focus sequence does NOT change genotype between partitions. The
  window-scan signal was not recombination. Look instead at sequencing error:
  check its private substitutions and whether they cluster near the shared
  missing-data blocks.""")
PY

echo
echo "Output in $OUT/  (partition_placement.tsv has the table)"
