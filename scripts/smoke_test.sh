#!/usr/bin/env bash
#
# smoke_test.sh
# -------------
# Fast check that everything is installed and the file formats are right, before
# committing to a long run. Takes under a minute.
#
# It does NOT touch results/. Everything happens in test_run/, which you can delete.
#
# This is the same idea as the Nextflow `test` profile and the GitHub Actions CI
# that come later in the project: a tiny dataset that exercises every step, so a
# broken tool or a malformed file is caught in seconds instead of after 20 minutes.
#
# Run from the project folder:  bash scripts/smoke_test.sh

set -uo pipefail          # deliberately NOT -e: we want every check to run and report

TEST="test_run"
PASS=0; FAIL=0
pass () { echo "  PASS  $1"; PASS=$((PASS+1)); }
fail () { echo "  FAIL  $1"; FAIL=$((FAIL+1)); }

echo "=============================================="
echo "  1. Are the tools installed?"
echo "=============================================="
for t in python nextclade mafft iqtree treetime; do
  if command -v "$t" >/dev/null 2>&1; then
    v=$("$t" --version 2>&1 | head -1)
    pass "$t  ($v)"
  else
    # iqtree may be installed as iqtree2
    if [ "$t" = "iqtree" ] && command -v iqtree2 >/dev/null 2>&1; then
      pass "iqtree2 ($(iqtree2 --version 2>&1 | head -1))"
    else
      fail "$t not found - check: conda activate denv3"
    fi
  fi
done

echo
echo "  Python packages:"
for m in pandas Bio matplotlib; do
  if python -c "import $m" 2>/dev/null; then pass "import $m"; else fail "import $m"; fi
done

echo
echo "  R packages (only needed for the figures, step 10):"
if command -v Rscript >/dev/null 2>&1; then
  for m in ggtree treeio ape ggplot2; do
    if Rscript -e "suppressMessages(library($m))" >/dev/null 2>&1; then
      pass "R: $m"
    else
      fail "R: $m  (run: mamba env update -f environment.yml)"
    fi
  done
else
  echo "  SKIP  Rscript not installed yet - fine until step 10"
fi

echo
echo "=============================================="
echo "  2. Do the scripts parse? (syntax only)"
echo "=============================================="
for f in scripts/*.py; do
  [ -e "$f" ] || continue
  if python -m py_compile "$f" 2>/dev/null; then pass "$f"; else fail "$f has a syntax error"; fi
done
for f in scripts/*.sh; do
  [ -e "$f" ] || continue
  if bash -n "$f" 2>/dev/null; then pass "$f"; else fail "$f has a syntax error"; fi
done
rm -rf scripts/__pycache__

echo
echo "=============================================="
echo "  3. End-to-end run on a tiny dataset"
echo "=============================================="

if [ ! -f results/datasets/wg.fasta ]; then
  echo "  SKIP  results/datasets/wg.fasta not found."
  echo "        Run scripts/04_select_sequences.py first, then this test."
else
  rm -rf "$TEST"; mkdir -p "$TEST"

  # Take a dozen taxa: every Bangladeshi 2025 and 2023 sequence we can find,
  # plus whatever else comes first, so the test exercises the sequences that
  # actually matter rather than a random dozen.
  python - <<'PY'
import pandas as pd
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser

D = Path("results/datasets"); T = Path("test_run")
traits = pd.read_csv(D / "wg_traits.tsv", sep="\t")
dates = pd.read_csv(D / "wg_dates.tsv", sep="\t")
m = traits.merge(dates, on="name", how="left")
m["year"] = m["date"].astype(str).str.slice(0, 4)

# priority: 2025 Bangladesh, the 2023 bridge genome, some 2019, one outgroup
pick = []
for cond in [(m.country == "Bangladesh") & (m.year == "2025"),
             (m.country == "Bangladesh") & (m.year == "2023"),
             (m.country == "Bangladesh") & (m.year == "2019"),
             (m.clade == "3I_A")]:
    for n in m[cond]["name"].tolist():
        if n not in pick:
            pick.append(n)
        if len(pick) >= 12:
            break
    if len(pick) >= 12:
        break
pick = pick[:12]
keep = set(pick)

with (T / "aln.fasta").open("w") as fh:
    for h, s in SimpleFastaParser((D / "wg.fasta").open()):
        acc = h.split()[0]
        if acc in keep:
            fh.write(f">{acc}\n{s}\n")
sel = m[m["name"].isin(keep)]
sel[["name", "date"]].to_csv(T / "dates.tsv", sep="\t", index=False)
print(f"  test set: {len(sel)} taxa")
print("  " + ", ".join(f"{r['name']}({r['year']})" for _, r in sel.iterrows()))
PY

  n=$(grep -c '^>' "$TEST/aln.fasta" 2>/dev/null || echo 0)
  if [ "$n" -ge 4 ]; then pass "built a test alignment with $n taxa"
  else fail "could not build a test alignment"; fi

  IQ=$(command -v iqtree2 2>/dev/null || command -v iqtree 2>/dev/null || true)
  if [ "$n" -ge 4 ] && [ -z "$IQ" ]; then
    echo "  SKIP  iqtree not installed, so the end-to-end run cannot be tested"
  elif [ "$n" -ge 4 ]; then
    echo
    echo "  running iqtree on $n taxa (fast settings) ..."
    if "$IQ" -s "$TEST/aln.fasta" --prefix "$TEST/test" --seqtype DNA \
         -m GTR+G -B 1000 --alrt 1000 -T 2 --redo > "$TEST/iqtree.out" 2>&1; then
      pass "iqtree produced $TEST/test.treefile"
      head -c 200 "$TEST/test.treefile"; echo " ..."
    else
      fail "iqtree failed - see $TEST/iqtree.out"
      tail -15 "$TEST/iqtree.out"
    fi

    if [ -f "$TEST/test.treefile" ] && command -v treetime >/dev/null 2>&1; then
      echo
      echo "  running treetime clock ..."
      if treetime clock --tree "$TEST/test.treefile" --aln "$TEST/aln.fasta" \
           --dates "$TEST/dates.tsv" --outdir "$TEST/clock" \
           --reroot least-squares > "$TEST/treetime.out" 2>&1; then
        pass "treetime ran; output in $TEST/clock/"
        grep -iE "rate|r\^?2" "$TEST/treetime.out" | head -5
      else
        fail "treetime failed - see $TEST/treetime.out"
        tail -15 "$TEST/treetime.out"
      fi
    fi
  fi
fi

echo
echo "=============================================="
echo "  SUMMARY:  $PASS passed, $FAIL failed"
echo "=============================================="
if [ "$FAIL" -eq 0 ]; then
  echo "  Everything works. Safe to run the full scripts."
else
  echo "  Fix the failures above before the long runs."
fi
echo
echo "  Delete the test output with:  rm -rf $TEST"
exit 0
