#!/usr/bin/env bash
#
# 13_robustness.sh
# ----------------
# The three robustness checks that still stand between the analysis and the paper.
#
#   outgroup   Re-root the tree properly. The main analysis is clock-rooted
#              because lineage 3I is a paraphyletic PARENT label and no single
#              branch can separate it. 3II_B is a defined sister lineage, so it
#              is a genuine outgroup. If the clock root and the outgroup root
#              agree, the root is trustworthy; if they disagree, that is a result
#              to report, not a choice to make.
#
#   seed2      Re-draw the context sequences with a different random seed.
#   seed3      And again. Context is capped at two sequences per country per
#              year, so which two get drawn is arbitrary - the conclusion must
#              not depend on it.
#
# Each run writes to its own directory and NOTHING touches the main results.
# Afterwards the key question is asked of every run: do the 2025 Bangladeshi
# sequences still sit inside a clade of Bangladeshi sequences, and with what
# support? If the answer is the same three times, the result is robust.
#
# Runtime: roughly 10 minutes per run for the whole-genome tree. The E-gene tree
# is skipped here (SETS="wg") because the whole-genome tree is what carries the
# argument and the E-gene tree takes 20 minutes on its own; pass SETS="wg eg" to
# include it.
#
# Run from the project folder:  bash scripts/13_robustness.sh

set -uo pipefail

SETS="${SETS:-wg}"
ROOT="results/robustness"
mkdir -p "$ROOT"

# label | extra args to 04_select_sequences.py
RUNS=(
  "outgroup|--seed 1 --outgroup-lineages 3II_B --outgroup-max 4"
  "seed2|--seed 2"
  "seed3|--seed 3"
)

for entry in "${RUNS[@]}"; do
  label="${entry%%|*}"
  extra="${entry#*|}"
  d="$ROOT/$label"
  echo
  echo "############################################################"
  echo "#  $label"
  echo "#  04_select_sequences.py $extra"
  echo "############################################################"
  mkdir -p "$d"

  python scripts/04_select_sequences.py --outdir "$d/datasets" $extra \
    > "$d/select.log" 2>&1
  if [ $? -ne 0 ]; then
    echo "  selection FAILED - see $d/select.log"; tail -15 "$d/select.log"; continue
  fi
  grep -E "Output directory|Outgroup lineages|EXCLUDED|RESCUED|^================" \
       "$d/select.log" | head -8

  DATA="$d/datasets" OUT="$d/trees" SETS="$SETS" \
    bash scripts/06_build_trees.sh > "$d/trees.log" 2>&1
  if [ $? -ne 0 ]; then
    echo "  tree building FAILED - see $d/trees.log"; tail -15 "$d/trees.log"; continue
  fi
  grep -E "Best-fit model according|Log-likelihood of the tree|outgroup:|Branch separating" \
       "$d/trees.log" | head -6
done

echo
echo "############################################################"
echo "#  Comparison: does the conclusion hold across all runs?"
echo "############################################################"
python scripts/13b_compare_runs.py --root "$ROOT" --main-data results/datasets \
       --main-trees results/trees --sets $SETS
