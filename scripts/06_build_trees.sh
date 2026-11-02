#!/usr/bin/env bash
#
# 06_build_trees.sh
# -----------------
# Maximum-likelihood trees for both datasets with IQ-TREE 2.
#
# Input : results/datasets/{wg,eg}.fasta      (+ *_traits.tsv for the outgroup)
# Output: results/trees/{wg,eg}.treefile      the ML tree
#         results/trees/{wg,eg}.iqtree        the full report: model chosen, likelihoods
#         results/trees/{wg,eg}.log
#
# Run from the project folder:  bash scripts/06_build_trees.sh
#
# Choices made here, and why:
#
#  -m MFP      ModelFinder Plus picks the substitution model by BIC rather than
#              assuming one. For a coding virus genome the winner is usually a
#              GTR variant with rate heterogeneity. The chosen model goes in the
#              methods section - read it out of the .iqtree file, do not guess.
#
#  -B 1000     Ultrafast bootstrap, 1000 replicates. Interpretation differs from
#              the classical bootstrap: UFBoot >= 95 is roughly the threshold for
#              a well-supported branch, not 70.
#
#  --alrt 1000 SH-aLRT, a second support measure from a different principle.
#              A branch supported by both (UFBoot >= 95 and SH-aLRT >= 80) is
#              solid; one supported by only one of them is not. Reporting two
#              independent measures is harder to argue with than reporting one.
#
#  -o          Root on the ancestral 3I sequences. Outgroup rooting is one of two
#              independent rootings used in this project; TreeTime's clock-based
#              rooting is the other, in the next step. If they agree, the root is
#              trustworthy. If they disagree, that disagreement is a result and
#              must be reported.
#
#  --seqtype DNA   stated explicitly so a gap-heavy alignment is never guessed
#                  as protein.
#
# IQ-TREE collapses identical sequences before inference and prints which ones.
# That is expected with closely related outbreak sequences; the dropped names are
# in the .log and must be mentioned when reporting the number of taxa.

set -euo pipefail

# Parameterised so each robustness run reads and writes its own directory, and
# so a workflow manager (Nextflow) can drive it. Defaults reproduce the main run.
DATA="${DATA:-results/datasets}"
OUT="${OUT:-results/trees}"
SETS="${SETS:-wg eg}"
# THREADS: AUTO lets IQ-TREE choose; a workflow manager should pass its own
# allocation instead, or IQ-TREE will grab every core on a shared runner.
THREADS="${THREADS:-AUTO}"
# MODEL: MFP runs ModelFinder, which is right for the real analysis but wasteful
# on a twelve-sequence test dataset.
MODEL="${MODEL:-MFP}"
BOOT="${BOOT:-1000}"
mkdir -p "$OUT"
echo "input datasets : $DATA"
echo "output trees   : $OUT"
echo "threads/model  : $THREADS / $MODEL, $BOOT bootstrap replicates"

echo "== IQ-TREE version =="
iqtree2 --version 2>/dev/null || iqtree --version
IQ=$(command -v iqtree2 || command -v iqtree)

for SET in $SETS; do
  FASTA="$DATA/$SET.fasta"
  TRAITS="$DATA/${SET}_traits.tsv"
  [ -f "$FASTA" ] || { echo "Missing $FASTA - run scripts/04_select_sequences.py first"; exit 1; }

  # Outgroup: whichever sequences script 04 marked, read from the named column
  # rather than hardcoding a lineage. Falls back to the old behaviour if the
  # column is absent, so older datasets still work.
  OUTGROUP=$(awk -F'\t' '
    NR==1 { for (i=1;i<=NF;i++) if ($i=="is_outgroup") c=i; next }
    c && $c=="yes" { printf "%s%s", sep, $1; sep="," }
  ' "$TRAITS")
  if [ -z "$OUTGROUP" ]; then
    OUTGROUP=$(awk -F'\t' 'NR>1 && $4=="3I" {printf "%s%s", sep, $1; sep=","}' "$TRAITS")
  fi

  echo
  echo "=============================================================="
  echo "  $SET : $(grep -c '^>' "$FASTA") sequences"
  echo "  alignment length: $(awk 'NR==2{print length($0)}' "$FASTA") bases"
  if [ -n "$OUTGROUP" ]; then
    echo "  outgroup: $OUTGROUP"
  else
    echo "  no 3I sequences found - tree will be left unrooted for TreeTime to root"
  fi
  echo "=============================================================="

  ARGS=(-s "$FASTA" --prefix "$OUT/$SET" --seqtype DNA
        -m "$MODEL" -B "$BOOT" --alrt "$BOOT" -T "$THREADS" --redo)
  [ -n "$OUTGROUP" ] && ARGS+=(-o "$OUTGROUP")

  "$IQ" "${ARGS[@]}"

  echo
  echo "---- $SET summary ----"
  grep -E "Best-fit model|Model of rate heterogeneity" "$OUT/$SET.iqtree" || true
  grep -E "^Log-likelihood of the tree" "$OUT/$SET.iqtree" || true
  awk '/identical/ {print "  " $0}' "$OUT/$SET.log" | head -20 || true
done

echo
echo "== Done =="
echo "Trees: $OUT/wg.treefile and $OUT/eg.treefile"
echo
echo "Record in the methods, read from the .iqtree files (do not guess):"
echo "  - the best-fit model for each dataset"
echo "  - the number of taxa actually analysed after identical sequences were merged"
echo
echo "Next: python scripts/07_clock_check.py"
