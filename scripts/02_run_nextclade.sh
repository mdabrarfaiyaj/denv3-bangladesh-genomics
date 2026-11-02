#!/usr/bin/env bash
#
# 02_run_nextclade.sh
# -------------------
# Assign DENV-3 lineages and run sequence QC with Nextclade.
#
# Input : data/raw/denv3_all/ncbi_dataset/data/genomic.fna   (all 9,815 sequences)
# Output: results/nextclade/nextclade.tsv                    lineage + QC per sequence
#         results/nextclade/nextclade.aligned.fasta          all sequences on reference coordinates
#         results/nextclade/dataset/                         the reference dataset, saved with the project
#
# Run from the project folder:  bash scripts/02_run_nextclade.sh

set -euo pipefail

# The dataset version is pinned on purpose. Lineage names change as the dataset is
# updated, so an unpinned run would give different answers in six months.
DATASET="community/v-gen-lab/dengue/denv3"
TAG="2026-10-06--17-37-18Z"

IN="data/raw/denv3_all/ncbi_dataset/data/genomic.fna"
OUTDIR="results/nextclade"

[ -f "$IN" ] || { echo "Input not found: $IN (download the data first, see README)"; exit 1; }
mkdir -p "$OUTDIR"

echo "== Nextclade version =="
nextclade --version

echo
echo "== Fetching dataset $DATASET @ $TAG =="
nextclade dataset get \
  --name "$DATASET" \
  --tag "$TAG" \
  --output-dir "$OUTDIR/dataset"

echo
echo "== Running Nextclade on all sequences =="
nextclade run \
  --input-dataset "$OUTDIR/dataset" \
  --output-tsv "$OUTDIR/nextclade.tsv" \
  --output-fasta "$OUTDIR/nextclade.aligned.fasta" \
  --jobs 4 \
  "$IN"

echo
echo "== Done =="
wc -l < "$OUTDIR/nextclade.tsv" | xargs echo "Rows in nextclade.tsv (incl. header):"
echo "Now run: python scripts/02b_summarise_nextclade.py"
