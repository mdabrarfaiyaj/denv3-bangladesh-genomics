#!/usr/bin/env bash
#
# 00_download.sh
# --------------
# Fetch the DENV-3 sequences from NCBI and record exactly what was fetched.
#
# This is deliberately NOT a step in the Nextflow pipeline. Three reasons:
#
#   1. GenBank changes daily, so a download step makes the pipeline
#      irreproducible by construction - the same command gives different data
#      next week. Pinning the data means recording a checksum, not re-fetching.
#   2. Re-downloading 9,815 sequences on every pipeline run is slow and impolite
#      to a public service.
#   3. Continuous integration must not depend on an external database being up.
#
# So the download happens once, here, and writes a manifest. The pipeline takes
# the downloaded files as input. This is how nf-core pipelines handle reference
# data, for the same reasons.
#
# Run from the project folder:  bash scripts/00_download.sh

set -euo pipefail

TAXID="${TAXID:-11069}"            # dengue virus type 3
RAW="data/raw"
PKG="$RAW/denv3_all.zip"
DIR="$RAW/denv3_all"

command -v datasets >/dev/null 2>&1 || {
  echo "ncbi-datasets-cli not found. conda activate denv3"; exit 1; }

mkdir -p "$RAW"

echo "== Downloading DENV-3 (taxid $TAXID) from NCBI Virus =="
datasets download virus genome taxon "$TAXID" --filename "$PKG"

echo
echo "== Unpacking =="
unzip -o -q "$PKG" -d "$DIR"

DATA="$DIR/ncbi_dataset/data"
FNA="$DATA/genomic.fna"
JSONL="$DATA/data_report.jsonl"
for f in "$FNA" "$JSONL"; do
  [ -s "$f" ] || { echo "Expected file missing or empty: $f"; exit 1; }
done

N=$(grep -c '^>' "$FNA")
RECORDS=$(wc -l < "$JSONL")

echo
echo "== Recording what was downloaded =="
MANIFEST="$RAW/download_manifest.txt"
{
  echo "# DENV-3 download manifest"
  echo "# This file pins the data the analysis was run on. If you re-download"
  echo "# later you will get MORE sequences, and the numbers in the paper will"
  echo "# not match. That is expected: cite this manifest, not a date alone."
  echo
  echo "downloaded_utc:    $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "taxid:             $TAXID"
  echo "tool:              $(datasets --version 2>&1 | head -1)"
  echo "sequences_fasta:   $N"
  echo "records_jsonl:     $RECORDS"
  echo
  echo "sha256:"
  sha256sum "$PKG" "$FNA" "$JSONL"
} > "$MANIFEST"

cat "$MANIFEST"

echo
if [ "$N" -ne "$RECORDS" ]; then
  echo "WARNING: $N sequences in the FASTA but $RECORDS records in the JSONL."
  echo "They should match. Check before going further."
else
  echo "FASTA and JSONL agree: $N sequences."
fi

echo
echo "Next:  nextflow run . -profile conda"
echo "       (or step by step, starting with scripts/01_parse_metadata.py)"
echo
echo "Commit $MANIFEST. Do NOT commit $RAW itself - it is in .gitignore,"
echo "because a 14 MB archive that changes weekly does not belong in git."
