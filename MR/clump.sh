#!/bin/bash

COHORT="CA"   # （CA / CT）

ROOT=$(cd "$(dirname "$0")/.." && pwd)

PLINK="plink"

BFILE="$ROOT/data/plink_ref/${COHORT}/${COHORT}ID"
GWAS_DIR="$ROOT/data/exposure_gwas/${COHORT}"
OUTDIR="$ROOT/results/${COHORT}_clump"

mkdir -p "$OUTDIR"

P1="1.78e-6"
R2="0.1"
KB="50"

echo "===== Clump: $COHORT ====="

for gwas in "$GWAS_DIR"/meta*.assoc.txt; do
    name=$(basename "$gwas" .assoc.txt)

    "$PLINK" \
      --bfile "$BFILE" \
      --allow-extra-chr \
      --clump "$gwas" \
      --clump-p1 "$P1" \
      --clump-r2 "$R2" \
      --clump-kb "$KB" \
      --clump-snp-field rs \
      --clump-field p_wald \
      --chr-set 35 \
      --out "$OUTDIR/${name}"

    echo "$name done"
done