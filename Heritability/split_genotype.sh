#!/usr/bin/env bash
set -euo pipefail

PLINK=plink

BFILE_PREFIX="data/raw/txx"
CA_KEEP="data/sample_lists/CAID.txt"
CT_KEEP="data/sample_lists/CTID.txt"

mkdir -p results/CA results/CT

${PLINK} \
  --bfile "${BFILE_PREFIX}" \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 40 \
  --keep "${CA_KEEP}" \
  --make-bed \
  --out results/CA/CAID

${PLINK} \
  --bfile "${BFILE_PREFIX}" \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 40 \
  --keep "${CT_KEEP}" \
  --make-bed \
  --out results/CT/CTID