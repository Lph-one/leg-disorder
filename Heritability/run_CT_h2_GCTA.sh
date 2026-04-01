#!/usr/bin/env bash
set -euo pipefail

python scripts/run_h2_by_group_GCTA.py \
  --group CT \
  --gcta gcta64 \
  --bfile-prefix results/CT/CTID \
  --phe-file data/pheno/CTID_phe.txt \
  --qcovar-file data/covar/CTID_Cov.txt \
  --output-dir results/CT/heritability