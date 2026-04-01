#!/usr/bin/env bash
set -euo pipefail

python scripts/run_h2_by_group_GCTA.py \
  --group CA \
  --gcta gcta64 \
  --bfile-prefix results/CA/CAID \
  --phe-file data/pheno/CAID_phe.txt \
  --qcovar-file data/covar/CAID_Cov.txt \
  --output-dir results/CA/heritability