#!/usr/bin/env bash
set -euo pipefail

bash scripts/run_hiblup_by_group_Hiblup.sh \
  --group CT \
  --pheno-file data/pheno/CTHiblupPhe.txt \
  --geno-file results/CT/CTID \
  --output-dir results/CT/hiblup \
  --hiblup hiblup \
  --qcovar 2,3,4 \
  --start-col 5 \
  --end-col 4495 \
  --parallel-jobs 49