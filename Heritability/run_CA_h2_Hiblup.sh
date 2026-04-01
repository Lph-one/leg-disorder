#!/usr/bin/env bash
set -euo pipefail

bash scripts/run_hiblup_by_group_Hiblup.sh \
  --group CA \
  --pheno-file data/pheno/CAHiblupPhe.txt \
  --geno-file results/CA/CAID \
  --output-dir results/CA/hiblup \
  --hiblup hiblup \
  --qcovar 2,3,4 \
  --start-col 5 \
  --end-col 4495 \
  --parallel-jobs 49