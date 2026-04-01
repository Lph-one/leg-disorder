#!/usr/bin/env bash
set -euo pipefail

PLINK=plink

${PLINK} \
  --bfile results/CA/CAID \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 40 \
  --pca 3 \
  --out results/CA/snp.pca3

awk '{print 1,$3,$4,$5}' results/CA/snp.pca3.eigenvec > results/CA/c.txt

${PLINK} \
  --bfile results/CT/CTID \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 40 \
  --pca 3 \
  --out results/CT/snp.pca3

awk '{print 1,$3,$4,$5}' results/CT/snp.pca3.eigenvec > results/CT/c.txt