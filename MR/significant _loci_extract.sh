#!/bin/bash

ROOT=$(cd "$(dirname "$0")/.." && pwd)

IN_DIR="$ROOT/data/gwas_raw"
OUT_DIR="$ROOT/results/gwas_sig"

mkdir -p "$OUT_DIR"

SUMMARY="$OUT_DIR/summary.txt"
> "$SUMMARY"

echo "开始筛选 GWAS 文件 (P < 1e-8)..."

for file in "$IN_DIR"/GWAS_SMR_YZ126_*.fastGWA; do
    name=$(basename "$file")
    out="$OUT_DIR/${name%.fastGWA}_P1e-8.txt"

    header=$(head -n 1 "$file")

    awk 'NR>1 && $10 < 1e-8' "$file" > tmp

    { echo "$header"; cat tmp; } > "$out"
    rm tmp

    count=$(($(wc -l < "$out") - 1))
    echo -e "$name\t$count" >> "$SUMMARY"

    echo "$name → $count SNPs"
done

echo "完成"