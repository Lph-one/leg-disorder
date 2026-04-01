#!/usr/bin/env bash
set -euo pipefail

GROUP=""
PHENO_FILE=""
GENO_FILE=""
OUTPUT_DIR=""
HIBLUP="hiblup"
QCOVAR="2,3,4"
START_COL=5
END_COL=4495
PARALLEL_JOBS=10

while [[ $# -gt 0 ]]; do
  case "$1" in
    --group)
      GROUP="$2"
      shift 2
      ;;
    --pheno-file)
      PHENO_FILE="$2"
      shift 2
      ;;
    --geno-file)
      GENO_FILE="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --hiblup)
      HIBLUP="$2"
      shift 2
      ;;
    --qcovar)
      QCOVAR="$2"
      shift 2
      ;;
    --start-col)
      START_COL="$2"
      shift 2
      ;;
    --end-col)
      END_COL="$2"
      shift 2
      ;;
    --parallel-jobs)
      PARALLEL_JOBS="$2"
      shift 2
      ;;
    *)
      echo "[ERROR] Unknown argument: $1"
      exit 1
      ;;
  esac
done

if [[ -z "$GROUP" || -z "$PHENO_FILE" || -z "$GENO_FILE" || -z "$OUTPUT_DIR" ]]; then
  echo "Usage:"
  echo "  bash scripts/run_hiblup_by_group.sh --group CA --pheno-file data/pheno/CAHiblupPhe.txt --geno-file results/CA/CAID --output-dir results/CA/hiblup"
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

TASK_FILE="$OUTPUT_DIR/tasks.txt"
: > "$TASK_FILE"

for ((i=START_COL; i<=END_COL; i++)); do
  echo "$i meta$((i - START_COL + 1))" >> "$TASK_FILE"
done

process_trait() {
  col="$1"
  trait="$2"
  output_prefix="$OUTPUT_DIR/${trait}"
  log_file="$OUTPUT_DIR/${trait}.log"

  "$HIBLUP" \
    --single-trait \
    --pheno "$PHENO_FILE" \
    --pheno-pos "$col" \
    --qcovar "$QCOVAR" \
    --bfile "$GENO_FILE" \
    --add \
    --out "$output_prefix" >> "$log_file" 2>&1
}

export -f process_trait
export PHENO_FILE GENO_FILE OUTPUT_DIR HIBLUP QCOVAR

cat "$TASK_FILE" | parallel --colsep ' ' --jobs "$PARALLEL_JOBS" \
  "process_trait {1} {2}"

SUMMARY_FILE="$OUTPUT_DIR/summary_h2.txt"
echo -e "Trait\tCol\tVar_GA\tVar_GA_SE\tVar_e\tVar_e_SE\th2\th2_SE\th2_p_value" > "$SUMMARY_FILE"

for ((i=START_COL; i<=END_COL; i++)); do
  trait="meta$((i - START_COL + 1))"
  result_file="$OUTPUT_DIR/${trait}.vars"

  if [[ -f "$result_file" ]]; then
    awk -v trait="$trait" -v col="$i" '
    BEGIN {OFS="\t"}
    /^GA/ {var_ga=$2; var_ga_se=$3; h2=$4; h2_se=$5; h2_p=$6}
    /^e/  {var_e=$2; var_e_se=$3;
             print trait, col, var_ga, var_ga_se, var_e, var_e_se, h2, h2_se, h2_p}
    ' "$result_file" >> "$SUMMARY_FILE"
  fi
done

echo "[INFO] HIBLUP analysis for ${GROUP} finished. Summary saved to: ${SUMMARY_FILE}"