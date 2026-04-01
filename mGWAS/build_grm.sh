#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage:
  bash scripts/build_grm.sh --group CA --bfile data/genotype/CA/CAID --gemma /path/to/gemma [--outdir results/CA/grm]

Required arguments:
  --group     Group name, e.g. CA or CT
  --bfile     PLINK file prefix (without .bed/.bim/.fam)
  --gemma     Path to GEMMA executable

Optional arguments:
  --outdir    Output directory (default: results/<group>/grm)
EOF
}

GROUP=""
BFILE=""
GEMMA=""
OUTDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --group)
      GROUP="$2"
      shift 2
      ;;
    --bfile)
      BFILE="$2"
      shift 2
      ;;
    --gemma)
      GEMMA="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown argument: $1"
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$GROUP" || -z "$BFILE" || -z "$GEMMA" ]]; then
  echo "[ERROR] Missing required arguments."
  usage
  exit 1
fi

if [[ -z "$OUTDIR" ]]; then
  OUTDIR="results/${GROUP}/grm"
fi

mkdir -p "$OUTDIR"

echo "[INFO] Group      : $GROUP"
echo "[INFO] BFILE      : $BFILE"
echo "[INFO] GEMMA      : $GEMMA"
echo "[INFO] Output dir : $OUTDIR"

"$GEMMA" -bfile "$BFILE" -gk 1 -o "${GROUP}.snp.all.matrix"

if [[ -d "output" ]]; then
  mv output/* "$OUTDIR"/
  rmdir output 2>/dev/null || true
fi

echo "[INFO] GRM construction finished."
echo "[INFO] Kinship matrix should be in: $OUTDIR"