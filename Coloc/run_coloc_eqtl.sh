#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage:
  bash run_coloc_eqtl.sh \\
    --gwas-file FILE \\
    --qtl-dir QTL_DIR \\
    --out-dir OUT_DIR \\
    --r-script coloc_eqtl.R \\
    --sample-size-file SAMPLE_SIZE_TXT \\
    [--docker-image myimage:tag] \\
    [--num-cores 8] \\
    [--gwas-n 100]

Description:
  Run coloc analysis between one GWAS regional file and eQTL datasets.
EOF
}

GWAS_FILE=""
QTL_DIR=""
OUT_DIR=""
R_SCRIPT=""
SAMPLE_SIZE_FILE=""
DOCKER_IMAGE="chgyi/binbash-r-base:4.2.2"
NUM_CORES="8"
GWAS_N="92"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --gwas-file) GWAS_FILE="$2"; shift 2 ;;
        --qtl-dir) QTL_DIR="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --r-script) R_SCRIPT="$2"; shift 2 ;;
        --sample-size-file) SAMPLE_SIZE_FILE="$2"; shift 2 ;;
        --docker-image) DOCKER_IMAGE="$2"; shift 2 ;;
        --num-cores) NUM_CORES="$2"; shift 2 ;;
        --gwas-n) GWAS_N="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1"; usage; exit 1 ;;
    esac
done

[[ -z "$GWAS_FILE" || -z "$QTL_DIR" || -z "$OUT_DIR" || -z "$R_SCRIPT" || -z "$SAMPLE_SIZE_FILE" ]] && {
    usage
    exit 1
}

mkdir -p "$OUT_DIR"

docker run --rm \
    -v /:/host \
    -w /host \
    "$DOCKER_IMAGE" \
    Rscript "/host/${R_SCRIPT#/}" \
        "/host/${GWAS_FILE#/}" \
        "/host/${QTL_DIR#/}" \
        "/host/${OUT_DIR#/}" \
        "/host/${SAMPLE_SIZE_FILE#/}" \
        "$NUM_CORES" \
        "$GWAS_N"