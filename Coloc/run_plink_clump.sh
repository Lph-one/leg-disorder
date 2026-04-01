#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage:
  bash run_plink_clump.sh \\
    --input-dir SIG_DIR \\
    --output-dir CLUMP_DIR \\
    --bfile PLINK_BFILE_PREFIX \\
    [--file-list FILE_LIST] \\
    [--pattern "*.assoc.txt"] \\
    [--plink plink] \\
    [--clump-p1 1] \\
    [--clump-r2 0.1] \\
    [--clump-kb 100] \\
    [--chr-set 35]

Description:
  Run PLINK clump for significant GWAS result files.
EOF
}

INPUT_DIR=""
OUTPUT_DIR=""
BFILE=""
FILE_LIST=""
PATTERN="*.assoc.txt"
PLINK_BIN="plink"
CLUMP_P1="1"
CLUMP_R2="0.1"
CLUMP_KB="100"
CHR_SET="35"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir) INPUT_DIR="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --bfile) BFILE="$2"; shift 2 ;;
        --file-list) FILE_LIST="$2"; shift 2 ;;
        --pattern) PATTERN="$2"; shift 2 ;;
        --plink) PLINK_BIN="$2"; shift 2 ;;
        --clump-p1) CLUMP_P1="$2"; shift 2 ;;
        --clump-r2) CLUMP_R2="$2"; shift 2 ;;
        --clump-kb) CLUMP_KB="$2"; shift 2 ;;
        --chr-set) CHR_SET="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1"; usage; exit 1 ;;
    esac
done

[[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$BFILE" ]] && { usage; exit 1; }

mkdir -p "$OUTPUT_DIR"

detect_fields() {
    local file="$1"
    local header
    header=$(head -n 1 "$file")

    local snp_field="SNP"
    local p_field="p"

    if [[ "$header" == *$'\trs\t'* || "$header" == rs$'\t'* || "$header" == *$'\trs' ]]; then
        snp_field="rs"
    elif [[ "$header" == *"SNP"* ]]; then
        snp_field="SNP"
    fi

    if [[ "$header" == *"p_wald"* ]]; then
        p_field="p_wald"
    elif [[ "$header" == *$'\tp\t'* || "$header" == p$'\t'* || "$header" == *$'\tp' ]]; then
        p_field="p"
    fi

    echo "${snp_field} ${p_field}"
}

run_one() {
    local file="$1"
    local in_file="${INPUT_DIR}/${file}"

    if [[ ! -f "$in_file" ]]; then
        echo "[WARN] Missing file: $in_file"
        return
    fi

    local prefix="$file"
    prefix="${prefix%.assoc.txt}"

    read -r snp_field p_field < <(detect_fields "$in_file")

    echo "[INFO] Processing: $file"
    echo "       SNP field: $snp_field"
    echo "       P field  : $p_field"

    "$PLINK_BIN" \
        --bfile "$BFILE" \
        --allow-extra-chr \
        --chr-set "$CHR_SET" \
        --clump-p1 "$CLUMP_P1" \
        --clump-r2 "$CLUMP_R2" \
        --clump-kb "$CLUMP_KB" \
        --clump "$in_file" \
        --clump-snp-field "$snp_field" \
        --clump-field "$p_field" \
        --out "${OUTPUT_DIR}/${prefix}"

    if [[ -f "${OUTPUT_DIR}/${prefix}.clumped" ]]; then
        local count
        count=$(tail -n +2 "${OUTPUT_DIR}/${prefix}.clumped" | wc -l | awk '{print $1}')
        echo "[OK] ${prefix}.clumped generated (${count} lead variants)"
    else
        echo "[WARN] No clumped result for: $file"
        [[ -f "${OUTPUT_DIR}/${prefix}.log" ]] && tail -20 "${OUTPUT_DIR}/${prefix}.log"
    fi

    echo "----------------------------------------"
}

if [[ -n "$FILE_LIST" ]]; then
    while read -r file; do
        file=$(echo "$file" | tr -d '\r\n ')
        [[ -z "$file" ]] && continue
        run_one "$file"
    done < "$FILE_LIST"
else
    shopt -s nullglob
    for path in "${INPUT_DIR}"/$PATTERN; do
        run_one "$(basename "$path")"
    done
fi

echo "All clump analyses completed."