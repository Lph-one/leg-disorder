#!/usr/bin/env python3
import os
import math
import argparse
import numpy as np
import pandas as pd


def calculate_lambda(p_values):
    p_values = np.array(p_values, dtype=float)
    p_values = p_values[~np.isnan(p_values)]

    if len(p_values) == 0:
        return None

    chisq_values = []
    for p in p_values:
        if p <= 0 or p >= 1:
            continue
        z = -4.91 * (math.pow(p, 0.14) - math.pow(1 - p, 0.14))
        chisq = z * z
        chisq_values.append(chisq)

    if len(chisq_values) == 0:
        return None

    median_chisq = np.median(chisq_values)
    lambda_gc = median_chisq / 0.4549364
    return lambda_gc


def detect_p_column(df):
    candidate_cols = ["p_wald", "p_lrt", "p_score", "p", "P", "P-value", "P_value"]
    for col in candidate_cols:
        if col in df.columns:
            return col

    for col in df.columns:
        if "p" in col.lower() and col.lower() != "pos":
            return col

    return None


def main():
    parser = argparse.ArgumentParser(description="Calculate genomic inflation factor (lambda) from GEMMA association files.")
    parser.add_argument("--input-dir", required=True, help="Directory containing *.assoc.txt files")
    parser.add_argument("--output-file", required=True, help="Output summary file")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_file = args.output_file

    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    files = [
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.endswith(".assoc.txt")
    ]

    if not files:
        raise FileNotFoundError(f"No .assoc.txt files found in: {input_dir}")

    files.sort()

    print(f"[INFO] Found {len(files)} association files")

    with open(output_file, "w") as fout:
        fout.write("Trait\tLambda\tEffectiveSNP\n")

        for input_file in files:
            trait_name = os.path.basename(input_file).replace(".assoc.txt", "")

            try:
                df = pd.read_csv(input_file, sep="\t")
                p_col = detect_p_column(df)

                if p_col is None:
                    fout.write(f"{trait_name}\tNA\t0\n")
                    print(f"[WARN] {trait_name}: no p-value column found")
                    continue

                p_values = pd.to_numeric(df[p_col], errors="coerce")
                valid_count = int(np.sum(~np.isnan(p_values)))

                lambda_gc = calculate_lambda(p_values)

                if lambda_gc is None:
                    fout.write(f"{trait_name}\tNA\t{valid_count}\n")
                    print(f"[WARN] {trait_name}: lambda could not be calculated")
                else:
                    fout.write(f"{trait_name}\t{lambda_gc:.6f}\t{valid_count}\n")
                    print(f"[INFO] {trait_name}: lambda={lambda_gc:.6f}, SNPs={valid_count}")

            except Exception as e:
                fout.write(f"{trait_name}\tERROR\t0\n")
                print(f"[ERROR] {trait_name}: {e}")

    print(f"[INFO] Results saved to: {output_file}")


if __name__ == "__main__":
    main()