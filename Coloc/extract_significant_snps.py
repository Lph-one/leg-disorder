#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import argparse
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed


def process_file(file_path, output_dir, p_col, threshold, sep):
    try:
        df = pd.read_csv(file_path, sep=sep)
        if p_col not in df.columns:
            return {
                "file": os.path.basename(file_path),
                "success": False,
                "n_sig": 0,
                "message": f"Missing column: {p_col}"
            }

        sig_df = df[df[p_col] < threshold].copy()

        if sig_df.empty:
            return {
                "file": os.path.basename(file_path),
                "success": True,
                "n_sig": 0,
                "message": "No significant variants"
            }

        out_path = os.path.join(output_dir, os.path.basename(file_path))
        sig_df.to_csv(out_path, sep="\t", index=False)

        return {
            "file": os.path.basename(file_path),
            "success": True,
            "n_sig": len(sig_df),
            "message": "Saved"
        }

    except Exception as e:
        return {
            "file": os.path.basename(file_path),
            "success": False,
            "n_sig": 0,
            "message": f"Error: {e}"
        }


def main():
    parser = argparse.ArgumentParser(
        description="Extract significant variants from GWAS summary statistic files."
    )
    parser.add_argument("--input-dir", required=True, help="Directory containing GWAS result files")
    parser.add_argument("--output-dir", required=True, help="Directory to save significant result files")
    parser.add_argument("--pattern", default="*.assoc.txt", help="File pattern, e.g. 'meta*.assoc.txt'")
    parser.add_argument("--p-col", default="p_wald", help="P-value column name")
    parser.add_argument("--threshold", type=float, default=1.78e-6, help="Significance threshold")
    parser.add_argument("--sep", default="\t", help="Input file separator")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel workers")
    parser.add_argument("--log-file", default=None, help="Optional log file path (tsv)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(args.input_dir, args.pattern)))

    print(f"Found {len(files)} files in: {args.input_dir}")

    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                process_file,
                file_path=f,
                output_dir=args.output_dir,
                p_col=args.p_col,
                threshold=args.threshold,
                sep=args.sep,
            ): f for f in files
        }

        for future in as_completed(futures):
            res = future.result()
            results.append(res)
            if res["success"] and res["message"] == "Saved":
                print(f"[OK] {res['file']} -> {res['n_sig']} significant variants")
            elif res["success"]:
                print(f"[INFO] {res['file']} -> {res['message']}")
            else:
                print(f"[ERR] {res['file']} -> {res['message']}")

    if args.log_file:
        log_df = pd.DataFrame(results)
        log_df.to_csv(args.log_file, sep="\t", index=False)
        print(f"Log saved to: {args.log_file}")

    print("Done.")


if __name__ == "__main__":
    main()