#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed


def infer_traits(clump_dir, assoc_dir, trait_list_file=None, assoc_suffix=".assoc.txt", clump_suffix=".clumped"):
    if trait_list_file:
        with open(trait_list_file) as f:
            traits = [x.strip() for x in f if x.strip()]
        return sorted(set(traits))

    clump_traits = {
        f[:-len(clump_suffix)]
        for f in os.listdir(clump_dir)
        if f.endswith(clump_suffix)
    }
    assoc_traits = {
        f[:-len(assoc_suffix)]
        for f in os.listdir(assoc_dir)
        if f.endswith(assoc_suffix)
    }
    return sorted(clump_traits & assoc_traits)


def process_single_trait(
    trait,
    clump_dir,
    assoc_dir,
    output_dir,
    window_bp,
    assoc_suffix,
    clump_suffix,
    chr_col,
    pos_col,
):
    clump_file = os.path.join(clump_dir, f"{trait}{clump_suffix}")
    assoc_file = os.path.join(assoc_dir, f"{trait}{assoc_suffix}")
    out_file = os.path.join(output_dir, f"{trait}_region.txt")

    if not os.path.exists(clump_file):
        return trait, False, f"Missing clump file: {clump_file}"
    if not os.path.exists(assoc_file):
        return trait, False, f"Missing assoc file: {assoc_file}"

    try:
        clump_df = pd.read_csv(clump_file, sep=r"\s+")
        if clump_df.empty:
            return trait, False, "Empty clump file"

        required_clump_cols = {"CHR", "BP", "SNP"}
        if not required_clump_cols.issubset(clump_df.columns):
            return trait, False, f"Clump file missing columns: {required_clump_cols - set(clump_df.columns)}"

        assoc_df = pd.read_csv(assoc_file, sep="\t")
        if chr_col not in assoc_df.columns or pos_col not in assoc_df.columns:
            return trait, False, f"Assoc file missing chr/pos columns: {chr_col}, {pos_col}"

        windows = []
        for _, row in clump_df.iterrows():
            bp = int(row["BP"])
            windows.append({
                "chr": str(row["CHR"]),
                "start": max(1, bp - window_bp),
                "end": bp + window_bp,
            })

        selected_parts = []
        for chr_value in sorted(set(w["chr"] for w in windows)):
            chr_windows = [w for w in windows if w["chr"] == chr_value]
            chr_df = assoc_df[assoc_df[chr_col].astype(str) == chr_value].copy()

            if chr_df.empty:
                continue

            mask = pd.Series(False, index=chr_df.index)
            for w in chr_windows:
                mask = mask | ((chr_df[pos_col] >= w["start"]) & (chr_df[pos_col] <= w["end"]))

            selected = chr_df[mask]
            if not selected.empty:
                selected_parts.append(selected)

        if not selected_parts:
            return trait, False, "No variants found in lead-SNP windows"

        result_df = pd.concat(selected_parts, ignore_index=True).drop_duplicates()
        result_df.to_csv(out_file, sep="\t", index=False)

        return trait, True, f"Saved {len(result_df)} variants from {len(clump_df)} lead SNP windows"

    except Exception as e:
        return trait, False, f"Error: {e}"


def main():
    parser = argparse.ArgumentParser(description="Extract loci regions around lead SNPs from clump results.")
    parser.add_argument("--clump-dir", required=True, help="Directory containing .clumped files")
    parser.add_argument("--assoc-dir", required=True, help="Directory containing GWAS assoc files")
    parser.add_argument("--output-dir", required=True, help="Directory for extracted region files")
    parser.add_argument("--trait-list", default=None, help="Optional file listing trait names")
    parser.add_argument("--window-bp", type=int, default=500000, help="Window size on each side of lead SNP")
    parser.add_argument("--assoc-suffix", default=".assoc.txt", help="Suffix of GWAS assoc files")
    parser.add_argument("--clump-suffix", default=".clumped", help="Suffix of clump files")
    parser.add_argument("--chr-col", default="chr", help="Chromosome column in assoc files")
    parser.add_argument("--pos-col", default="ps", help="Position column in assoc files")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel workers")
    parser.add_argument("--log-file", default=None, help="Optional output log file")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    traits = infer_traits(
        clump_dir=args.clump_dir,
        assoc_dir=args.assoc_dir,
        trait_list_file=args.trait_list,
        assoc_suffix=args.assoc_suffix,
        clump_suffix=args.clump_suffix,
    )

    print(f"Traits to process: {len(traits)}")
    print(f"Window size: ±{args.window_bp} bp")

    results = []
    success_count = 0
    fail_count = 0

    with ProcessPoolExecutor(max_workers=min(args.workers, max(1, len(traits)))) as executor:
        futures = {
            executor.submit(
                process_single_trait,
                trait=t,
                clump_dir=args.clump_dir,
                assoc_dir=args.assoc_dir,
                output_dir=args.output_dir,
                window_bp=args.window_bp,
                assoc_suffix=args.assoc_suffix,
                clump_suffix=args.clump_suffix,
                chr_col=args.chr_col,
                pos_col=args.pos_col,
            ): t
            for t in traits
        }

        for i, future in enumerate(as_completed(futures), 1):
            trait, success, msg = future.result()
            results.append((trait, success, msg))
            if success:
                success_count += 1
                print(f"{i:4d}/{len(traits)} [OK]   {trait}: {msg}")
            else:
                fail_count += 1
                print(f"{i:4d}/{len(traits)} [FAIL] {trait}: {msg}")

    print("=" * 60)
    print(f"Success: {success_count}")
    print(f"Fail   : {fail_count}")
    print(f"Output : {args.output_dir}")

    if args.log_file:
        pd.DataFrame(results, columns=["trait", "success", "message"]).to_csv(
            args.log_file, sep="\t", index=False
        )
        print(f"Log saved to: {args.log_file}")


if __name__ == "__main__":
    main()