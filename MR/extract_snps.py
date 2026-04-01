#!/usr/bin/env python3

import os
import pandas as pd

COHORT = "CA"   # (CA/CT)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

clump_dir = os.path.join(ROOT, "results", f"{COHORT}_clump")
assoc_dir = os.path.join(ROOT, "data", "exposure_gwas", COHORT)
out_dir   = os.path.join(ROOT, "results", f"{COHORT}_extract")

os.makedirs(out_dir, exist_ok=True)

for fname in os.listdir(clump_dir):
    if not fname.endswith(".clumped"):
        continue

    prefix = fname.replace(".clumped", "")

    clump_path = os.path.join(clump_dir, fname)
    assoc_path = os.path.join(assoc_dir, f"{prefix}.assoc.txt")

    if not os.path.exists(assoc_path):
        continue

    print(f"处理 {prefix}")

    clump = pd.read_csv(clump_path, delim_whitespace=True)

    if "SNP" not in clump.columns:
        continue

    snps = set(clump["SNP"].dropna())

    assoc = pd.read_csv(assoc_path, delim_whitespace=True)

    res = assoc[assoc["rs"].isin(snps)]

    res.to_csv(os.path.join(out_dir, f"{prefix}_extract.txt"), sep="\t", index=False)