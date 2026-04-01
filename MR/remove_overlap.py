#!/usr/bin/env python3

import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

COHORT = "CA"   # (CA/CT)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

gwas_dir = os.path.join(ROOT, "results", "gwas_sig")
meta_dir = os.path.join(ROOT, "results", f"{COHORT}_extract")
out_dir  = os.path.join(ROOT, "results", f"{COHORT}_exp")

os.makedirs(out_dir, exist_ok=True)

gwas_files = [f for f in os.listdir(gwas_dir)]
meta_files = [f for f in os.listdir(meta_dir)]

def run(meta_file, gwas_file):
    meta = pd.read_csv(os.path.join(meta_dir, meta_file), sep="\t")
    gwas = pd.read_csv(os.path.join(gwas_dir, gwas_file), sep="\t")

    if "rs" not in meta.columns or "SNP" not in gwas.columns:
        return

    res = meta[~meta["rs"].isin(gwas["SNP"])]

    meta_id = meta_file.split("_")[0]
    pheno = gwas_file.split("_")[3]

    out = f"{meta_id}_P{pheno}exp.txt"
    res.to_csv(os.path.join(out_dir, out), sep="\t", index=False)

tasks = [(m, g) for m in meta_files for g in gwas_files]

with ProcessPoolExecutor(16) as pool:
    pool.map(lambda x: run(*x), tasks)