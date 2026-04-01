#  Metabolite–Complex Trait Mendelian Randomization

This repository contains a set of standalone scripts for performing Mendelian randomization (MR) analyses between metabolite GWAS and complex trait GWAS summary statistics.


---

## Overview

The pipeline includes four main steps:

1. extract significant loci from the outcome GWAS (complex trait).
2. perform clumping on the exposure GWAS (metabolite) to obtain independent instrumental variables.
3. extract the clumped SNPs from the exposure GWAS
4. remove SNPs that overlap between the exposure instrumental variables and the outcome significant loci
5. conduct two‑sample Mendelian randomization analysis using `TwoSampleMR`


---

## Scope of this repository
These scripts do not implement new statistical methods. Instead, they integrate existing tools into a consistent and reproducible workflow.

---

## Requirements
Software
- Python ≥ 3.7
- R ≥ 4.0

Python packages
- pandas

R packages
- data.table
- TwoSampleMR
- dplyr

---

## Workflow
### 1. Extract significant loci from complex trait GWAS

The script significant_loci_extract.sh scans all complex trait GWAS files and keeps loci with:

P < 1e-8

```bash
bash scripts/significant_loci_extract.sh
```


### 2. Clump metabolite GWAS
The script clump.sh uses PLINK to perform LD clumping on metabolite GWAS summary statistics, producing approximately independent instrumental variables.

Key parameters:

- -clump-p1 = 1.78e-6
- clump-r2 = 0.1
- clump-kb = 50

```bash
bash scripts/clump.sh
```


### 3. Extract clumped SNPs from exposure GWAS
The script extract_snps.py reads PLINK .clumped files and extracts the corresponding SNPs from the original exposure GWAS summary statistics.

This produces one extracted exposure file per metabolite.

```python
python scripts/extract_snps.py
```


### 4. Remove overlapping SNPs
The script remove_overlap.py compares each extracted exposure file against each significant complex trait GWAS file and removes overlapping SNPs.

This generates exposure datasets suitable for downstream MR analyses.

```python
python scripts/remove_overlap.py
```


### 5. Run Mendelian randomization
The scripts run_mr.sh and R/run_mr.R perform two-sample MR analysis using the TwoSampleMR R package.

For each exposure–outcome pair, the code:

- harmonizes exposure and outcome alleles
- runs MR analysis
- calculates odds ratios
- computes heterogeneity statistics
- computes pleiotropy statistics
- estimates instrument strength using F statistics

```bash
bash scripts/run_mr.sh
```

---






## License
For academic and research use only.