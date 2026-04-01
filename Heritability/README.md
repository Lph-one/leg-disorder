# Genotype-Based PCA and Heritability Analysis (CA vs CT)

A reproducible pipeline for genotype-based PCA and large-scale heritability estimation using GCTA and HIBLUP.

This repository provides a reproducible workflow for genotype-based analysis, including group-wise data splitting, principal component analysis (PCA), and heritability estimation using both GCTA and HIBLUP. The pipeline is designed for comparative analysis between two groups (CA and CT) and supports large-scale multi-trait analysis in a parallelized manner.

---

## Overview

The pipeline consists of four main steps:

1. Split genotype data into CA and CT groups  
2. Perform PCA for each group  
3. Estimate heritability using GCTA  
4. Estimate heritability using HIBLUP  


---

## Scope of this repository
These scripts do not implement new statistical methods. Instead, they integrate existing tools into a consistent and reproducible workflow.

---

## Requirements

The following software is required:

- PLINK (v1.9 or later)
- GCTA (v1.93 or later)
- HIBLUP
- Python ≥ 3.7

Make sure all tools are available in your environment (`$PATH`).  

---

## Workflow
### 1. Split genotype data

```bash
bash split_genotype.sh
```

### 2. Run PCA

```bash
bash run_pca.sh
```

### 3. GCTA heritability analysis

```bash
bash run_CA_h2_GCTA.sh
bash run_CT_h2_GCTA.sh
```

### 4. HIBLUP heritability analysis

```bash
bash run_CA_h2_Hiblup.sh
bash run_CT_h2_Hiblup.sh
```


---

## Output
For each group (CA and CT), the pipeline generates:

- PCA results
- GRM files (GCTA)
- Trait-level analysis outputs
- Summary heritability files




---

## License
For academic and research use only.