#  GWAS Loci Extraction and Colocalization

This repository contains a set of standalone scripts for performing colocalization (coloc) analyses between mQTL and molecular QTL (molQTL) summary statistics.


---

## Overview

The pipeline includes four main steps:

1. extracting significant variants from GWAS summary statistics
2. identifying independent lead variants using PLINK clumping
3. extracting regional loci around lead variants
4. performing colocalization analysis with eQTL or sQTL datasets


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
- coloc
- dplyr
- tools
- parallel

---

## Workflow
### 1. Extract significant variants

The script extract_significant_snps.py extracts significant variants from a directory of GWAS summary statistic files.

```python
python extract_significant_snps.py \
  --input-dir path/to/gwas_results \
  --output-dir path/to/significant_results \
  --pattern "*.assoc.txt" \
  --p-col p_wald \
  --threshold 5e-8 \
  --workers 16 \
  --log-file path/to/significant_extraction_log.tsv
```


### 2. Run PLINK clumping on all files in a directory
The script run_plink_clump.sh uses PLINK to clump on significant GWAS result files.

```bash
bash run_plink_clump.sh \
  --input-dir path/to/significant_results \
  --output-dir path/to/clump_results \
  --bfile path/to/plink_reference_prefix \
  --pattern "*.assoc.txt" \
  --plink plink \
  --clump-p1 1 \
  --clump-r2 0.1 \
  --clump-kb 100 \
  --chr-set 35
```


### 3. Extract loci regions around lead variants
The script extract_loci_regions.py extracts all variants located within a specified genomic window around lead variants from clumping results.

```python
python extract_loci_regions.py \
  --clump-dir path/to/clump_results \
  --assoc-dir path/to/gwas_results \
  --output-dir path/to/locus_regions \
  --window-bp 500000 \
  --chr-col chr \
  --pos-col ps \
  --workers 16
```


### 4. Run eQTL colocalization
The scripts run_coloc_eqtl.sh and R/rcoloc_eqtl.R perform colocalization analysis using the coloc R package.

```bash
bash run_coloc_eqtl.sh \
  --gwas-file path/to/locus_regions/trait1_region.txt \
  --qtl-dir path/to/eqtl_root \
  --out-dir path/to/eqtl_coloc_results \
  --r-script path/to/coloc_eqtl.R \
  --sample-size-file path/to/sample_size.txt \
  --num-cores 16 \
  --gwas-n SampleSize
```


### 5. Run sQTL colocalization
The scripts run_coloc_sqtl.sh and R/coloc_sqtl.R perform colocalization analysis using the coloc R package.

```bash
bash run_coloc_sqtl.sh \
  --gwas-file path/to/locus_regions/trait1_region.txt \
  --qtl-dir path/to/sqtl_root \
  --out-dir path/to/sqtl_coloc_results \
  --r-script path/to/coloc_sqtl.R \
  --sample-size-file path/to/sample_size.txt \
  --num-cores 16 \
  --gwas-n SampleSize
```

---






## License
For academic and research use only.