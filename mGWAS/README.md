#  Metabolomics GWAS

A reproducible pipeline for metabolite genome-wide association analysis (mGWAS), including genomic relationship matrix construction, GWAS, genomic inflation factor calculation, and locus annotation.


---

## Overview

This repository provides a modular workflow for metabolite genome-wide association analysis (mGWAS). The pipeline includes four main steps:

1. Build the genomic relationship matrix (GRM) using GEMMA  
2. Perform GWAS for multiple metabolite traits  
3. Calculate the genomic inflation factor (\u03bb) for quality assessment  
4. Annotate significant or selected loci based on genome annotation  


---

## Scope of this repository
These scripts do not implement new statistical methods. Instead, they integrate existing tools into a consistent and reproducible workflow.

---
## Requirements
Software
- GEMMA (v0.98.4 or later)
- Python ≥ 3.7
- R ≥ 4.0

Python packages
- numpy
- pandas

R packages
- ChIPpeakAnno
- GenomicRanges

---

## Workflow
### 1. Build genomic relationship matrix (GRM)

Use build_grm.sh to generate the relatedness matrix with GEMMA.

```bash
bash scripts/build_grm.sh \
  --group sample \
  --bfile data/genotype/sample_prefix \
  --gemma /path/to/gemma \
  --outdir results/grm
```

This step generates GEMMA kinship matrix files, including:

- *.cXX.txt
- log files

### 2. Run GWAS
Use run_gwas.sh to perform GWAS for all traits in the phenotype file.

```bash
bash scripts/run_gwas.sh \
  --group sample \
  --gemma /path/to/gemma \
  --bfile data/genotype/sample_prefix \
  --kfile results/grm/sample.snp.all.matrix.cXX.txt \
  --cfile data/covariates/c.txt \
  --pheno-file data/phenotype/phenotype.txt \
  --outdir results/GWAS_results \
  --start-col 3
```

For each trait, GEMMA produces association result files such as:

- <trait_name>.assoc.txt
- <trait_name>.log.txt

### 3. Calculate genomic inflation factor
Use calculate_lambda.py to summarize genomic inflation across all GWAS traits.

```python
python scripts/calculate_lambda.py \
  --input-dir results/GWAS_results \
  --output-file results/lambda/lambda_results.txt
```


### 4. Annotate loci
Use annotate_peaks.R to annotate genomic intervals by overlapping them with genome annotation.

```R
Rscript scripts/annotate_peaks.R \
  data/annotation/genome_annotation.txt \
  data/peaks/peaks.txt \
  results/annotation/peaks_annotation.txt \
  CHR BIN_START BIN_END
```

---

## Notes
The GRM generated in Step 1 is required for GWAS in Step 2.
The order of phenotype columns determines the trait analysis order.
The annotation file and peak file must use consistent chromosome naming.

---

## License
For academic and research use only.