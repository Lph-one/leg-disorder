This repository contains the code used for integrative analysis of serum metabolomics, genomic variation, and phenotypic traits in chickens with and without leg disorders.

The study is based on a balanced cohort comprising chickens with natural leg disease (n = 102) and healthy controls (n = 92). For each individual, whole-genome resequencing data (10× coverage), serum metabolite profiles, growth traits, and measurements of nine bone metabolism regulators were collected.

The analytical workflow includes:

- Weighted correlation network analysis（WGCNA）
- Estimation of genetic components (heritability, h²)
- Metabolite GWAS (mGWAS) to identify candidate loci
- Colocalization (coloc) analysis between mQTLs and molecular QTLs (molQTLs) from Chicken Sex GTEx
- Mendelian randomization (MR) analysis integrating mQTLs with complex traits from Chicken GTEx

This repository provides scripts to reproduce the main computational analyses described in the manuscript.
