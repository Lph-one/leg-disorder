# WGCNA Analysis and Visualization for CA and CT Groups

A reproducible workflow for weighted correlation network analysis (WGCNA) and cross-group module comparison using metabolomics data.

---

## Overview

This repository provides a reproducible workflow for performing WGCNA on metabolomics data from two groups (CA and CT), together with downstream visualization including clustered heatmaps and Sankey plots for cross-group module comparison.

The analysis starts from two input datasets:

- `CAmeta.csv`
- `CTmeta.csv`

In these files:
- Rows represent metabolites  
- Columns represent samples  

---

## Scope of this repository
These scripts do not implement new statistical methods. Instead, they integrate existing tools into a consistent and reproducible workflow.

---

## Requirements
The following R packages are required:

- WGCNA
- pheatmap
- dplyr
- grid
- tidyr
- ggplot2
- ggalluvial

---

## Workflow
### 1. Perform WGCNA separately for CA and CT 
- run_CA_WGCNA.r
- run_CT_WGCNA.r

### 2. Generate TOM-based heatmaps  
heatmap_by_group.r

### 3. Compare modules across groups using Sankey plots  
sankey_plot.r

---

## Output
- TOM matrices
- Module color assignments
- Module eigengenes
- ytoscape network files




---


## License
For academic and research use only.

