#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript coloc_eqtl.R <gwas_file> <qtl_root> <out_root> <sample_size_file> <num_cores> <gwas_n>")
}

gwas_path     <- args[1]
qtl_root      <- args[2]
out_root      <- args[3]
sample_size_f <- args[4]
num_cores     <- as.integer(args[5])
gwas_N        <- as.integer(args[6])

library(data.table)
library(dplyr)
library(coloc)
library(tools)
library(parallel)

pheno_name <- file_path_sans_ext(basename(gwas_path))
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
out_dir_pheno <- file.path(out_root, pheno_name)
dir.create(out_dir_pheno, showWarnings = FALSE, recursive = TRUE)

cat("===== coloc analysis started =====\n")
cat("GWAS file:", gwas_path, "\n")
cat("QTL root :", qtl_root, "\n")
cat("Out dir  :", out_dir_pheno, "\n")
cat("Cores    :", num_cores, "\n")
cat("GWAS N   :", gwas_N, "\n")

safe_name <- function(x) {
  x <- gsub("[/\\\\:*?\"<>|]", "_", as.character(x))
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

# ---- Read GWAS ----
gwas <- fread(gwas_path)
required_gwas_cols <- c("rs", "af", "beta", "se", "p_wald")
if (!all(required_gwas_cols %in% names(gwas))) {
  stop("GWAS file must contain columns: ", paste(required_gwas_cols, collapse = ", "))
}

gwas <- gwas %>%
  mutate(maf = ifelse(af <= 0.5, af, 1 - af)) %>%
  rename(SNP = rs) %>%
  select(SNP, beta, se, p_wald, maf)

setDT(gwas)
setkey(gwas, SNP)

# ---- Read sample sizes ----
qtl_ss <- fread(sample_size_f, header = FALSE, sep = "\t")
names(qtl_ss) <- c("Tissue", "SampleSize")
qtl_ss <- qtl_ss %>%
  mutate(Tissue = trimws(Tissue), SampleSize = as.integer(SampleSize)) %>%
  filter(!is.na(SampleSize), SampleSize > 0)

# ---- Find eQTL files ----
qtl_files <- list.files(
  path = qtl_root,
  pattern = "\\.cis_qtl_pairs\\.significant\\.txt\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(qtl_files) == 0) {
  stop("No QTL files found under: ", qtl_root)
}

qtl_paths <- data.table(
  Tissue = basename(dirname(dirname(qtl_files))),
  FilePath = qtl_files
)
qtl_paths <- merge(qtl_paths, qtl_ss, by = "Tissue", all.x = TRUE)
qtl_paths <- qtl_paths[!is.na(SampleSize)]

run_coloc_pair <- function(tissue, qtl_path, qtl_N) {
  qtl <- fread(qtl_path,
               select = c("pheno_id", "variant_id", "af", "beta_g1", "beta_se_g1", "pval_g1"))
  names(qtl) <- make.names(names(qtl), unique = TRUE)

  qtl <- qtl %>%
    mutate(af = ifelse(af <= 0.5, af, 1 - af)) %>%
    rename(
      phenotype_id = pheno_id,
      slope = beta_g1,
      slope_se = beta_se_g1,
      pval_nominal = pval_g1
    )

  setDT(qtl)
  setkey(qtl, variant_id)

  kept_genes <- 0L
  out_summary <- list()

  for (gene in unique(qtl$phenotype_id)) {
    qtl_target <- qtl[phenotype_id == gene]
    common_snps <- intersect(qtl_target$variant_id, gwas$SNP)
    if (!length(common_snps)) next

    df_target <- merge(
      gwas[J(common_snps)],
      qtl_target[J(common_snps)],
      by.x = "SNP",
      by.y = "variant_id"
    )

    if (!nrow(df_target)) next
    if (anyNA(df_target[, .(beta, se, maf, slope, slope_se, af)])) next

    df1 <- list(
      beta = df_target$beta,
      varbeta = df_target$se^2,
      snp = df_target$SNP,
      type = "quant",
      N = gwas_N,
      MAF = df_target$maf
    )

    df2 <- list(
      beta = df_target$slope,
      varbeta = df_target$slope_se^2,
      snp = df_target$SNP,
      type = "quant",
      N = qtl_N,
      MAF = df_target$af
    )

    res <- tryCatch(
      coloc.abf(dataset1 = df1, dataset2 = df2),
      error = function(e) {
        message("[ERROR] coloc failed: Tissue=", tissue, " Gene=", gene, " -> ", e$message)
        NULL
      }
    )
    if (is.null(res)) next
    if (as.numeric(res$summary["PP.H4.abf"]) <= 0.5) next

    kept_genes <- kept_genes + 1

    gene_tag <- safe_name(gene)
    tissue_tag <- safe_name(tissue)
    pheno_tag <- safe_name(pheno_name)
    lead_snp <- df_target$SNP[which.min(df_target$p_wald)]

    sum_df <- as.data.frame(t(unlist(res$summary)))
    sum_df$lead_SNP <- lead_snp
    sum_df$gene <- gene
    sum_df$tissue <- tissue
    sum_df$phenotype <- pheno_name

    det_df <- as.data.frame(res$results)
    det_df$gene <- gene
    det_df$tissue <- tissue
    det_df$phenotype <- pheno_name

    fwrite(sum_df, file.path(out_dir_pheno,
                             paste0(pheno_tag, "__", gene_tag, "__", tissue_tag, "_summary.txt")),
           sep = "\t")
    fwrite(det_df, file.path(out_dir_pheno,
                             paste0(pheno_tag, "__", gene_tag, "__", tissue_tag, "_details.txt")),
           sep = "\t")

    out_summary[[length(out_summary) + 1]] <- sum_df
  }

  if (length(out_summary) > 0) {
    summary_df <- rbindlist(out_summary, fill = TRUE)
    fwrite(summary_df, file.path(out_dir_pheno, paste0(safe_name(tissue), "_summary_all.txt")), sep = "\t")
  }

  data.table(Tissue = tissue, Genes_kept = kept_genes)
}

results_all <- mclapply(
  seq_len(nrow(qtl_paths)),
  function(i) {
    run_coloc_pair(qtl_paths$Tissue[i], qtl_paths$FilePath[i], qtl_paths$SampleSize[i])
  },
  mc.cores = num_cores
)

results_all <- rbindlist(results_all, fill = TRUE)
fwrite(results_all, file.path(out_dir_pheno, "coloc_run_summary.txt"), sep = "\t")

cat("===== coloc analysis finished =====\n")
cat("GWAS =", pheno_name, "\n")