library(data.table)
library(TwoSampleMR)

args <- commandArgs(TRUE)
COHORT <- args[1]

ROOT <- normalizePath(".")

exp_dir <- file.path(ROOT, "results", paste0(COHORT, "_exp"))
outcome_dir <- file.path(ROOT, "data", "gwas_raw")
output_dir <- file.path(ROOT, "results", "MR")

dir.create(output_dir, showWarnings = FALSE)

exposure_n <- xxx ##(xxx is Sample Size)
outcome_n <- xxx ##(xxx is Sample Size)

exp_files <- list.files(exp_dir, full.names = TRUE)
out_files <- list.files(outcome_dir, full.names = TRUE)

res_all <- data.frame()

for (exp_file in exp_files) {

  exp <- fread(exp_file)
  if (nrow(exp) == 0) next

  exp_dat <- data.frame(
    SNP = exp$rs,
    beta.exposure = exp$beta,
    se.exposure = exp$se,
    effect_allele.exposure = exp$allele1,
    other_allele.exposure = exp$allele0,
    eaf.exposure = exp$af,
    pval.exposure = exp$p_wald,
    samplesize.exposure = exposure_n
  )

  for (of in out_files) {

    out <- fread(of)

    sub <- out[out$SNP %in% exp_dat$SNP,]
    if (nrow(sub) == 0) next

    out_dat <- data.frame(
      SNP = sub$SNP,
      beta.outcome = sub$BETA,
      se.outcome = sub$SE,
      effect_allele.outcome = sub$A1,
      other_allele.outcome = sub$A2,
      eaf.outcome = sub$AF1,
      pval.outcome = sub$P,
      samplesize.outcome = outcome_n
    )

    dat <- harmonise_data(exp_dat, out_dat)
    if (nrow(dat) == 0) next

    res <- mr(dat)
    res_all <- rbind(res_all, res)
  }
}

fwrite(res_all, file.path(output_dir, paste0(COHORT, "_MR_results.csv")))