#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ChIPpeakAnno)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  cat("
Usage:
  Rscript scripts/annotate_peaks.R <annotation_file> <peak_file> <output_file> <chr_col> <start_col> <end_col>

Example:
  Rscript scripts/annotate_peaks.R data/annotation/Gallus_7gff.txt data/peaks/CAspecific.txt results/CA/annotation/CAspecific_annotation.txt CHR BIN_START BIN_END
")
  quit(status = 1)
}

annotation_file <- args[1]
peak_file <- args[2]
output_file <- args[3]
chr_col <- args[4]
start_col <- args[5]
end_col <- args[6]

if (!file.exists(annotation_file)) {
  stop(paste("Annotation file not found:", annotation_file))
}
if (!file.exists(peak_file)) {
  stop(paste("Peak file not found:", peak_file))
}

options(scipen = 100)

cat("[INFO] Reading annotation file...\n")
ann_df <- read.table(annotation_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

required_ann_cols <- c("Chr", "Gene_start", "Gene_end", "Gene_ID")
missing_ann <- setdiff(required_ann_cols, colnames(ann_df))
if (length(missing_ann) > 0) {
  stop(paste("Missing required columns in annotation file:", paste(missing_ann, collapse = ", ")))
}

ann <- GRanges(
  seqnames = as.character(ann_df$Chr),
  ranges = IRanges(start = as.numeric(ann_df$Gene_start), end = as.numeric(ann_df$Gene_end)),
  id = as.character(ann_df$Gene_ID)
)

cat("[INFO] Reading peak file...\n")
peak_df <- read.table(peak_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

required_peak_cols <- c(chr_col, start_col, end_col)
missing_peak <- setdiff(required_peak_cols, colnames(peak_df))
if (length(missing_peak) > 0) {
  stop(paste("Missing required columns in peak file:", paste(missing_peak, collapse = ", ")))
}

if (!"peak_ID" %in% colnames(peak_df)) {
  peak_df$peak_ID <- paste0("peak_", seq_len(nrow(peak_df)))
}

peak <- GRanges(
  seqnames = as.character(peak_df[[chr_col]]),
  ranges = IRanges(start = as.numeric(peak_df[[start_col]]), end = as.numeric(peak_df[[end_col]])),
  id = as.character(peak_df$peak_ID)
)

cat("[INFO] Finding overlaps...\n")
hits <- findOverlaps(ann, peak)

if (length(hits) == 0) {
  cat("[WARN] No overlaps found.\n")
  out_df <- data.frame(
    peak_ID = character(),
    CHR = character(),
    BIN_START = numeric(),
    BIN_END = numeric(),
    Gene_ID = character(),
    stringsAsFactors = FALSE
  )
  write.table(out_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  quit(status = 0)
}

ann_hits <- ann_df[queryHits(hits), , drop = FALSE]
peak_hits <- peak_df[subjectHits(hits), , drop = FALSE]

out_df <- data.frame(
  peak_ID = peak_hits$peak_ID,
  CHR = peak_hits[[chr_col]],
  BIN_START = peak_hits[[start_col]],
  BIN_END = peak_hits[[end_col]],
  Gene_ID = ann_hits$Gene_ID,
  stringsAsFactors = FALSE
)

out_dir <- dirname(output_file)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

write.table(out_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("[INFO] Annotation results saved to:", output_file, "\n")