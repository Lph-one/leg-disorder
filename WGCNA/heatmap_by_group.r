library(pheatmap)
library(dplyr)
library(grid)

plot_heatmap_by_group <- function(group_name) {
  ## Read TOM matrix and module color data
  TOM_mat <- read.csv(file.path("output", paste0(group_name, "_WGCNA"), paste0(group_name, "_Metabolite_TOM_matrix.csv")), row.names = 1, check.names = FALSE)
  TOM_mat <- as.matrix(TOM_mat)

  mod_colors <- read.csv(file.path("output", paste0(group_name, "_WGCNA"), paste0(group_name, "_Metabolite_colors.csv")), sep = ",", stringsAsFactors = FALSE)

  ## Match metabolites
  common_ids <- intersect(rownames(TOM_mat), mod_colors$Metabolite)
  TOM_mat <- TOM_mat[common_ids, common_ids]
  mod_colors <- mod_colors[match(common_ids, mod_colors$Metabolite), ]
  order_ids <- mod_colors %>% arrange(ModuleColor) %>% pull(Metabolite)
  TOM_mat <- TOM_mat[order_ids, order_ids]
  mod_colors <- mod_colors[match(order_ids, mod_colors$Metabolite), ]

  ## Plot heatmap
  col_fun <- colorRampPalette(c("#F6EDED", "#e36e8C"))(100)

  p1 <- pheatmap(TOM_mat, color = col_fun, cluster_rows = FALSE, cluster_cols = FALSE,
                 show_rownames = FALSE, show_colnames = FALSE, main = "", border_color = NA)
  p1

  grid.force()
  vp_names <- grep("matrix", grid.ls(viewports = TRUE, print = FALSE)$name, value = TRUE)
  seekViewport(vp_names[1])

  module_sizes <- table(mod_colors$ModuleColor)
  module_edges <- cumsum(module_sizes)
  module_start <- c(0, module_edges[-length(module_edges)])
  n <- nrow(TOM_mat)
  m <- ncol(TOM_mat)

  for (i in seq_along(module_edges)) {
    start <- module_start[i]
    end <- module_edges[i]
    grid.rect(
      x = (start + end) / 2 / m,
      y = 1 - (start + end) / 2 / n,
      width = (end - start) / m,
      height = (end - start) / n,
      gp = gpar(lwd = 2, col = "black", fill = NA))
  }

  grid.rect(x = 0.5, y = 0.5, width = 1, height = 1, gp = gpar(lwd = 2, col = "black", fill = NA))

  ## Export row annotation
  row_info <- data.frame(
    RowIndex = seq_len(nrow(TOM_mat)),
    Metabolite = rownames(TOM_mat),
    ModuleColor = mod_colors$ModuleColor)

  write.csv(row_info, file.path("output", paste0(group_name, "_WGCNA"), paste0(group_name, "_Heatmap_row_info.csv")), row.names = FALSE)
}

plot_heatmap_by_group("CA")
plot_heatmap_by_group("CT")