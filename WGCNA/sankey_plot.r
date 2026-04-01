library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)

## Read heatmap row annotation from two groups
group1 <- read.csv(file.path("output", "CT_WGCNA", "CT_Heatmap_row_info.csv"), stringsAsFactors = FALSE)
group2 <- read.csv(file.path("output", "CA_WGCNA", "CA_Heatmap_row_info.csv"), stringsAsFactors = FALSE)

## Merge matched metabolites
merged <- merge(group1, group2, by = "Metabolite", suffixes = c("_G1", "_G2"))
merged$ModuleColor_G1 <- factor(merged$ModuleColor_G1, levels = unique(group1$ModuleColor))
merged$ModuleColor_G2 <- factor(merged$ModuleColor_G2, levels = unique(group2$ModuleColor))

data_l <- merged %>% select(ModuleColor_G1, ModuleColor_G2) %>% rename(Source = ModuleColor_G1, Target = ModuleColor_G2)
data_lodes <- to_lodes_form(data_l, key = "x", value = "stratum", id = "alluvium", axes = 1:2)

color_map <- c("black" = "#3B3A3D", "blue" = "#195FB4", "brown" = "#5F3F33", "green" = "#055E1D", "greenyellow" = "#B9E0A2", "grey" = "#888B8A", "magenta" = "#C63187", "pink" = "#F6C4CF", "purple" = "#AC6BB1", "red" = "#74000B", "turquoise" = "#4CC4D4", "yellow" = "#E4AA21")

p1 <- ggplot(data_lodes, aes(x = x, stratum = stratum, alluvium = alluvium, fill = stratum)) +
  geom_flow(alpha = 0.4, width = 0, knot.pos = 0.3) +
  geom_stratum(width = 0.075, aes(color = stratum)) +
  scale_x_discrete(labels = c("", "")) +
  scale_fill_manual(values = color_map) +
  scale_color_manual(values = color_map) +
  guides(fill = FALSE) +
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

p1