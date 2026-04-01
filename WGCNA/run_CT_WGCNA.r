library(WGCNA)
options(stringsAsFactors = FALSE)

group_name <- "CT"

## WGCNA for selected group
# Read metabolite data
metaData <- read.csv(paste0(group_name, "meta.csv"), stringsAsFactors = FALSE)
datExpr0 <- as.data.frame(t(metaData[, -1]))
names(datExpr0) <- metaData$ID
rownames(datExpr0) <- names(metaData)[-1]

# Check sample and metabolite quality
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  if(sum(!gsg$goodGenes) > 0)
    print(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ",")))
  if(sum(!gsg$goodSamples) > 0)
    print(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ",")))
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Recheck data
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) stop("Data still not clean")

# Sample clustering to detect outliers
sampleTree <- hclust(dist(datExpr0), method = "average")
plot(sampleTree, main = paste(group_name, "sample clustering to detect outliers"), xlab = "", sub = "")

# Remove outlier samples by cut height
clust <- cutreeStatic(sampleTree, cutHeight = 8e6, minSize = 10)
table(clust)
keepSamples <- (clust == 1)
datExpr <- datExpr0[keepSamples, ]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
print(paste("Samples:", nSamples, "Genes(metabolites):", nGenes))

## 2. Select soft-thresholding power
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot soft-threshold results
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (R^2)",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.9, col = "red")

plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

# Choose soft-threshold power
softPower <- 6

## 3. Construct network and detect modules
# Adjacency matrix and TOM similarity
adjacency <- adjacency(datExpr, power = softPower)
TOM <- TOMsimilarity(adjacency)
rownames(TOM) <- colnames(TOM) <- names(datExpr)
dissTOM <- 1 - TOM

# Hierarchical clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = paste(group_name, "metabolite clustering on TOM-based dissimilarity"),
     xlab = "", sub = "", labels = FALSE, hang = 0.04)

## Dynamic tree cut for module detection
minModuleSize <- 100
dynamicMods <- cutreeDynamic(dendro = geneTree,
                             distM = dissTOM,
                             deepSplit = 2,
                             pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

# Visualize module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = paste(group_name, "metabolite dendrogram and module colors"))

## 4. Calculate module eigengenes and merge modules
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres <- 0.2
abline(h = MEDissThres, col = "red")

merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# Visualize merged modules
plotDendroAndColors(geneTree,
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged modules"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = paste(group_name, "metabolite dendrogram and merged module colors"))

# Update final results
moduleColors <- mergedColors
MEs <- mergedMEs

## Set output directory
outputDir <- file.path("output", paste0(group_name, "_WGCNA"))
if (!dir.exists(outputDir)) {dir.create(outputDir, recursive = TRUE)}

## Save module results
save(datExpr, moduleColors, MEs, geneTree, TOM,
     file = file.path(outputDir, paste0(group_name, "_WGCNA_metabolite_modules.RData")))

## Export network edges and nodes for each module
threshold <- 0.05
metaboliteNames <- colnames(datExpr)

uniqueModules <- unique(moduleColors)
for (module in uniqueModules) {
  if (module == "grey") next
  probes <- metaboliteNames[moduleColors == module]
  moduleTOM <- TOM[probes, probes]
  moduleTOM[moduleTOM < threshold] <- 0
  edgeFile <- file.path(outputDir, paste0(group_name, "_CytoscapeInput_edges_", module, ".txt"))
  nodeFile <- file.path(outputDir, paste0(group_name, "_CytoscapeInput_nodes_", module, ".txt"))
  cyt <- exportNetworkToCytoscape(
    moduleTOM,
    edgeFile = edgeFile,
    nodeFile = nodeFile,
    weighted = TRUE,
    threshold = threshold,
    nodeNames = probes,
    altNodeNames = probes,
    nodeAttr = module)
  cat("Module", module, "exported to", outputDir, "with", length(probes), "metabolites\n")
}

## Export TOM matrix and module annotation
subsetTOM <- TOM
diag(subsetTOM) <- NA

write.csv(subsetTOM, file = file.path(outputDir, paste0(group_name, "_Metabolite_TOM_matrix.csv")), quote = FALSE, row.names = TRUE)

module_info <- data.frame(
  Metabolite = colnames(subsetTOM),
  ModuleColor = moduleColors)
write.csv(module_info, file = file.path(outputDir, paste0(group_name, "_Metabolite_colors.csv")), row.names = FALSE)