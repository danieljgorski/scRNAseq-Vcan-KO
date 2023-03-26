# Re-clustering cleaned object, using reference based integration with RPCA

# Load libraries
library(Seurat)
library(clustree)

# Load object
load("results/objects/obj_clean.Rdata")

# Set up output dirs
output_dirs <- c("results",
                 "results/integrated-clustering")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Split object on samples
obj_list <- SplitObject(obj, split.by = "sample")

# SCTransform objects individually with sped up glmGamPoi method
obj_list <- lapply(X = obj_list, FUN = SCTransform, method = "glmGamPoi")

# Find integration features
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)

# Prep for SCT integration, run PCA on individual samples
obj_list <- PrepSCTIntegration(object.list = obj_list,
                               anchor.features = features)
obj_list <- lapply(X = obj_list, FUN = RunPCA, features = features)

# Find integration anchors, using WT samples as reference
controls <- c(2, 3) # indices of control samples
obj_anchors <- FindIntegrationAnchors(object.list = obj_list,
                                      normalization.method = "SCT",
                                      reference = controls,
                                      anchor.features = features,
                                      dims = 1:25,
                                      reduction = "rpca",
                                      k.anchor = 5)

# Save anchors, due to the large memory requirements of integration
save(obj_anchors, file = "results/objects/obj_anchors_clean.Rdata")
remove(obj)
remove(obj_list)

# Integrate data
obj <- IntegrateData(anchorset = obj_anchors,
                     normalization.method = "SCT",
                     dims = 1:25)

# Run PCA, UMAP and cluster integrated object
obj <- RunPCA(obj, npcs = 50, verbose = T)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:25, verbose = T)
obj <- FindNeighbors(obj, dims = 1:25, verbose = T)
obj <- FindClusters(obj, resolution = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2), verbose = T)

# Find the optimum clustering resolution with clustree
# https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html
clust <- clustree(obj, prefix = "integrated_snn_res.", layout = "sugiyama")
pdf(file ="results/integrated-clustering/clustree.pdf",
    useDingbats = F)
print(clust)
dev.off()

# Normalizing and scaling RNA assay
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000,
                     verbose = T)
obj <- ScaleData(obj, features = rownames(obj), verbose = T)

# Factor genotype level
obj@meta.data$genotype <- factor(obj@meta.data$genotype, levels = c("Vcan-WT", "Vcan-KO"))

# Save object
save(obj, file = "results/objects/obj_integrated_clean.Rdata")

# Clear memory
rm(list = ls())
gc()
