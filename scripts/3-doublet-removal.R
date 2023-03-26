# Doublet removal

# Load libraries
library(Seurat)
library(dplyr)

# Load object
load("results/objects/obj_db_classified.Rdata")

# Set up output dirs
output_dirs <- c("results",
                 "results/doublet-removal")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Cluster merged object to visualize doublets
obj <- SCTransform(obj, verbose = T)
obj <- RunPCA(obj, verbose = T)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:25, verbose = T)
obj <- FindNeighbors(obj, dims = 1:25, verbose = T)
obj <- FindClusters(obj, resolution = 0.8, verbose = T)
obj@meta.data$doublet_classification <- factor(obj@meta.data$doublet_classification,
                                               levels = c("Singlet", "Doublet"))
p <- DimPlot(obj, group.by = "doublet_classification", raster = F)
pdf(file = "results/doublet-removal/DimPlot_doublet_classification.pdf",
    useDingbats = F)
print(p)
dev.off()

# Subset for singlets
obj <- subset(x = obj, subset = doublet_classification == "Singlet")

# Remove extra meta data columns
obj@meta.data <- select(obj@meta.data,
                        -contains("pANN"),
                        -contains("SCT"),
                        -c(seurat_clusters, doublet_classification))

# Keep only necessary counts and data slots
DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj,
                  counts = T,
                  data = T,
                  scale.data = F,
                  features = NULL,
                  assays = "RNA",
                  dimreducs = NULL,
                  graphs = NULL)

# Save slimmed down object with no doublets
save(obj, file = "results/objects/obj_db_removed.Rdata")

# Clear memory
rm(list = ls())
gc()
