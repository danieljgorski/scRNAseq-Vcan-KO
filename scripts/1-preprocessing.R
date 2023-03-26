# Preprocessing count matrices

# Load libraries
library(Seurat)
library(dplyr)

# Set up output dirs
output_dirs <- c("results",
                 "results/objects",
                 "results/preprocessing")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}
  
# Load in cellranger aggr data
obj <- Read10X("data/300_NGS_111119_new_aggr/outs/filtered_feature_bc_matrix")
obj <- CreateSeuratObject(counts = obj,
                          project = "300_NGS_111119",
                          min.cells = 3,
                          min.features = 200,
                          names.field = 2,
                          names.delim = "-")

# Add sample metadata
obj@meta.data <- obj@meta.data %>%
  mutate(sample = case_when(orig.ident == "1" ~ "Blau4KO",
                            orig.ident == "2" ~ "Blau5Control",
                            orig.ident == "3" ~ "Rot2Control",
                            orig.ident == "4" ~ "Rot4KO",
                            orig.ident == "5" ~ "Rot5KO"))
obj@meta.data <- obj@meta.data %>%
  mutate(genotype = case_when(orig.ident == "1" ~ "Vcan-KO",
                              orig.ident == "2" ~ "Vcan-WT",
                              orig.ident == "3" ~ "Vcan-WT",
                              orig.ident == "4" ~ "Vcan-KO",
                              orig.ident == "5" ~ "Vcan-KO"))

# Add percent mitochondrial reads
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")


# Pre-filter 
cells_pre_qc <- length(colnames(obj))
p <- VlnPlot(obj,
             features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             group.by = "sample",
             pt.size = 0)

# Plot pre-filer QC metrics
pdf(file = "results/preprocessing/VlnPlot_QC_metrics_pre-filter.pdf",
    width = 12,
    height = 6,
    pointsize = 12,
    useDingbats = F)
print(p)
dev.off()

# Filter
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 10)

# Post-filter
cells_post_qc <- length(colnames(obj))
percent_passed <- (cells_post_qc / cells_pre_qc) * 100
qc_filter <- as.data.frame(cells_pre_qc)
qc_filter$cells_post_qc <- cells_post_qc
qc_filter$percent_passed <- percent_passed

# Write filter percentage
write.csv(qc_filter,
          file = "results/preprocessing/qc_filter.csv",
          row.names = F)
p <- VlnPlot(obj,
             features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             group.by = "sample",
             pt.size = 0)

# Plot post-filer QC metrics
pdf(file = "results/preprocessing/VlnPlot_QC_metrics_post-filter.pdf",
    width = 12,
    height = 6,
    pointsize = 12,
    useDingbats = F)
print(p)
dev.off()

# Save object
save(obj, file = "results/objects/obj_preprocessed.Rdata")

# Clear memory
rm(list = ls())
gc()