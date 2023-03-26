# Basic cluster annotation based on canonical cell-type markers

# Load libraries
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ape)

# Load object
load("results/objects/obj_integrated_clean.Rdata")

# Set up output dirs
output_dirs <- c("results",
                 "results/basic-annotation")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Set clustering idents to resolution 0.6, based on the clustree, this seems
# to be stable
Idents(obj) <- "integrated_snn_res.0.6"

# Calculate seurat_cluster markers, with high stringency to return very
# specific markers, quickly.
seurat_cluster_markers <- FindAllMarkers(obj,
                                  assay = "RNA",
                                  logfc.threshold = 1.5,
                                  min.pct = 0.5,
                                  only.pos = T,
                                  return.thresh = 0.0001,
                                  densify = T,
                                  verbose = T)
write.csv(seurat_cluster_markers,
          file = "results/basic-annotation/seurat_cluster_markers.csv",
          row.names = F)
top5 <- seurat_cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Export a basic DimPlot of clusters
p <- DimPlot(obj, label = T, raster = F)
pdf(file = "results/basic-annotation/DimPlot.pdf",
    useDingbats = F)
print(p)
dev.off()

# Loop through top5 cluster markers and generate Feature and VlnPlots of expression
for (i in top5$gene) {
  p1 <- FeaturePlot(obj, features = i, label = T, raster = F)
  p2 <- VlnPlot(obj, features = i, pt.size = 0, sort = T) + NoLegend() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  pdf(file = paste0("results/basic-annotation/Feat_VlnPlot_", i, ".pdf"),
      width = 16,
      height = 5,
      useDingbats = F)
  print(p1 + p2 + plot_layout(ncol = 2, widths = c(1,2)))
  dev.off()
}

# Read in markers
markers <- read.csv(file = "data/canonical_markers.csv")
markers <- markers$All

# Loop through canonical markers and generate Feature and VlnPlots of expression
for (i in markers) {
  p1 <- FeaturePlot(obj, features = i, label = T, raster = F)
  p2 <- VlnPlot(obj, features = i, pt.size = 0, sort = T) + NoLegend() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  pdf(file = paste0("results/basic-annotation/Feat_VlnPlot_", i, ".pdf"),
      width = 16,
      height = 5,
      useDingbats = F)
  print(p1 + p2 + plot_layout(ncol = 2, widths = c(1,2)))
  dev.off()
}

# Dendrogram
obj <- BuildClusterTree(obj, dims = 1:25)
PlotClusterTreeDJG <- function(object, ...) {
  if (is.null(x = Tool(object = object, slot = "BuildClusterTree"))) {
    stop("Phylogenetic tree does not exist, build using BuildClusterTree")
  }
  data.tree <- Tool(object = object, slot = "BuildClusterTree")
  plot.phylo(x = data.tree, font = 2, direction = "rightwards", label.offset = 2, edge.width = 1)
}
pdf(file = "results/basic-annotation/Dendrogram.pdf",
    useDingbats = F)
PlotClusterTreeDJG(obj)
dev.off()


###############################################################################
# Annotation notes
###############################################################################

#-0-Myofibroblast (Fibro-Myo-1): Tcf21+, Cthrc1+, Postn+, Col1a1-high

#-3-Myofibroblast (Fibro-Myo-2): Tcf21+, Cthrc1+, Postn+, Col1a1-high

#-8-Myofibroblast (Fibro-Myo-3): Tcf21+, Cthrc1+, Postn+, Col1a1-high

#-11-Interferon Stimulated Fibroblasts (Fibro-IFN): Tcf21+, Ifit3+, Ifti1+, Ifit3b+, Igs15+

#-1-Activating Fibroblasts (Fibro-Act-1): Tcf21+, Postn+, but no Cthrc1 yet, lower Col1a1

#-2-Activating Fibroblasts (Fibro-Act-2): Tcf21+, Postn+, but no Cthrc1 yet, lower Col1a1

#-5-Resting Fibroblasts (Fibro-Rest): Tc21+ (highest), Gsn-High, Dcn-High, 
#     Dpep1+ (Farbehi et al 2019, Forte et al. 2020)

#-12-Cycling Fibroblasts (Fibro-Cyc): Tcf21+, Mki67+, Ccnb2+

#-4-Activated Epicardial cells (Epi-Act-1): Wt1+, Dmnk1+, Clu+, Lcn2+, Saa3+, Msln-

#-7-Activated Epicardial cells (Epi-Act-2): Wt1+, Dmnk1-low, Clu+, Saa3+, Msln-

#-9-Activated Epicardial cells (Epi-Act-3): Myofibroblastic, but likely epicardial origin 
#     Cthrc1+, Postn+, Acta2-High, Tcf21-, Col1a1-high, Wt1-low, Dmkn-low,
#     Clu-low

#-6-Resting Epicardial cells (Epi-Rest): Wt1+, Dmkn+, Clu+, Lcn2+, Saa3+, Krt19+, Gpm6a+
#     Msln+ Epicardial marker confined to epicardium even after MI, Zhou et al

#-10-Cycling Epicardial cells (Epi-Cyc): Wt1+, Dmkn+, Clu+, Saa3+, Serpinb2+, Krt19+, Gpm6a+
#     Mki67+, Ccnb2+

###############################################################################

# Rename Idents to annotations
obj <- RenameIdents(obj,
                    "0" = "Fibro-Myo-1",
                    "1" = "Fibro-Act-1",
                    "2" = "Fibro-Act-2",
                    "3" = "Fibro-Myo-2",
                    "4" = "Epi-Act-1",
                    "5" = "Fibro-Rest",
                    "6" = "Epi-Rest",
                    "7" = "Epi-Act-2",
                    "8" = "Fibro-Myo-3",
                    "9" = "Epi-Act-3",
                    "10" = "Epi-Cyc",
                    "11" = "Fibro-IFN",
                    "12" = "Fibro-Cyc")

# Store renamed idents as a new meta data column, set as Idents
obj@meta.data$basic_annotation <- Idents(obj)

# Refactor annotation levels
source("scripts/etc/dimplotlevels.R")
obj@meta.data$basic_annotation <- factor(obj@meta.data$basic_annotation,
                                         levels = dimplotlevels)
DimPlot(obj,
        group.by = "basic_annotation",
        label = T,
        repel = T)

# Set Idents as re-factored basic_annotation identities
Idents(obj) <- "basic_annotation"

# Save object with basic annotations
save(obj, file = "results/objects/obj_annotated.Rdata")

# Saved basic annotation, barcodes and UMAP embeddings etc. for consistency in
# external usage, de-comment to overwrite
# barcodes <- rownames(obj@meta.data)
# annotation <- obj@meta.data$basic_annotation
# genotype <- obj@meta.data$genotype
# sample <- obj@meta.data$sample
# umap_1 <- Embeddings(obj[["umap"]])[,1]
# umap_2 <- Embeddings(obj[["umap"]])[,2]
# basic_annotation <- data.frame(barcodes,
#                                annotation,
#                                genotype,
#                                sample,
#                                umap_1,
#                                umap_2)
# write.csv(basic_annotation, file = "data/basic_annotation.csv", row.names = F)

# Quality control summary of basic annotated object

# Counting cluster markers
marker_count <- seurat_cluster_markers %>%
  group_by(cluster) %>%
  summarise(nMarkers = n())
marker_count$cluster <- as.character(marker_count$cluster)

# Finding most common phase designation
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

# Summarizing QC metrics
qc_summary <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(mean(percent.mt),
            mean(nFeature_RNA),
            mean(nCount_RNA),
            calculate_mode(basic_annotation))
qc_summary <- qc_summary %>% left_join(marker_count,
                                       by = c("seurat_clusters" = "cluster"))

colnames(qc_summary) <- c("cluster",
                          "mean_percent_mt",
                          "mean_nFeatures",
                          "mean_nCounts",
                          "basic_annotation",
                          "nMarkers")
write.csv(qc_summary,
          file = "results/basic-annotation/final_qc_summary.csv",
          row.names = F)

# Plotting QC summary
p <- ggplot(qc_summary, aes(x = mean_percent_mt,
                            y = mean_nFeatures,
                            colour = nMarkers,
                            label = basic_annotation)) +
  scale_colour_viridis_c() +
  geom_point(size = 3) +
  geom_text_repel() +
  ggtitle("Final clustering QC")
pdf(file = "results/basic-annotation/final_qc_summary_scatter.pdf",
    width = 8,
    height = 6,
    useDingbats = F)
print(p)
dev.off()

# Plotting QC overview summary + nFeature + Clusters
p1 <- DimPlot(obj, label = T, raster = F) + NoLegend()
p2 <- FeaturePlot(obj,
                  features = "nFeature_RNA",
                  label = T,
                  raster = F) + NoLegend()
p3 <- ggplot(qc_summary, aes(x = mean_percent_mt,
                             y = mean_nFeatures,
                             colour = nMarkers,
                             label = basic_annotation)) +
  scale_colour_viridis_c() +
  geom_point(size = 3) +
  geom_text_repel() +
  ggtitle("Final clustering QC")
pdf(file = "results/basic-annotation/final_qc_overview.pdf",
    width = 13,
    height = 10,
    useDingbats = F)
print((p1 / p2) | p3 )
dev.off()
