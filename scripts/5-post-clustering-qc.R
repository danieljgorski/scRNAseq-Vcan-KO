# Quality control after initial clustering

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(ggrepel)

# Load object
load("results/objects/obj_integrated.Rdata")

# Set up output dirs
output_dirs <- c("results",
                 "results/post-clustering-qc")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Switch default assay to RNA
DefaultAssay(obj) <- "RNA"

# Dissociation related genes
diss_genes <- unique(c("Atf3", "Btg2", "Cebpb", "Cebpb",
                       "Cxcl3", "Cxcl2", "Cxcl1",
                       "Dnaja1", "Dnajb1", "Dusp1",
                       "Egr1", "Fos", "Fosb", "Hsp90aa1",
                       "Hsp90ab1", "Hspa1a", "Hspa1b",
                       "Hspa1a", "Hspa1b", "Hspa8",
                       "Hspb1", "Hspe1", "Hsph1", "Id3",
                       "Ier2", "Jun", "Junb", "Jund",
                       "Mt1", "Nfkbia", "Nr4a1", "Ppp1r15a",
                       "Socs3", "Zfp36"))
obj <- AddModuleScore(obj,
                      features = list(diss_genes),
                      ctrl = 50,
                      name = "diss_genes")
p1 <- FeaturePlot(obj, features = "diss_genes1", label = T, raster = F) +
  ggtitle("Dissociation gene score")
p2 <- VlnPlot(obj, features = "diss_genes1", pt.size = 0, sort = T) +
  ggtitle("Dissociation gene score") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NoLegend()
pdf(file = "results/post-clustering-qc/dissocation_genes.pdf",
    width = 14,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(ncol = 2)
dev.off()

# Batch effect exploration
p <- DimPlot(obj, group.by = "sample", raster = F)
pdf(file = "results/post-clustering-qc/DimPlot_sample.pdf")
print(p)
dev.off()

# nFeature_RNA
p1 <- FeaturePlot(obj, features = "nFeature_RNA", label = T, raster = F)
p2 <- VlnPlot(obj, features = "nFeature_RNA", pt.size = 0, sort = T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NoLegend()
pdf(file = "results/post-clustering-qc/nFeature_RNA.pdf",
    width = 14,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(ncol = 2)
dev.off()

# nCount_RNA
p1 <- FeaturePlot(obj, features = "nCount_RNA", label = T, raster = F)
p2 <- VlnPlot(obj, features = "nCount_RNA", pt.size = 0, sort = T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NoLegend()
pdf(file = "results/post-clustering-qc/nCount_RNA.pdf",
    width = 14,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(ncol = 2)
dev.off()

# percent.mt
p1 <- FeaturePlot(obj, features = "percent.mt", label = T, raster = F)
p2 <- VlnPlot(obj, features = "percent.mt", pt.size = 0, sort = T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NoLegend()
pdf(file = "results/post-clustering-qc/percent_mt.pdf",
    width = 14,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(ncol = 2)
dev.off()

# Cell-cycle phase
cc_genes <- read_csv("data/seurat_cell_cycle.csv")
s_genes <- as.character(cc_genes$mouse.s.genes)
g2m_genes <- as.character(cc_genes$mouse.g2m.genes)
remove(cc_genes)
obj <- CellCycleScoring(obj,
                        s.features = s_genes,
                        g2m.features = g2m_genes)
obj@meta.data$Phase <- factor(obj@meta.data$Phase,
                              levels = c("G1", "S", "G2M"))
phase_membership <- (prop.table(table(obj$Phase, Idents(obj)),
                                margin = 2) * 100)
phase_membership <- as.data.frame(phase_membership)
colnames(phase_membership) <- c("Phase", "Cluster", "Percent")
phase_membership$percent_round <- round(phase_membership$Percent)
write.csv(phase_membership,
          file = "results/post-clustering-qc/phase_membership.csv",
          row.names = F)
phase_membership$Phase <- factor(phase_membership$Phase,
                                 levels = c("G1", "S", "G2M"))
p1 <- DimPlot(obj,
              group.by = "Phase",
              cols = c("#e5e5e5", "#3a86ff", "#ffaa00"),
              raster = F)
p2 <- ggplot(phase_membership,
             aes(fill = Phase,
                 y = Percent,
                 x = Cluster)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Percent of total") +
  xlab("Identity") +
  geom_text(aes(label = paste0(percent_round, "%", sep = "")),
            position = position_stack(vjust = 0.5), size = 2.5) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "Black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = c("#e5e5e5", "#3a86ff", "#ffaa00"))
pdf(file = "results/post-clustering-qc/cell_cycle_phase.pdf",
    width = 17,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(guides = "collect")
dev.off()

# Find markers of initial clustering
cluster_markers <- FindAllMarkers(obj,
                                  assay = "RNA",
                                  logfc.threshold = 0.5,
                                  min.pct = 0.5,
                                  only.pos = T,
                                  return.thresh = 0.001,
                                  densify = T)
write.csv(cluster_markers,
          file = "results/post-clustering-qc/cluster_markers.csv",
          row.names = F)
top5 <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(top5,
          file = "results/post-clustering-qc/cluster_markers_top_5.csv",
          row.names = F)
top10 <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(top10,
          file = "results/post-clustering-qc/cluster_markers_top_10.csv",
          row.names = F)

# Counting cluster markers
marker_count <- cluster_markers %>%
  group_by(cluster) %>%
  summarise(nMarkers = n())
marker_count$cluster <- as.character(marker_count$cluster)

# Summarizing quality control
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
qc_summary <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(mean(percent.mt),
            mean(diss_genes1),
            mean(nFeature_RNA),
            mean(nCount_RNA),
            calculate_mode(Phase))
qc_summary <- qc_summary %>% left_join(marker_count,
                                      by = c("seurat_clusters" = "cluster"))
colnames(qc_summary) <- c("cluster",
                          "mean_percent_mt",
                          "mean_dissociation_score",
                          "mean_nFeatures",
                          "mean_nCounts",
                          "cell_cycle",
                          "nMarkers")
write.csv(qc_summary,
          file = "results/post-clustering-qc/qc_summary.csv",
          row.names = F)

# Plotting QC summary
p <- ggplot(qc_summary, aes(x = mean_percent_mt,
                       y = mean_nFeatures,
                       colour = nMarkers,
                       label = cluster)) +
  scale_colour_viridis_c() +
  geom_point(aes(shape = cell_cycle), size = 3) +
  geom_text_repel() +
  ggtitle("Initial clustering QC")
pdf(file = "results/post-clustering-qc/qc_summary_scatter.pdf",
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
                            label = cluster)) +
  scale_colour_viridis_c() +
  geom_point(aes(shape = cell_cycle), size = 3) +
  geom_text_repel() +
  ggtitle("Initial clustering QC")
pdf(file = "results/post-clustering-qc/qc_overview.pdf",
    width = 13,
    height = 10,
    useDingbats = F)
print((p1 / p2) | p3 )
dev.off()

# Simple annotation, to identify remaining possible heterotypic multiplets

# Read in markers
markers <- read.csv(file = "data/canonical_markers.csv")
markers <- markers$All

# Loop through markers and generate Feature and VlnPlots of expression
for (i in markers) {
  p1 <- FeaturePlot(obj, features = i, label = T, raster = F)
  p2 <- VlnPlot(obj, features = i, pt.size = 0, sort = T) + NoLegend() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  pdf(file = paste0("results/post-clustering-qc/Feat_VlnPlot_", i, ".pdf"),
      width = 16,
      height = 5,
      useDingbats = F)
  print(p1 + p2 + plot_layout(ncol = 2, widths = c(1,2)))
  dev.off()
}

# Typical cell-types with top canonical markers
# Fibroblast markers, Col1a1, Tcf21, Postn, Gsn, Cthrc1
# Macrophage markers, Fcgr1, Adgre1,
# EC markers, Cdh5, Kdr, Pecam1
# Granulocytes, S100a8, S100a9,
# Proliferating, Mki67, Ccnb2
# DC_like Cd209a, H2-Ab1, Cd74
# Cardiomyocyte, Nppa, Actc1
# Bcell Ms4a1, Cd79a
# Tcell Cd3e, Cd3d
# NKcell Ncr1, Klrk1, Ccl5
# Mural Myh11, Vtn
# Epicardial Wt1, Dmkn, Krt19, Krt8
# Glial Plp1, Kcna1

# Notes on canonical marker expression across clusters
# 0 - Fibroblast - Tcf21, Col1a1, Postn, Gsn (Fibroblast)
# 1 - Fibroblast - Tcf21, Col1a1, Postn, Gsn, Cthrc1
# 2 - Fibroblast - Tcf21, Col1a1, Postn, Gsn, Cthrc1
# 3 - Fibroblast - Tcf21, Col1a1, Postn, Gsn
# 4 - Fibroblast - Tcf21, Col1a1, Gsn
# 5 - Epicardial - Wt1, Dmkn, Clu
# 6 - Fibroblast - Tcf21, Col1a1, Postn, Gsn, Cthrc1, Clu
# 7 - Epicardial - Wt1, Dmkn, Cthrc1
# 8 - Epicardial - Wt1, Dmkn, Krt19, Krt8, Clu
# 9 - Fibroblast - Tcf21, Col1a1, Postn, Gsn, Cthrc1
# 10 - Epicardial - Cthrc1, Clu 
# 11 - Epicardial - Wt1, Cthrc1, Clu
# 12 - Epicardial, Proliferating - Wt1, Dmkn, Krt19, Krt8, Clu, Ccnb2, Stmn1, Mik67
# 13 - Endothelial - Kdr, Pecam1, Cd36, Cdh5
# 14 - Fibroblast - Tcf21, Col1a1, Postn, Gsn, Cthrc1
# 15 - Mural Cells - Rgs5, Myh11, Vtn
# 16 - Cardiomyocytes - Actc1, Tnnt2, Tnnc1
# 17 - Fibroblast, Proliferating - Tcf21, Mik67, Stmn1, Ccnb2
# 18 - Granulocytes - S100a9, S100a8,
# 19 - Glial - Kcna1, Plp1
# 20 - B-Cell - Cd79a, Cd74, H2-Aa

# Excluding non-fibroblast, non-epicardial cells that were captured despite
# CD31-CD45- sorting, as well as low quality fibroblast (9) and epicardial
# cells (10) with very low nFeature, nCount
obj <- subset(x = obj,
              idents = c("9",
                         "10",
                         "13",
                         "15",
                         "16",
                         "18",
                         "19",
                         "20"),
              invert = TRUE)

# Remove extra meta data columns created during integrated clustering
obj@meta.data <- select(obj@meta.data,
                        -c(nCount_SCT,
                           nFeature_SCT,
                           integrated_snn_res.0.8,
                           seurat_clusters))

# Keep only necessary counts and data slots
obj <- DietSeurat(obj,
                  counts = T,
                  data = T,
                  scale.data = F,
                  features = NULL,
                  assays = "RNA",
                  dimreducs = NULL,
                  graphs = NULL)

# Save clean object
save(obj, file = "results/objects/obj_clean.Rdata")

# Clear memory
rm(list = ls())
gc()
