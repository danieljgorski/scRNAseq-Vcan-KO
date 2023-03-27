# Miscellaneous figures for the manuscript

# Load libraries
library(Seurat) #>=4.0.1
library(scater)
library(stringr)
library(readr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggrepel)
library(viridis)
source("scripts/etc/colors.R")
source("scripts/etc/geno_colors.R")
source("scripts/etc/dimplotlevels.R")

# Load object, slim down and save
load("results/objects/obj_annotated.Rdata")
obj <- DietSeurat(obj, counts = F, dimreducs = "umap")
save(obj, file = "results/objects/obj_annotated_diet.Rdata")

# Metadata to df for ggplot2 QC plots
df <- obj@meta.data
df$timepoint <- factor(df$timepoint, levels = c("day_3", "day_7"),
                       labels = c("3 d", "7 d"))
df$sample <- factor(df$sample, levels = c("S47_B5_WT",
                                          "S47_B6_WT",
                                          "S48_R5_WT",
                                          "S47_B1_OE",
                                          "S47_B3_OE",
                                          "S48_R1_OE",
                                          "S48_R3_OE",
                                          "S18_G2_WT",
                                          "S18_G3_WT",
                                          "S18_G8_WT",
                                          "S18_G1_OE",
                                          "S18_G6_OE",
                                          "S18_G7_OE",
                                          "S18_R1_OE"))
# nFeature Plot
p <- ggplot(data = df, aes(x=sample, y=nFeature_RNA, fill = timepoint)) +
  geom_violin() +
  labs(x="Sample", y="nFeature") +
  scale_fill_manual(values = c("#F8BC24", "#F58800")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 10))
p
pdf(file = "results/misc-figures/nFeature.pdf",
    height = 3.5,
    width = 7,
    useDingbats = F)
print(p)
dev.off()

# nCount Plot
p <- ggplot(data = df, aes(x=sample, y=nCount_RNA, fill = timepoint)) +
  geom_violin() +
  labs(x="Sample", y="nCount") +
  scale_fill_manual(values = c("#F8BC24", "#F58800")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 10))
p
pdf(file = "results/misc-figures/nCount.pdf",
    height = 3.5,
    width = 7,
    useDingbats = F)
print(p)
dev.off()

# percent.mt Plot
p <- ggplot(data = df, aes(x=sample, y=percent.mt, fill = timepoint)) +
  geom_violin() +
  labs(x="Sample", y="% mitochondrial reads") +
  scale_fill_manual(values = c("#F8BC24", "#F58800")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 10))
p
pdf(file = "results/misc-figures/percent_mt.pdf",
    height = 3.5,
    width = 7,
    useDingbats = F)
print(p)
dev.off()

# DimPlot with correct sizing
p <- DimPlot(obj,
             reduction = "umap",
             pt.size = .3,
             raster = F,
             label = F,
             cols = colors,
             group.by = "basic_annotation") +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  labs(color = "Identity") +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        legend.justification = "top",
        legend.key.size = unit(3, "point"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 32),
        plot.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4.25),
                              nrow = 37))
q <- LabelClusters(plot = p,
                   id = "basic_annotation",
                   repel = T,
                   force = 0.25,
                   box = T,
                   fill = alpha("white", 0.45),
                   size = 6,
                   label.r = unit(0.25, "lines"),
                   label.size = NA)
pdf(file = "results/misc-figures/DimPlot_basic_annotation.pdf",
    useDingbats = F,
    width = 8)
print(q)
dev.off()


# DimPlot with correct sizing, split by time point
obj@meta.data$timepoint <- factor(obj@meta.data$timepoint,
                                  levels = c("day_3", "day_7"),
                                  labels = c("3 d", "7 d"))
p <- DimPlot(obj,
             reduction = "umap",
             pt.size = .3,
             raster = F,
             label = T,
             label.size = 5,
             cols = colors,
             group.by = "basic_annotation",
             split.by = "timepoint") +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  labs(color = "Identity") +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 24),
        legend.justification = "top",
        legend.key.size = unit(3, "point"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 32),
        plot.title = element_blank(),
        strip.text.x = element_text(size = 22)) +
  guides(color = guide_legend(override.aes = list(size = 3.25),
                              nrow = 37))
p
pdf(file = "results/misc-figures/DimPlot_basic_annotation_timepoint_split.pdf",
    useDingbats = F,
    width = 16)
print(p)
dev.off()

# DimPlot with correct sizing, split by genotype
obj@meta.data$genotype <- factor(obj@meta.data$genotype,
                                  levels = c("WT", "OE"),
                                  labels = c("WT", "Hmmr-OE"))
p <- DimPlot(obj,
             reduction = "umap",
             pt.size = .3,
             raster = F,
             label = T,
             label.size = 5,
             cols = colors,
             group.by = "basic_annotation",
             split.by = "genotype") +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  labs(color = "Identity") +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 24),
        legend.justification = "top",
        legend.key.size = unit(3, "point"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 32),
        plot.title = element_blank(),
        strip.text.x = element_text(size = 22)) +
  guides(color = guide_legend(override.aes = list(size = 3.25),
                              nrow = 37))
p
pdf(file = "results/misc-figures/DimPlot_basic_annotation_genotype_split.pdf",
    useDingbats = F,
    width = 16)
print(p)
dev.off()


# OE Hmmr feature plot
f <- FeaturePlot(subset(obj, genotype == "OE"),
            features = "Hmmr",
            reduction = "umap",
            pt.size = .3,
            raster = F) & scale_color_viridis_c(limits = c(0, 5)) 
f <- f + labs(x = "UMAP-1", y = "UMAP-2", title = expression(paste(italic("Hmmr"), "-OE")), col = "Hmmr") +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 20, face = "italic"),
        legend.key.size = unit(.75, "cm"),
        legend.position = "right",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 32),
        plot.title = element_text(size = 50, hjust = 0))
f
pdf(file = "results/misc-figures/featureplot_hmmr-oe.pdf",
    useDingbats = F,
    width = 8)
print(f)
dev.off()

# WT Hmmr feature plot
f <- FeaturePlot(subset(obj, genotype == "WT"),
                 features = "Hmmr",
                 reduction = "umap",
                 pt.size = .3,
                 raster = F) & scale_color_viridis_c(limits = c(0, 5)) 
f <- f + labs(x = "UMAP-1", y = "UMAP-2", title = expression(paste("WT")), col = "Hmmr") +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 20, face = "italic"),
        legend.key.size = unit(.75, "cm"),
        legend.position = "right",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 32),
        plot.title = element_text(size = 50, hjust = 0))
f
pdf(file = "results/misc-figures/featureplot_hmmr-wt.pdf",
    useDingbats = F,
    width = 8)
print(f)
dev.off()

# Hmmr split vln plot
v <- VlnPlot(obj, features = "Hmmr",
        split.by = "genotype",
        cols = c("#999999", "#008000"),
        pt.size = 0,
        raster = F,
        split.plot = T) +
  theme(plot.title = element_text(size = 50, face = "italic", hjust = 0),
        legend.text = element_text(size = 24),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 24))
pdf(file = "results/misc-figures/vlnplot_hmmr-split.pdf",
    useDingbats = F,
    width = 16)
print(v)
dev.off()

# Feature plots of canonical markers, sce style
sce <- as.SingleCellExperiment(obj, assay = "RNA")
markers <- read.csv(file = "data/canonical_markers.csv")
markers <- markers$All
for (i in markers) {
  p <- plotReducedDim(sce,
                      dimred = "UMAP",
                      colour_by = i,
                      by_exprs_values = "logcounts") +
    scale_fill_viridis_b() +
    labs(x = "UMAP-1", y = "UMAP-2", title = i) +
    theme(axis.text = element_blank(),
          legend.key.size = unit(.75, "cm"),
          legend.position = c(0.01,0.95),
          legend.direction = "horizontal",
          legend.title = element_blank(),
          plot.title = element_text(size = 50, face = "italic"),
          legend.spacing.y = unit(0.15, "cm"),
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 32))
  pdf(file = paste0("results/misc-figures/featureplot_", i, ".pdf"),
      useDingbats = F)
  print(p)
  dev.off()
}
