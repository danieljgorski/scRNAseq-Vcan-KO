# Miscellaneous figures for the manuscript

# Load libraries
library(Seurat)
library(stringr)
library(readr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggrepel)
source("scripts/etc/colors.R")
source("scripts/etc/geno_colors.R")
source("scripts/etc/dimplotlevels.R")

# Set up output dirs
output_dirs <- c("results",
                 "results/misc-figures")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Load object, slim down and save
load("results/objects/obj_annotated.Rdata")
obj <- DietSeurat(obj, counts = F, dimreducs = "umap")
save(obj, file = "results/objects/obj_annotated_diet.Rdata")

# Metadata to df for ggplot2 QC plots
df <- obj@meta.data
df$genotype <- factor(df$genotype, levels = c("Vcan-WT", "Vcan-KO"),
                       labels = c("WT", "KO"))
df$sample <- factor(df$sample, levels = c("Blau5Control",
                                          "Rot2Control",
                                          "Blau4KO",
                                          "Rot4KO",
                                          "Rot5KO"))
# nFeature Plot
p <- ggplot(data = df, aes(x=sample, y=nFeature_RNA, fill = genotype)) +
  geom_violin() +
  labs(x="Sample", y="nFeature") +
  scale_fill_manual(values = c("#999999", "#55A0FB")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(colour = "black", size = 20),
        axis.title = element_text(size = 24),
        legend.title = element_blank(),
        legend.text = element_text(size = 20))
p
pdf(file = "results/misc-figures/nFeature.pdf",
    useDingbats = F)
print(p)
dev.off()

# nCount Plot
p <- ggplot(data = df, aes(x=sample, y=nCount_RNA, fill = genotype)) +
  geom_violin() +
  labs(x="Sample", y="nCount") +
  scale_fill_manual(values = c("#999999", "#55A0FB")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(colour = "black", size = 20),
        axis.title = element_text(size = 24),
        legend.title = element_blank(),
        legend.text = element_text(size = 20))
p
pdf(file = "results/misc-figures/nCount.pdf",
    useDingbats = F)
print(p)
dev.off()

# percent.mt Plot
p <- ggplot(data = df, aes(x=sample, y=percent.mt, fill = genotype)) +
  geom_violin() +
  labs(x="Sample", y="% mitochondrial reads") +
  scale_fill_manual(values = c("#999999", "#55A0FB")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(colour = "black", size = 20),
        axis.title = element_text(size = 24),
        legend.title = element_blank(),
        legend.text = element_text(size = 20))
p
pdf(file = "results/misc-figures/percent_mt.pdf",
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
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.justification = "top",
        legend.key.size = unit(3, "point"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 18),
        plot.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4.25),
                              nrow = 37))
q <- LabelClusters(plot = p,
                   id = "basic_annotation",
                   repel = T,
                   force = 0.25,
                   box = T,
                   fill = alpha("white", 0.45),
                   size = 5,
                   label.r = unit(0.25, "lines"),
                   label.size = NA)
pdf(file = "results/misc-figures/DimPlot_basic_annotation.pdf",
    useDingbats = F,
    width = 8)
print(q)
dev.off()

# Feature plots of canonical markers
gois <- c("Col1a1", "Wt1", "Dmkn", "Tcf21", "Pdgfra", "Gsn",
          "Vcan", "Postn", "Cthrc1", "Mki67")
for (i in gois) {
  p <- FeaturePlot(obj, features = i, cols = c("lightgrey", "#fe3108")) +
    labs(x = "UMAP-1", y = "UMAP-2", title = i) +
    theme(axis.text = element_blank(),
          legend.key.size = unit(.75, "cm"),
          legend.position = c(0.01,0.95),
          legend.direction = "horizontal",
          legend.title = element_blank(),
          plot.title = element_text(size = 40, face = "italic", hjust = 0),
          legend.spacing.y = unit(0.15, "cm"),
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 24))
  pdf(file = paste0("results/misc-figures/featureplot_", i, ".pdf"),
      useDingbats = F)
  print(p)
  dev.off()
}

# VlnPlot of Vcan expression
v <- VlnPlot(obj, features = "Vcan", cols = colors, pt.size = 0) +
  labs(y = expression(paste(italic("Vcan"), " expression level"))) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_blank())
pdf(file = "results/misc-figures/VlnPlot_Vcan.pdf",
    useDingbats = F,
    height = 3.5)
print(v)
dev.off()
