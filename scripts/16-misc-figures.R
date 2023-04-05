# Miscellaneous figures for the manuscript

# Load libraries
library(Seurat)
library(stringr)
library(readr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggrepel)
library(dplyr)
library(org.Mm.eg.db)
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


## Quality control plots----
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

## DimPlot with correct sizing----
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

## Feature plots of canonical markers----
gois <- c("Col1a1", "Wt1", "Dmkn", "Tcf21", "Pdgfra", "Gsn",
          "Vcan", "Postn", "Cthrc1", "Mki67", "Ctgf")
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

## Multi-panel volcano plots (certain order)----
# Retrieve GO terms of interest
res_to_ERstress <- AnnotationDbi::select(org.Mm.eg.db, keytype="GOALL",
                                         keys="GO:0034976", columns="SYMBOL")
res_to_ERstress <- unique(res_to_ERstress$SYMBOL)
translation <- AnnotationDbi::select(org.Mm.eg.db, keytype="GOALL",
                                     keys="GO:0006412", columns="SYMBOL")
translation <- unique(translation$SYMBOL)

# Load unfiltered DE
df <- read.csv(file = "results/differential-gene-expression/dge_no_threshold.csv")

# Calculate -log10(Padj)
df <- df %>% mutate(neglog10p = -(log10(df$p_val_adj)))

# Factor order
df$cluster <- factor(df$cluster, levels = c("Epi-Rest",
                                            "Fibro-Rest",
                                            "Epi-Act-1",
                                            "Fibro-Act-1",
                                            "Epi-Act-2",
                                            "Fibro-Act-2",
                                            "Epi-Act-3",
                                            "Fibro-Myo-1",
                                            "Epi-Cyc",
                                            "Fibro-Myo-2",
                                            "Fibro-IFN",
                                            "Fibro-Myo-3",
                                            "Fibro-Cyc"))


# Indicate if genes are statically significant and have abs log2 FC > 0.25
df <- df %>% mutate(significance = case_when(p_val_adj < 0.01 &
                                               abs(avg_log2FC) > 0.25
                                             ~ "DE",
                                             .default = "Not-DE"))

# Indicate if DE genes are in GO term translation or response to ER stress
df <- df %>% mutate(significance = case_when(significance == "DE" &
                                               gene %in% translation ~
                                               "DE: Translation - GO:0006412",
                                             TRUE ~ as.character(significance)))
df <- df %>% mutate(significance = case_when(significance == "DE" &
                                               gene %in% res_to_ERstress ~
                                               "DE: Response to ER Stress - GO:0034976",
                                             TRUE ~ as.character(significance)))
df$significance <- factor(df$significance, levels = c("Not-DE",
                                                      "DE",
                                                      "DE: Response to ER Stress - GO:0034976",
                                                      "DE: Translation - GO:0006412"))

# Plot
v <- ggplot(df, aes(x = avg_log2FC, y = neglog10p)) +
  facet_wrap(~ cluster, scales = "free", ncol = 2) +
  geom_point(aes(color = significance)) +
  scale_color_manual(values = c("#EBEBEB", "#B8B8B8",  "#f94144", "#858ae3")) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed",
             color = "#B8B8B8", linewidth = 0.2) +
  geom_vline(xintercept = 0.25, linetype = "dashed",
             color = "#B8B8B8", linewidth = 0.2) +
  geom_vline(xintercept = -0.25, linetype = "dashed",
             color = "#B8B8B8", linewidth = 0.2) +
  labs(x = expression("avg log"[2] * "(fold change)"),
       y = expression("-log"[10] * "(P"[adj] * ")"),
       title = expression(paste("Differential gene expression: ", italic("Vcan"), "-KO/WT"))) +
  theme_bw() +
  theme(plot.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = c(0.99, 0.015),
        legend.justification = c(1, 0),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
pdf(file = "results/misc-figures/VolcanoPlot_ERstress_Translation.pdf",
    useDingbats = F,
    height = 14,
    width = 6.5)
print(v)
dev.off()

## Split Violin of DE UPR genes----
DE_UPR <- c("Creld2", "Manf", "Hspa5", "Dnajc3")
for (i in DE_UPR) {
  v <- VlnPlot(obj,
               features = i,
               pt.size = 0,
               split.by = "genotype",
               split.plot = TRUE,
               cols = geno_colors) +
    labs(title = i,
         y = "Expression level") +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          plot.title = element_text(face = "italic", hjust = 0, size = 18))
  pdf(file = paste0("results/misc-figures/VlnPlot_", i, "_split.pdf"),
      useDingbats = F,
      height = 3,
      width = 5.5)
      print(v)
      dev.off()
}


## Violin of Vcan----
v <- VlnPlot(obj,
             features = "Vcan",
             pt.size = 0,
             cols = colors,
             sort = T) +
  labs(y = expression(paste(italic("Vcan"), " expression level"))) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 17),
        plot.title = element_blank())
pdf(file = "results/misc-figures/VlnPlot_Vcan.pdf",
    height = 4.75,
    width = 8.25,
    useDingbats = F)
print(v)
dev.off()
