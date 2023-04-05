# TI and differential progression analysis on fibroblasts in monocle3

# Load libraries
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(monocle3)
source("scripts/etc/geno_colors.R")
source("scripts/etc/colors.R")

# Set up output dirs
output_dirs <- c("results",
                 "results/fibro-differential-progression")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Load the fibroblast subset object
load("results/objects/obj_fibro_subset.Rdata")

# Select root cell IDs 
root <- CellSelector(plot = DimPlot(obj, reduction = "umap"))
save(root, file = "results/fibro-differential-progression/root.Rdata")

# Fibroblast TI analysis with monocle3

# Convert Seurat object to cds
cds <- as.cell_data_set(obj, reduction = "umap")
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(obj)
cds <- cluster_cells(cds, reduction_method = "UMAP")

# Select out main fibroblast population
cds <- choose_cells(cds)
cds <- cluster_cells(cds, reduction_method = "UMAP")

# Infer trajectory
cds <- learn_graph(cds, use_partition = F)
plot_cells(cds, label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE)

# Order cells along pseudotime, using the root cells as start
cds <- order_cells(cds, root_cells = root)
plot_cells(cds, color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE, 
           label_branch_points = F)
save(cds, file = "results/objects/obj_fibro_subset_cds.Rdata")

# Subset seurat obj match manual subset selected by cluster_cells()
obj[["cell_ids"]] <- rownames(obj@meta.data)
obj <- subset(obj, subset = cell_ids %in% colnames(cds))

# Save pseudotime values inside seurat object metadata
obj[["monocle3_pseudotime"]] <- pseudotime(cds, reduction_method = "UMAP")

# Statistical testing of pseudotime values
data <- data.frame(obj$monocle3_pseudotime)
colnames(data)[1] <- "pseudotime"
data$genotype <- obj$genotype
data$genotype <- factor(obj$genotype, levels = c("Vcan-WT", "Vcan-KO"),
                        labels = c("WT", "Vcan-KO"))

# Kolmogorov-Smirnov Test of pseudotime distributions
ks <- ks.test(data[data$genotype=='Vcan-KO',1], data[data$genotype=='WT',1])

# Plot the pseudotime distributions of each genotype
p <- ggplot(data, aes(x = pseudotime,
                     fill = genotype,
                     color = genotype)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values = geno_colors) +
  scale_color_manual(values = geno_colors) +
  labs(x = "Pseudotime", fill = "Genotype", y = "Density", color = "Genotype") +
  annotate("text", x = 4, y = .0725, label = "P < 2.2e-16", size = 9) +
  theme_classic() +
  theme(legend.position = "top",
        legend.key.height = unit(1, "cm"),
        legend.key.width =  unit(1, "cm"),
        legend.justification = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24))
pdf(file = "results/fibro-differential-progression/differential-progression.pdf",
    width = 8,
    useDingbats = F)
print(p)
dev.off()

# Plot pseudotime UMAP
p <- plot_cells(cds,
                color_cells_by = "pseudotime",
                label_cell_groups = FALSE,
                label_leaves = FALSE, 
                label_branch_points = F,
                trajectory_graph_color = "black",
                trajectory_graph_segment_size = 1.5,
                label_roots = F,
                cell_size = .1,
                cell_stroke = 1) +
  labs(title="Pseudotime", x = "UMAP-1", y = "UMAP-2") +
  scale_color_viridis_c() +
  theme(axis.text = element_blank(),
        legend.key.size = unit(.75, "cm"),
        legend.position = c(0.275,0.94),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        plot.title = element_text(size = 40, hjust = 0),
        legend.spacing.y = unit(0.15, "cm"),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 24))
p
pdf(file = "results/fibro-differential-progression/pseudotime.pdf",
    useDingbats = F,
    width = 3.5)
print(p)
dev.off()

# Plot basic_annotation UMAP
p <- plot_cells(cds,
                color_cells_by = "basic_annotation",
                label_leaves = FALSE, 
                label_branch_points = F,
                trajectory_graph_color = "black",
                trajectory_graph_segment_size = 1.5,
                label_cell_groups = F,
                label_roots = F,
                group_label_size = 12,
                cell_size = .05,
                cell_stroke = .5) +
  labs(title = "Fibroblast subset", x = "UMAP-1", y = "UMAP-2") +
  scale_color_manual(values = c("#F5705B",
                              "#F1AC78",
                              "#F3A29B",
                              "#32CEF0",
                              "#FF843F",
                              "#7F7F9D",
                              "#FF5D4E")) +
  theme(axis.text = element_blank(),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 40, hjust = 0),
        legend.title = element_blank(),
        legend.justification = "top",
        legend.text = element_text(size = 20))
p
pdf(file = "results/fibro-differential-progression/fibro-subset.pdf",
    useDingbats = F,
    width = 5.75)
print(p)
dev.off()

# GOI through pseudotime
Genesofinterest <- c("Gsn", "Vcan", "Postn", "Cthrc1", "Col1a1")
PScds <- cds[rowData(cds)$gene_short_name %in% Genesofinterest]
p <- plot_genes_in_pseudotime(PScds,
                              trend_formula = "~ splines::ns(pseudotime, df=3)", 
                              vertical_jitter = T,
                              horizontal_jitter = T,
                              cell_size = 0.75,
                              label_by_short_name = F,
                              ncol = 5,
                              nrow = 1,
                              panel_order = Genesofinterest) + 
  scale_color_viridis_c() +
  labs(color = "Pseudotime")+
  xlab("Pseudotime") +
  geom_line(aes(x = pseudotime, y = expectation), color = "black", size = 3) +
  theme(strip.text.x = element_text(size = 42, face = "italic", hjust = 0),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 26),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 26))
pdf(file = "results/fibro-differential-progression/GOI-ps.pdf", height = 7, width = 21)
print(p)
dev.off()

