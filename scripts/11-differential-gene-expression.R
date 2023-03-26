# Differential gene expression analysis

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
source("scripts/etc/dimplotlevels.R")

# Load object
load("results/objects/obj_annotated.Rdata")

# Set up output dirs
output_dirs <- c("results",
                 "results/differential-gene-expression")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Seurat based Wilcox method, no threshold, loop through each cluster
genes <- list()
for (i in levels(Idents(obj))) {
  results <- FindMarkers(obj,
                         subset.ident = i,
                         group.by = "genotype",
                         ident.1 = "Vcan-KO",
                         base = 2,
                         logfc.threshold = 0,
                         densify = T)
  results$cluster <- i
  results$gene <- row.names(results)
  row.names(results) <- NULL
  results$regulation <- ifelse(results$avg_log2FC > 0, "Up", "Down")
  genes[[i]] <- results
}
dge_no_threshold <- do.call(rbind, genes)
rownames(dge_no_threshold) <- NULL
write.csv(dge_no_threshold,
          file = "results/differential-gene-expression/dge_no_threshold.csv",
          row.names = F)

# Filter out non-significant genes
dge <- dge_no_threshold[dge_no_threshold$p_val_adj < 0.01 &
                          abs(dge_no_threshold$avg_log2FC) > 0.25, ]

# Save significant deg
write.csv(dge,
          file = "results/differential-gene-expression/dge.csv",
          row.names = F)

# Count and plot deg per cluster
dge$regulation <- factor(dge$regulation, levels = c("Up", "Down"))
dge$cluster <- factor(dge$cluster, levels = rev(dimplotlevels))
p <- dge %>%
  group_by(cluster, regulation) %>%
  ggplot(aes(x = cluster)) +
  geom_bar(aes(fill = regulation)) +
  coord_flip() +
  scale_fill_manual(name = "Regulation in Vcan-KO",
    values = c("#D75438", "#4878CD")) +
  theme_bw() +
  ggtitle("Number of differentially expressed genes")
pdf(file = "results/differential-gene-expression/dge_count.pdf",
    useDingbats = F)
print(p)
dev.off()
