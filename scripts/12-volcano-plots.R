# Volcano plots of differential gene expression analysis

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
source("scripts/etc/VolcanoPlot.R")

# Load dge data
dge_no_threshold <-
  read.csv(file = "results/differential-gene-expression/dge_no_threshold.csv")

# Set up output dirs
output_dirs <- c("results",
                 "results/volcano-plots")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Loop through clusters and output volcano plots
for (i in unique(dge_no_threshold$cluster)) {
  pdf(file = paste0("results/volcano-plots/VolcanoPlot_", i, ".pdf"),
    height = 6,
    width = 8,
    useDingbats = F)
  VolcanoPlot(df = dge_no_threshold,
              identity = i,
              top_n_stat = 10,
              top_n_fc = 0,
              title = paste0(i, " differential gene expression: ", "Vcan-KO/WT"))
  dev.off()
}
