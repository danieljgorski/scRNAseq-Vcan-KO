# Doublet classification with DoubletFinder

# Load libraries
library(Seurat)
library(DoubletFinder) #v2.0.3
library(dplyr)
library(ggplot2)

# Load object
load("results/objects/obj_preprocessed.Rdata")

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

# Loop through each sample and classify doublets
obj_list <- list()
for (i in unique(obj@meta.data$sample)) {
  # Cluster and sweep
  obj_sub <- subset(x = obj, subset = sample == i)
  obj_sub <- SCTransform(obj_sub)
  obj_sub <- RunPCA(obj_sub)
  obj_sub <- RunUMAP(obj_sub, reduction = "pca", dims = 1:18, verbose = T)
  obj_sub <- FindNeighbors(obj_sub, dims = 1:18, verbose = T)
  obj_sub <- FindClusters(obj_sub, resolution = 0.8, verbose = T)
  sweep_res_list_obj <- paramSweep_v3(obj_sub, PCs = 1:18, sct = T)
  sweep_stats_obj <- summarizeSweep(sweep_res_list_obj, GT = FALSE)
  bcmvn_obj <- find.pK(sweep_stats_obj)

  # Model homotypic doublets
  # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/54
  homotypic_prop <- modelHomotypic(as.character(obj_sub@active.ident))
  # assuming ~ 0.8% per 1000 cells recovered
  doublet_rate <- (ncol(obj_sub) / 1000) * 0.008
  nexp_poi <- round(doublet_rate * length(colnames(obj_sub)))
  nexp_poi_adj <- round(nexp_poi * (1 - homotypic_prop))

  # Choose pK
  # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/62
  pK = as.numeric(as.character(bcmvn_obj$pK))
  BCmetric = bcmvn_obj$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  par(mar = c(5, 4, 4, 8) + 1, cex.main = 1.2, font.main = 2)

  # BCmvn distributions
  pdf(file = paste0("results/doublet-removal/BCmvn_distributions_",
                    i,
                    ".pdf"),
      useDingbats = F)
  plot(x = pK,
       y = BCmetric,
       pch = 16,
       type = "b",
       col = "blue",
       lty = 1)
  abline(v = pK_choose, lwd = 2, col = "red", lty = 2)
  title(paste0("BCmvn_distributions_", i))
  text(pK_choose, max(BCmetric), as.character(pK_choose), pos = 4, col = "red")
  dev.off()

  # Classification
  obj_sub <- doubletFinder_v3(obj_sub,
                              PCs = 1:18,
                              pN = 0.25,
                              pK = pK_choose,
                              nExp = nexp_poi_adj,
                              reuse.pANN = FALSE,
                              sct = T)
  colnames(obj_sub@meta.data)[ncol(obj_sub@meta.data)] <- "doublet_classification"
  obj_sub@meta.data$doublet_classification <- factor(obj_sub@meta.data$doublet_classification,
                                                     levels = c("Singlet", "Doublet"))

  # DimPlot of classifications
  p <- DimPlot(obj_sub, group.by = "doublet_classification") + ggtitle(i)
  pdf(file = paste0("results/doublet-removal/DimPlot_classification_",
                    i,
                    ".pdf"),
      useDingbats = F)
  print(p)
  dev.off()

  # VlnPlot of classification
  v <- VlnPlot(obj_sub,
               features = "nFeature_RNA",
               group.by = "doublet_classification") +
    ylab("nFeature_RNA") +
    ggtitle(i) +
    NoLegend()
  pdf(file = paste0("results/doublet-removal/VlnPlot_nFeature_RNA_classification_",
                    i,
                    ".pdf"),
      useDingbats = F)
  print(v)
  dev.off()

  # Append classified objects to list
  obj_list[i] <- obj_sub

  # Status
  print(paste0(i, " done.---------------------------------------------------"))
}

# Merge classified objects
obj <- merge(x = obj_list[[1]], y = c(obj_list[2:length(obj_list)]))

# Save barcode classifications and summary
barcode <- rownames(obj@meta.data)
classification <- obj@meta.data$doublet_classification
sample <- obj@meta.data$sample
barcode_classification <- data.frame(barcode, classification, sample)
write.csv(barcode_classification,
          file = "results/doublet-removal/barcode_classification.csv",
          row.names = F)

# Calculate total doublet frequency
doublet_summary <- barcode_classification %>%
  group_by(classification) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))
write.csv(doublet_summary,
          file = "results/doublet-removal/doublet_summary.csv",
          row.names = F)

# Calculate doublet frequency by sample
doublet_summary_by_sample <- barcode_classification %>%
  group_by(sample, classification) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))
write.csv(doublet_summary_by_sample,
          file = "results/doublet-removal/doublet_summary_by_sample.csv",
          row.names = F)

# Save object
save(obj, file = "results/objects/obj_db_classified.Rdata")

# Clear memory
rm(list = ls())
gc()
