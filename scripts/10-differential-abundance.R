# Differential abundance analysis

# Load libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(edgeR)
source("scripts/etc/colors.R")
source("scripts/etc/dimplotlevels.R")

# Load object
load("results/objects/obj_annotated.Rdata")

# Set up output dirs
output_dirs <- c("results",
                 "results/differential-abundance")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Proportion calculation----
prop_by_sample <- (prop.table(table(obj@meta.data$basic_annotation,
                                          obj@meta.data$sample), margin = 2))
prop_by_sample <- as.data.frame(prop_by_sample)
colnames(prop_by_sample) <- c("cluster", "sample", "fraction")
prop_by_sample$genotype[str_detect(prop_by_sample$sample, "OE")] <- "OE"
prop_by_sample$genotype[str_detect(prop_by_sample$sample, "WT")] <- "WT"
write.csv(prop_by_sample,
          file = "results/differential-abundance/prop_by_sample.csv",
          row.names = F)

# miloR----
# Convert to SingleCellExperiment
obj_sce <- as.SingleCellExperiment(obj, assay = "RNA")

# Convert to milo object
obj_milo <- Milo(obj_sce)
head(colData(obj_milo))

# Build kNN graph, calculate neighborhood counts
obj_milo <- buildGraph(obj_milo, k = 40, d = 25, reduced.dim = "PCA")
obj_milo <- makeNhoods(obj_milo, prop = 0.1, k = 40, d = 25, refined = TRUE)

# Plot neighborhood size, peak should be between 50-100
p <- plotNhoodSizeHist(obj_milo)
pdf(file = "results/differential-abundance/milo_Nhood_size.pdf",
    useDingbats = F)
print(p)
dev.off()

# Count cells in neighborhoods
obj_milo <- countCells(obj_milo,
                       meta.data = data.frame(colData(obj_milo)),
                       sample = "sample")

# Compute neighborhood connectivity
obj_milo <- calcNhoodDistance(obj_milo, d = 25)

# Define experimental design
exp_design <- data.frame(colData(obj_milo))[, c("sample", "genotype")]
exp_design <- distinct(exp_design)
rownames(exp_design) <- exp_design$sample
exp_design

# Testing
milo_res <- testNhoods(obj_milo, design = ~ genotype, design.df = exp_design)
milo_res %>%  arrange(SpatialFDR) %>%  head()

# Inspecting and plotting results
p1 <- ggplot(milo_res, aes(PValue)) + geom_histogram(bins = 50)
p2 <- ggplot(milo_res, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
pdf(file = "results/differential-abundance/milo_pvalue_spatialfdr.pdf",
    useDingbats = F)
print(p1 + p2)
dev.off()

# Build Nhood graph
obj_milo <- buildNhoodGraph(obj_milo)

# Plot DA results next to UMAP
p1 <- plotReducedDim(obj_milo,
                     dimred = "UMAP",
                     colour_by = "basic_annotation",
                     text_by = "basic_annotation",
                     text_size = 3) +
  scale_color_manual(values = colors) +
  NoLegend()
p2 <- plotNhoodGraphDA(obj_milo, milo_res, alpha = 0.05)
                    # alpha here is statistical sig threshold, not transparency
pdf(file = "results/differential-abundance/milo_UMAP_NhoodGraph.pdf",
    height = 6.5,
    width = 12,
    useDingbats = F
    )
print(p1 + p2 + plot_layout(guides = "collect"))
dev.off()

# Annotate the Nhoods based on basic_annotation
milo_res <- annotateNhoods(obj_milo, milo_res, coldata_col = "basic_annotation")
unique(milo_res$basic_annotation)
ggplot(milo_res, aes(basic_annotation_fraction)) + geom_histogram(bins = 50)
milo_res$basic_annotation <- ifelse(milo_res$basic_annotation_fraction < 0.6,
                                    "Mixed",
                                    milo_res$basic_annotation)
milo_res$basic_annotation <- factor(milo_res$basic_annotation,
                                    levels = dimplotlevels)

# plot DAbeeswarm...ignored because none are significant
plotDAbeeswarm(milo_res,
               group.by = "basic_annotation",
               alpha = 0.05) +
  ggtitle("Cluster assignment of DA neighborhoods") +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# Save DA results and objects
write.csv(milo_res,
          file = "results/differential-abundance/milo_res.csv",
          row.names = F)

# Save milo object
save(obj_milo, file = "results/objects/obj_milo.Rdata")


# OCSA-DA Analysis with edgeR----
#http://bioconductor.org/books/3.14/OSCA.multisample/differential-abundance.html

# Abundances
abundances <- table(obj_sce$basic_annotation, obj_sce$sample)
abundances <- unclass(abundances)
head(abundances)

# Attaching some column metadata and making DGEList object
extra.info <- colData(obj_sce)[match(colnames(abundances), obj_sce$sample), ]
y.ab <- DGEList(abundances, samples = extra.info)
y.ab

# Filter low abundance
keep <- filterByExpr(y.ab, group = y.ab$samples$genotype)
y.ab <- y.ab[keep, ]
summary(keep)

# Design matrix
design <- model.matrix(~genotype, y.ab$samples)

# Estimate dispersion
y.ab <- estimateDisp(y.ab, design, trend = "none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex = 1)

# QL fit
fit.ab <- glmQLFit(y.ab, design, robust = TRUE, abundance.trend = FALSE)
summary(fit.ab$var.prior)
plotQLDisp(fit.ab, cex = 1)

# Test
ocsa_da_res <- glmQLFTest(fit.ab, coef = ncol(design))
summary(decideTests(ocsa_da_res))
topTags(ocsa_da_res, n = 34)
ocsa_da_res$table

# Workflow with normalization (assuming most labels do not change in abundance)
y.ab2 <- calcNormFactors(y.ab)
y.ab2$samples$norm.factors
y.ab2 <- estimateDisp(y.ab2, design, trend = "none")
fit.ab2 <- glmQLFit(y.ab2, design, robust = TRUE, abundance.trend = FALSE)
ocsa_da_res_norm <- glmQLFTest(fit.ab2, coef = ncol(design))
topTags(ocsa_da_res_norm)

# Save results
write.csv(topTags(ocsa_da_res, n = 34),
          file = "results/differential-abundance/ocsa_da_res.csv")
write.csv(topTags(ocsa_da_res_norm, n = 34),
          file = "results/differential-abundance/ocsa_da_res_norm.csv")

# Save sce object
save(obj_sce, file = "results/objects/obj_sce.Rdata")
