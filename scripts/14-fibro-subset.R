# Subsetting fibroblast clusters for re-clustering and targeted analysis

# Load libraries
library(Seurat)

# Load full object
load("results/objects/obj_annotated.Rdata")

# Subset fibroblast clusters, without cycling fibroblasts
fibro <- subset(x = obj, idents = c("Fibro-Myo-1",
                                    "Fibro-Myo-2",
                                    "Fibro-Myo-3",
                                    "Fibro-Act-1",
                                    "Fibro-Act-2",
                                    "Fibro-IFN",
                                    "Fibro-Rest"))

# Save fibroblast object
save(fibro, file = "results/objects/fibro.Rdata")

