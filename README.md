# scRNAseq-Vcan-KO
This repository contains single-cell RNA sequencing analysis of CD31^-^CD45^-^ cardiac interstitial cells from WT and Versican(*Vcan*) knockout (KO) mice, 7 days (n = 3, 4) after acute myocardial infarction. Data were processed with the Seurat toolkit, using SCTransform normalization and reference-based integration with reciprocal PCA. WT samples were used as reference, KO samples as query.

## Sequencing data
Sequencing data, including fastq files and count matrices will be available upon publication.

## Analysis
To recreate the full analysis you can follow the steps below. If you would like to process the data with your own custom workflow, a final list of cells (barcodes + metadata) after preprocessing, doublet removal and low-quality cluster removal can be found in:

* `data/basic_annotation.csv`

### Libraries
The following scRNA-seq-specific libraries were used:

* [Seurat](https://satijalab.org/seurat/index.html) v4.1.1
* [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) v2.0.3
* [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) v2.13.1
* [miloR](https://marionilab.github.io/miloR/) v.1.2.0
* [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) v1.16.0
* [scater](https://bioconductor.org/packages/release/bioc/html/scater.html) v1.22.0
* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) v3.36.0
* [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) v4.2.1
* [enrichplot](https://bioconductor.org/packages/release/bioc/html/enrichplot.html) v1.14.1
* [org.Mm.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html) v3.14.0
* [yulab.utils](https://cran.r-project.org/package=yulab.utils) v0.0.4

They can be installed with the following commands:
```R
# Seurat
remotes::install_version("Seurat", version = "4.0.1")

# DoubletFinder
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

# ComplexHeatmap
library(devtools)
install_github("jokergoo/ComplexHeatmap")

# miloR
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("miloR")

# SingleCellExperiment
BiocManager::install("SingleCellExperiment")

# scater
BiocManager::install("scater")

# edgeR
BiocManager::install("edgeR")

# clusterProfiler etc.
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("org.Mm.eg.db")
install.packages("yulab.utils")
```

Additionally, the following standard R libraries were used:
* dplyr
* ggplot2
* patchwork
* readr
* ggrepel
* stringr
* knitr
* kableExtra
* forcats

### Instructions
To reproduce the analysis, clone this repository and place the count matrices inside the `data` folder. It should then contain the following:

## Examples

## To-do
* Add abstract at submission
* Add full author list at submission
* Add mechanism schematic after publication
* Updated "About" short description after publication
