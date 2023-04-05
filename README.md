# scRNAseq-Vcan-KO
This repository contains single-cell RNA sequencing analysis of CD31-CD45- cardiac interstitial cells from WT and Versican knockout (*Vcan*-KO) mice, 7 days (n = 2, 3) after acute myocardial infarction. Data were processed with the Seurat toolkit, using SCTransform normalization and reference-based integration with reciprocal PCA. WT samples were used as reference, *Vcan*-KO samples as query.

## Sequencing data
Sequencing data is available at ArrayExpress accession E-MTAB-12880, a public link will be available upon publication or request.

## Analysis
To recreate the full analysis you can run each script in order. If you would like to process the data with your own custom workflow, a final list of cells (barcodes, annotation, genotype and embeddings) after preprocessing, doublet removal and low-quality cluster removal can be found in `data/basic_annotation.csv`.

## Libraries
The following scRNA-seq-specific libraries were used with **R version 4.1.2**:

* [Seurat](https://satijalab.org/seurat/index.html) v4.3.0
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
* [monocle3](https://cole-trapnell-lab.github.io/monocle3/) v1.0.0

## Examples
<p align="center">
  <img src="/examples/DimPlot_basic_annotation.png" width="1000">
</p>

<p align="center">
  <img src="/examples/Heatmap.png" width="1000">
</p>
