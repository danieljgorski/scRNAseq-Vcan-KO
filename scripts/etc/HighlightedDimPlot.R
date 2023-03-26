HighlightedDimPlot <- function (seurat_object, identity) {
  dd <- DimPlot(seurat_object, 
                reduction = "umap", 
                cells.highlight = WhichCells(seurat_object,
                                             idents = as.character(identity)),
                sizes.highlight = 0.3,
                cols.highlight = "#F06745",
                raster = F,
                pt.size = .3) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank()) +
    NoLegend()
  dd <- LabelClusters(plot = dd, 
                      id = "ident",
                      clusters = identity,
                      repel = F, 
                      box = T,
                      fill = alpha("white", 0.45),
                      size = 4,
                      label.r = unit(0.25, "lines"),
                      label.size = NA)
  print(dd)
}


