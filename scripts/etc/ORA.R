# Performs a gene ontology biological processes over-representation analysis
# using the ClusterProfiler package. Calculates and plots a barchart of
# fold enrichment instead of GeneRatio.

ORA <- function(genes_of_interest, background_genes, plot_title) {
  GO_BP_over_rep <- enrichGO(
    gene = genes_of_interest,
    universe = background_genes,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    keyType = "SYMBOL",
    pAdjustMethod = "bonferroni",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = F
  )
  df <- GO_BP_over_rep@result
  df$mapped <- sub(".*/", "", df$GeneRatio)
  df$mapped <- as.numeric(df$mapped)
  df$GO_in_background <- sub("/.*", "", df$BgRatio)
  df$GO_in_background <- as.numeric(df$GO_in_background)
  df$background_mapped <- sub(".*/", "", df$BgRatio)
  df$background_mapped <- as.numeric(df$background_mapped)
  df$Fold_enrichment <- df$Count / (df$mapped * (df$GO_in_background / df$background_mapped))
  GO_BP_over_rep@result <- df
  p <- barplot(GO_BP_over_rep,
    showCategory = 10,
    x = "Fold_enrichment",
    color = "p.adjust",
    order = T,
    title = plot_title,
    font.size = 10
  ) +
    theme(plot.title = element_text(size = 10)) +
    xlab("Fold enrichment")
  print(p)
}
