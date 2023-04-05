VolcanoPlot <- function(df,
                        identity,
                        title = "DE Genes",
                        top_n_stat = 5,
                        top_n_fc = 5,
                        colors = c("black", "lightgrey")) {
  # Filter by identity
  df <- df %>% filter(cluster == identity)
  
  # Calculate -log10(Padj)
  df <- df %>% mutate(neglog10p = -(log10(df$p_val_adj)))
  
  # Indicate if genes are statically significant and have abs log2 FC > 0.25
  df <- df %>% mutate(significance = case_when(p_val_adj < 0.01 &
                                                 abs(avg_log2FC) > 0.25
                                               ~ "DE",
                                               .default = "Not-DE"))
  
  # Select top n DE genes most statistically significant
  df_top_n_stat <- df %>%
    filter(significance == "DE") %>%
    group_by(regulation) %>%
    slice_min(p_val_adj, n = top_n_stat)
  
  # Select top n DE genes with largest FC
  df_top_n_fc <- df %>%
    filter(significance == "DE") %>%
    group_by(regulation) %>%
    slice_max(abs(avg_log2FC), n = top_n_fc)
  
  # Plot
  v <- ggplot(df, aes(x = avg_log2FC, y = neglog10p)) +
    geom_point(aes(color = significance)) +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed",
               color = "black", linewidth = 0.2) +
    geom_vline(xintercept = 0.25, linetype = "dashed",
               color = "black", linewidth = 0.2) +
    geom_vline(xintercept = -0.25, linetype = "dashed",
               color = "black", linewidth = 0.2) +
    geom_label_repel(data = df_top_n_stat, aes(label = gene)) +
    geom_label_repel(data = df_top_n_fc, aes(label = gene)) +
    labs(x = expression("avg log"[2] * "(fold change)"),
         y = expression("-log"[10] * "(P"[adj] * ")"),
         title = title) +
    theme_bw() +
    theme(plot.title = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  print(v)
}

