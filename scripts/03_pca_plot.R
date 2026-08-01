## =============================================================
## Stage 4.3 — PCA Plot: Tumor vs Normal Clustering
## NSCLC GSE40419
## =============================================================

library(DESeq2)
library(ggplot2)
library(ggrepel)

# ---- 1. Reuse the vsd object from script 02, or rebuild ----
# vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 4, alpha = 0.85) +
    geom_text_repel(size = 3, show.legend = FALSE) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    scale_color_manual(values = c("Tumor" = "firebrick3", "Normal" = "steelblue")) +
    ggtitle("PCA — NSCLC Tumor vs Adjacent Normal (GSE40419)") +
    theme_bw(base_size = 14) +
    theme(legend.title = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))

ggsave("figures/PCA_plot_NSCLC.pdf", plot = p, width = 8, height = 6)
ggsave("figures/PCA_plot_NSCLC.png", plot = p, width = 8, height = 6, dpi = 300)

print(p)
