## =============================================================
## Stage 4.5 — Novel Target Discovery Funnel Visualization
## 2387 DEGs -> 1189 upregulated -> 63 no approved drug -> 4 oncogene/TSG confirmed
## =============================================================

library(ggplot2)
library(dplyr)

# ---- 1. Funnel counts (your actual filtering results) ----
funnel_df <- data.frame(
    stage = factor(
        c("All DEGs", "Upregulated", "No approved\ndrug (DGIdb)", "Oncogene/TSG\n(OncoKB) - Final"),
        levels = c("All DEGs", "Upregulated", "No approved\ndrug (DGIdb)", "Oncogene/TSG\n(OncoKB) - Final")
    ),
    n_genes = c(2387, 1189, 63, 4)
)

# ---- 2. Bar chart funnel ----
p_funnel <- ggplot(funnel_df, aes(x = stage, y = n_genes, fill = stage)) +
    geom_col(width = 0.65, show.legend = FALSE) +
    geom_text(aes(label = n_genes), vjust = -0.5, size = 5, fontface = "bold") +
    scale_fill_manual(values = colorRampPalette(c("steelblue", "firebrick3"))(nrow(funnel_df))) +
    labs(
        title = "Novel Target Discovery Funnel — NSCLC (GSE40419)",
        subtitle = "Upregulated DEGs with no approved drug, confirmed as oncogene/tumor suppressor",
        x = NULL, y = "Number of Genes"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 11, color = "gray30"),
        axis.text.x = element_text(size = 11)
    )

ggsave("figures/novel_target_funnel.pdf", plot = p_funnel, width = 9, height = 6)
ggsave("figures/novel_target_funnel.png", plot = p_funnel, width = 9, height = 6, dpi = 300)

# ---- 3. Final target summary table/chart (4 novel targets) ----
targets_df <- data.frame(
    gene = c("CDC25A", "E2F2", "E2F5", "PAX5"),
    role = c("Oncogene", "Tumor Suppressor", "Tumor Suppressor", "Oncogene")
)

p_targets <- ggplot(targets_df, aes(x = reorder(gene, gene), fill = role)) +
    geom_bar(width = 0.5) +
    coord_flip() +
    scale_fill_manual(values = c("Oncogene" = "firebrick3", "Tumor Suppressor" = "steelblue")) +
    labs(
        title = "Final Novel Targets — No Approved Drug, Confirmed Oncogene/TSG",
        x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.title = element_blank())

ggsave("figures/final_novel_targets.pdf", plot = p_targets, width = 8, height = 5)
ggsave("figures/final_novel_targets.png", plot = p_targets, width = 8, height = 5, dpi = 300)

cat("Saved: figures/novel_target_funnel.pdf/png\n")
cat("Saved: figures/final_novel_targets.pdf/png\n")
