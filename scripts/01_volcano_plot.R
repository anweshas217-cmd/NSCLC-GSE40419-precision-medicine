## =============================================================
## Stage 4.1 — Volcano Plot
## NSCLC GSE40419: Tumor vs Normal DEGs
## =============================================================

# install.packages("BiocManager")
# BiocManager::install(c("EnhancedVolcano"))
library(EnhancedVolcano)
library(dplyr)

# ---- 1. Load data ----
# Actual columns confirmed from your file: SYMBOL, log2FoldChange, padj, etc.
res <- read.csv("data/merged_annotated_DESeq2.csv", stringsAsFactors = FALSE)

# Remove rows with NA padj (DESeq2 sets these for low-count/outlier genes)
res <- res %>% filter(!is.na(padj))

# ---- 2. Define your final novel targets to label ----
# (No approved drug + confirmed oncogene/tumor suppressor via OncoKB)
top_targets <- c("CDC25A", "E2F2", "E2F5", "PAX5")

# ---- 3. Volcano plot ----
p <- EnhancedVolcano(res,
    lab = res$SYMBOL,
    x = "log2FoldChange",
    y = "padj",
    selectLab = top_targets,
    xlim = c(min(res$log2FoldChange, na.rm = TRUE) - 1,
             max(res$log2FoldChange, na.rm = TRUE) + 1),
    title = "NSCLC Tumor vs Normal — Differential Expression",
    subtitle = "GSE40419 (LC_S1-5 vs LC_C2/3/5/7/9)",
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 2.0,
    labSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    boxedLabels = TRUE,
    colAlpha = 0.6,
    legendPosition = "right",
    col = c("grey70", "steelblue", "goldenrod", "firebrick")
)

# ---- 4. Save ----
ggplot2::ggsave("volcano_plot_NSCLC.pdf", plot = p, width = 10, height = 8)
ggplot2::ggsave("volcano_plot_NSCLC.png", plot = p, width = 10, height = 8, dpi = 300)

print(p)
