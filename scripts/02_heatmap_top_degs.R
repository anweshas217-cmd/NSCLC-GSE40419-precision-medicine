## =============================================================
## Stage 4.2 — Heatmap of Top DEGs Across All 10 Samples
## NSCLC GSE40419
## =============================================================

library(DESeq2)
library(pheatmap)
library(dplyr)

# ---- 1. Load raw counts ----
counts <- read.csv("data/featurecounts_matrix.csv", row.names = 1, check.names = FALSE)
# rownames(counts) = Geneid (Entrez IDs); columns = LC_S1..LC_S5, LC_C2,C3,C5,C7,C9

coldata <- data.frame(
    sample = colnames(counts),
    condition = ifelse(grepl("^LC_S", colnames(counts)), "Tumor", "Normal")
)
rownames(coldata) <- coldata$sample

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
dds <- DESeq(dds)

# Variance-stabilizing transform (better than raw log2 for heatmaps)
vsd <- vst(dds, blind = FALSE)

# ---- 2. Load DEG results and pick top genes ----
# IMPORTANT: match on ENTREZID, since that's what rownames(counts) actually are.
# SYMBOL is only used afterward, to make the heatmap labels readable.
res <- read.csv("data/merged_annotated_DESeq2.csv", stringsAsFactors = FALSE) %>%
    filter(!is.na(padj)) %>%
    arrange(padj)

top_n <- 50
top_res <- head(res, top_n)
top_entrez <- as.character(top_res$ENTREZID)

# Match to expression matrix by Entrez ID
mat <- assay(vsd)[rownames(assay(vsd)) %in% top_entrez, ]

# Relabel rows with gene SYMBOL instead of Entrez ID, for a readable heatmap
symbol_lookup <- setNames(top_res$SYMBOL, as.character(top_res$ENTREZID))
rownames(mat) <- symbol_lookup[rownames(mat)]

# ---- 3. Scale by row (z-score) for visualization ----
mat_scaled <- t(scale(t(mat)))

# ---- 4. Annotation for columns ----
annotation_col <- data.frame(Condition = coldata$condition)
rownames(annotation_col) <- coldata$sample

# ---- 5. Plot heatmap ----
pheatmap(mat_scaled,
    annotation_col = annotation_col,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 6,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    main = paste("Top", top_n, "DEGs — NSCLC Tumor vs Normal (GSE40419)"),
    filename = "figures/heatmap_top50_DEGs.pdf",
    width = 8, height = 12
)

cat("Heatmap saved to figures/heatmap_top50_DEGs.pdf\n")
cat("Genes actually plotted:", nrow(mat), "out of top", top_n, "requested\n")
