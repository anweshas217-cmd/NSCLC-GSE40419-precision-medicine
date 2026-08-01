## =============================================================
## Stage 4.4 — GO / KEGG Pathway Enrichment
## NSCLC GSE40419 — 2,387 significant DEGs
## =============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(ggplot2)

# ---- 1. Load DEGs ----
res <- read.csv("data/merged_annotated_DESeq2.csv", stringsAsFactors = FALSE) %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

cat("Total significant DEGs used for enrichment:", nrow(res), "\n")

# ---- 2. Use ENTREZID directly (already in your file — no SYMBOL->ENTREZ conversion needed) ----
entrez_ids <- as.character(res$ENTREZID)
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

# ---- 3. GO enrichment (Biological Process) ----
go_bp <- enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
)

# Dotplot of top GO terms
p_go <- dotplot(go_bp, showCategory = 20, title = "GO Biological Process — NSCLC DEGs")
ggsave("figures/GO_BP_dotplot.pdf", plot = p_go, width = 10, height = 9)
ggsave("figures/GO_BP_dotplot.png", plot = p_go, width = 10, height = 9, dpi = 300)

# ---- 4. KEGG pathway enrichment ----
kegg <- enrichKEGG(
    gene          = entrez_ids,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
)
kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

p_kegg <- dotplot(kegg, showCategory = 20, title = "KEGG Pathways — NSCLC DEGs")
ggsave("figures/KEGG_dotplot.pdf", plot = p_kegg, width = 10, height = 9)
ggsave("figures/KEGG_dotplot.png", plot = p_kegg, width = 10, height = 9, dpi = 300)

# ---- 5. Save full enrichment tables ----
write.csv(as.data.frame(go_bp), "results_GO_BP_enrichment.csv", row.names = FALSE)
write.csv(as.data.frame(kegg), "results_KEGG_enrichment.csv", row.names = FALSE)

cat("Saved: figures/GO_BP_dotplot.pdf/png\n")
cat("Saved: figures/KEGG_dotplot.pdf/png\n")
cat("Saved: results_GO_BP_enrichment.csv\n")
cat("Saved: results_KEGG_enrichment.csv\n")
