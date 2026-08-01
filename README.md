[README.md](https://github.com/user-attachments/files/30612220/README.md)
# Identification of Differentially Expressed Genes and Novel Therapeutic Targets in Non-Small Cell Lung Cancer (NSCLC) using RNA-seq and Systematic Druggability Screening

## Overview

This project uses RNA-seq data from lung adenocarcinoma patients to identify genes that are significantly dysregulated in tumor tissue, then applies a systematic druggability screen to shortlist genes that (a) drive cancer biology and (b) currently have **no approved drug** — i.e. genuinely novel therapeutic targets, rather than candidates for drug repurposing.

The pipeline moves from raw sequencing reads through differential expression analysis, druggability filtering, and biological validation via pathway enrichment.

---

## Dataset

**GSE40419** (Gene Expression Omnibus)
- 5 lung adenocarcinoma tumor samples: LC_S1, LC_S2, LC_S3, LC_S4, LC_S5
- 5 adjacent normal lung tissue samples: LC_C2, LC_C3, LC_C5, LC_C7, LC_C9
- Paired-end RNA-seq

---

## Pipeline

### Stage 1 — Data Collection
- Downloaded GSE40419 from GEO using `fasterq-dump` in Galaxy EU
- 10 paired-end FASTQ files imported (5 tumor + 5 normal)

### Stage 2 — RNA-seq Processing (Galaxy EU)
| Step | Tool |
|---|---|
| Quality control | FastQC |
| Trimming | Trimmomatic (ILLUMINACLIP + SLIDINGWINDOW 4:20 + MINLEN 36) |
| Alignment | HISAT2 (hg38 reference genome) |
| Read counting | featureCounts (paired-end, counted as 1 fragment) |
| Differential expression | DESeq2 (tumor vs. normal) |
| Annotation | annotateMyIDs (Entrez ID → gene symbol) |

**Result:** 2,387 significant DEGs (padj < 0.05, \|log2FC\| > 1)

### Stage 3 — Systematic Druggability Screening

Rather than identifying drug-repurposing candidates (genes with existing approved drugs), this project specifically filters for genes **without** an approved drug, to surface novel target opportunities.

**Filtering funnel:**

| Stage | Gene Count |
|---|---|
| All significant DEGs | 2,387 |
| Upregulated in tumor | 1,189 |
| No approved drug (DGIdb) | 63 |
| Confirmed oncogene/tumor suppressor (OncoKB) | **4** |

**Final novel therapeutic targets:**

| Gene | Role (OncoKB) |
|---|---|
| CDC25A | Oncogene |
| E2F2 | Tumor Suppressor |
| E2F5 | Tumor Suppressor |
| PAX5 | Oncogene |

### Stage 4 — Visualization (R)

All figures generated in R 4.6.1 / RStudio from the DESeq2 output and count matrix.

| Figure | Description |
|---|---|
| `volcano_plot_NSCLC.pdf/png` | All DEGs by log2FC vs. padj, with the 4 final novel targets labeled |
| `heatmap_top50_DEGs.pdf` | Expression pattern of top 50 DEGs across all 10 samples — clean separation between tumor and normal groups |
| `PCA_plot_NSCLC.pdf/png` | Sample clustering by condition (PC1: 39% variance, PC2: 21% variance) |
| `GO_BP_dotplot.pdf/png` | GO Biological Process enrichment — dominated by cell cycle, mitotic division, and chromosome segregation terms |
| `KEGG_dotplot.pdf/png` | KEGG pathway enrichment — Cell Cycle, FoxO signaling, IL-17 signaling |
| `novel_target_funnel.pdf/png` | Visual funnel of the druggability screening (2387 → 1189 → 63 → 4) |
| `final_novel_targets.pdf/png` | Final 4 targets by oncogene/tumor suppressor role |

**Notable finding:** GO/KEGG enrichment independently supports the target selection — cell cycle and mitotic pathways are among the most significantly enriched terms, directly reinforcing the biological relevance of CDC25A and the E2F family (both core cell-cycle regulators) as shortlisted targets.

### Stage 5 — Documentation
This README, repository organization, and reproducibility notes.

---

## Repository Structure

```
NSCLC_GSE40419_project/
├── data/
│   ├── merged_annotated_DESeq2.csv      # Full DEG list with gene annotations
│   ├── druggability_matrix.csv          # DGIdb + OncoKB screening results
│   └── featurecounts_matrix.csv         # Raw count matrix (28,395 genes × 10 samples)
├── scripts/
│   ├── 00_merge_featurecounts.R
│   ├── 01_volcano_plot.R
│   ├── 02_heatmap_top_degs.R
│   ├── 03_pca_plot.R
│   ├── 04_go_kegg_enrichment.R
│   └── 05_druggability_funnel.R
├── figures/
│   └── (all output PDFs/PNGs)
├── results/
│   ├── results_GO_BP_enrichment.csv
│   └── results_KEGG_enrichment.csv
└── README.md
```

---

## Tools & Packages

**Upstream pipeline (Galaxy EU):** FastQC, Trimmomatic, HISAT2, featureCounts, DESeq2, annotateMyIDs

**Druggability screening:** DGIdb (drug-gene interaction database), OncoKB (oncogene/tumor suppressor annotation)

**Visualization (R 4.6.1):** DESeq2, ggplot2, pheatmap, EnhancedVolcano, clusterProfiler, org.Hs.eg.db, dplyr, ggrepel

---

## How to Reproduce

1. Install R and RStudio
2. Install required packages:
   ```r
   install.packages(c("ggplot2", "pheatmap", "dplyr", "ggrepel"))
   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   BiocManager::install(c("DESeq2", "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db", "enrichplot"))
   ```
3. Place the three data files above into `data/`
4. Run scripts in `scripts/` in numeric order (00 → 05)
5. Outputs are saved to `figures/` and `results/`

---

## Limitations & Future Work

- Small sample size (n = 5 per group), consistent with the original GEO dataset — results would benefit from validation in a larger cohort
- Target discovery is transcriptomic and database-driven (DGIdb, OncoKB); the 4 shortlisted genes are computational candidates and would need experimental validation (e.g. knockdown studies, expanded cohort validation) before any therapeutic development claim
- No structural or molecular docking analysis was performed in this phase — the current scope is target identification, not drug design
- Future work could extend the screen to downregulated genes, or explore combination/polypharmacology angles for the 4 targets identified

---

## Author

*Anwesha Sarkar, Final Year, KIIT School of Biotechnology* 

