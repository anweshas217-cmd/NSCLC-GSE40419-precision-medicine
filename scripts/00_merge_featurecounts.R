## =============================================================
## Stage 4.0 — Merge 10 featureCounts files into one count matrix
## NSCLC GSE40419
## =============================================================

# ---- 1. List your 10 files ----
# Make sure all 10 are in data/ and named clearly, e.g.:
# LC_S1.tabular, LC_S2.tabular, LC_S3.tabular, LC_S4.tabular, LC_S5.tabular,
# LC_C2.tabular, LC_C3.tabular, LC_C5.tabular, LC_C7.tabular, LC_C9.tabular

files <- list.files("data/", pattern = "\\.tabular$", full.names = TRUE)
print(files)   # <- CHECK this prints exactly 10 files, in a sensible order

# Extract sample names from filenames (strips folder path and .tabular extension)
sample_names <- gsub("\\.tabular$", "", basename(files))
print(sample_names)   # <- CHECK these match LC_S1, LC_S2, ..., LC_C9 exactly

# ---- 2. Read each file, keep only Geneid + count column, rename count column ----
count_list <- lapply(seq_along(files), function(i) {
    df <- read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # First column = Geneid, second column = count (regardless of its messy header name)
    df <- df[, 1:2]
    colnames(df) <- c("Geneid", sample_names[i])
    df
})

# ---- 3. Merge all 10 on Geneid ----
merged <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), count_list)

# ---- 4. Check for NAs (would mean a gene was missing in some file — shouldn't happen, but verify) ----
cat("Any NAs after merge?", any(is.na(merged)), "\n")
cat("Total genes:", nrow(merged), " | Total samples:", ncol(merged) - 1, "\n")

# ---- 5. Save ----
write.csv(merged, "data/featurecounts_matrix.csv", row.names = FALSE)
cat("Saved: data/featurecounts_matrix.csv\n")
