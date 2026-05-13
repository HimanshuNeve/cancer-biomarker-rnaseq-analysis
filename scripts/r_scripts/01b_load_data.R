# ============================================================
# Script: 01b_load_data.R
# Purpose: Load and prepare BRCA expression matrix
# UPDATED: Using 200 tumor + 100 normal samples
# ============================================================

library(data.table)
library(dplyr)

setwd("~/Desktop/Mini Project")

# ============================================================
# Load tumor TPM matrix
# ============================================================

cat("Loading tumor matrix...\n")

tumor_counts <- fread(
  "data/raw_data/GSE62944/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_TPM.txt.gz",
  sep = "\t",
  header = TRUE
)

# ============================================================
# Load normal TPM matrix
# ============================================================

cat("Loading normal matrix...\n")

normal_counts <- fread(
  "data/raw_data/GSE62944/GSM1697009_06_01_15_TCGA_24.normal_Rsubread_TPM.txt.gz",
  sep = "\t",
  header = TRUE
)

cat("Files loaded successfully!\n")

# ============================================================
# Check dimensions
# ============================================================

cat("\nTumor matrix dimensions:\n")
cat("Genes:", nrow(tumor_counts), "\n")
cat("Samples:", ncol(tumor_counts)-1, "\n")

cat("\nNormal matrix dimensions:\n")
cat("Genes:", nrow(normal_counts), "\n")
cat("Samples:", ncol(normal_counts)-1, "\n")

# ============================================================
# Keep gene column
# ============================================================

gene_col <- colnames(tumor_counts)[1]

# ============================================================
# Select 200 tumor samples
# ============================================================

tumor_subset <- tumor_counts[
  ,
  c(gene_col, colnames(tumor_counts)[2:201]),
  with = FALSE
]

# ============================================================
# Select 100 normal samples
# ============================================================

normal_subset <- normal_counts[
  ,
  c(gene_col, colnames(normal_counts)[2:101]),
  with = FALSE
]

# ============================================================
# Rename gene column
# ============================================================

colnames(tumor_subset)[1]  <- "Gene"
colnames(normal_subset)[1] <- "Gene"

# ============================================================
# Merge matrices
# ============================================================

count_matrix <- merge(
  tumor_subset,
  normal_subset,
  by = "Gene"
)

cat("\nSelected samples:\n")
cat("Tumor:", ncol(tumor_subset)-1, "\n")
cat("Normal:", ncol(normal_subset)-1, "\n")
cat("Total:", ncol(count_matrix)-1, "\n")

# ============================================================
# Convert to data frame
# ============================================================

count_df <- as.data.frame(count_matrix)

rownames(count_df) <- count_df$Gene

count_df$Gene <- NULL

# ============================================================
# Final matrix dimensions
# ============================================================

cat("\nFinal expression matrix:\n")
cat("Genes:", nrow(count_df), "\n")
cat("Samples:", ncol(count_df), "\n")

# ============================================================
# Create metadata
# ============================================================

tumor_names  <- colnames(tumor_subset)[-1]
normal_names <- colnames(normal_subset)[-1]

metadata <- data.frame(
  sample = c(tumor_names, normal_names),
  condition = c(
    rep("Tumor", length(tumor_names)),
    rep("Normal", length(normal_names))
  )
)

rownames(metadata) <- metadata$sample

cat("\nMetadata summary:\n")
print(table(metadata$condition))

# ============================================================
# Create folders if missing
# ============================================================

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Save files
# ============================================================

cat("\nSaving processed files...\n")

saveRDS(count_df,
        "data/processed/count_matrix_BRCA.rds")

saveRDS(metadata,
        "data/processed/metadata_BRCA.rds")

write.csv(count_df,
          "data/processed/count_matrix_BRCA.csv")

write.csv(metadata,
          "data/metadata/sample_metadata.csv")

cat("\nFiles saved successfully!\n")

cat("\nGenerated files:\n")
cat("✔ count_matrix_BRCA.rds\n")
cat("✔ metadata_BRCA.rds\n")
cat("✔ count_matrix_BRCA.csv\n")
cat("✔ sample_metadata.csv\n")

cat("\nReady for 02_preprocess.R\n")