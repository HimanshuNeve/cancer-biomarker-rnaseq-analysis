# ============================================================
# Script: 01b_load_data.R
# Purpose: Load, subset, and save BRCA expression data
# ============================================================

# Load libraries
library(data.table)
library(dplyr)

# Set working directory
setwd("~/Desktop/Mini Project")

# ============================================================
# Load tumor TPM matrix
# ============================================================

tumor_counts <- fread(
  "data/raw_data/GSE62944/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_TPM.txt.gz",
  sep = "\t",
  header = TRUE
)

# ============================================================
# Load normal TPM matrix
# ============================================================

normal_counts <- fread(
  "data/raw_data/GSE62944/GSM1697009_06_01_15_TCGA_24.normal_Rsubread_TPM.txt.gz",
  sep = "\t",
  header = TRUE
)

# ============================================================
# Check dimensions
# ============================================================

cat("Tumor matrix dimensions:\n")
print(dim(tumor_counts))

cat("\nNormal matrix dimensions:\n")
print(dim(normal_counts))

# ============================================================
# Select subset of samples
# ============================================================

# Keep first column (gene names)
gene_col <- colnames(tumor_counts)[1]

# Select 50 tumor samples
tumor_subset <- tumor_counts[
  ,
  c(gene_col, colnames(tumor_counts)[2:51]),
  with = FALSE
]

# Select 20 normal samples
normal_subset <- normal_counts[
  ,
  c(gene_col, colnames(normal_counts)[2:21]),
  with = FALSE
]

# Rename gene column
colnames(tumor_subset)[1] <- "Gene"
colnames(normal_subset)[1] <- "Gene"

# Merge matrices
combined_counts <- merge(
  tumor_subset,
  normal_subset,
  by = "Gene"
)

# Convert to dataframe
count_df <- as.data.frame(combined_counts)

rownames(count_df) <- count_df$Gene
count_df$Gene <- NULL

# ============================================================
# Create metadata
# ============================================================

sample_names <- colnames(count_df)

condition <- c(
  rep("Tumor", 50),
  rep("Normal", 20)
)

metadata <- data.frame(
  sample = sample_names,
  condition = condition
)

rownames(metadata) <- metadata$sample

# ============================================================
# Save files
# ============================================================

saveRDS(
  count_df,
  "data/processed/count_matrix_BRCA.rds"
)

saveRDS(
  metadata,
  "data/processed/metadata_BRCA.rds"
)

cat("Reduced dataset saved successfully!\n")