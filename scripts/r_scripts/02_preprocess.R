# ============================================================
# Script: 02_preprocess.R
# Purpose: Quality control and preprocessing of count matrix
# ============================================================
library(DESeq2)
library(ggplot2)
library(dplyr)
setwd("~/Desktop/Mini Project")

# ── Load data ──
count_df <- readRDS('data/processed/count_matrix_BRCA.rds')
metadata <- readRDS('data/processed/metadata_BRCA.rds')
# ── Step 1: Ensure all count values are integers ──
# DESeq2 requires integer counts (not TPM or normalized values)
count_matrix <- round(as.matrix(count_df)) # Round to nearest integer
count_matrix[count_matrix < 0] <- 0 # Remove negative values if any
cat('Count matrix dimensions:', dim(count_matrix), '\n')


# ── Step 2: Check for missing values ──
na_count <- sum(is.na(count_matrix))
cat('Missing values:', na_count, '\n')
if (na_count > 0) {
  count_matrix[is.na(count_matrix)] <- 0 # Replace NA with 0
  cat('Replaced NAs with 0\n')
}


# ── Step 3: Filter low-expression genes ──
# Keep only genes with at least 10 counts in at least 10 samples
# This removes noise from genes that are essentially not expressed
keep <- rowSums(count_matrix >= 10) >= 10
count_matrix_filtered <- count_matrix[keep, ]
cat('Genes before filtering:', nrow(count_matrix), '\n')
cat('Genes after filtering:', nrow(count_matrix_filtered), '\n')


# ── Step 4: Check library sizes (sequencing depth per sample) ──
lib_sizes <- colSums(count_matrix_filtered)
cat('\nLibrary size summary:\n')
print(summary(lib_sizes))


# Plot library sizes
pdf('plots/deseq2_plots/library_sizes.pdf', width=10, height=5)
barplot(lib_sizes / 1e6,
        col = ifelse(metadata$condition == 'Tumor', '#E74C3C', '#3498DB'),
        main = 'Library Sizes per Sample',
        ylab = 'Total Counts (millions)',
        xlab = 'Samples',
        las = 2, cex.names = 0.5)
legend('topright', legend = c('Tumor','Normal'),
       fill = c('#E74C3C','#3498DB'))
dev.off()


# ── Step 5: Create DESeqDataSet ──
# Ensure metadata rows match count matrix columns
metadata_ordered <- metadata[colnames(count_matrix_filtered), ]
metadata_ordered$condition <- factor(metadata_ordered$condition,
                                     levels = c('Normal', 'Tumor'))
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = metadata_ordered,
  design = ~ condition # Compare tumor vs normal
)
cat('\nDESeqDataSet created successfully!\n')
print(dds)


# ── Step 6: Variance Stabilizing Transformation (VST) for visualization ──
# VST makes the variance roughly constant across the dynamic range
# This is NOT used for DESeq2 — only for visualization and ML
vst_data <- vst(dds, blind = TRUE)

# Extract the transformed matrix
vst_matrix <- assay(vst_data)
cat('\nVST transformation complete!\n')
cat('VST matrix dimensions:', dim(vst_matrix), '\n')

# ── Save objects ──
saveRDS(dds, 'data/processed/dds_BRCA.rds')
saveRDS(vst_data, 'data/processed/vst_BRCA.rds')
saveRDS(vst_matrix, 'data/processed/vst_matrix_BRCA.rds')
cat('\nPreprocessing complete! Objects saved.\n')

