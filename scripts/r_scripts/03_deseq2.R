# ============================================================
# Script: 03_deseq2.R
# Purpose: Differential expression analysis with DESeq2
# Person 1 Task — Day 2
# ============================================================
library(DESeq2)
library(dplyr)
library(readr)
library(apeglm)
setwd("~/Desktop/Mini Project")

# ── Load preprocessed DESeqDataSeat ──
dds <- readRDS('data/processed/dds_BRCA.rds')

# ── Run DESeq2 ──
# This performs:
# size factor estimation,
# dispersion estimation,
# and hypothesis testing
# It may take 2-5 minutes depending on dataset size
cat('Running DESeq2... This may take a few minutes.\n')
dds <- DESeq(dds)


# ── Extract Results ──
# contrast: compare Tumor vs Normal
# reference level is 'Normal' (set in preprocessing script)
res <- results(dds,
               contrast = c('condition', 'Tumor', 'Normal'),
               alpha = 0.05) # FDR threshold for independent filtering

# Save raw DESeq2 results for MA plot
saveRDS(res, 'data/processed/deseq2_raw_results.rds')

# Print summary of results
cat('\n── DESeq2 Results Summary ──\n')
summary(res)

# ── Apply LFC Shrinkage (for better fold-change estimates) ──
# apeglm shrinkage reduces noise in low-count gene estimates
# This is recommended for visualization but NOT required
res_shrunk <- lfcShrink(dds,
                        coef = 'condition_Tumor_vs_Normal',
                        type = 'apeglm')

# ── Convert to data frame ──
res_df <- as.data.frame(res_shrunk)
res_df$gene <- rownames(res_df)


# ── Filter significant DEGs ──
# Criteria: adjusted p-value < 0.05 AND |log2FC| > 1
sig_degs <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) > 1) %>%
  arrange(padj)

# Separate upregulated and downregulated
up_genes <- sig_degs %>% filter(log2FoldChange > 1)
down_genes <- sig_degs %>% filter(log2FoldChange < -1)
cat('\n── Significant DEGs ──\n')
cat('Total significant DEGs:', nrow(sig_degs), '\n')
cat('Upregulated in tumor:', nrow(up_genes), '\n')
cat('Downregulated in tumor:', nrow(down_genes), '\n')

# ── Preview top 10 DEGs ──
cat('\nTop 10 most significant DEGs:\n')
print(head(sig_degs[, c('gene','log2FoldChange','padj','baseMean')], 10))


# ── Save Results ──
write_csv(res_df, 'results/deseq2/all_genes_results.csv')
write_csv(sig_degs, 'results/deseq2/significant_DEGs.csv')
write_csv(up_genes, 'results/deseq2/upregulated_genes.csv')
write_csv(down_genes, 'results/deseq2/downregulated_genes.csv')

# Save top 50 DEGs (used for ML feature selection)
top50_degs <- head(sig_degs, 50)
write_csv(top50_degs, 'results/deseq2/top50_DEGs.csv')

# Save DESeq object
saveRDS(dds, 'data/processed/dds_BRCA_analyzed.rds')
saveRDS(res_df, 'data/processed/deseq2_results.rds')
cat('\nDESeq2 analysis complete!\n')
cat('Results saved to results/deseq2/\n')