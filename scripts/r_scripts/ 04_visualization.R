# ============================================================
# Script: 04_visualization.R
# Purpose: Generate all publication-quality plots
# ============================================================
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(DESeq2)
library(RColorBrewer)
library(dplyr)
setwd("~/Desktop/Mini Project")

# ── Load data ──
res_df <- readRDS('data/processed/deseq2_results.rds')
res <- readRDS('data/processed/deseq2_raw_results.rds')
dds <- readRDS('data/processed/dds_BRCA_analyzed.rds')
vst_data <- readRDS('data/processed/vst_BRCA.rds')
metadata <- readRDS('data/processed/metadata_BRCA.rds')

# =============================================================
# PLOT 1: VOLCANO PLOT
# =============================================================
# The EnhancedVolcano package creates beautiful, publication-ready volcano plots
png('plots/deseq2_plots/volcano_plot.png', width=2400, height=2000, res=200)
EnhancedVolcano(res_df,
                lab = res_df$gene, # Gene labels
                x = 'log2FoldChange', # X-axis: fold change
                y = 'padj', # Y-axis: adjusted p-value
                title = 'Volcano Plot: Breast Cancer Tumor vs Normal',
                subtitle = 'DESeq2 | GSE62944 | Threshold: padj<0.05, |log2FC|>1',
                pCutoff = 0.05, # Significance threshold
                FCcutoff = 1.0, # Fold change threshold
                pointSize = 2.0,
                labSize = 3.5,
                col = c('grey70','grey70','#3498DB','#E74C3C'),
                colAlpha = 0.7,
                legendLabels = c('NS','|log2FC|>1','padj<0.05','padj<0.05 & |log2FC|>1'),
                drawConnectors = TRUE,
                widthConnectors = 0.3,
                max.overlaps = 20
)
dev.off()
cat('Volcano plot saved!\n')

# =============================================================
# PLOT 2: HEATMAP OF TOP 50 DEGs
# =============================================================
# Get top 50 most significant DEGs
sig_degs <- read.csv('results/deseq2/significant_DEGs.csv')
top50_genes <- head(sig_degs$gene, 50)
# Extract VST expression values for top 50 genes
vst_matrix <- assay(vst_data)
heatmap_mat <- vst_matrix[top50_genes, ]
# Scale rows (z-score per gene) for better visualization
# Z-score: (value - mean) / SD — centers each gene around 0
heatmap_scaled <- t(scale(t(heatmap_mat)))
# Create annotation for samples (tumor vs normal)
annotation_col <- data.frame(
  Condition = metadata[colnames(heatmap_scaled), 'condition']
)
rownames(annotation_col) <- colnames(heatmap_scaled)
# Color annotations
ann_colors <- list(Condition = c(Tumor = '#E74C3C', Normal = '#3498DB'))
# Draw heatmap
png('plots/deseq2_plots/heatmap_top50DEGs.png', width=2800, height=3200, res=200)
pheatmap(heatmap_scaled,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = FALSE, # Sample names too long to show
         show_rownames = TRUE, # Show gene names
         cluster_rows = TRUE, # Cluster genes by expression pattern
         cluster_cols = TRUE, # Cluster samples
         color = colorRampPalette(c('#3498DB','white','#E74C3C'))(100),
         main = 'Top 50 DEGs - Breast Cancer Tumor vs Normal',
         fontsize_row = 8,
         fontsize = 10
)
dev.off()
cat('Heatmap saved!\n')

# =============================================================
# PLOT 3: PCA PLOT
# =============================================================
# PCA on VST-transformed data
# plotPCA is a DESeq2 built-in function for PCA visualization
pca_data <- plotPCA(vst_data,
                    intgroup = 'condition',
                    returnData = TRUE) # Return data for custom ggplot
# Get variance explained by each PC
percentVar <- round(100 * attr(pca_data, 'percentVar'))
# Create publication-quality PCA plot using ggplot2
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2,
                                 color = condition, shape = condition)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c(Normal = '#3498DB', Tumor = '#E74C3C')) +
  labs(
    title = 'PCA: Breast Cancer Samples (GSE62944)',
    subtitle = 'VST-transformed RNA-Seq data',
    x = paste0('PC1: ', percentVar[1], '% variance'),
    y = paste0('PC2: ', percentVar[2], '% variance'),
    color = 'Sample Type',
    shape = 'Sample Type'
  ) +
  stat_ellipse(level = 0.95, linetype = 2) + # 95% confidence ellipses
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face='bold'),
    legend.position = 'bottom'
  )
ggsave('plots/deseq2_plots/PCA_plot.png', pca_plot,
       width = 8, height = 7, dpi = 200)
cat('PCA plot saved!\n')
cat('\nAll visualizations complete!\n')
cat('Check plots/deseq2_plots/ for all figures\n')

# =============================================================
# PLOT 1B: MA PLOT
# =============================================================

png(
  'plots/deseq2_plots/MA_plot.png',
  width = 2200,
  height = 1800,
  res = 200
)

plotMA(
  res,
  ylim = c(-5, 5),
  main = 'MA Plot: Breast Cancer Tumor vs Normal'
)

dev.off()

cat('MA plot saved!\n')