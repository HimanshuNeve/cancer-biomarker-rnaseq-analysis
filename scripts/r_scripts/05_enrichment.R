# ============================================================
# Script: 05_enrichment.R
# Purpose: GO and KEGG pathway enrichment analysis
# ============================================================
library(clusterProfiler) # Core enrichment package
library(org.Hs.eg.db) # Human gene annotation
library(AnnotationDbi) # Annotation interface
library(ggplot2)
library(dplyr)
setwd("~/Desktop/Mini Project")

# ── Load significant DEGs ──
sig_degs <- read.csv('results/deseq2/significant_DEGs.csv')
gene_list <- sig_degs$gene
cat('Input genes for enrichment:', length(gene_list), '\n')

# ── Convert gene symbols to Entrez IDs ──
# clusterProfiler requires Entrez Gene IDs (not HGNC symbols)
gene_entrez <- bitr(gene_list,
                    fromType = 'SYMBOL',
                    toType = 'ENTREZID',
                    OrgDb = org.Hs.eg.db)
cat('Genes converted to Entrez IDs:', nrow(gene_entrez), '\n')


# =============================================================
# PART A: GO ENRICHMENT ANALYSIS
# =============================================================
# Biological Process GO enrichment
ego_BP <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = 'BP', # Biological Process
  pAdjustMethod = 'BH', # Benjamini-Hochberg correction
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE # Show gene symbols instead of IDs
)
cat('\nTop GO Biological Process terms:\n')
print(head(as.data.frame(ego_BP)[, c('Description','GeneRatio','p.adjust')], 10))

# ── GO Dot Plot (Top 20 Terms) ──
png('plots/deseq2_plots/GO_BP_dotplot.png', width=2400, height=2400, res=200)
dotplot(ego_BP, showCategory=20,
        title='GO Biological Process Enrichment (Breast Cancer DEGs)')
dev.off()

# ── GO Bar Plot ──
png('plots/deseq2_plots/GO_BP_barplot.png', width=2400, height=2000, res=200)
barplot(ego_BP, showCategory=15,
        title='Top 15 GO Biological Processes')
dev.off()


# =============================================================
# PART B: KEGG PATHWAY ANALYSIS
# =============================================================
ekegg <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = 'hsa', # hsa = Homo sapiens
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH'
)
cat('\nTop KEGG pathways:\n')
print(head(as.data.frame(ekegg)[, c('Description','GeneRatio','p.adjust')], 10))

# ── KEGG Dot Plot ──
png('plots/deseq2_plots/KEGG_dotplot.png', width=2400, height=2000, res=200)
dotplot(ekegg, showCategory=20,
        title='KEGG Pathway Enrichment (Breast Cancer DEGs)')
dev.off()

# ── Save Enrichment Results ──
write.csv(as.data.frame(ego_BP), 'results/enrichment/GO_BP_results.csv')
write.csv(as.data.frame(ekegg), 'results/enrichment/KEGG_results.csv')
cat('\nEnrichment analysis complete!\n')