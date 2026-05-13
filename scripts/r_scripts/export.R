library(DESeq2)
setwd("~/Desktop/Mini Project")

vst_matrix <- assay(readRDS('data/processed/vst_BRCA.rds'))
write.csv(vst_matrix, 'data/processed/vst_matrix_BRCA.csv')
file.exists('data/processed/vst_matrix_BRCA.csv')
