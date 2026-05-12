# ============================================================
# Script: 01_download.R
# Purpose: Download GSE62944 dataset from NCBI GEO
# ============================================================

# Load required libraries
library(GEOquery)
library(data.table)

# Set working directory to your project folder
setwd("~/Desktop/Mini Project")

# Download GEO dataset
gse <- getGEO("GSE62944", GSEMatrix = TRUE, getGPL = FALSE)

# Save GEO object
saveRDS(gse, file = "data/raw_data/gse62944_object.rds")

# Download supplementary files
getGEOSuppFiles(
  "GSE62944",
  makeDirectory = TRUE,
  baseDir = "data/raw_data/"
)

cat("Dataset download complete!\n")
setwd("~/Desktop/Mini Project")
# ── Download supplementary count files ──
# The count matrix is provided as a supplementary file

getGEOSuppFiles(
  "GSE62944",
  makeDirectory = TRUE,
  baseDir = "data/raw_data/"
)

# ── List downloaded files ──
list.files("data/raw_data/GSE62944/")

cat("Dataset download complete!\n")
cat("Check data/raw_data/GSE62944/ for downloaded files\n")