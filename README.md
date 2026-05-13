# Cancer Biomarker Discovery using RNA-Seq and Machine Learning

This is our short bioinformatics project where we built a pipeline to identify potential cancer biomarkers from publicly available RNA-Seq data. The idea was simple — take real breast cancer gene expression data, find which genes behave differently in tumors vs normal tissue, and then see if a machine learning model can use those genes to predict cancer.

We used the TCGA breast cancer dataset (GSE62944) from NCBI GEO. It has RNA-Seq count data from real patient samples which made the results feel meaningful rather than just an academic exercise.

---
## Live Web Application

🔗 Streamlit Dashboard:
https://cancer-biomarker-rnaseq.streamlit.app

This interactive dashboard enables:
- Breast cancer prediction from RNA-Seq expression profiles
- Biomarker importance visualization
- Machine learning–based tumor classification
- Exploration of top cancer-associated genes
## What this project does

Starting from raw count data, the pipeline goes through:

- Quality filtering and normalization of raw RNA-Seq counts
- Differential expression analysis using DESeq2 to find genes that are significantly up or downregulated in tumor samples
- Volcano plots, heatmaps, and PCA to visualize the expression patterns
- GO and KEGG pathway enrichment to understand what biological processes are affected
- Machine learning models (Logistic Regression, Random Forest, SVM) trained on the top differentially expressed genes
- Feature importance analysis to extract the most predictive biomarker genes
- A small Streamlit web app where you can upload expression data and get a prediction

---

## Dataset

**GEO Accession:** GSE62944  
**Cancer type:** Breast Invasive Carcinoma (BRCA)  
**Samples used:** 50 tumor + 20 normal  
**Source:** TCGA via NCBI GEO  
**Paper:** Rahman et al., Scientific Data, 2015  
**Link:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944

---

## Folder structure

```
cancer-biomarker-rnaseq/
├── data/
│   ├── metadata/            sample group labels
│   └── processed/           normalized CSVs used in analysis
├── scripts/
│   ├── r_scripts/           01 to 05 — bioinformatics side
│   └── python_scripts/      06 to 09 — machine learning side
├── results/
│   ├── deseq2/              DEG tables
│   ├── enrichment/          GO and KEGG outputs
│   └── ml_models/           saved models and evaluation results
├── plots/
│   ├── deseq2_plots/        volcano, heatmap, PCA
│   └── ml_plots/            ROC curves, confusion matrix, feature importance
├── dashboard/
│   └── app.py               Streamlit app
├── docs/
│   └── session_info.txt     R package versions
├── requirements.txt
└── .gitignore
```

---

## How to run it

### Requirements
- R 4.2 or higher
- Python 3.9 or higher
- Around 8GB RAM
- Ubuntu or macOS

### 1. Clone the repo

```bash
git clone https://github.com/YOUR_USERNAME/cancer-biomarker-rnaseq.git
cd cancer-biomarker-rnaseq
```

### 2. Set up Python

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### 3. Set up R packages

Open R or RStudio and run:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "GEOquery", "clusterProfiler",
                       "org.Hs.eg.db", "EnhancedVolcano", "pheatmap"))
install.packages(c("ggplot2", "dplyr", "RColorBrewer", "data.table"))
```

### 4. Run the bioinformatics pipeline (R)

```bash
Rscript scripts/r_scripts/01_download.R
Rscript scripts/r_scripts/01b_load_data.R
Rscript scripts/r_scripts/02_preprocess.R
Rscript scripts/r_scripts/03_deseq2.R
Rscript scripts/r_scripts/04_visualization.R
Rscript scripts/r_scripts/05_enrichment.R
```

### 5. Run the ML pipeline (Python)

```bash
source venv/bin/activate
python3 scripts/python_scripts/06_ml_preprocess.py
python3 scripts/python_scripts/07_ml_models.py
python3 scripts/python_scripts/08_evaluation.py
python3 scripts/python_scripts/09_biomarkers.py
```

### 6. Launch the dashboard (optional)

```bash
streamlit run dashboard/app.py
```

Then open http://localhost:8501 in your browser.

---

## Results summary

All three models performed well on the test set. Random Forest gave the best overall results:

| Model | Accuracy | ROC-AUC |
|---|---|---|
| Logistic Regression | ~0.92 | ~0.95 |
| Random Forest | ~0.95 | ~0.97 |
| SVM | ~0.93 | ~0.96 |

The top biomarker genes identified by feature importance included well-known cancer-associated genes like MKI67, PCNA, and ERBB2, which gave us confidence the model was picking up on real biology rather than noise.

---

## Reproducibility

We made an effort to keep this reproducible:

- Random seeds are fixed (`set.seed(42)` in R, `random_state=42` in Python)
- Processed data CSVs are committed so you can skip the raw download step
- R package versions are saved in `docs/session_info.txt`
- Python package versions are pinned in `requirements.txt`

---

## Team

This was a 2-person project. We split the work by domain:

- **Person 1** handled the bioinformatics side — data download, DESeq2, enrichment analysis
- **Person 2** handled the ML side — feature engineering, model training, dashboard

---

## Tools used

R, DESeq2, clusterProfiler, ggplot2, pheatmap, EnhancedVolcano, Python, scikit-learn, pandas, matplotlib, seaborn, Streamlit
