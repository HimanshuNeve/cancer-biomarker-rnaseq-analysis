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

## Repository Structure

```text
cancer-biomarker-rnaseq-analysis/
│
├── data/
│   └── metadata/
│       └── sample_metadata.csv
│
├── dashboard/
│   └── app.py
│
├── notebooks/
│   └── Section11_ML_Pipeline.ipynb
│
├── plots/
│   ├── ml_plots/
│   │   ├── ROC_curves.png
│   │   ├── confusion_matrices.png
│   │   └── feature_importance.png
│   │
│   ├── shap/
│   │   ├── shap_summary.png
│   │   ├── shap_bar.png
│   │
│   │
│   └── workflow/
│       ├── workflow_overview.png
│       └── pipeline_architecture.png
│
├── results/
│   └── ml_models/
│       ├── random_forest.pkl
│       ├── scaler.pkl
│       ├── feature_names.pkl
│       ├── label_encoder.pkl
│       ├── biomarker_importance.csv
│       └── model_comparison.csv
│
├── scripts/
│   ├── python_scripts/
│   │   └── shap_analysis.py
│   │
│   └── r_scripts/
│       ├── 01_download.R
│       ├── 01b_load_data.R
│       ├── 02_preprocess.R
│       ├── 03_deseq2.R
│       ├── 04_visualization.R
│       ├── 05_enrichment.R
│       └── export.R
│
├── LICENSE
├── Manual.pdf
├── README.md
├── requirements.txt
├── requirements_R.txt
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

## 🧠 Explainable AI (SHAP)

This project integrates SHAP (SHapley Additive Explanations)
to interpret machine learning predictions and identify
gene-level contributions toward breast cancer classification.

### SHAP Features

#### Global Explainability
- Global biomarker contribution analysis
- Interpretable Random Forest predictions
- Gene importance visualization
- Global feature contribution ranking

#### Dynamic Sample-Level Explainability
- Real-time SHAP interpretation for uploaded samples
- Interactive sample selection
- Dynamic waterfall-based prediction explanation
- Personalized gene contribution visualization

### Generated SHAP Outputs

| Plot | Description |
|---|---|
| `shap_summary.png` | Displays overall impact of genes on model predictions |
| `shap_bar.png` | Ranks genes by global contribution importance |
| Dynamic Waterfall Plot | Explains prediction contribution for uploaded samples |

### Biological Relevance

SHAP enables identification of:
- tumor-driving biomarkers
- protective gene signatures
- feature contribution patterns
- interpretable cancer prediction mechanisms

This improves the transparency and biological interpretability
of the RNA-Seq machine learning pipeline.

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
