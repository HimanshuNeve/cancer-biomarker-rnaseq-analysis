# ============================================================
# SHAP Explainable AI Analysis
# ============================================================

import os
import joblib
import shap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

BASE = os.path.expanduser("~/Desktop/Mini Project")

MODEL_DIR = f"{BASE}/results/ml_models"
PLOT_DIR  = f"{BASE}/plots/shap"

os.makedirs(PLOT_DIR, exist_ok=True)

# ------------------------------------------------------------
# Load models
# ------------------------------------------------------------

rf_model = joblib.load(f"{MODEL_DIR}/random_forest.pkl")
scaler   = joblib.load(f"{MODEL_DIR}/scaler.pkl")
genes    = joblib.load(f"{MODEL_DIR}/feature_names.pkl")

# ------------------------------------------------------------
# Load expression matrix
# ------------------------------------------------------------

df = pd.read_csv(
    f"{BASE}/data/processed/vst_matrix_BRCA.csv",
    index_col=0
)
df = df.T
X = df[genes]

# Scale
X_scaled = scaler.transform(X)

# ------------------------------------------------------------
# SHAP explainer
# ------------------------------------------------------------

explainer = shap.TreeExplainer(rf_model)

shap_values = explainer.shap_values(X_scaled)

# ============================================================
# SHAP Summary Plot
# ============================================================

plt.figure(figsize=(10, 6))

shap.summary_plot(
    shap_values,
    X,
    feature_names=genes,
    show=False
)

plt.savefig(
    f"{PLOT_DIR}/shap_summary.png",
    dpi=300,
    bbox_inches='tight'
)

plt.close()

# ============================================================
# SHAP Feature Importance Bar Plot
# ============================================================

plt.figure(figsize=(10, 6))

shap.plots.bar(
    shap.Explanation(
        values=np.abs(shap_values).mean(0),
        feature_names=genes
    ),
    show=False
)

plt.savefig(
    f"{PLOT_DIR}/shap_bar.png",
    dpi=300,
    bbox_inches='tight'
)

plt.close()


print("SHAP analysis completed successfully")
print(f"Plots saved in: {PLOT_DIR}")