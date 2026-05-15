# ============================================================
# File: dashboard/app.py
# Purpose: Cancer prediction web dashboard
# Run locally: streamlit run dashboard/app.py
# ============================================================

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import joblib
import os
import shap

# ============================================================
# Page Configuration
# ============================================================

st.set_page_config(
    page_title='Cancer Biomarker Predictor',
    layout='wide'
)

# ============================================================
# Global Paths
# ============================================================

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

MODEL_DIR = os.path.abspath(
    os.path.join(BASE_DIR, '..', 'results', 'ml_models')
)

# ============================================================
# Load Models
# ============================================================

@st.cache_resource
def load_models():

    rf = joblib.load(f'{MODEL_DIR}/random_forest.pkl')

    scaler = joblib.load(f'{MODEL_DIR}/scaler.pkl')

    genes = joblib.load(f'{MODEL_DIR}/feature_names.pkl')

    le = joblib.load(f'{MODEL_DIR}/label_encoder.pkl')

    return rf, scaler, genes, le


rf_model, scaler, feature_genes, le = load_models()

# ============================================================
# Header
# ============================================================

st.title('Cancer Biomarker Discovery Dashboard')

st.markdown(
    '### RNA-Seq + Machine Learning | Breast Cancer Classification'
)

st.markdown(
    '**Dataset:** GSE62944 (TCGA Breast Cancer)  |  '
    '**Best Model:** Random Forest'
)

st.divider()

# ============================================================
# Tabs
# ============================================================

tab1, tab2, tab3, tab4 = st.tabs([
    'Prediction',
    'Biomarker Importance',
    'Global SHAP Analysis',
    'Project Info'
])

# ============================================================
# TAB 1 — Prediction
# ============================================================

with tab1:

    st.header('Sample Prediction')

    st.markdown(
        '''
        Upload a CSV file containing gene expression values.

        Expected format:
        - Rows = Samples
        - Columns = Genes
        '''
    )

    uploaded_file = st.file_uploader(
        'Upload expression matrix (.csv)',
        type=['csv']
    )

    if uploaded_file is not None:

        # ====================================================
        # Read Uploaded File
        # ====================================================

        input_df = pd.read_csv(
            uploaded_file,
            index_col=0
        )

        # ====================================================
        # Auto Fix Orientation
        # ====================================================

        overlap = len(
            set(feature_genes).intersection(
                input_df.columns
            )
        )

        if overlap < 5:

            st.warning(
                'Detected gene names in rows. '
                'Automatically transposing matrix...'
            )

            input_df = input_df.T

        st.success(
            f'Data loaded successfully: '
            f'{input_df.shape[0]} samples × '
            f'{input_df.shape[1]} genes'
        )

        # ====================================================
        # Check Required Genes
        # ====================================================

        missing = [
            g for g in feature_genes
            if g not in input_df.columns
        ]

        if missing:

            st.error(
                f'Missing {len(missing)} required biomarker genes.'
            )

            with st.expander('View Missing Genes'):

                st.write(missing)

        else:

            # =================================================
            # Prediction
            # =================================================

            X = input_df[feature_genes]

            X_scaled = scaler.transform(X)

            predictions = rf_model.predict(X_scaled)

            probabilities = rf_model.predict_proba(
                X_scaled
            )[:, 1]

            results = pd.DataFrame({
                'Sample': input_df.index,
                'Prediction': le.inverse_transform(
                    predictions
                ),
                'Tumor Probability': probabilities.round(4)
            })

            # =================================================
            # Show Results
            # =================================================

            st.success('Prediction completed successfully!')

            st.subheader('Prediction Results')

            st.dataframe(results)
            
            # =================================================
            # Dynamic SHAP Explanation
            # =================================================

            (st.subheader("🧬 Sample-Specific SHAP Explanation"))

            try:

                explainer = shap.TreeExplainer(
                    rf_model
                )

                shap_values = explainer.shap_values(
                    X_scaled
                )

                selected_sample = st.selectbox(
                    "Select sample for SHAP explanation",
                    input_df.index
                )

                sample_index = list(
                    input_df.index
                ).index(selected_sample)

                st.markdown(
                    '''
                    SHAP explains how individual genes
                    contributed toward the prediction.
                    '''
                )

                shap_exp = shap.Explanation(
                    values=shap_values[:, :, 1][sample_index],
                    base_values=explainer.expected_value[1],
                    data=X.iloc[sample_index],
                    feature_names=feature_genes
                )

                fig_shap = plt.figure(
                    figsize=(10, 6)
                )

                shap.plots.waterfall(
                    shap_exp,
                    max_display=15,
                    show=False
                )

                st.pyplot(fig_shap)

            except Exception as e:

                st.warning(
                    f'SHAP explanation could not '
                    f'be generated: {e}'
                )
                
            # =================================================
            # Download Button
            # =================================================

            csv = results.to_csv(
                index=False
            ).encode('utf-8')

            st.download_button(
                label='📥 Download Prediction Results',
                data=csv,
                file_name='prediction_results.csv',
                mime='text/csv'
            )

            # =================================================
            # Probability Plot
            # =================================================

            st.subheader(
                'Tumor Probability Distribution'
            )

            fig, ax = plt.subplots(
                figsize=(8, 4)
            )

            ax.hist(
                probabilities,
                bins=20,
                color='#E74C3C',
                alpha=0.75,
                edgecolor='black'
            )

            ax.axvline(
                0.5,
                color='black',
                linestyle='--',
                label='Decision Boundary'
            )

            ax.set_xlabel(
                'Tumor Probability'
            )

            ax.set_ylabel(
                'Number of Samples'
            )

            ax.set_title(
                'Distribution of Tumor Prediction Probabilities'
            )

            ax.legend()

            st.pyplot(fig)

    else:

        st.info(
            'Please upload a CSV file to generate predictions.'
        )

# ============================================================
# TAB 2 — Biomarker Importance
# ============================================================

with tab2:

    st.header('Top Cancer Biomarker Genes')

    try:

        imp_df = pd.read_csv(
            f'{MODEL_DIR}/biomarker_importance.csv'
        )

        top20 = imp_df.head(20)

        fig2, ax2 = plt.subplots(
            figsize=(10, 6)
        )

        ax2.barh(
            top20['Gene'],
            top20['Importance'],
            color='#3498DB'
        )

        ax2.invert_yaxis()

        ax2.set_xlabel(
            'Feature Importance'
        )

        ax2.set_title(
            'Top 20 Biomarker Genes '
            '(Random Forest Model)'
        )

        st.pyplot(fig2)

        st.subheader('Top Biomarkers')

        st.dataframe(top20)

    except FileNotFoundError:

        st.warning(
            'Biomarker importance file not found. '
            'Please run the ML pipeline first.'
        )
# ── Tab 3: Explainable AI ──
with tab3:

    st.header("Explainable AI using SHAP")

    st.markdown("""
    SHAP (SHapley Additive exPlanations) identifies
    how individual genes contribute to tumor prediction.
    """)

    shap_summary = os.path.join(
        BASE_DIR,
        '..',
        'plots',
        'shap',
        'shap_summary.png'
    )

    shap_bar = os.path.join(
        BASE_DIR,
        '..',
        'plots',
        'shap',
        'shap_bar.png'
    )

    if os.path.exists(shap_summary):
        st.subheader("SHAP Summary Plot")
        st.image(shap_summary)

    if os.path.exists(shap_bar):
        st.subheader("Global Biomarker Contributions")
        st.image(shap_bar)

# ============================================================
# TAB 4 — Project Information
# ============================================================

with tab4:

    st.header('About This Project')

    st.markdown('''
    | Detail | Information |
    |---|---|
    | **Project Title** | Cancer Biomarker Discovery using RNA-Seq and Machine Learning |
    | **Dataset** | GSE62944 (TCGA Breast Cancer) |
    | **Samples** | 200 Tumor + 100 Normal |
    | **Pipeline** | RNA-Seq → DEGs → Machine Learning |
    | **Tools** | DESeq2, clusterProfiler, scikit-learn, Streamlit |
    | **ML Models** | Logistic Regression, Random Forest, SVM |
    | **Best Model** | Random Forest |
    | **Best ROC-AUC** | 1.000 |
    | **Explainable AI** | SHAP (SHapley Additive exPlanations) |
    ''')

    st.markdown('---')

    # ========================================================
    # Contributors
    # ========================================================

    st.subheader('👨‍💻 Contributors')

    st.markdown('''
    ### Himanshu Neve
    - Lead developer and principal researcher
    - Designed and implemented the complete RNA-Seq analysis workflow
    - Performed differential gene expression analysis using DESeq2
    - Developed machine learning classification pipelines
    - Implemented Random Forest, Logistic Regression, and SVM models
    - Conducted biomarker ranking and feature importance analysis
    - Managed GitHub repository, deployment workflow, and documentation
    - Integrated R and Python pipelines for reproducible bioinformatics analysis

    ### Sara
    - Assisted in dashboard testing and workflow validation
    - Contributed to UI refinement and usability improvements
    - Developed Streamlit-based interactive prediction dashboard
    - Supported debugging, project organization, and implementation review
    - Assisted in workflow development discussions and result verification
    - Contributed to presentation structure and project coordination
    ''')

    st.markdown('---')

    # ========================================================
    # Analysis Pipeline
    # ========================================================

    st.subheader('🧪 Analysis Pipeline')

    st.code(
        '''
Raw RNA-Seq Data (GSE62944)
            ↓
Data Preprocessing and Quality Control
            ↓
DESeq2 Differential Expression Analysis
            ↓
Top Differentially Expressed Genes
            ↓
Machine Learning Classification
            ↓
Random Forest Biomarker Ranking
            ↓
Interactive Prediction Dashboard
        ''',
        language='text'
    )

    # ========================================================
    # Project Highlights
    # ========================================================

    st.subheader('📌 Project Highlights')

    st.markdown('''
    - End-to-end RNA-Seq biomarker discovery workflow
    - Machine learning–based breast cancer classification
    - Automated biomarker ranking using Random Forest
    - Interactive prediction dashboard using Streamlit
    - Integrated R + Python hybrid bioinformatics pipeline
    - Reproducible GitHub-based research workflow
    - Deployment-ready bioinformatics web application
    ''')


# ============================================================
# Footer
# ============================================================

st.divider()

st.caption(
    'Cancer Biomarker Discovery Dashboard | '
    'RNA-Seq + Machine Learning | Developed by Himanshu Neve and Sara'
)
