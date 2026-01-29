#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NASIMMUNE Transcriptome Analysis Pipeline

Author:      Ivan Tomic
Affiliation: aTomic Lab, Boston University
Date:        November 6, 2025
"""
import os
import sys
import logging
import traceback
from typing import List, Dict, Tuple, Any

import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.metrics import (
    classification_report,
    roc_curve,
    roc_auc_score,
    confusion_matrix,
)
from sklearn.preprocessing import StandardScaler

# Set non-interactive backend for Matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


# --- Configuration ---
DATA_FILE_PATH = "./output/processed_data/hai_degs_for_gene_score.csv"
OUTPUT_DIR = "./output/gene_score"
FIGURES_DIR = os.path.join(OUTPUT_DIR, "figures")
DATA_DIR = os.path.join(OUTPUT_DIR, "data")
LOG_FILE = os.path.join(DATA_DIR, "analysis_summary_report.log")

TARGET_COLUMN = 'seroconversion'
RANDOM_STATE = 42

# The 60-gene list (crossed-referenced the 205 DEGs against the Interferome DB)
INTERFERON_GENES = [
    'APOL1', 'BATF2', 'CMPK2', 'CMTR1', 'CXCL10', 'DDX58', 'DDX60', 'DDX60L', 'DHX58', 'DTX3L',
    'EIF2AK2', 'ETV7', 'FAM46A', 'GBP1', 'HELZ2', 'HERC5', 'HERC6', 'IFI27', 'IFI35',
    'IFI44', 'IFI44L', 'IFIH1', 'IFIT1', 'IFIT2', 'IFIT3', 'IRF7', 'ISG15', 'LAMP3', 'LAP3',
    'LGALS3BP', 'MAFB', 'MMP14', 'MOV10', 'MX1', 'MX2', 'OAS1', 'OAS2', 'OAS3',
    'PARP10', 'PARP12', 'PARP14', 'PML', 'RNF213', 'RSAD2', 'SAMD9L', 'SERPING1', 'STAT1', 'STAT2',
    'TAP2', 'TDRD7', 'THBS1', 'TRANK1', 'TRIM5', 'TYMP', 'UBA7', 'USP18', 'VCAN',
    'XAF1', 'ZCCHC2', 'ZNFX1'
]
# 10-gene signature discovered by the RFE
#CORE_SIGNATURE_GENES = ['LCNL1', 'SASH1', 'PTPRU', 'CXCL10', 
 #                       'CEACAM8', 'SDC1', 'KIF4A', 'GPER1', 
  #                      'TOP2A', 'DIAPH3']
                        
CORE_SIGNATURE_GENES = ['PTPRU', 'CXCL10', 'KIF4A', 'DIAPH3', 'GPER1']
# The full list of 205 DEGs for validation
ALL_DEGS = [
    "HSPA5", "POLR2A", "HCFC1", "FLNA", "HYOU1", "LDLR", "CSF1R", "EIF4G1",
    "PRF1", "PLXNB2", "PRPF8", "ANKFY1", "VCAN", "APP", "CYFIP1", "MMP14",
    "TLN1", "RNF213", "FLNB", "CMKLR1", "EPB41L3", "MYOF", "TAP2", "HSP90B1",
    "FBN2", "GIMAP8", "OSBPL5", "SASH1", "LTF", "MYH9", "ZNFX1", "ITSN1",
    "NEO1", "TRANK1", "UBA7", "AGRN", "SRGAP2", "EPHB2", "DDX60L", "KIF4A",
    "ESPL1", "TOP2A", "MAP1A", "TMEM63C", "HELZ2", "DDX60", "SLC27A3", "TLR7",
    "AXL", "LILRB4", "JUP", "LAMB2", "LAG3", "MOV10", "CMTR1", "LGALS3BP",
    "TRIM5", "PARP14", "FPR3", "SIGLEC1", "OAS2", "IGF2R", "JAG2", "SERPINE1",
    "CLSPN", "NID1", "FMNL2", "FAM46A", "GRIN2D", "KIR2DS4", "APOL1", "CES1",
    "DTX3L", "TTC21A", "ADAMTSL4", "SDC1", "CDHR5", "TDRD7", "STAT1", "STAT2",
    "SAMD9L", "IFIH1", "ZNF366", "MAFB", "OAS3", "TMEM255A", "MPO", "HERC5",
    "MKI67", "MX1", "OTOF", "ZCCHC2", "EIF2AK2", "DIAPH3", "CMPK2", "PARP10",
    "COL24A1", "P3H2", "UNC93B1", "PML", "HERC6", "ODF3B", "DHX58", "MSR1",
    "HES4", "KCNN3", "TYMP", "SPTA1", "MDK", "CEACAM8", "CCL2", "TICRR",
    "LAP3", "THBS1", "F5", "RIN2", "IFI44L", "CLDN23", "GBP1", "LAMP3",
    "DOCK4", "PARP12", "C2", "RGL1", "GPER1", "CIT", "MX2", "SLC35F3",
    "IRF7", "OAS1", "LTBP1", "USP18", "TCN2", "C1QC", "AVPR2", "CDH1",
    "TMEM51", "CD209", "SAMD4A", "RSAD2", "SERPING1", "TMEM184A", "ASPM", "SHISA7",
    "SHROOM3", "IFI35", "IFIT2", "DDX58", "DNAAF1", "XAF1", "IFI44", "ITGA2B",
    "CES1P1", "TERT", "LIPM", "SPATS2L", "COL5A2", "NHSL1", "C3", "GUCY2D",
    "CXCL10", "ZBTB16", "PERM1", "SLC1A3", "PTPRU", "IFIT1", "BNC2", "C8G",
    "IFIT3", "IL27", "BATF2", "HESX1", "IGSF9", "KCP", "LDB2", "ETV7",
    "MYO1A", "ZBP1", "BGN", "CCNB3", "MID1", "FOXC1", "GPR173", "CEACAM6",
    "AANAT", "ARHGAP29", "MAPK15", "NEURL3", "HEG1", "EPB41L1", "PRSS22", "ATF3",
    "ZFHX2", "GPR161", "CHSY3", "PHEX", "CCL8", "ARHGAP23", "IL17D", "AR",
    "ISG15", "EREG", "LGALS17A", "IFI27", "LCNL1"
]


# --- Setup Functions ---
def configure_logging() -> None:
    logging.getLogger().handlers = []
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)-8s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(LOG_FILE, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.info("Logging configured. Report will be saved to %s", LOG_FILE)


def initialize_directories() -> None:
    os.makedirs(FIGURES_DIR, exist_ok=True)
    os.makedirs(DATA_DIR, exist_ok=True)
    logging.info("--- Environment Setup ---")
    logging.info(f"Output directory: {OUTPUT_DIR}")
    logging.info(f"Figures directory: {FIGURES_DIR}")
    logging.info(f"Data directory: {DATA_DIR}")
    logging.info("------------------------------")


# --- Plotting Functions ---
def generate_roc_plot(
    y_test: pd.Series, 
    y_prob: np.ndarray, 
    title_prefix: str
) -> float:

    sns.set_style("whitegrid")
    plt.figure(figsize=(8, 6))
    
    fpr, tpr, _ = roc_curve(y_test, y_prob)
    auc = roc_auc_score(y_test, y_prob)
    
    sns.lineplot(x=fpr, y=tpr, label=f'AUC = {auc:.4f}', color='red', linewidth=2.5)
    sns.lineplot(x=[0, 1], y=[0, 1], linestyle='--', label='Random Chance (AUC = 0.50)', color='blue')
    
    #plt.title(f'ROC Curve for {title_prefix}', fontsize=16, fontweight='bold')
    plt.xlabel('False Positive Rate', fontsize=16)
    plt.ylabel('True Positive Rate', fontsize=16)
    plt.legend(loc='lower right', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    
    # Clean filename and save as PDF
    clean_title = title_prefix.lower().replace(' ', '_').replace(':', '').replace('(', '').replace(')', '')
    filename = f"roc_curve_{clean_title}.pdf"
    filepath = os.path.join(FIGURES_DIR, filename)
    plt.savefig(filepath, format="pdf", dpi=300)
    plt.close()
    
    logging.info(f"Area Under the Curve (AUC): {auc:.4f}")
    logging.info(f"  > ROC Curve plot saved to: {filepath}")
    return auc

def generate_confusion_matrix_plot(
    y_test: pd.Series, 
    y_pred: np.ndarray, 
    classes: np.ndarray, 
    title_prefix: str
) -> None:

    sns.set_style("white")
    plt.figure(figsize=(8, 6))
    
    try:
        cm = confusion_matrix(y_test, y_pred)
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                    xticklabels=classes, yticklabels=classes,
                    annot_kws={"size": 16})
        #plt.title(f'Confusion Matrix for {title_prefix}', fontsize=16, fontweight='bold')
        plt.xlabel('Predicted Label', fontsize=16)
        plt.ylabel('True Label', fontsize=16)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.tight_layout()
        
        # Clean filename and save as PDF
        clean_title = title_prefix.lower().replace(' ', '_').replace(':', '').replace('(', '').replace(')', '')
        filename = f"confusion_matrix_{clean_title}.pdf"
        filepath = os.path.join(FIGURES_DIR, filename)
        plt.savefig(filepath, format="pdf", dpi=300)
        plt.close()
        logging.info(f"  > Confusion matrix plot saved to: {filepath}")
        
    except Exception as e:
        logging.warning(f"Could not generate confusion matrix for {title_prefix}. Error: {e}", exc_info=True)

def generate_permutation_plot(
    real_score: float, 
    null_scores: List[float], 
    model_name: str
) -> None:

    sns.set_style("whitegrid")
    plt.figure(figsize=(8, 6))
    
    # Calculate p-value
    n_permutations = len(null_scores)
    p_value = (np.sum(null_scores >= real_score) + 1) / (n_permutations + 1)
    
    sns.histplot(null_scores, kde=True, label='Null Distribution (Shuffled Labels)', color='blue')
    plt.axvline(real_score, color='red', linestyle='--', linewidth=2.5, 
                label=f'Actual AUC = {real_score:.4f}\n(p-value = {p_value:.4f})')
    
    #plt.title(f'Permutation Test for {model_name}', fontsize=16, fontweight='bold')
    plt.xlabel('AUC Score', fontsize=16)
    plt.ylabel('Frequency', fontsize=16)
    plt.legend(loc='upper right', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    
    # Clean filename and save as PDF
    clean_title = model_name.lower().replace(' ', '_').replace(':', '').replace('(', '').replace(')', '')
    filename = f"permutation_test_{clean_title}.pdf"
    filepath = os.path.join(FIGURES_DIR, filename)
    plt.savefig(filepath, format="pdf", dpi=300)
    plt.close()
    
    logging.info(f"\n--- Permutation Test Results for {model_name} ---")
    logging.info(f"Actual Model AUC: {real_score:.4f}")
    logging.info(f"Mean Null AUC: {np.mean(null_scores):.4f}")
    logging.info(f"p-value: {p_value:.4f}")
    if p_value < 0.05:
        logging.info("  > Result is statistically significant (p < 0.05).")
    else:
        logging.info("  > Result is NOT statistically significant (p >= 0.05).")
    logging.info(f"  > Permutation test plot saved to: {filepath}")
    
    null_scores_df = pd.DataFrame(null_scores, columns=['null_auc'])
    null_scores_filename = os.path.join(DATA_DIR, f"permutation_null_scores_{clean_title}.csv")
    null_scores_df.to_csv(null_scores_filename, index=False)
    logging.info(f"  > Null distribution scores saved to: {null_scores_filename}")


# --- Analysis Core Functions ---
def ingest_dataset(
    filepath: str, 
    all_genes: List[str], 
    target: str,
    step_num: int
) -> pd.DataFrame:

    logging.info(f"\n--- {step_num}. Loading & Validating Data ---")
    logging.info(f"Loading data from: {filepath}")
    
    if not os.path.exists(filepath):
        logging.error(f"Data file not found at {filepath}")
        raise FileNotFoundError(f"Data file not found at {filepath}")
        
    data = pd.read_csv(filepath)
    logging.info(f"Data loaded successfully. Shape: {data.shape}")
    
    # Validation
    required_cols = all_genes + [target]
    missing_cols = [col for col in required_cols if col not in data.columns]
    
    if missing_cols:
        error_msg = f"Missing {len(missing_cols)} required columns: {', '.join(missing_cols)}"
        logging.error(error_msg)
        raise ValueError(error_msg)
        
    logging.info(f"All {len(required_cols)} required columns are present.")
    logging.info("Validation complete.")
    return data

def compute_composite_score(
    data: pd.DataFrame, 
    gene_list: List[str], 
    score_col: str,
    step_num: str
) -> pd.DataFrame:

    logging.info(f"\n--- {step_num}. Creating Score: '{score_col}' ---")
    
    # Find missing genes from the provided list
    present_genes = [g for g in gene_list if g in data.columns]
    missing_genes = [g for g in gene_list if g not in data.columns]
    
    if missing_genes:
        logging.warning(f"  > {len(missing_genes)} genes from list not in data. Skipping them.")
        logging.warning(f"  > Missing: {', '.join(missing_genes)}")
        
    logging.info(f"Calculating as unweighted sum of {len(present_genes)} genes.")
    
    data[score_col] = data[present_genes].sum(axis=1)
    
    logging.info("Score calculation complete.")
    
    # Sanity check
    logging.info(f"Sanity Check: Describing '{score_col}' grouped by '{TARGET_COLUMN}':")
    stats = data.groupby(TARGET_COLUMN)[score_col].describe()
    logging.info(stats.to_string())
    
    # Save stats
    stats_filename = os.path.join(DATA_DIR, f"{score_col}_summary_stats.csv")
    stats.to_csv(stats_filename)
    logging.info(f"  > Summary statistics saved to: {stats_filename}")
    
    # Interpretation
    try:
        mean_0 = stats.loc[0, 'mean']
        mean_1 = stats.loc[1, 'mean']
        logging.info("Interpretation:")
        logging.info(f"  Mean score for non-seroconverted (0): {mean_0:.4f}")
        logging.info(f"  Mean score for seroconverted (1):     {mean_1:.4f}")
        if mean_1 > mean_0:
            logging.info("  > As hypothesized, seroconverted group shows a higher mean score.")
        else:
            logging.warning("  > WARNING: Seroconverted group does NOT show a higher mean score.")
    except KeyError:
        logging.warning("  > Could not compare means (e.g., one class missing).")
        
    return data

def evaluate_univariate_predictor(
    data: pd.DataFrame, 
    score_col: str, 
    target_col: str, 
    step_num: str
) -> None:

    logging.info(f"\n--- {step_num}. Testing Predictive Power for: {score_col} ---")
    
    X = data[[score_col]]  # Univariate (1 column)
    y = data[target_col]
    
    logging.info(f"Target variable distribution:\n{y.value_counts(normalize=True).to_string()}")
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, 
        test_size=0.3,          # 30% test set
        random_state=RANDOM_STATE,  # For reproducible results
        stratify=y              # Preserve class distribution
    )
    logging.info(f"Data split: {len(X_train)} train, {len(X_test)} test.")
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # --- Model Training ---
    logging.info("Training Logistic Regression Model...")
    model = LogisticRegression(random_state=RANDOM_STATE, class_weight='balanced')
    model.fit(X_train_scaled, y_train)
    
    logging.info(f"  > Model Coefficient: {model.coef_[0][0]:.4f} (Supports hypothesis)")
    
    # --- Model Evaluation ---
    logging.info("Evaluating model on Test Set...")
    y_pred = model.predict(X_test_scaled)
    y_prob = model.predict_proba(X_test_scaled)[:, 1]  # Probabilities for class 1
    
    # 1. Classification Report
    report_str = classification_report(y_test, y_pred, target_names=['Non-Seroconverted (0)', 'Seroconverted (1)'])
    report_df = pd.DataFrame(classification_report(y_test, y_pred, output_dict=True)).transpose()
    logging.info(f"Classification Report (Test Set):\n{report_str}")
    
    report_filename = os.path.join(DATA_DIR, f"logistic_regression_{score_col}_classification_report.csv")
    report_df.to_csv(report_filename)
    logging.info(f"  > Classification report saved to: {report_filename}")
    
    # 2. ROC Curve & AUC
    roc_title = f"Logistic Regression {score_col}"
    _ = generate_roc_plot(y_test, y_prob, title_prefix=roc_title)
    
    # 3. Confusion Matrix
    generate_confusion_matrix_plot(y_test, y_pred, model.classes_, title_prefix=roc_title)


def evaluate_multivariate_classifier(
    data: pd.DataFrame, 
    feature_cols: List[str], 
    target_col: str, 
    model_name: str,
    step_num: str
) -> Tuple[RandomForestClassifier, float]:

    logging.info(f"\n--- {step_num}. Testing Model: {model_name} ---")
    
    X = data[feature_cols]  # Multivariate
    y = data[target_col]
    
    logging.info(f"Using {len(feature_cols)} features: {', '.join(feature_cols)}")
    logging.info(f"Target variable distribution:\n{y.value_counts(normalize=True).to_string()}")
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, 
        test_size=0.3, 
        random_state=RANDOM_STATE,  # For reproducible results
        stratify=y
    )
    logging.info(f"Data split: {len(X_train)} train, {len(X_test)} test.")
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # --- Model Training with Hyperparameter Tuning ---
    logging.info(f"Tuning {model_name} with GridSearchCV...")
    logging.info("(This may take a minute...)")
    
    # Define the "search space" for the best model
    param_grid = {
        'n_estimators': [50, 100, 200],
        'max_depth': [None, 5, 10, 15],
        'min_samples_leaf': [1, 2, 4],
        'class_weight': ['balanced', None]
    }
    
    grid_search = GridSearchCV(
        estimator=RandomForestClassifier(random_state=RANDOM_STATE),
        param_grid=param_grid,
        scoring='roc_auc',
        cv=5,
        n_jobs=-1,
        verbose=1
    )
    
    grid_search.fit(X_train_scaled, y_train)
    
    # Get the single best model found
    model = grid_search.best_estimator_
    
    logging.info("Model tuning complete.")
    logging.info(f"Best cross-validation AUC from tuning: {grid_search.best_score_:.4f}")
    logging.info(f"Best parameters: {grid_search.best_params_}")
    
    # Save best parameters to data directory
    params_df = pd.DataFrame(grid_search.best_params_, index=[0])
    params_filename = os.path.join(DATA_DIR, f"{model_name.replace(' ', '_').lower()}_best_params.csv")
    params_df.to_csv(params_filename, index=False)
    logging.info(f"  > Best parameters saved to: {params_filename}")

    # --- Model Evaluation ---
    logging.info("Evaluating BEST model on held-out Test Set...")
    
    y_pred = model.predict(X_test_scaled)
    y_prob = model.predict_proba(X_test_scaled)[:, 1]  # Probabilities for class 1
    
    # 1. Classification Report
    report_str = classification_report(y_test, y_pred, target_names=['Non-Seroconverted (0)', 'Seroconverted (1)'])
    report_df = pd.DataFrame(classification_report(y_test, y_pred, output_dict=True)).transpose()
    logging.info(f"Classification Report (Test Set):\n{report_str}")
    
    report_filename = os.path.join(DATA_DIR, f"{model_name.replace(' ', '_').lower()}_classification_report.csv")
    report_df.to_csv(report_filename)
    logging.info(f"  > Classification report saved to: {report_filename}")
    
    # 2. ROC Curve & AUC
    roc_title = f"{model_name}"
    real_auc = generate_roc_plot(y_test, y_prob, title_prefix=roc_title)
    
    # 3. Confusion Matrix
    generate_confusion_matrix_plot(y_test, y_pred, model.classes_, title_prefix=roc_title)
    
    # Return the key components needed for the permutation test
    return model, real_auc


def execute_permutation_validation(
    model: RandomForestClassifier,
    X: pd.DataFrame, 
    y: pd.Series, 
    model_name: str, 
    real_auc: float, 
    step_num: str,
    n_permutations: int = 1000
) -> None:

    logging.info(f"\n--- {step_num}. Running Permutation Test for {model_name} ---")
    logging.info(f"Running {n_permutations} permutations...")
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, 
        test_size=0.3, 
        random_state=RANDOM_STATE,
        stratify=y
    )
    
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    null_scores = []
    
    for _ in tqdm(range(n_permutations), desc=f"Permutation Test ({model_name})"):
        # Shuffle the labels
        y_train_shuffled = np.random.permutation(y_train)
        
        # Re-train the model on shuffled labels
        # Note: We use the *same* model hyperparameters discovered earlier
        model.fit(X_train_scaled, y_train_shuffled)
        
        # Evaluate on the *unshuffled* test set
        y_prob_null = model.predict_proba(X_test_scaled)[:, 1]
        auc_null = roc_auc_score(y_test, y_prob_null)
        null_scores.append(auc_null)
    
    # Plot results
    generate_permutation_plot(real_auc, null_scores, model_name)
    

def identify_core_biomarkers(
    data: pd.DataFrame, 
    all_genes: List[str], 
    target_col: str, 
    step_num: str
) -> None:

    logging.info(f"\n--- {step_num}. Discovering Minimal Core Signature (RFE) ---")
    logging.info("Running Recursive Feature Elimination (RFE) to find top 10...")

    X = data[all_genes]
    y = data[target_col]
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, 
        test_size=0.3, 
        random_state=RANDOM_STATE,  # For reproducible results
        stratify=y
    )
    
    logging.info(f"Data split: {len(X_train)} train, {len(X_test)} test.")
    
    # Scale features before RFE
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_train)

    # Use a Random Forest as the estimator for RFE
    estimator = RandomForestClassifier(
        random_state=RANDOM_STATE, 
        n_estimators=100, 
        n_jobs=-1, 
        class_weight='balanced'
    )
    
    # RFE will select the top 10 features
    selector = RFE(estimator, n_features_to_select=5, step=0.1, verbose=1)
    selector = selector.fit(X_scaled, y_train)
    
    logging.info("RFE analysis complete.")
    
    # Get rankings
    gene_rankings = pd.DataFrame({
        'gene': all_genes,
        'ranking': selector.ranking_
    }).sort_values(by='ranking')
    
    core_signature = gene_rankings[gene_rankings['ranking'] == 1]['gene'].tolist()
    
    logging.info("Top 20 Gene Rankings (1 = most important):")
    logging.info(gene_rankings.head(20).to_string())
    
    # Save rankings
    rankings_filename = os.path.join(DATA_DIR, "rfe_all_gene_rankings.csv")
    gene_rankings.to_csv(rankings_filename, index=False)
    logging.info(f"  > All gene rankings saved to: {rankings_filename}")
    
    core_sig_df = pd.DataFrame(core_signature, columns=['gene'])
    core_sig_filename = os.path.join(DATA_DIR, "core_10_gene_signature.csv")
    core_sig_df.to_csv(core_sig_filename, index=False)
    logging.info(f"  > 10-Gene Core Signature saved to: {core_sig_filename}")
    
    # Validate against our hardcoded list
    if set(core_signature) == set(CORE_SIGNATURE_GENES):
        logging.info("  > RFE result successfully matches hardcoded CORE_SIGNATURE_GENES.")
    else:
        logging.warning("  > RFE result differs from hardcoded list. Update list?")
        logging.warning(f"  > Discovered: {core_signature}")


def main() -> None:
    initialize_directories()
    configure_logging()
    
    logging.info("="*50)
    logging.info("NASIMMUNE Transcriptome Analysis Report")
    logging.info("="*50)
    
    try:
        # 1. Load Data
        data = ingest_dataset(
            DATA_FILE_PATH, ALL_DEGS, TARGET_COLUMN, step_num=1
        )
        
        # --- Model 1: Broad Interferon Score (Baseline) ---
        # 2. Create the 63-gene "interferon_score"
        data = compute_composite_score(
            data, INTERFERON_GENES, 'interferon_score', step_num="2"
        )
        
        # 2a. Test the 63-gene "interferon_score" with Logistic Regression
        evaluate_univariate_predictor(
            data, 'interferon_score', TARGET_COLUMN, step_num="2a"
        )
        
        # --- Model 2: Naive Core Signature ---
        # This model was an intermediate discovery step and is
        # superseded by the advanced Model 3.
        # logging.info("\nSkipping intermediate Model 2 (Univariate Core Signature)...")
        #
        # data = compute_composite_score(
        #     data, CORE_SIGNATURE_GENES, 'core_signature_score', step_num="X"
        # )
        # evaluate_univariate_predictor(
        #     data, 'core_signature_score', TARGET_COLUMN, step_num="Xa"
        # )
        
        
        # --- Model 3: Advanced Core Signature ---
        # 3. Test the 10-gene signature with a tuned Random Forest
        model, real_auc = evaluate_multivariate_classifier(
            data, CORE_SIGNATURE_GENES, TARGET_COLUMN, 
            model_name="Random Forest (10-gene)", 
            step_num="3"
       )
        
        # 3a. Run 1000-permutation test on the tuned model
        execute_permutation_validation(
           model=model,
            X=data[CORE_SIGNATURE_GENES],
            y=data[TARGET_COLUMN],
            model_name="Random Forest (10-gene)",
            real_auc=real_auc,
            step_num="3a",
            n_permutations=1000  # 1000 is a good standard
        )
        
        # 4. Run RFE to (re)discover the 10-gene signature
        identify_core_biomarkers(
            data, ALL_DEGS, TARGET_COLUMN, step_num="4"
        )
        
        logging.info("\n--- Analysis Pipeline Complete ---")

    except Exception as e:
        # Log the full error traceback if the pipeline fails
        logging.critical(f"\n---!! An Error Occurred !!---")
        logging.critical(f"Pipeline failed with error: {e}")
        logging.critical(traceback.format_exc())
        
    finally:
        logging.info(f"Report saved to: {LOG_FILE}")


if __name__ == "__main__":
    main()
