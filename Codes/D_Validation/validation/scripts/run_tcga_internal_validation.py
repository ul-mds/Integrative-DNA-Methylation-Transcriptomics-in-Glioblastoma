#!/usr/bin/env python3
"""
TCGA Internal Validation: Train/Test Split
Validates 6 biomarkers on held-out TCGA samples.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

BIOMARKERS = ['B4GALT1', 'EGFR', 'VASP', 'GFPT2', 'PHF11', 'TSPYL5']
BASE = Path(__file__).resolve().parent.parent.parent.parent
METH_FILE = BASE / "methylation_epic.txt"
CLINICAL_FILE = BASE / "clinical_patient.txt"
VALIDATION_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = VALIDATION_DIR / "results"
PLOTS_DIR = VALIDATION_DIR / "plots"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)


def load_tcga_data():
    """Load methylation and clinical from TCGA."""
    clinical = pd.read_csv(CLINICAL_FILE, sep="\t", skiprows=4)
    clinical.columns = clinical.columns.str.strip()
    clinical['PATIENT_ID'] = clinical['PATIENT_ID'].str.strip()
    clinical['OS_EVENT'] = (clinical['VITAL_STATUS'].str.strip().str.lower() == 'deceased').astype(int)
    death = pd.to_numeric(clinical['PATH_DIAG_TO_DEATH_DAYS'], errors='coerce')
    contact = pd.to_numeric(clinical['PATH_DIAG_TO_LAST_CONTACT_DAYS'], errors='coerce')
    clinical['OS_DAYS'] = np.where(clinical['OS_EVENT'] == 1, death, contact)
    clinical = clinical.dropna(subset=['OS_DAYS', 'OS_EVENT'])
    clinical = clinical[clinical['OS_DAYS'] > 0]
    
    # Load methylation for biomarker genes only
    load_cols = ['Hugo_Symbol'] + [c for c in pd.read_csv(METH_FILE, sep="\t", nrows=0).columns if c.startswith('C3')]
    chunks = []
    for chunk in pd.read_csv(METH_FILE, sep="\t", usecols=lambda c: c in load_cols, chunksize=50000):
        chunks.append(chunk)
    meth = pd.concat(chunks, ignore_index=True)
    gene_means = meth.groupby("Hugo_Symbol")[load_cols[1:]].mean()
    
    return clinical, gene_means


def run_survival(meth_df, clinical, gene, threshold=0.7):
    """Kaplan-Meier by methylation high/low. meth_df: genes x patients."""
    if gene not in meth_df.index:
        return None
    vals = meth_df.loc[gene]  # Series: patient_id -> methylation value
    common = [p for p in vals.index if p in clinical['PATIENT_ID'].values]
    if len(common) < 20:
        return None
    cl = clinical[clinical['PATIENT_ID'].isin(common)].copy()
    cl['meth'] = cl['PATIENT_ID'].map(vals)
    cl = cl.dropna(subset=['meth'])
    high = cl[cl['meth'] > threshold]
    low = cl[cl['meth'] <= threshold]
    if len(high) < 5 or len(low) < 5:
        return None
    try:
        lr = logrank_test(high['OS_DAYS'], low['OS_DAYS'], high['OS_EVENT'], low['OS_EVENT'])
        pval = lr.p_value
    except:
        pval = np.nan
    return {
        'n_high': len(high), 'n_low': len(low),
        'median_high': high['OS_DAYS'].median(), 'median_low': low['OS_DAYS'].median(),
        'survival_diff': high['OS_DAYS'].median() - low['OS_DAYS'].median(),
        'logrank_p': pval, 'significant': pval < 0.05 if not np.isnan(pval) else False
    }


def main():
    print("TCGA Internal Validation (Train/Test Split)")
    clinical, meth = load_tcga_data()
    patients = [p for p in meth.columns if p in clinical['PATIENT_ID'].values]
    np.random.seed(42)
    np.random.shuffle(patients)
    n_train = int(0.7 * len(patients))
    train_patients, test_patients = set(patients[:n_train]), set(patients[n_train:])
    
    cl_train = clinical[clinical['PATIENT_ID'].isin(train_patients)]
    cl_test = clinical[clinical['PATIENT_ID'].isin(test_patients)]
    meth_train = meth[[c for c in meth.columns if c in train_patients]]
    meth_test = meth[[c for c in meth.columns if c in test_patients]]
    
    print(f"Train: {len(train_patients)}, Test: {len(test_patients)}")
    
    results = []
    for gene in BIOMARKERS:
        r_train = run_survival(meth_train, cl_train, gene)
        r_test = run_survival(meth_test, cl_test, gene)
        if r_train and r_test:
            results.append({
                'Gene': gene,
                'Train_P': r_train['logrank_p'], 'Train_Sig': r_train['significant'],
                'Train_Diff': r_train['survival_diff'],
                'Test_P': r_test['logrank_p'], 'Test_Sig': r_test['significant'],
                'Test_Diff': r_test['survival_diff'],
                'Concordant': (r_train['survival_diff'] > 0) == (r_test['survival_diff'] > 0)
            })
    
    df = pd.DataFrame(results)
    if df.empty:
        print("No results - insufficient data in test set")
        return None
    df.to_csv(RESULTS_DIR / "tcga_internal_validation.csv", index=False)
    
    # Plot
    genes = df['Gene'].tolist()
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(genes))
    w = 0.35
    ax.bar(x - w/2, df['Train_Diff'], w, label='Train (70%)')
    ax.bar(x + w/2, df['Test_Diff'], w, label='Test (30%)')
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(genes)
    ax.set_ylabel('Survival Difference (High - Low Methylation) days')
    ax.set_title('TCGA Internal Validation: Train vs Test Set')
    ax.legend()
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "tcga_internal_validation.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    print("\nResults:")
    print(df.to_string(index=False))
    print("\nSaved to:", RESULTS_DIR, PLOTS_DIR)
    return df


if __name__ == "__main__":
    main()
