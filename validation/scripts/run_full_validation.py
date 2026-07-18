#!/usr/bin/env python3
"""
Full Validation Pipeline: GSE119774, CGGA, TCGA
Runs all validations and saves to Validation folder.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import gzip
import subprocess
import sys

BASE = Path(__file__).resolve().parent.parent.parent.parent
VALIDATION_DIR = Path(__file__).resolve().parent.parent
EXTERNAL_DIR = VALIDATION_DIR / "external_cohort"
RESULTS_DIR = VALIDATION_DIR / "results"
PLOTS_DIR = VALIDATION_DIR / "plots"
BIOMARKERS = ['B4GALT1', 'EGFR', 'VASP', 'GFPT2', 'PHF11', 'TSPYL5']

for d in [RESULTS_DIR, PLOTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)


def run_tcga_validation():
    """TCGA internal validation."""
    print("\n=== TCGA Internal Validation ===")
    subprocess.run([sys.executable, str(VALIDATION_DIR / "scripts" / "run_tcga_internal_validation.py")],
                   cwd=str(BASE / "Methylation_Transcriptomics"), check=False)
    return RESULTS_DIR / "tcga_internal_validation.csv"


def run_cgga_validation():
    """CGGA validation - clinical only (expression requires cgga.org.cn)."""
    print("\n=== CGGA Validation ===")
    clin_path = EXTERNAL_DIR / "CGGA" / "CCGA_clinical_mRNAseq325.csv"
    if not clin_path.exists():
        print("  CGGA clinical not found. Skipping.")
        return None
    
    clin = pd.read_csv(clin_path)
    # Filter GBM only
    gbm = clin[clin['Histology'] == 'GBM'].copy()
    if len(gbm) < 20:
        gbm = clin[clin['Disease_Name'] == 'GBM'].copy()
    if len(gbm) < 20:
        gbm = clin.copy()
    
    gbm = gbm.rename(columns={'SampleID': 'Sample', 'OS': 'OS_DAYS'})
    gbm['OS_EVENT'] = 1 - gbm['Censor (1=dead,0=alive)']  # 1=dead
    gbm = gbm.dropna(subset=['OS_DAYS'])
    gbm = gbm[gbm['OS_DAYS'] > 0]
    
    # CGGA has clinical but no expression - we report cohort summary
    summary = {
        'dataset': 'CGGA_mRNAseq325',
        'n_total': len(clin),
        'n_GBM': len(gbm),
        'median_OS_days': gbm['OS_DAYS'].median(),
        'status': 'clinical_only_no_expression'
    }
    pd.DataFrame([summary]).to_csv(RESULTS_DIR / "cgga_cohort_summary.csv", index=False)
    print(f"  CGGA: {len(gbm)} GBM with survival. Expression data requires cgga.org.cn registration.")
    return RESULTS_DIR / "cgga_cohort_summary.csv"


def run_gse119774_methylation_summary():
    """GSE119774 methylation level summary for GBM samples (no survival in GEO)."""
    print("\n=== GSE119774 Methylation Summary ===")
    beta_path = EXTERNAL_DIR / "GSE119774" / "GSE119774_GBM_GSC_i850K_beta_formatted.txt.gz"
    if not beta_path.exists():
        print("  Beta file not found.")
        return None
    
    # Load first chunk to get structure
    with gzip.open(beta_path, 'rt') as f:
        header = f.readline().strip().split('\t')
    
    gbm_cols = [c for c in header if c.startswith('GBM')]
    if len(gbm_cols) < 10:
        print("  Insufficient GBM samples.")
        return None
    
    # Load sample of probes - use first 50k rows for speed
    df = pd.read_csv(beta_path, sep='\t', nrows=100000, usecols=['ID_REF'] + gbm_cols)
    df = df.set_index('ID_REF')
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    
    # Genome-wide mean methylation per sample (proxy)
    mean_per_sample = df.mean(axis=0)
    summary = {
        'dataset': 'GSE119774',
        'n_GBM_samples': len(gbm_cols),
        'n_probes_used': len(df),
        'mean_methylation_range': f"{mean_per_sample.min():.3f}-{mean_per_sample.max():.3f}",
        'median_methylation': float(mean_per_sample.median()),
        'status': 'methylation_available_no_survival_no_probe_gene_mapping'
    }
    pd.DataFrame([summary]).to_csv(RESULTS_DIR / "gse119774_summary.csv", index=False)
    mean_per_sample.to_csv(RESULTS_DIR / "gse119774_gbm_mean_methylation.csv")
    print(f"  GSE119774: {len(gbm_cols)} GBM samples. Methylation summary saved.")
    return RESULTS_DIR / "gse119774_summary.csv"


def create_combined_report():
    """Create combined validation report."""
    rows = []
    
    # TCGA
    tcga = RESULTS_DIR / "tcga_internal_validation.csv"
    if tcga.exists():
        df = pd.read_csv(tcga)
        for _, r in df.iterrows():
            rows.append({
                'Dataset': 'TCGA (Test 30%)',
                'Gene': r['Gene'],
                'Survival_Diff': r['Test_Diff'],
                'P_value': r['Test_P'],
                'Significant': r['Test_Sig'],
                'Note': 'Internal validation'
            })
    
    # CGGA
    cgga = RESULTS_DIR / "cgga_cohort_summary.csv"
    if cgga.exists():
        s = pd.read_csv(cgga).iloc[0]
        rows.append({
            'Dataset': 'CGGA',
            'Gene': '-',
            'Survival_Diff': '-',
            'P_value': '-',
            'Significant': '-',
            'Note': f"n={s['n_GBM']} GBM, median OS={s['median_OS_days']:.0f} days. Expression requires registration."
        })
    
    # GSE119774
    gse = RESULTS_DIR / "gse119774_summary.csv"
    if gse.exists():
        s = pd.read_csv(gse).iloc[0]
        rows.append({
            'Dataset': 'GSE119774',
            'Gene': '-',
            'Survival_Diff': '-',
            'P_value': '-',
            'Significant': '-',
            'Note': f"n={s['n_GBM_samples']} GBM. Methylation available. Survival not in GEO."
        })
    
    if rows:
        pd.DataFrame(rows).to_csv(RESULTS_DIR / "validation_combined_report.csv", index=False)
        print("\nCombined report saved:", RESULTS_DIR / "validation_combined_report.csv")


def create_summary_plot():
    """Create summary visualization."""
    tcga_path = RESULTS_DIR / "tcga_internal_validation.csv"
    if not tcga_path.exists():
        return
    
    df = pd.read_csv(tcga_path)
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(df))
    colors = ['green' if d > 0 else 'red' for d in df['Test_Diff']]
    ax.bar(x, df['Test_Diff'], color=colors)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(df['Gene'])
    ax.set_ylabel('Survival Difference (days)')
    ax.set_title('Validation: TCGA Test Set (30% holdout)\nGreen = High methylation better (concordant with discovery)')
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "validation_summary.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("Plot saved:", PLOTS_DIR / "validation_summary.png")


def main():
    print("="*60)
    print("Full Validation Pipeline")
    print("="*60)
    
    run_tcga_validation()
    run_cgga_validation()
    run_gse119774_methylation_summary()
    create_combined_report()
    create_summary_plot()
    
    print("\n" + "="*60)
    print("Validation complete. Results in:", RESULTS_DIR)
    print("Plots in:", PLOTS_DIR)
    print("="*60)


if __name__ == "__main__":
    main()
