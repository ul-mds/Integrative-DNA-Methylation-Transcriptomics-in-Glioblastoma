#!/usr/bin/env python3
"""
External Cohort Validation for GBM Biomarkers
Validates 6 biomarkers (B4GALT1, EGFR, VASP, GFPT2, PHF11, TSPYL5) across:
- GSE119774 (methylation)
- GSE4271 (expression)
- CGGA mRNAseq-325 (expression)
- CGGA Methylation

Compares with TCGA-GBM results and saves outputs to Validation folder.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Add parent for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

# Biomarkers to validate
BIOMARKERS = ['B4GALT1', 'EGFR', 'VASP', 'GFPT2', 'PHF11', 'TSPYL5']

# Paths
BASE = Path(__file__).resolve().parent.parent.parent.parent  # GBM root
VALIDATION_DIR = Path(__file__).resolve().parent.parent
EXTERNAL_DIR = VALIDATION_DIR / "external_cohort"
EXTERNAL_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR = VALIDATION_DIR / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR = VALIDATION_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# TCGA reference results (from biomarker_survival_analysis.csv)
TCGA_RESULTS = {
    'B4GALT1': {'median_high': 496, 'median_low': 340, 'diff': 156, 'p': 0.000245, 'direction': 'high_better'},
    'EGFR': {'median_high': 401, 'median_low': 192, 'diff': 209, 'p': 0.000551, 'direction': 'high_better'},
    'VASP': {'median_high': 627, 'median_low': 344.5, 'diff': 282.5, 'p': 0.001687, 'direction': 'high_better'},
    'GFPT2': {'median_high': 402, 'median_low': 218.5, 'diff': 183.5, 'p': 0.002374, 'direction': 'high_better'},
    'PHF11': {'median_high': 577, 'median_low': 344.5, 'diff': 232.5, 'p': 0.006932, 'direction': 'high_better'},
    'TSPYL5': {'median_high': 568, 'median_low': 344.5, 'diff': 223.5, 'p': 0.008472, 'direction': 'high_better'},
}


def run_survival_analysis(gene_values, clinical_df, time_col, event_col, 
                          split_by='median', high_better_expected=True):
    """
    Run Kaplan-Meier survival analysis.
    gene_values: Series of patient_id -> value (methylation or expression)
    clinical_df: DataFrame with time_col, event_col, and patient IDs as index
    split_by: 'median' or 'threshold' (0.7 for methylation)
    Returns dict with results.
    """
    df = clinical_df.copy()
    # Map gene values to samples
    if isinstance(clinical_df.index, pd.MultiIndex) or 'Sample' not in df.columns:
        df['value'] = df.index.map(gene_values)
    else:
        df['value'] = df['Sample'].map(gene_values)
    
    df = df.dropna(subset=['value', time_col, event_col])
    df = df[df[time_col] > 0]
    
    if len(df) < 15:
        return None
    
    if split_by == 'median':
        thresh = df['value'].median()
    else:
        thresh = 0.7  # methylation high/low
    
    df['group'] = np.where(df['value'] > thresh, 'high', 'low')
    
    high = df[df['group'] == 'high']
    low = df[df['group'] == 'low']
    
    if len(high) < 5 or len(low) < 5:
        return None
    
    try:
        lr = logrank_test(high[time_col], low[time_col], high[event_col], low[event_col])
        pval = lr.p_value
    except:
        pval = np.nan
    
    med_high = high[time_col].median()
    med_low = low[time_col].median()
    diff = med_high - med_low
    
    # For methylation: high methylation = better survival (high_better=True)
    # For expression: high methylation = low expression, so low expression = better
    # So for expression we expect low = better (high_better_expected=False)
    concordant = (high_better_expected and diff > 0) or (not high_better_expected and diff < 0)
    
    return {
        'n_high': len(high), 'n_low': len(low),
        'median_high': med_high, 'median_low': med_low,
        'survival_diff': diff, 'logrank_p': pval,
        'concordant_with_TCGA': concordant,
        'significant': pval < 0.05 if not np.isnan(pval) else False
    }


def load_tcga_biomarker_results():
    """Load TCGA biomarker survival results."""
    p = BASE / "Methylation_Transcriptomics" / "results" / "biomarker_survival_analysis.csv"
    if p.exists():
        return pd.read_csv(p)
    return None


def validate_gse119774():
    """Validate using GSE119774 (methylation 850K)."""
    print("\n" + "="*60)
    print("GSE119774 Validation (Methylation 850K)")
    print("="*60)
    
    try:
        import GEOparse
    except ImportError:
        print("GEOparse not installed. Run: pip install GEOparse")
        return None
    
    dest = EXTERNAL_DIR / "GSE119774"
    dest.mkdir(exist_ok=True)
    
    try:
        gse = GEOparse.get_GEO(geo="GSE119774", destdir=str(dest))
    except Exception as e:
        print(f"Could not download GSE119774: {e}")
        # Try loading from cache
        soft = dest / "GSE119774_family.soft.gz"
        if not soft.exists():
            soft = dest / "GSE119774.soft.gz"
        if soft.exists():
            gse = GEOparse.get_GEO(filepath=str(soft))
        else:
            return None
    
    # Get methylation matrix and sample metadata
    all_gsms = list(gse.gsms.keys())
    samples = []
    for gsm_id in all_gsms:
        gsm = gse.gsms[gsm_id]
        # Extract metadata for survival
        meta = gsm.metadata
        samples.append({
            'Sample': gsm_id,
            'title': meta.get('title', [''])[0] if meta.get('title') else '',
            'characteristics': str(meta.get('characteristics_ch1', []))
        })
    
    # Get table data - methylation values
    gpl = None
    if gse.gpls:
        gpl = list(gse.gpls.values())[0]
    
    # Build methylation matrix from GSMs
    tables = []
    for gsm_id in all_gsms[:5]:  # Sample first 5 to get structure
        gsm = gse.gsms[gsm_id]
        if gsm.table is not None and not gsm.table.empty:
            tables.append(gsm.table)
            break
    
    if not tables:
        print("No methylation table data in GSE119774")
        return None
    
    # Get full data - iterate all samples
    probe_col = 'ID_REF' if 'ID_REF' in tables[0].columns else tables[0].columns[0]
    value_col = 'VALUE' if 'VALUE' in tables[0].columns else [c for c in tables[0].columns if c != probe_col][0]
    
    meth_data = {}
    for gsm_id in all_gsms:
        gsm = gse.gsms[gsm_id]
        if gsm.table is not None and not gsm.table.empty:
            if value_col in gsm.table.columns:
                meth_data[gsm_id] = dict(zip(gsm.table[probe_col], pd.to_numeric(gsm.table[value_col], errors='coerce')))
    
    if not meth_data:
        print("No methylation values extracted")
        return None
    
    # Map probes to genes - need GPL
    probe_to_gene = {}
    if gpl and gpl.table is not None:
        if 'UCSC_RefGene_Name' in gpl.table.columns:
            for _, row in gpl.table.iterrows():
                pid = row.get('ID', row.iloc[0])
                genes = str(row.get('UCSC_RefGene_Name', '')).split(';')
                for g in genes:
                    g = g.strip()
                    if g and g in BIOMARKERS:
                        probe_to_gene[pid] = probe_to_gene.get(pid, []) + [g]
    
    # Aggregate to gene level (mean of probes per gene)
    gene_meth = {g: {} for g in BIOMARKERS}
    for sample, probes in meth_data.items():
        for probe, val in probes.items():
            if pd.isna(val):
                continue
            genes = probe_to_gene.get(probe, [])
            for g in genes:
                if g not in gene_meth:
                    continue
                gene_meth[g].setdefault(sample, []).append(val)
    
    for g in BIOMARKERS:
        for s in gene_meth[g]:
            gene_meth[g][s] = np.mean(gene_meth[g][s])
    
    # GSE119774 may not have survival - check metadata
    # Try to extract OS from characteristics
    clinical_list = []
    for gsm_id in all_gsms:
        gsm = gse.gsms[gsm_id]
        ch = gsm.metadata.get('characteristics_ch1', [])
        os_days = np.nan
        os_event = 0
        for c in ch:
            if isinstance(c, str):
                if 'survival' in c.lower() or 'os' in c.lower() or 'days' in c.lower():
                    try:
                        os_days = float(''.join(filter(lambda x: x.isdigit() or x=='.', c)))
                    except:
                        pass
                if 'deceased' in c.lower() or 'dead' in c.lower():
                    os_event = 1
        clinical_list.append({'Sample': gsm_id, 'OS_DAYS': os_days, 'OS_EVENT': os_event})
    
    clinical = pd.DataFrame(clinical_list)
    valid_clinical = clinical.dropna(subset=['OS_DAYS'])
    if len(valid_clinical) < 10:
        print("GSE119774: Insufficient survival data in metadata. Skipping survival validation.")
        return {'dataset': 'GSE119774', 'n_samples': len(all_gsms), 'survival_n': len(valid_clinical), 
                'status': 'no_survival_data', 'biomarker_results': {}}
    
    clinical = valid_clinical.set_index('Sample')
    results = {}
    for gene in BIOMARKERS:
        if gene not in gene_meth or not gene_meth[gene]:
            continue
        vals = pd.Series(gene_meth[gene])
        common = vals.index.intersection(clinical.index)
        if len(common) < 10:
            continue
        res = run_survival_analysis(vals, clinical.loc[common], 'OS_DAYS', 'OS_EVENT', 
                                    split_by='median', high_better_expected=True)
        if res:
            results[gene] = res
    
    return {'dataset': 'GSE119774', 'n_samples': len(all_gsms), 'survival_n': len(clinical),
            'biomarker_results': results}


def validate_gse4271():
    """Validate using GSE4271 (gene expression)."""
    print("\n" + "="*60)
    print("GSE4271 Validation (Gene Expression)")
    print("="*60)
    
    try:
        import GEOparse
    except ImportError:
        return None
    
    dest = EXTERNAL_DIR / "GSE4271"
    dest.mkdir(exist_ok=True)
    
    try:
        gse = GEOparse.get_GEO(geo="GSE4271", destdir=str(dest))
    except Exception as e:
        print(f"Could not download GSE4271: {e}")
        for f in dest.glob("*.soft*"):
            try:
                gse = GEOparse.get_GEO(filepath=str(f))
                break
            except:
                continue
        else:
            return None
    
    # Get expression matrix
    if gse.gsms:
        gsm = list(gse.gsms.values())[0]
        if gsm.table is not None:
            probe_col = gsm.table.columns[0]
        else:
            return None
    else:
        return None
    
    # Build expression matrix
    expr = {}
    for gsm_id, gsm in gse.gsms.items():
        if gsm.table is not None:
            for _, row in gsm.table.iterrows():
                probe = row[probe_col]
                val = row.get('VALUE', row.get('AVG_SIGNAL', row.iloc[-1]))
                try:
                    val = float(val)
                except:
                    continue
                expr.setdefault(probe, {})[gsm_id] = val
    
    # Map probes to genes via GPL
    gpl = list(gse.gpls.values())[0] if gse.gpls else None
    probe_to_gene = {}
    if gpl and gpl.table is not None:
        gene_col = None
        for c in ['Gene Symbol', 'GENE_SYMBOL', 'GeneSymbol', 'SYMBOL']:
            if c in gpl.table.columns:
                gene_col = c
                break
        if gene_col:
            for _, row in gpl.table.iterrows():
                pid = row.iloc[0]
                g = str(row[gene_col]).strip()
                if g and g in BIOMARKERS:
                    probe_to_gene[pid] = g
    
    # Gene-level expression
    gene_expr = {g: {} for g in BIOMARKERS}
    for probe, samples in expr.items():
        g = probe_to_gene.get(probe)
        if not g:
            continue
        for s, v in samples.items():
            gene_expr[g].setdefault(s, []).append(v)
    
    for g in BIOMARKERS:
        for s in list(gene_expr[g].keys()):
            gene_expr[g][s] = np.median(gene_expr[g][s])
    
    # Clinical - extract from characteristics
    clinical_list = []
    for gsm_id in gse.gsms:
        gsm = gse.gsms[gsm_id]
        ch = gsm.metadata.get('characteristics_ch1', [])
        os_days, os_event = np.nan, 0
        for c in ch:
            if isinstance(c, str):
                if 'survival' in c.lower() or 'os' in c.lower() or 'days' in c.lower() or 'follow' in c.lower():
                    try:
                        os_days = float(''.join(filter(lambda x: x.isdigit() or x=='.', c)))
                    except:
                        pass
                if 'deceased' in c.lower() or 'dead' in c.lower():
                    os_event = 1
        clinical_list.append({'Sample': gsm_id, 'OS_DAYS': os_days, 'OS_EVENT': os_event})
    
    clinical = pd.DataFrame(clinical_list).dropna(subset=['OS_DAYS'])
    if len(clinical) < 10:
        return {'dataset': 'GSE4271', 'n_samples': len(gse.gsms), 'survival_n': 0, 'status': 'no_survival', 'biomarker_results': {}}
    
    clinical = clinical.set_index('Sample')
    results = {}
    for gene in BIOMARKERS:
        if not gene_expr.get(gene):
            continue
        vals = pd.Series(gene_expr[gene])
        common = vals.index.intersection(clinical.index)
        if len(common) < 10:
            continue
        # High methylation = low expression, so LOW expression = better survival
        res = run_survival_analysis(vals, clinical.loc[common], 'OS_DAYS', 'OS_EVENT',
                                    split_by='median', high_better_expected=False)
        if res:
            results[gene] = res
    
    return {'dataset': 'GSE4271', 'n_samples': len(gse.gsms), 'survival_n': len(clinical),
            'biomarker_results': results}


def validate_cgga_mrna():
    """Validate using CGGA mRNAseq-325 (if data available)."""
    print("\n" + "="*60)
    print("CGGA mRNAseq-325 Validation")
    print("="*60)
    
    # CGGA data typically in external_cohort/CGGA/
    cgga_dir = EXTERNAL_DIR / "CGGA"
    mrna_file = cgga_dir / "mRNAseq_325" / "mRNAseq_325_expression.csv"
    if not mrna_file.exists():
        mrna_file = cgga_dir / "mRNAseq_325_expression.csv"
    clinical_file = cgga_dir / "mRNAseq_325_clinical.csv"
    if not clinical_file.exists():
        clinical_file = cgga_dir / "CGGA_clinical_mRNAseq325.csv"
    
    if not mrna_file.exists() or not clinical_file.exists():
        print("CGGA mRNAseq-325 data not found. Place files in validation/external_cohort/CGGA/")
        return {'dataset': 'CGGA_mRNAseq325', 'status': 'data_not_found'}
    
    expr = pd.read_csv(mrna_file)
    clinical = pd.read_csv(clinical_file)
    
    # Normalize column names
    if 'Gene_name' in expr.columns or 'gene_name' in expr.columns:
        gene_col = 'Gene_name' if 'Gene_name' in expr.columns else 'gene_name'
    else:
        gene_col = expr.columns[0]
    
    sample_cols = [c for c in expr.columns if c != gene_col and not c.startswith('Unnamed')]
    
    # Match clinical - CGGA uses CGGA_ID or similar
    id_col = None
    for c in ['CGGA_ID', 'ID', 'Sample', 'sample']:
        if c in clinical.columns:
            id_col = c
            break
    if not id_col:
        id_col = clinical.columns[0]
    
    time_col = None
    for c in ['OS', 'Overall_survival', 'survival_days', 'OS_days']:
        if c in clinical.columns:
            time_col = c
            break
    if not time_col:
        time_col = clinical.columns[1]  # Often second column
    
    event_col = None
    for c in ['Censor', 'Status', 'Event', 'Vital_status']:
        if c in clinical.columns:
            event_col = c
            break
    
    clinical = clinical.rename(columns={id_col: 'Sample'})
    if time_col and time_col != 'OS_DAYS':
        clinical = clinical.rename(columns={time_col: 'OS_DAYS'})
    if event_col and event_col != 'OS_EVENT':
        # Censor: 0=event, 1=censored typically
        clinical['OS_EVENT'] = 1 - clinical[event_col] if 'Censor' in event_col else clinical[event_col]
    
    clinical = clinical.dropna(subset=['OS_DAYS']).copy()
    clinical = clinical[clinical['OS_DAYS'] > 0]
    clinical = clinical.set_index('Sample')
    
    expr_long = expr.set_index(gene_col)
    results = {}
    for gene in BIOMARKERS:
        if gene not in expr_long.index:
            continue
        vals = expr_long.loc[gene, sample_cols]
        vals = pd.Series(vals).astype(float)
        common = vals.index.intersection(clinical.index)
        if len(common) < 15:
            continue
        res = run_survival_analysis(vals, clinical.loc[common], 'OS_DAYS', 'OS_EVENT',
                                    split_by='median', high_better_expected=False)
        if res:
            results[gene] = res
    
    return {'dataset': 'CGGA_mRNAseq325', 'n_samples': len(sample_cols), 'survival_n': len(clinical),
            'biomarker_results': results}


def validate_cgga_methylation():
    """Validate using CGGA Methylation (if data available)."""
    print("\n" + "="*60)
    print("CGGA Methylation Validation")
    print("="*60)
    
    cgga_dir = EXTERNAL_DIR / "CGGA"
    meth_file = cgga_dir / "Methylation" / "methylation_beta.csv"
    if not meth_file.exists():
        meth_file = cgga_dir / "methylation_beta.csv"
    clinical_file = cgga_dir / "Methylation" / "clinical.csv"
    if not clinical_file.exists():
        clinical_file = cgga_dir / "methylation_clinical.csv"
    
    if not meth_file.exists() or not clinical_file.exists():
        print("CGGA Methylation data not found.")
        return {'dataset': 'CGGA_Methylation', 'status': 'data_not_found'}
    
    meth = pd.read_csv(meth_file)
    clinical = pd.read_csv(clinical_file)
    
    # Similar processing as mRNA
    probe_col = meth.columns[0]
    sample_cols = [c for c in meth.columns if c != probe_col]
    
    # Map probes to genes - CGGA may have gene names
    if 'Gene' in meth.columns or 'gene' in meth.columns:
        gene_col = 'Gene' if 'Gene' in meth.columns else 'gene'
        meth = meth.rename(columns={gene_col: 'Gene'})
        gene_meth = meth.groupby('Gene')[sample_cols].mean()
    else:
        # Use probe as proxy - would need mapping
        gene_meth = meth.set_index(probe_col)[sample_cols]
    
    id_col = [c for c in clinical.columns if 'ID' in c or 'Sample' in c][0]
    time_col = [c for c in clinical.columns if 'OS' in c or 'survival' in c.lower()][0]
    event_col = [c for c in clinical.columns if 'Censor' in c or 'Status' in c][0]
    
    clinical = clinical.rename(columns={id_col: 'Sample', time_col: 'OS_DAYS'})
    clinical['OS_EVENT'] = 1 - clinical[event_col] if 'Censor' in event_col else clinical[event_col]
    clinical = clinical.dropna(subset=['OS_DAYS'])[clinical['OS_DAYS'] > 0].set_index('Sample')
    
    results = {}
    for gene in BIOMARKERS:
        if gene not in gene_meth.index:
            continue
        vals = gene_meth.loc[gene]
        common = vals.index.intersection(clinical.index)
        if len(common) < 15:
            continue
        res = run_survival_analysis(vals, clinical.loc[common], 'OS_DAYS', 'OS_EVENT',
                                    split_by='threshold', high_better_expected=True)
        if res:
            results[gene] = res
    
    return {'dataset': 'CGGA_Methylation', 'n_samples': len(sample_cols), 'survival_n': len(clinical),
            'biomarker_results': results}


def create_comparison_report(all_results):
    """Create comparison report and plots."""
    rows = []
    for res in all_results:
        if res is None or 'biomarker_results' not in res or not res['biomarker_results']:
            continue
        ds = res['dataset']
        for gene, br in res['biomarker_results'].items():
            tcga = TCGA_RESULTS.get(gene, {})
            rows.append({
                'Dataset': ds,
                'Gene': gene,
                'N_High': br['n_high'], 'N_Low': br['n_low'],
                'Median_High': br['median_high'], 'Median_Low': br['median_low'],
                'Survival_Diff': br['survival_diff'],
                'LogRank_P': br['logrank_p'],
                'Significant': br['significant'],
                'Concordant_TCGA': br['concordant_with_TCGA'],
                'TCGA_P': tcga.get('p'), 'TCGA_Diff': tcga.get('diff')
            })
    
    df = pd.DataFrame(rows)
    if df.empty:
        return None
    
    df.to_csv(RESULTS_DIR / "validation_survival_comparison.csv", index=False)
    
    # Summary table
    summary = []
    for gene in BIOMARKERS:
        sub = df[df['Gene'] == gene]
        n_val = len(sub)
        n_sig = sub['Significant'].sum()
        n_concordant = sub['Concordant_TCGA'].sum()
        summary.append({
            'Gene': gene,
            'TCGA_P': TCGA_RESULTS.get(gene, {}).get('p'),
            'TCGA_Survival_Diff_Days': TCGA_RESULTS.get(gene, {}).get('diff'),
            'N_Validated_Datasets': n_val,
            'N_Significant': int(n_sig),
            'N_Concordant': int(n_concordant)
        })
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(RESULTS_DIR / "validation_summary.csv", index=False)
    
    return df, summary_df


def plot_validation_results(df, summary_df):
    """Create validation comparison plots."""
    if df is None or df.empty:
        return
    
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    
    # Plot 1: Survival difference by dataset
    ax1 = axes[0]
    pivot = df.pivot_table(index='Gene', columns='Dataset', values='Survival_Diff')
    pivot.plot(kind='bar', ax=ax1, width=0.8)
    ax1.axhline(0, color='black', linewidth=0.5)
    ax1.set_ylabel('Survival Difference (High - Low) days')
    ax1.set_title('Biomarker Validation: Survival Difference Across Datasets')
    ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    ax1.tick_params(axis='x', rotation=45)
    
    # Plot 2: -log10(p) by dataset
    ax2 = axes[1]
    df['neg_log_p'] = -np.log10(df['LogRank_P'].replace(0, 1e-300))
    pivot2 = df.pivot_table(index='Gene', columns='Dataset', values='neg_log_p')
    pivot2.plot(kind='bar', ax=ax2, width=0.8)
    ax2.axhline(-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
    ax2.set_ylabel('-log10(p-value)')
    ax2.set_title('Statistical Significance Across Datasets')
    ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "validation_comparison.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    # Summary heatmap
    fig2, ax = plt.subplots(figsize=(8, 5))
    sig_matrix = df.pivot_table(index='Gene', columns='Dataset', values='Significant').fillna(False).astype(int)
    conc_matrix = df.pivot_table(index='Gene', columns='Dataset', values='Concordant_TCGA').fillna(False).astype(int)
    combined = sig_matrix.astype(str) + '/' + conc_matrix.astype(str)
    sns.heatmap(sig_matrix, annot=combined, fmt='', cmap='RdYlGn', ax=ax, cbar_kws={'label': 'Significant'})
    ax.set_title('Validation: Significant (1/0) / Concordant (1/0)')
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "validation_heatmap.png", dpi=150, bbox_inches='tight')
    plt.close()


def main():
    print("="*60)
    print("External Cohort Validation for GBM Biomarkers")
    print("="*60)
    print("Biomarkers:", BIOMARKERS)
    print("Output:", VALIDATION_DIR)
    
    all_results = []
    
    # Run validations
    r1 = validate_gse119774()
    if r1:
        all_results.append(r1)
    
    r2 = validate_gse4271()
    if r2:
        all_results.append(r2)
    
    r3 = validate_cgga_mrna()
    if r3 and r3.get('status') != 'data_not_found':
        all_results.append(r3)
    
    r4 = validate_cgga_methylation()
    if r4 and r4.get('status') != 'data_not_found':
        all_results.append(r4)
    
    # Create report
    out = create_comparison_report(all_results)
    if out:
        df, summary_df = out
        plot_validation_results(df, summary_df)
        print("\n" + "="*60)
        print("VALIDATION SUMMARY")
        print("="*60)
        print(summary_df.to_string(index=False))
        print("\nResults saved to:", RESULTS_DIR)
        print("Plots saved to:", PLOTS_DIR)
    
    return all_results


if __name__ == "__main__":
    main()
