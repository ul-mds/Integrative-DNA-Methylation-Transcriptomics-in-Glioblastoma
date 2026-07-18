#!/usr/bin/env python3
"""
External Cohort Validation for GBM Biomarkers
Validates 6 biomarkers (B4GALT1, EGFR, VASP, GFPT2, PHF11, TSPYL5) across:
- GSE119774 (methylation)
- GSE4271 (expression)
- CGGA mRNAseq-325 (expression) - when data provided
- CGGA Methylation - when data provided

Compares with TCGA-GBM results. Saves to Validation folder.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import gzip
import re
import warnings
warnings.filterwarnings('ignore')

try:
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test
except ImportError:
    logrank_test = None

BIOMARKERS = ['B4GALT1', 'EGFR', 'VASP', 'GFPT2', 'PHF11', 'TSPYL5']

BASE = Path(__file__).resolve().parent.parent.parent.parent
VALIDATION_DIR = Path(__file__).resolve().parent.parent
EXTERNAL_DIR = VALIDATION_DIR / "external_cohort"
PLOTS_DIR = VALIDATION_DIR / "plots"
RESULTS_DIR = VALIDATION_DIR / "results"

for d in [EXTERNAL_DIR, PLOTS_DIR, RESULTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# TCGA reference
TCGA_RESULTS = {
    'B4GALT1': {'median_high': 496, 'median_low': 340, 'diff': 156, 'p': 0.000245},
    'EGFR': {'median_high': 401, 'median_low': 192, 'diff': 209, 'p': 0.000551},
    'VASP': {'median_high': 627, 'median_low': 344.5, 'diff': 282.5, 'p': 0.001687},
    'GFPT2': {'median_high': 402, 'median_low': 218.5, 'diff': 183.5, 'p': 0.002374},
    'PHF11': {'median_high': 577, 'median_low': 344.5, 'diff': 232.5, 'p': 0.006932},
    'TSPYL5': {'median_high': 568, 'median_low': 344.5, 'diff': 223.5, 'p': 0.008472},
}


def parse_beta_file(filepath):
    """Parse GSE119774 beta formatted file (probes x samples)."""
    opener = gzip.open if str(filepath).endswith('.gz') else open
    df = pd.read_csv(opener(filepath), sep='\t', nrows=100000, low_memory=False)
    df = df.set_index(df.columns[0])
    for c in df.columns:
        try:
            df[c] = pd.to_numeric(df[c], errors='coerce')
        except:
            pass
    return df, {}


def parse_geo_series_matrix(filepath):
    """Parse GEO series matrix file. Returns (data_df, sample_metadata)."""
    opener = gzip.open if str(filepath).endswith('.gz') else open
    with opener(filepath, 'rt', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    in_table = False
    headers = []
    table_rows = []
    sample_meta = {}
    
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('!series_matrix_table_begin'):
            in_table = True
            continue
        if line.startswith('!series_matrix_table_end'):
            break
        if in_table:
            parts = line.split('\t')
            if not table_rows:
                headers = [p.strip('"') for p in parts]
            table_rows.append([p.strip('"') for p in parts])
        else:
            if line.startswith('!Sample_') and '=' in line:
                m = re.match(r'!(Sample_\w+)\s*=\s*(.+)', line)
                if m:
                    key, val = m.group(1), m.group(2).strip()
                    sample_meta.setdefault(key, []).append(val)
    
    if not table_rows:
        return None, {}
    
    df = pd.DataFrame(table_rows[1:], columns=headers)
    df = df.set_index(df.columns[0])
    for c in df.columns:
        try:
            df[c] = pd.to_numeric(df[c], errors='coerce')
        except:
            pass
    
    return df, sample_meta


def extract_survival_from_geo_metadata(sample_meta, gsm_ids):
    """Extract survival from GEO sample characteristics."""
    # Sample_characteristics_ch1 often has "survival: 365" or "os: 400"
    char_key = 'Sample_characteristics_ch1'
    if char_key not in sample_meta:
        for k in sample_meta:
            if 'char' in k.lower():
                char_key = k
                break
    
    rows = []
    for i, gsm in enumerate(gsm_ids):
        chars = sample_meta.get(char_key, [[]])[i] if i < len(sample_meta.get(char_key, [])) else ''
        if isinstance(chars, list):
            chars = ';'.join(chars)
        os_days, os_event = np.nan, 0
        for part in str(chars).split(';'):
            part = part.strip().lower()
            if 'survival' in part or 'os' in part or 'days' in part or 'follow' in part:
                nums = re.findall(r'[\d.]+', part)
                if nums:
                    os_days = float(nums[0])
            if 'deceased' in part or 'dead' in part or 'death' in part:
                os_event = 1
            if 'alive' in part or 'living' in part or 'censored' in part:
                os_event = 0
        rows.append({'Sample': gsm, 'OS_DAYS': os_days, 'OS_EVENT': os_event})
    
    return pd.DataFrame(rows).dropna(subset=['OS_DAYS'])


def run_survival(gene_vals, clinical, time_col='OS_DAYS', event_col='OS_EVENT', 
                 high_better=True):
    """Kaplan-Meier + log-rank. gene_vals: Series index=sample."""
    if logrank_test is None:
        return None
    df = clinical.copy()
    df['val'] = df['Sample'].map(gene_vals)
    df = df.dropna(subset=['val', time_col, event_col])
    df = df[df[time_col] > 0]
    if len(df) < 15:
        return None
    thresh = df['val'].median()
    high = df[df['val'] > thresh]
    low = df[df['val'] <= thresh]
    if len(high) < 5 or len(low) < 5:
        return None
    try:
        lr = logrank_test(high[time_col], low[time_col], high[event_col], low[event_col])
        pval = lr.p_value
    except:
        pval = np.nan
    diff = high[time_col].median() - low[time_col].median()
    concordant = (high_better and diff > 0) or (not high_better and diff < 0)
    return {
        'n_high': len(high), 'n_low': len(low),
        'median_high': high[time_col].median(), 'median_low': low[time_col].median(),
        'survival_diff': diff, 'logrank_p': pval,
        'concordant_with_TCGA': concordant,
        'significant': pval < 0.05 if not np.isnan(pval) else False
    }


def load_probe_to_gene_mapping(gpl_file):
    """Load probe->gene from GPL annotation."""
    mapping = {}
    if not gpl_file or not Path(gpl_file).exists():
        return mapping
    try:
        df = pd.read_csv(gpl_file, sep='\t', low_memory=False, nrows=50000)
        gene_col = None
        for c in ['Gene Symbol', 'GENE_SYMBOL', 'GeneSymbol', 'SYMBOL', 'UCSC_RefGene_Name']:
            if c in df.columns:
                gene_col = c
                break
        if gene_col:
            for _, row in df.iterrows():
                pid = str(row.iloc[0])
                genes = str(row.get(gene_col, '')).replace(';', ',').split(',')
                for g in genes:
                    g = g.strip()
                    if g in BIOMARKERS:
                        mapping.setdefault(pid, []).append(g)
    except Exception as e:
        print(f"GPL parse warning: {e}")
    return mapping


def validate_gse119774():
    """GSE119774 - Methylation 850K."""
    print("\n--- GSE119774 (Methylation) ---")
    # Prefer beta file (has actual methylation values)
    beta_f = EXTERNAL_DIR / "GSE119774" / "GSE119774_GBM_GSC_i850K_beta_formatted.txt.gz"
    series_f = EXTERNAL_DIR / "GSE119774" / "GSE119774_series_matrix.txt.gz"
    f = beta_f if beta_f.exists() else series_f
    if not f.exists():
        f = EXTERNAL_DIR / "GSE119774" / "GSE119774_series_matrix.txt"
    if not f.exists():
        print("  File not found. Run download_geo_data.py or place manually.")
        return None
    
    if 'beta' in str(f):
        df, meta = parse_beta_file(f)
    else:
        df, meta = parse_geo_series_matrix(f)
    if df is None or df.empty:
        return None
    
    # GEO methylation: rows=probes, columns=samples. Probe IDs like cg00000029
    # Need GPL for probe->gene. GPL23976 for EPIC 850K
    gpl = EXTERNAL_DIR / "GSE119774" / "GPL23976.txt"
    probe_to_gene = load_probe_to_gene_mapping(gpl)
    
    samples = [c for c in df.columns if c.startswith('GSM')]
    clinical = extract_survival_from_geo_metadata(meta, samples)
    if len(clinical) < 10:
        return {'dataset': 'GSE119774', 'n': len(samples), 'survival_n': len(clinical),
                'status': 'insufficient_survival', 'biomarker_results': {}}
    
    gene_data = {g: {} for g in BIOMARKERS}
    for probe in df.index:
        genes = probe_to_gene.get(probe, [])
        for g in genes:
            for s in samples:
                v = df.loc[probe, s]
                if pd.notna(v) and np.isfinite(float(v) if isinstance(v, (int, float)) else 0):
                    gene_data[g].setdefault(s, []).append(float(v))
    
    for g in BIOMARKERS:
        for s in gene_data[g]:
            gene_data[g][s] = np.mean(gene_data[g][s])
    
    results = {}
    for g in BIOMARKERS:
        if not gene_data[g]:
            continue
        vals = pd.Series(gene_data[g])
        common = list(set(vals.index) & set(clinical['Sample']))
        if len(common) < 10:
            continue
        cl = clinical[clinical['Sample'].isin(common)]
        res = run_survival(vals, cl, high_better=True)
        if res:
            results[g] = res
    
    return {'dataset': 'GSE119774', 'n': len(samples), 'survival_n': len(clinical),
            'biomarker_results': results}


def validate_gse4271():
    """GSE4271 - Gene expression."""
    print("\n--- GSE4271 (Expression) ---")
    f = EXTERNAL_DIR / "GSE4271" / "GSE4271_series_matrix.txt.gz"
    if not f.exists():
        f = EXTERNAL_DIR / "GSE4271" / "GSE4271_series_matrix.txt"
    if not f.exists():
        print("  File not found.")
        return None
    
    df, meta = parse_geo_series_matrix(f)
    if df is None or df.empty:
        return None
    
    samples = [c for c in df.columns if c.startswith('GSM')]
    clinical = extract_survival_from_geo_metadata(meta, samples)
    if len(clinical) < 10:
        return {'dataset': 'GSE4271', 'n': len(samples), 'survival_n': len(clinical),
                'status': 'insufficient_survival', 'biomarker_results': {}}
    
    # GSE4271 may use probe IDs. Try GPL
    gpl = EXTERNAL_DIR / "GSE4271" / "GPL96.txt"
    if not gpl.exists():
        gpl = EXTERNAL_DIR / "GSE4271" / "GPL97.txt"
    probe_to_gene = load_probe_to_gene_mapping(gpl)
    
    results = {}
    for g in BIOMARKERS:
        if probe_to_gene:
            probes = [p for p, genes in probe_to_gene.items() if g in genes]
        else:
            probes = [g] if g in df.index else []
        if not probes:
            continue
        sub = df.loc[df.index.isin(probes), samples].mean(axis=0)
        vals = sub.astype(float)
        common = list(set(vals.index) & set(clinical['Sample']))
        if len(common) < 10:
            continue
        cl = clinical[clinical['Sample'].isin(common)]
        res = run_survival(vals, cl, high_better=False)  # low expr = better (high meth)
        if res:
            results[g] = res
    
    return {'dataset': 'GSE4271', 'n': len(samples), 'survival_n': len(clinical),
            'biomarker_results': results}


def validate_cgga():
    """CGGA mRNAseq-325 and Methylation."""
    print("\n--- CGGA ---")
    cgga = EXTERNAL_DIR / "CGGA"
    
    # mRNAseq-325
    mrna_files = list((cgga / "mRNAseq_325").glob("*.csv")) if (cgga / "mRNAseq_325").exists() else []
    mrna_files += list(cgga.glob("*mRNA*325*.csv"))
    clinical_files = list(cgga.glob("*clinical*.csv")) + list(cgga.glob("*Clinical*.csv"))
    
    all_res = []
    if mrna_files and clinical_files:
        expr = pd.read_csv(mrna_files[0])
        clin = pd.read_csv(clinical_files[0])
        gene_col = [c for c in expr.columns if 'gene' in c.lower() or c == 'Gene_name'][0]
        sample_cols = [c for c in expr.columns if c != gene_col and not str(c).startswith('Unnamed')]
        id_col = [c for c in clin.columns if 'ID' in c or 'Sample' in c][0]
        time_col = [c for c in clin.columns if 'OS' in c or 'survival' in c.lower()][0]
        event_col = [c for c in clin.columns if 'Censor' in c or 'Status' in c][0]
        clin = clin.rename(columns={id_col: 'Sample', time_col: 'OS_DAYS'})
        clin['OS_EVENT'] = 1 - clin[event_col] if 'Censor' in event_col else clin[event_col]
        clin = clin.dropna(subset=['OS_DAYS'])[clin['OS_DAYS'] > 0]
        expr = expr.set_index(gene_col)
        results = {}
        for g in BIOMARKERS:
            if g not in expr.index:
                continue
            vals = expr.loc[g, sample_cols].astype(float)
            common = list(set(vals.index) & set(clin['Sample']))
            if len(common) < 15:
                continue
            cl = clin[clin['Sample'].isin(common)]
            res = run_survival(vals, cl, high_better=False)
            if res:
                results[g] = res
        all_res.append({'dataset': 'CGGA_mRNAseq325', 'n': len(sample_cols), 'survival_n': len(clin),
                        'biomarker_results': results})
    
    if not all_res:
        return {'dataset': 'CGGA', 'status': 'data_not_found'}
    return all_res[0]


def build_report(all_results):
    rows = []
    for res in all_results:
        if res is None or 'biomarker_results' not in res or not res.get('biomarker_results'):
            continue
        ds = res['dataset']
        for gene, br in res['biomarker_results'].items():
            tcga = TCGA_RESULTS.get(gene, {})
            rows.append({
                'Dataset': ds, 'Gene': gene,
                'N_High': br['n_high'], 'N_Low': br['n_low'],
                'Median_High': br['median_high'], 'Median_Low': br['median_low'],
                'Survival_Diff': br['survival_diff'], 'LogRank_P': br['logrank_p'],
                'Significant': br['significant'], 'Concordant_TCGA': br['concordant_with_TCGA'],
                'TCGA_P': tcga.get('p'), 'TCGA_Diff': tcga.get('diff')
            })
    
    if not rows:
        return None
    
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "validation_survival_comparison.csv", index=False)
    
    summary = []
    for gene in BIOMARKERS:
        sub = df[df['Gene'] == gene]
        summary.append({
            'Gene': gene, 'TCGA_P': TCGA_RESULTS.get(gene, {}).get('p'),
            'TCGA_Survival_Diff_Days': TCGA_RESULTS.get(gene, {}).get('diff'),
            'N_Validated_Datasets': len(sub), 'N_Significant': int(sub['Significant'].sum()),
            'N_Concordant': int(sub['Concordant_TCGA'].sum())
        })
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(RESULTS_DIR / "validation_summary.csv", index=False)
    
    return df, summary_df


def plot_results(df, summary_df):
    if df is None or df.empty:
        return
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    
    pivot = df.pivot_table(index='Gene', columns='Dataset', values='Survival_Diff')
    pivot.plot(kind='bar', ax=axes[0], width=0.8)
    axes[0].axhline(0, color='black', linewidth=0.5)
    axes[0].set_ylabel('Survival Difference (days)')
    axes[0].set_title('Biomarker Validation: Survival Difference Across Datasets')
    axes[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    axes[0].tick_params(axis='x', rotation=45)
    
    df['neg_log_p'] = -np.log10(df['LogRank_P'].replace(0, 1e-300))
    pivot2 = df.pivot_table(index='Gene', columns='Dataset', values='neg_log_p')
    pivot2.plot(kind='bar', ax=axes[1], width=0.8)
    axes[1].axhline(-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
    axes[1].set_ylabel('-log10(p-value)')
    axes[1].set_title('Statistical Significance')
    axes[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    axes[1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "validation_comparison.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    # Heatmap
    fig2, ax = plt.subplots(figsize=(8, 5))
    sig = df.pivot_table(index='Gene', columns='Dataset', values='Significant').fillna(False).astype(int)
    conc = df.pivot_table(index='Gene', columns='Dataset', values='Concordant_TCGA').fillna(False).astype(int)
    annot = sig.astype(str) + '/' + conc.astype(str)
    sns.heatmap(sig, annot=annot, fmt='', cmap='RdYlGn', ax=ax)
    ax.set_title('Validation: Significant/Concordant with TCGA')
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "validation_heatmap.png", dpi=150, bbox_inches='tight')
    plt.close()


def main():
    print("="*60)
    print("External Cohort Validation for GBM Biomarkers")
    print("="*60)
    
    all_results = []
    for fn in [validate_gse119774, validate_gse4271, validate_cgga]:
        try:
            r = fn()
            if r:
                all_results.append(r)
        except Exception as e:
            print(f"  Error: {e}")
    
    out = build_report(all_results)
    if out:
        df, summary_df = out
        plot_results(df, summary_df)
        print("\n" + "="*60)
        print("VALIDATION SUMMARY")
        print("="*60)
        print(summary_df.to_string(index=False))
        print("\nResults:", RESULTS_DIR)
        print("Plots:", PLOTS_DIR)
    
    return all_results


if __name__ == "__main__":
    main()
