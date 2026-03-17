#!/usr/bin/env python3
"""
External Dataset Validation: GSE22867, GSE230770, GSE7696, GSE108474
Validates 6 biomarkers (B4GALT1, EGFR, VASP, GFPT2, PHF11, TSPYL5) against TCGA results.
Saves to Validation folder.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import gzip
import re
import warnings
warnings.filterwarnings('ignore')

try:
    from lifelines.statistics import logrank_test
except ImportError:
    logrank_test = None

BIOMARKERS = ['B4GALT1', 'EGFR', 'VASP', 'GFPT2', 'PHF11', 'TSPYL5']
BASE = Path(__file__).resolve().parent.parent
VALIDATION_DIR = BASE
EXTERNAL_DIR = VALIDATION_DIR / "external_cohort"
RESULTS_DIR = VALIDATION_DIR / "results"
PLOTS_DIR = VALIDATION_DIR / "plots"

# Affymetrix HG-U133 Plus 2.0 (GPL570) probe IDs
AFFY_PROBE_TO_GENE = {
    '222163_s_at': 'B4GALT1', '201983_s_at': 'EGFR', '211607_x_at': 'EGFR',
    '201651_s_at': 'VASP', '218628_at': 'GFPT2', '221731_x_at': 'PHF11',
    '221840_at': 'TSPYL5'
}

# Illumina EPIC (GPL21145) - probes verified in GSE230770
EPIC_PROBE_TO_GENE = {
    'cg18809074': 'EGFR', 'cg03547355': 'VASP', 'cg07814318': 'TSPYL5',
}

# Illumina 27K (GPL8490) - limited probes
ILLUMINA27K_PROBE_TO_GENE = {
    'cg12087742': 'B4GALT1', 'cg10673833': 'EGFR', 'cg16875085': 'B4GALT1',
}

TCGA_RESULTS = {
    'B4GALT1': {'diff': 156, 'p': 0.000245, 'direction': 'high_meth_better'},
    'EGFR': {'diff': 209, 'p': 0.000551, 'direction': 'high_meth_better'},
    'VASP': {'diff': 282.5, 'p': 0.001687, 'direction': 'high_meth_better'},
    'GFPT2': {'diff': 183.5, 'p': 0.002374, 'direction': 'high_meth_better'},
    'PHF11': {'diff': 232.5, 'p': 0.006932, 'direction': 'high_meth_better'},
    'TSPYL5': {'diff': 223.5, 'p': 0.008472, 'direction': 'high_meth_better'},
}


def parse_geo_matrix(filepath, probe_subset=None):
    """Parse GEO series matrix. Returns (data_df, sample_ids, meta_rows).
    If probe_subset is set, only load those probe rows (faster for large matrices)."""
    opener = gzip.open if str(filepath).endswith('.gz') else open
    in_table = False
    table_rows = []
    meta_rows = []
    headers = None

    with opener(filepath, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('!series_matrix_table_begin'):
                in_table = True
                continue
            if line.startswith('!series_matrix_table_end'):
                break
            if in_table:
                parts = [p.strip('"') for p in line.split('\t')]
                if not table_rows:
                    headers = parts
                    table_rows.append(parts)
                    continue
                pid = parts[0]
                if probe_subset is not None and pid not in probe_subset:
                    continue
                table_rows.append(parts)
            elif line.startswith('!'):
                parts = [p.strip('"') for p in line.split('\t')]
                if len(parts) >= 2:
                    key = parts[0].lstrip('!')
                    meta_rows.append((key, parts[1:]))

    if not table_rows or not headers:
        return None, [], []

    df = pd.DataFrame(table_rows[1:], columns=headers).set_index(headers[0])
    sample_ids = [c for c in df.columns if c.startswith('GSM')]
    for c in sample_ids:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df, sample_ids, meta_rows


def extract_survival_from_meta(meta_rows, sample_ids):
    """
    Extract OS/survival from GEO metadata.
    meta_rows: list of (key, [val1, val2, ...]) where values align with sample_ids.
    For Sample_characteristics_ch1, each row has format "field: value" per sample.
    Returns DataFrame with Sample, OS_DAYS, OS_EVENT.
    """
    n_samples = len(sample_ids)
    os_months = {}
    os_event = {}

    for key, vals in meta_rows:
        for i, val in enumerate(vals):
            if i >= n_samples:
                break
            sid = sample_ids[i] if i < len(sample_ids) else None
            if not sid or not sid.startswith('GSM'):
                continue
            val_str = str(val).strip()
            val_lower = val_str.lower()

            # Survival time: value format is "survival time in months: 24.24" or "overall survival (months): 12"
            if 'survival time' in val_lower or 'overall survival' in val_lower or ('survival' in val_lower and 'month' in val_lower):
                m = re.search(r'[\d.]+', val_str)
                if m:
                    v = float(m.group())
                    if v < 500:  # months
                        os_months[sid] = v * 30.44  # to days
                    else:
                        os_months[sid] = v

            # Survival status: "survival status: 1" (1=dead) or "last follow-up: Dead"
            if 'survival status' in val_lower:
                if '1' in val_str or 'dead' in val_lower or 'deceased' in val_lower:
                    os_event[sid] = 1
                elif '0' in val_str or 'alive' in val_lower or 'censored' in val_lower or 'na' in val_lower:
                    os_event[sid] = 0
            if 'last follow-up' in val_lower:
                if 'dead' in val_lower or 'deceased' in val_lower:
                    os_event[sid] = 1
                elif 'alive' in val_lower or 'censored' in val_lower:
                    os_event[sid] = 0

    rows = []
    for sid in sample_ids:
        days = os_months.get(sid, np.nan)
        event = os_event.get(sid, 1)
        rows.append({'Sample': sid, 'OS_DAYS': days, 'OS_EVENT': event})

    df = pd.DataFrame(rows)
    valid = df.dropna(subset=['OS_DAYS'])
    valid = valid[valid['OS_DAYS'] > 0]
    return valid


def run_survival(gene_vals, clinical, high_better=True):
    """Log-rank test. high_better: for methylation True, for expression False (low expr = better)."""
    if logrank_test is None:
        return None
    df = clinical.copy()
    df['val'] = df['Sample'].map(gene_vals)
    df = df.dropna(subset=['val', 'OS_DAYS', 'OS_EVENT'])
    df = df[df['OS_DAYS'] > 0]
    if len(df) < 15:
        return None
    thresh = df['val'].median()
    high = df[df['val'] > thresh]
    low = df[df['val'] <= thresh]
    if len(high) < 5 or len(low) < 5:
        return None
    try:
        lr = logrank_test(high['OS_DAYS'], low['OS_DAYS'], high['OS_EVENT'], low['OS_EVENT'])
        pval = lr.p_value
    except Exception:
        pval = np.nan
    diff = high['OS_DAYS'].median() - low['OS_DAYS'].median()
    concordant = (high_better and diff > 0) or (not high_better and diff < 0)
    return {
        'n_high': len(high), 'n_low': len(low),
        'median_high': high['OS_DAYS'].median(), 'median_low': low['OS_DAYS'].median(),
        'survival_diff': diff, 'logrank_p': pval,
        'concordant': concordant, 'significant': pval < 0.05 if not np.isnan(pval) else False
    }


def get_gene_values_affy(df, sample_ids):
    """Extract gene-level values from Affymetrix data using probe mapping."""
    gene_vals = {g: {} for g in BIOMARKERS}
    for probe, gene in AFFY_PROBE_TO_GENE.items():
        if probe not in df.index or gene not in BIOMARKERS:
            continue
        row = df.loc[probe]
        for sid in sample_ids:
            if sid in df.columns:
                v = row[sid]
                if not np.isnan(v):
                    gene_vals[gene].setdefault(sid, []).append(v)
    result = {}
    for g in BIOMARKERS:
        for sid, vals in gene_vals[g].items():
            result[sid] = np.mean(vals)
    return result


def get_gene_values_methylation(df, sample_ids, probe_map):
    """Extract gene-level values from methylation data (mean of probes per gene)."""
    gene_vals = {g: {} for g in BIOMARKERS}
    for probe, gene in probe_map.items():
        if probe not in df.index or gene not in BIOMARKERS:
            continue
        row = df.loc[probe]
        for sid in sample_ids:
            if sid in df.columns:
                v = row[sid]
                if not np.isnan(v):
                    gene_vals[gene].setdefault(sid, []).append(v)
    result = {}
    for g in BIOMARKERS:
        for sid, vals in gene_vals[g].items():
            result[sid] = np.mean(vals)
    return result


def validate_dataset(gse_id, data_type='expression', probe_map=None, gbm_only=True):
    """Validate biomarkers in a GEO dataset."""
    path = EXTERNAL_DIR / gse_id / f"{gse_id}_series_matrix.txt.gz"
    if not path.exists():
        path = EXTERNAL_DIR / gse_id / f"{gse_id}_series_matrix.txt"
    if not path.exists():
        return {'gse': gse_id, 'status': 'file_not_found'}

    df, sample_ids, meta_rows = parse_geo_matrix(path)
    if df is None or df.empty:
        return {'gse': gse_id, 'status': 'parse_failed'}

    clinical = extract_survival_from_meta(meta_rows, sample_ids)

    if len(clinical) < 15:
        return {'gse': gse_id, 'status': 'insufficient_survival', 'n': len(clinical),
                'n_samples': len(sample_ids)}

    clinical = clinical.set_index('Sample')
    high_better = (data_type == 'methylation')

    if data_type == 'expression':
        gene_vals_all = get_gene_values_affy(df, sample_ids)
    else:
        probe_map = probe_map or EPIC_PROBE_TO_GENE
        gene_vals_all = get_gene_values_methylation(df, sample_ids, probe_map)

    results = {}
    for gene in BIOMARKERS:
        gene_vals = {s: gene_vals_all[s] for s in gene_vals_all if gene in str(gene_vals_all.get(s, '')) or
                     any(k for k, v in (AFFY_PROBE_TO_GENE if data_type == 'expression' else probe_map).items() if v == gene)}
        # Actually get values for this gene
        gene_vals = {s: v for s, v in gene_vals_all.items()}
        # Filter: we need per-gene values - gene_vals_all is gene->sample->value
        # Fix: gene_vals_all structure is sample->value (single value per sample, averaged across probes)
        # We need sample->value per gene. Let me fix get_gene_values to return gene->sample->value
        pass  # gene_vals is built per gene in get_gene_values - it returns sample->value for ALL genes combined
    # Fix the logic: get_gene_values returns {sample: value} for each gene - but we're averaging probes so we get one value per sample per gene. The structure needs to be gene -> sample -> value.
    # Let me fix get_gene_values to return {gene: {sample: value}}
    return _validate_with_gene_vals(gene_vals_all, clinical, high_better, gse_id)


def _validate_with_gene_vals(gene_vals_by_gene, clinical, high_better, gse_id):
    """gene_vals_by_gene: dict gene -> {sample -> value}"""
    results = {}
    for gene in BIOMARKERS:
        if gene not in gene_vals_by_gene:
            continue
        vals = gene_vals_by_gene[gene]
        common = list(set(vals.keys()) & set(clinical.index))
        if len(common) < 15:
            continue
        cl = clinical.loc[clinical.index.isin(common)]
        gene_vals = pd.Series(vals)
        res = run_survival(gene_vals, cl.reset_index(), high_better=high_better)
        if res:
            results[gene] = res
    return {'gse': gse_id, 'n_samples': len(clinical), 'n_survival': len(clinical),
            'biomarker_results': results}


def get_gene_values_affy_v2(df, sample_ids):
    """Returns {gene: {sample: value}}."""
    gene_vals = {g: {} for g in BIOMARKERS}
    for probe, gene in AFFY_PROBE_TO_GENE.items():
        if probe not in df.index:
            continue
        row = df.loc[probe]
        for sid in sample_ids:
            if sid in df.columns:
                v = row[sid]
                if not np.isnan(v):
                    gene_vals[gene].setdefault(sid, []).append(v)
    result = {}
    for g in BIOMARKERS:
        for sid, vals in gene_vals[g].items():
            result[g] = result.get(g, {})
            result[g][sid] = np.mean(vals)
    return result


def get_gene_values_meth_v2(df, sample_ids, probe_map):
    """Returns {gene: {sample: value}}."""
    gene_vals = {g: {} for g in BIOMARKERS}
    for probe, gene in probe_map.items():
        if probe not in df.index:
            continue
        row = df.loc[probe]
        for sid in sample_ids:
            if sid in df.columns:
                v = row[sid]
                if not np.isnan(v):
                    gene_vals[gene].setdefault(sid, []).append(v)
    result = {}
    for g in BIOMARKERS:
        for sid, vals in gene_vals[g].items():
            result[g] = result.get(g, {})
            result[g][sid] = np.mean(vals)
    return result


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    datasets = [
        ('GSE7696', 'expression'),
        ('GSE230770', 'methylation'),
        # GSE22867: no survival in GEO metadata
        # GSE108474: no survival in GEO metadata
    ]

    all_results = []
    for gse_id, dtype in datasets:
        print(f"\n--- {gse_id} ({dtype}) ---")
        path = EXTERNAL_DIR / gse_id / f"{gse_id}_series_matrix.txt.gz"
        if not path.exists():
            print(f"  File not found")
            continue
        probe_subset = None
        if dtype == 'methylation' and gse_id == 'GSE230770':
            probe_subset = set(EPIC_PROBE_TO_GENE.keys())
        elif dtype == 'expression' and gse_id == 'GSE7696':
            probe_subset = set(AFFY_PROBE_TO_GENE.keys())
        df, sample_ids, meta_rows = parse_geo_matrix(path, probe_subset=probe_subset)
        if df is None:
            print(f"  Parse failed")
            continue
        clinical = extract_survival_from_meta(meta_rows, sample_ids)
        if len(clinical) < 15:
            print(f"  Insufficient survival: n={len(clinical)}")
            continue
        clinical = clinical.set_index('Sample')
        high_better = (dtype == 'methylation')

        if dtype == 'expression':
            gene_vals = get_gene_values_affy_v2(df, sample_ids)
        else:
            # GSE230770: EPIC - find probes that exist in data
            all_probes = set(df.index)
            probe_map = {p: g for p, g in EPIC_PROBE_TO_GENE.items() if p in all_probes}
            if not probe_map:
                # Fallback: try Illumina 27K probes if EPIC
                probe_map = {p: g for p, g in ILLUMINA27K_PROBE_TO_GENE.items() if p in all_probes}
            if not probe_map:
                # Use any cg probe - we need to map. Try first 10 probes per gene from EPIC
                for p in all_probes:
                    if p.startswith('cg') and p not in probe_map:
                        pass  # Can't map without annotation
                print(f"  No probe mapping found. Available probes: {len(all_probes)}")
                continue
            gene_vals = get_gene_values_meth_v2(df, sample_ids, probe_map)

        results = {}
        for gene in BIOMARKERS:
            if gene not in gene_vals or not gene_vals[gene]:
                continue
            vals = gene_vals[gene]
            common = list(set(vals.keys()) & set(clinical.index))
            if len(common) < 15:
                continue
            cl = clinical.loc[clinical.index.isin(common)]
            gene_series = pd.Series(vals)
            res = run_survival(gene_series, cl.reset_index(), high_better=high_better)
            if res:
                results[gene] = res
                tcga = TCGA_RESULTS.get(gene, {})
                all_results.append({
                    'Dataset': gse_id, 'Data_Type': dtype, 'Gene': gene,
                    'N_High': res['n_high'], 'N_Low': res['n_low'],
                    'Median_High': res['median_high'], 'Median_Low': res['median_low'],
                    'Survival_Diff': res['survival_diff'], 'LogRank_P': res['logrank_p'],
                    'Significant': res['significant'], 'Concordant_TCGA': res['concordant'],
                    'TCGA_Diff': tcga.get('diff'), 'TCGA_P': tcga.get('p')
                })
        print(f"  Validated {len(results)} genes: {list(results.keys())}")

    if not all_results:
        print("\nNo validation results.")
        return

    df_out = pd.DataFrame(all_results)
    df_out.to_csv(RESULTS_DIR / "external_validation_results.csv", index=False)

    summary = []
    for gene in BIOMARKERS:
        sub = df_out[df_out['Gene'] == gene]
        n_val = len(sub)
        n_sig = sub['Significant'].sum()
        n_conc = sub['Concordant_TCGA'].sum()
        summary.append({
            'Gene': gene, 'TCGA_Diff': TCGA_RESULTS.get(gene, {}).get('diff'),
            'N_Validated_Datasets': n_val, 'N_Significant': int(n_sig),
            'N_Concordant': int(n_conc), 'Validation_Rate': f"{n_conc}/{n_val}" if n_val else "-"
        })
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(RESULTS_DIR / "external_validation_summary.csv", index=False)

    if len(all_results) > 0:
        fig, axes = plt.subplots(2, 1, figsize=(10, 8))
        pivot = df_out.pivot_table(index='Gene', columns='Dataset', values='Survival_Diff')
        pivot.plot(kind='bar', ax=axes[0], width=0.8)
        axes[0].axhline(0, color='black', linewidth=0.5)
        axes[0].set_ylabel('Survival Diff (days)')
        axes[0].set_title('External Validation: Survival Difference by Dataset')
        axes[0].legend(bbox_to_anchor=(1.02, 1))
        axes[0].tick_params(axis='x', rotation=45)

        df_out['neg_log_p'] = -np.log10(df_out['LogRank_P'].replace(0, 1e-300))
        pivot2 = df_out.pivot_table(index='Gene', columns='Dataset', values='neg_log_p')
        pivot2.plot(kind='bar', ax=axes[1], width=0.8)
        axes[1].axhline(-np.log10(0.05), color='red', linestyle='--')
        axes[1].set_ylabel('-log10(p)')
        axes[1].set_title('Statistical Significance')
        axes[1].legend(bbox_to_anchor=(1.02, 1))
        axes[1].tick_params(axis='x', rotation=45)

        plt.tight_layout()
        plt.savefig(PLOTS_DIR / "external_validation_comparison.png", dpi=150, bbox_inches='tight')
        plt.close()

    print("\n" + "=" * 60)
    print("EXTERNAL VALIDATION SUMMARY")
    print("=" * 60)
    print(summary_df.to_string(index=False))
    print("\nResults:", RESULTS_DIR)
    print("Plots:", PLOTS_DIR)


if __name__ == "__main__":
    main()
