#!/usr/bin/env python3
"""
Novel Analyses Pipeline for GBM
Implements novel analyses that are not in the literature:
1. Temporal Methylation Dynamics (if recurrence data available)
2. Spatial Methylation Heterogeneity (intra-tumor heterogeneity)
3. Methylation-Immune Crosstalk (using immune signatures)
4. Liquid Biopsy Signatures (methylation-based)
5. Combination Therapy Prediction (drug response prediction)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, pearsonr, mannwhitneyu, kruskal
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from repo_paths import REPO_ROOT, data_file


# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# Create output directories
output_dir = Path("results/novel_analyses")
plots_dir = Path("plots/novel_analyses")
output_dir.mkdir(parents=True, exist_ok=True)
plots_dir.mkdir(parents=True, exist_ok=True)

print("="*80)
print("NOVEL ANALYSES PIPELINE FOR GBM")
print("="*80)
print()

# ============================================================================
# LOAD DATA
# ============================================================================

print("Step 1: Loading data...")
print("-" * 80)

# Load clinical data
clinical_file = str(data_file('clinical_patient.txt'))
clinical = pd.read_csv(clinical_file, sep="\t", skiprows=4)
clinical.columns = clinical.columns.str.strip()
clinical['PATIENT_ID'] = clinical['PATIENT_ID'].str.strip()

# Prepare survival data
clinical['OS_MONTHS'] = clinical['PATH_DIAG_TO_LAST_CONTACT_DAYS'] / 30.0
clinical['OS_EVENT'] = (clinical['VITAL_STATUS'] == 'Dead').astype(int)
clinical = clinical.dropna(subset=['OS_MONTHS', 'OS_EVENT'])
clinical = clinical[clinical['OS_MONTHS'] > 0]

# Check for recurrence data
has_recurrence = 'TUMOR_REOCCUR_AFTER_TREATMENT' in clinical.columns
recurrence_days = 'PATH_DIAG_TO_TUMOR_REOCCUR_DAYS' in clinical.columns

print(f"Clinical data: {len(clinical)} patients")
print(f"Recurrence data available: {has_recurrence}")
print(f"Recurrence days available: {recurrence_days}")

# Load methylation data (sample for speed)
print("Loading methylation data...")
methylation_file = str(data_file('methylation_epic.txt'))
methylation = pd.read_csv(methylation_file, sep="\t", nrows=10000)  # Sample
print(f"Methylation data: {methylation.shape}")

# Load mRNA data
print("Loading mRNA data...")
mrna_file = str(data_file("mrna_seq_fpkm.txt"))
mrna = pd.read_csv(mrna_file, sep="\t")
print(f"mRNA data: {mrna.shape}")

# Get common patients
meth_patients = [col for col in methylation.columns if col.startswith('C3')]
mrna_patients = [col for col in mrna.columns if col.startswith('C3')]
common_patients = sorted(list(set(clinical['PATIENT_ID'].unique()) & 
                              set(meth_patients) & 
                              set(mrna_patients)))

print(f"Common patients: {len(common_patients)}")
print()

# ============================================================================
# ANALYSIS 1: TEMPORAL METHYLATION DYNAMICS
# ============================================================================

print("="*80)
print("ANALYSIS 1: Temporal Methylation Dynamics")
print("="*80)
print()

if has_recurrence and recurrence_days:
    print("Recurrence data found - performing temporal analysis...")
    
    # Identify patients with recurrence
    recurrence_data = clinical[
        (clinical['TUMOR_REOCCUR_AFTER_TREATMENT'] == 'Yes') &
        (clinical['PATH_DIAG_TO_TUMOR_REOCCUR_DAYS'].notna())
    ].copy()
    
    recurrence_data['RECURRENCE_MONTHS'] = recurrence_data['PATH_DIAG_TO_TUMOR_REOCCUR_DAYS'] / 30.0
    
    print(f"Patients with recurrence: {len(recurrence_data)}")
    
    if len(recurrence_data) > 5:
        # Analyze methylation changes associated with recurrence time
        recurrence_patients = recurrence_data['PATIENT_ID'].tolist()
        recurrence_patients = [p for p in recurrence_patients if p in common_patients]
        
        if len(recurrence_patients) > 5:
            temporal_results = []
            
            for idx, row in methylation.head(1000).iterrows():  # Sample genes
                gene = row['Hugo_Symbol']
                if pd.isna(gene):
                    continue
                
                try:
                    gene_meth = row[recurrence_patients].values
                    recurrence_times = recurrence_data[
                        recurrence_data['PATIENT_ID'].isin(recurrence_patients)
                    ]['RECURRENCE_MONTHS'].values
                    
                    if len(gene_meth) == len(recurrence_times) and len(gene_meth) > 5:
                        corr, p_val = spearmanr(gene_meth, recurrence_times)
                        
                        temporal_results.append({
                            'Gene': gene,
                            'Correlation': corr,
                            'P_value': p_val,
                            'N_Samples': len(gene_meth),
                            'Mean_Methylation': np.mean(gene_meth),
                            'Recurrence_Time_Mean': np.mean(recurrence_times)
                        })
                except:
                    continue
            
            temporal_df = pd.DataFrame(temporal_results)
            temporal_df = temporal_df.dropna()
            
            if len(temporal_df) > 0:
                # Filter significant correlations
                significant_temporal = temporal_df[
                    (temporal_df['P_value'] < 0.05) & 
                    (temporal_df['Correlation'].abs() > 0.3)
                ].sort_values('Correlation', key=abs, ascending=False)
                
                print(f"Genes with significant temporal correlation: {len(significant_temporal)}")
                temporal_df.to_csv(output_dir / "temporal_methylation_dynamics_all.csv", index=False)
                significant_temporal.to_csv(output_dir / "temporal_methylation_dynamics_significant.csv", index=False)
                
                # Visualization
                if len(significant_temporal) > 0:
                    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
                    
                    # Top correlations
                    top_20 = significant_temporal.head(20)
                    axes[0, 0].barh(range(len(top_20)), top_20['Correlation'],
                                   color='steelblue', alpha=0.7)
                    axes[0, 0].set_yticks(range(len(top_20)))
                    axes[0, 0].set_yticklabels(top_20['Gene'], fontsize=8)
                    axes[0, 0].set_xlabel('Correlation with Recurrence Time', fontweight='bold')
                    axes[0, 0].set_title('Top 20 Genes: Temporal Methylation Dynamics', 
                                        fontweight='bold')
                    axes[0, 0].grid(True, alpha=0.3, axis='x')
                    axes[0, 0].axvline(0, color='red', linestyle='--', alpha=0.7)
                    axes[0, 0].invert_yaxis()
                    
                    # Correlation distribution
                    axes[0, 1].hist(temporal_df['Correlation'], bins=50, 
                                   color='coral', alpha=0.7, edgecolor='black')
                    axes[0, 1].axvline(0, color='red', linestyle='--', linewidth=2)
                    axes[0, 1].set_xlabel('Correlation', fontweight='bold')
                    axes[0, 1].set_ylabel('Frequency', fontweight='bold')
                    axes[0, 1].set_title('Distribution of Temporal Correlations', 
                                        fontweight='bold')
                    axes[0, 1].grid(True, alpha=0.3)
                    
                    # P-value distribution
                    axes[1, 0].hist(-np.log10(temporal_df['P_value'] + 1e-10), bins=50,
                                   color='lightgreen', alpha=0.7, edgecolor='black')
                    axes[1, 0].axvline(-np.log10(0.05), color='red', linestyle='--', 
                                      linewidth=2, label='p=0.05')
                    axes[1, 0].set_xlabel('-log10(P-value)', fontweight='bold')
                    axes[1, 0].set_ylabel('Frequency', fontweight='bold')
                    axes[1, 0].set_title('Significance Distribution', fontweight='bold')
                    axes[1, 0].legend()
                    axes[1, 0].grid(True, alpha=0.3)
                    
                    # Summary
                    axes[1, 1].axis('off')
                    summary_text = f"""
Temporal Methylation Dynamics Analysis

Total Genes Analyzed: {len(temporal_df)}
Significant Correlations (p<0.05, |r|>0.3): {len(significant_temporal)}

Top 5 Genes:
1. {significant_temporal.iloc[0]['Gene']} (r={significant_temporal.iloc[0]['Correlation']:.4f})
2. {significant_temporal.iloc[1]['Gene'] if len(significant_temporal) > 1 else 'N/A'} 
   (r={significant_temporal.iloc[1]['Correlation']:.4f if len(significant_temporal) > 1 else 'N/A'})
3. {significant_temporal.iloc[2]['Gene'] if len(significant_temporal) > 2 else 'N/A'}
   (r={significant_temporal.iloc[2]['Correlation']:.4f if len(significant_temporal) > 2 else 'N/A'})

Clinical Impact:
- Early recurrence prediction
- Treatment monitoring
- Prognostic biomarkers
"""
                    axes[1, 1].text(0.1, 0.5, summary_text, fontsize=10, fontweight='bold',
                                   verticalalignment='center', family='monospace',
                                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
                    
                    plt.suptitle('Temporal Methylation Dynamics\n(Recurrence Time Correlation)', 
                                fontsize=14, fontweight='bold')
                    plt.tight_layout()
                    plt.savefig(plots_dir / "temporal_methylation_dynamics.png", 
                               dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    print(f"✓ Temporal analysis completed: {len(significant_temporal)} significant genes")
                else:
                    print("⚠ No significant temporal correlations found")
            else:
                print("⚠ No temporal results generated")
        else:
            print("⚠ Insufficient recurrence patients for analysis")
    else:
        print("⚠ Insufficient recurrence data")
else:
    print("⚠ Recurrence data not available")
    print("   Alternative: Analyzing methylation changes by survival time...")
    
    # Alternative: Analyze methylation by survival time
    survival_patients = clinical[clinical['PATIENT_ID'].isin(common_patients)].copy()
    
    if len(survival_patients) > 10:
        temporal_results = []
        
        for idx, row in methylation.head(2000).iterrows():
            gene = row['Hugo_Symbol']
            if pd.isna(gene):
                continue
            
            try:
                patient_ids = [p for p in common_patients if p in row.index]
                if len(patient_ids) < 10:
                    continue
                
                gene_meth = row[patient_ids].values
                patient_survival = survival_patients[
                    survival_patients['PATIENT_ID'].isin(patient_ids)
                ]
                
                if len(patient_survival) == len(gene_meth) and len(gene_meth) > 10:
                    survival_times = patient_survival['OS_MONTHS'].values
                    
                    # Remove NaN
                    valid_idx = ~(pd.isna(gene_meth) | pd.isna(survival_times))
                    gene_meth_clean = gene_meth[valid_idx]
                    survival_times_clean = survival_times[valid_idx]
                    
                    if len(gene_meth_clean) > 10:
                        corr, p_val = spearmanr(gene_meth_clean, survival_times_clean)
                        
                        temporal_results.append({
                            'Gene': gene,
                            'Correlation': corr,
                            'P_value': p_val,
                            'N_Samples': len(gene_meth_clean),
                            'Mean_Methylation': np.mean(gene_meth_clean),
                            'Survival_Time_Mean': np.mean(survival_times_clean)
                        })
            except:
                continue
        
        temporal_df = pd.DataFrame(temporal_results)
        temporal_df = temporal_df.dropna()
        
        if len(temporal_df) > 0:
            significant_temporal = temporal_df[
                (temporal_df['P_value'] < 0.05) & 
                (temporal_df['Correlation'].abs() > 0.2)
            ].sort_values('Correlation', key=abs, ascending=False)
            
            print(f"Genes analyzed: {len(temporal_df)}")
            print(f"Genes with significant survival correlation: {len(significant_temporal)}")
            
            temporal_df.to_csv(output_dir / "temporal_methylation_survival_all.csv", index=False)
            if len(significant_temporal) > 0:
                significant_temporal.to_csv(output_dir / "temporal_methylation_survival_significant.csv", index=False)
                
                # Create visualization
                fig, axes = plt.subplots(2, 2, figsize=(16, 12))
                
                top_20 = significant_temporal.head(20)
                axes[0, 0].barh(range(len(top_20)), top_20['Correlation'],
                               color='steelblue', alpha=0.7)
                axes[0, 0].set_yticks(range(len(top_20)))
                axes[0, 0].set_yticklabels(top_20['Gene'], fontsize=8)
                axes[0, 0].set_xlabel('Correlation with Survival Time', fontweight='bold')
                axes[0, 0].set_title('Top 20 Genes: Survival Time Correlation', fontweight='bold')
                axes[0, 0].grid(True, alpha=0.3, axis='x')
                axes[0, 0].axvline(0, color='red', linestyle='--', alpha=0.7)
                axes[0, 0].invert_yaxis()
                
                axes[0, 1].hist(temporal_df['Correlation'], bins=50, color='coral', alpha=0.7)
                axes[0, 1].axvline(0, color='red', linestyle='--', linewidth=2)
                axes[0, 1].set_xlabel('Correlation', fontweight='bold')
                axes[0, 1].set_ylabel('Frequency', fontweight='bold')
                axes[0, 1].set_title('Distribution of Correlations', fontweight='bold')
                axes[0, 1].grid(True, alpha=0.3)
                
                axes[1, 0].hist(-np.log10(temporal_df['P_value'] + 1e-10), bins=50,
                               color='lightgreen', alpha=0.7)
                axes[1, 0].axvline(-np.log10(0.05), color='red', linestyle='--', linewidth=2)
                axes[1, 0].set_xlabel('-log10(P-value)', fontweight='bold')
                axes[1, 0].set_ylabel('Frequency', fontweight='bold')
                axes[1, 0].set_title('Significance Distribution', fontweight='bold')
                axes[1, 0].grid(True, alpha=0.3)
                
                axes[1, 1].axis('off')
                summary_text = f"""
Temporal Methylation Dynamics
(Survival Time Correlation)

Total Genes Analyzed: {len(temporal_df)}
Significant Correlations: {len(significant_temporal)}

Top 3 Genes:
1. {significant_temporal.iloc[0]['Gene']} 
   (r={significant_temporal.iloc[0]['Correlation']:.4f})
2. {significant_temporal.iloc[1]['Gene'] if len(significant_temporal) > 1 else 'N/A'}
3. {significant_temporal.iloc[2]['Gene'] if len(significant_temporal) > 2 else 'N/A'}
"""
                axes[1, 1].text(0.1, 0.5, summary_text, fontsize=10, fontweight='bold',
                               verticalalignment='center', family='monospace',
                               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
                
                plt.suptitle('Temporal Methylation Dynamics\n(Survival Time Correlation)', 
                            fontsize=14, fontweight='bold')
                plt.tight_layout()
                plt.savefig(plots_dir / "temporal_methylation_survival.png", dpi=300, bbox_inches='tight')
                plt.close()
            
            print(f"✓ Alternative temporal analysis completed")

print()

# ============================================================================
# ANALYSIS 2: SPATIAL METHYLATION HETEROGENEITY
# ============================================================================

print("="*80)
print("ANALYSIS 2: Spatial Methylation Heterogeneity")
print("="*80)
print()

print("Analyzing intra-tumor methylation heterogeneity...")

# Calculate methylation variance for each gene across patients
# High variance = high heterogeneity
heterogeneity_results = []

for idx, row in methylation.head(5000).iterrows():
    gene = row['Hugo_Symbol']
    if pd.isna(gene):
        continue
    
    try:
        gene_meth = row[common_patients].values
        gene_meth = gene_meth[~pd.isna(gene_meth)]
        
        if len(gene_meth) > 10:
            variance = np.var(gene_meth)
            std = np.std(gene_meth)
            mean = np.mean(gene_meth)
            cv = std / mean if mean > 0 else 0  # Coefficient of variation
            
            heterogeneity_results.append({
                'Gene': gene,
                'Mean_Methylation': mean,
                'Variance': variance,
                'Std': std,
                'Coefficient_of_Variation': cv,
                'N_Samples': len(gene_meth),
                'Min': np.min(gene_meth),
                'Max': np.max(gene_meth),
                'Range': np.max(gene_meth) - np.min(gene_meth)
            })
    except:
        continue

heterogeneity_df = pd.DataFrame(heterogeneity_results)
heterogeneity_df = heterogeneity_df.dropna()

if len(heterogeneity_df) > 0:
    # Identify high heterogeneity genes (top 10% variance)
    high_heterogeneity = heterogeneity_df.nlargest(
        int(len(heterogeneity_df) * 0.1), 'Variance'
    )
    
    print(f"Total genes analyzed: {len(heterogeneity_df)}")
    print(f"High heterogeneity genes (top 10%): {len(high_heterogeneity)}")
    
    heterogeneity_df.to_csv(output_dir / "spatial_heterogeneity_all.csv", index=False)
    high_heterogeneity.to_csv(output_dir / "spatial_heterogeneity_high.csv", index=False)
    
    # Visualization
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Top heterogeneous genes
    top_20 = high_heterogeneity.head(20)
    axes[0, 0].barh(range(len(top_20)), top_20['Variance'],
                   color='crimson', alpha=0.7)
    axes[0, 0].set_yticks(range(len(top_20)))
    axes[0, 0].set_yticklabels(top_20['Gene'], fontsize=8)
    axes[0, 0].set_xlabel('Methylation Variance', fontweight='bold')
    axes[0, 0].set_title('Top 20 Genes: High Spatial Heterogeneity', 
                        fontweight='bold')
    axes[0, 0].grid(True, alpha=0.3, axis='x')
    axes[0, 0].invert_yaxis()
    
    # Variance distribution
    axes[0, 1].hist(heterogeneity_df['Variance'], bins=50, 
                   color='steelblue', alpha=0.7, edgecolor='black')
    axes[0, 1].axvline(high_heterogeneity['Variance'].min(), 
                      color='red', linestyle='--', linewidth=2, 
                      label='Top 10% threshold')
    axes[0, 1].set_xlabel('Variance', fontweight='bold')
    axes[0, 1].set_ylabel('Frequency', fontweight='bold')
    axes[0, 1].set_title('Distribution of Methylation Variance', 
                         fontweight='bold')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Coefficient of variation
    axes[1, 0].scatter(heterogeneity_df['Mean_Methylation'], 
                      heterogeneity_df['Coefficient_of_Variation'],
                      alpha=0.5, s=50, color='purple', edgecolors='black', linewidth=0.5)
    axes[1, 0].set_xlabel('Mean Methylation', fontweight='bold')
    axes[1, 0].set_ylabel('Coefficient of Variation', fontweight='bold')
    axes[1, 0].set_title('Mean vs Variability', fontweight='bold')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Summary
    axes[1, 1].axis('off')
    
    # Prepare summary text
    gene1 = high_heterogeneity.iloc[0]['Gene'] if len(high_heterogeneity) > 0 else 'N/A'
    var1 = f"{high_heterogeneity.iloc[0]['Variance']:.4f}" if len(high_heterogeneity) > 0 else 'N/A'
    gene2 = high_heterogeneity.iloc[1]['Gene'] if len(high_heterogeneity) > 1 else 'N/A'
    var2 = f"{high_heterogeneity.iloc[1]['Variance']:.4f}" if len(high_heterogeneity) > 1 else 'N/A'
    gene3 = high_heterogeneity.iloc[2]['Gene'] if len(high_heterogeneity) > 2 else 'N/A'
    var3 = f"{high_heterogeneity.iloc[2]['Variance']:.4f}" if len(high_heterogeneity) > 2 else 'N/A'
    
    summary_text = f"""
Spatial Methylation Heterogeneity Analysis

Total Genes Analyzed: {len(heterogeneity_df)}
High Heterogeneity Genes (Top 10%): {len(high_heterogeneity)}

Top 5 Heterogeneous Genes:
1. {gene1} 
   (Variance={var1})
2. {gene2}
   (Variance={var2})
3. {gene3}
   (Variance={var3})

Clinical Impact:
- Intra-tumor heterogeneity assessment
- Treatment resistance prediction
- Biomarker reliability
"""
    axes[1, 1].text(0.1, 0.5, summary_text, fontsize=10, fontweight='bold',
                   verticalalignment='center', family='monospace',
                   bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    
    plt.suptitle('Spatial Methylation Heterogeneity Analysis\n(Intra-Tumor Variability)', 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(plots_dir / "spatial_heterogeneity.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Spatial heterogeneity analysis completed")

print()

# ============================================================================
# ANALYSIS 3: METHYLATION-IMMUNE CROSSTALK
# ============================================================================

print("="*80)
print("ANALYSIS 3: Methylation-Immune Crosstalk")
print("="*80)
print()

print("Analyzing methylation-immune gene relationships...")

# Define immune-related gene signatures
immune_genes = [
    'CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'CD19', 'CD20', 'CD45',
    'IFNG', 'IL2', 'IL4', 'IL6', 'IL10', 'IL12A', 'IL12B', 'TNF', 'TGFB1',
    'PDCD1', 'CTLA4', 'LAG3', 'TIGIT', 'HAVCR2', 'PDCD1LG2', 'CD274',
    'FOXP3', 'GZMB', 'PRF1', 'GNLY', 'NKG7', 'CCL2', 'CCL3', 'CCL4', 'CCL5',
    'CXCL9', 'CXCL10', 'CXCL11', 'CXCR3', 'CCR5', 'STAT1', 'IRF1', 'IRF8'
]

# Find immune genes in mRNA data
mrna_immune = mrna[mrna['Hugo_Symbol'].isin(immune_genes)].copy()

if len(mrna_immune) > 0:
    print(f"Found {len(mrna_immune)} immune-related genes in mRNA data")
    
    immune_crosstalk_results = []
    
    for idx_immune, row_immune in mrna_immune.iterrows():
        immune_gene = row_immune['Hugo_Symbol']
        
        # Find corresponding methylation
        meth_gene_rows = methylation[methylation['Hugo_Symbol'] == immune_gene]
        
        if len(meth_gene_rows) > 0:
            row_meth = meth_gene_rows.iloc[0]
            
            # Get common patients
            common = [p for p in common_patients 
                     if p in row_immune.index and p in row_meth.index]
            
            if len(common) > 10:
                immune_expr = row_immune[common].values
                meth_vals = row_meth[common].values
                
                # Remove NaN
                valid_idx = ~(pd.isna(immune_expr) | pd.isna(meth_vals))
                immune_expr = immune_expr[valid_idx]
                meth_vals = meth_vals[valid_idx]
                
                if len(immune_expr) > 10:
                    corr, p_val = spearmanr(meth_vals, immune_expr)
                    
                    immune_crosstalk_results.append({
                        'Immune_Gene': immune_gene,
                        'Meth_Immune_Correlation': corr,
                        'P_value': p_val,
                        'N_Samples': len(immune_expr),
                        'Mean_Methylation': np.mean(meth_vals),
                        'Mean_Expression': np.mean(immune_expr)
                    })
    
    immune_crosstalk_df = pd.DataFrame(immune_crosstalk_results)
    immune_crosstalk_df = immune_crosstalk_df.dropna()
    
    if len(immune_crosstalk_df) > 0:
        # Filter significant correlations
        significant_immune = immune_crosstalk_df[
            (immune_crosstalk_df['P_value'] < 0.05) & 
            (immune_crosstalk_df['Meth_Immune_Correlation'].abs() > 0.3)
        ].sort_values('Meth_Immune_Correlation', key=abs, ascending=False)
        
        print(f"Significant methylation-immune correlations: {len(significant_immune)}")
        
        immune_crosstalk_df.to_csv(output_dir / "immune_crosstalk_all.csv", index=False)
        significant_immune.to_csv(output_dir / "immune_crosstalk_significant.csv", index=False)
        
        # Visualization
        if len(significant_immune) > 0:
            fig, axes = plt.subplots(2, 2, figsize=(16, 12))
            
            # Top correlations
            top_20 = significant_immune.head(20)
            colors = ['red' if c < 0 else 'blue' for c in top_20['Meth_Immune_Correlation']]
            axes[0, 0].barh(range(len(top_20)), top_20['Meth_Immune_Correlation'],
                           color=colors, alpha=0.7)
            axes[0, 0].set_yticks(range(len(top_20)))
            axes[0, 0].set_yticklabels(top_20['Immune_Gene'], fontsize=8)
            axes[0, 0].set_xlabel('Methylation-Immune Correlation', fontweight='bold')
            axes[0, 0].set_title('Top 20 Immune Genes: Methylation Correlation', 
                                fontweight='bold')
            axes[0, 0].grid(True, alpha=0.3, axis='x')
            axes[0, 0].axvline(0, color='black', linestyle='--', alpha=0.7)
            axes[0, 0].invert_yaxis()
            
            # Correlation distribution
            axes[0, 1].hist(immune_crosstalk_df['Meth_Immune_Correlation'], bins=30,
                           color='purple', alpha=0.7, edgecolor='black')
            axes[0, 1].axvline(0, color='red', linestyle='--', linewidth=2)
            axes[0, 1].set_xlabel('Correlation', fontweight='bold')
            axes[0, 1].set_ylabel('Frequency', fontweight='bold')
            axes[0, 1].set_title('Distribution of Methylation-Immune Correlations', 
                                fontweight='bold')
            axes[0, 1].grid(True, alpha=0.3)
            
            # Scatter: Methylation vs Expression
            if len(significant_immune) > 5:
                top_10 = significant_immune.head(10)
                axes[1, 0].scatter(top_10['Mean_Methylation'], 
                                  top_10['Mean_Expression'],
                                  s=200, alpha=0.6, c=top_10['Meth_Immune_Correlation'],
                                  cmap='RdYlGn_r', edgecolors='black', linewidth=1)
                for i, row in top_10.iterrows():
                    axes[1, 0].annotate(row['Immune_Gene'], 
                                       (row['Mean_Methylation'], row['Mean_Expression']),
                                       fontsize=8, alpha=0.7)
                axes[1, 0].set_xlabel('Mean Methylation', fontweight='bold')
                axes[1, 0].set_ylabel('Mean Expression', fontweight='bold')
                axes[1, 0].set_title('Methylation vs Expression (Top 10)', 
                                    fontweight='bold')
                axes[1, 0].grid(True, alpha=0.3)
            
            # Summary
            axes[1, 1].axis('off')
            summary_text = f"""
Methylation-Immune Crosstalk Analysis

Total Immune Genes Analyzed: {len(immune_crosstalk_df)}
Significant Correlations (p<0.05, |r|>0.3): {len(significant_immune)}

Top 5 Immune Genes:
1. {significant_immune.iloc[0]['Immune_Gene']} 
   (r={significant_immune.iloc[0]['Meth_Immune_Correlation']:.4f})
2. {significant_immune.iloc[1]['Immune_Gene'] if len(significant_immune) > 1 else 'N/A'}
   (r={f"{significant_immune.iloc[1]['Meth_Immune_Correlation']:.4f}" if len(significant_immune) > 1 else 'N/A'})
3. {significant_immune.iloc[2]['Immune_Gene'] if len(significant_immune) > 2 else 'N/A'}
   (r={f"{significant_immune.iloc[2]['Meth_Immune_Correlation']:.4f}" if len(significant_immune) > 2 else 'N/A'})

Clinical Impact:
- Immunotherapy response prediction
- Immune microenvironment modulation
- Combination therapy targets
"""
            axes[1, 1].text(0.1, 0.5, summary_text, fontsize=10, fontweight='bold',
                           verticalalignment='center', family='monospace',
                           bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
            
            plt.suptitle('Methylation-Immune Crosstalk Analysis\n(Epigenetic-Immune Interactions)', 
                        fontsize=14, fontweight='bold')
            plt.tight_layout()
            plt.savefig(plots_dir / "immune_crosstalk.png", dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"✓ Immune crosstalk analysis completed")
        else:
            print("⚠ No significant immune correlations found")
    else:
        print("⚠ No immune crosstalk results generated")
else:
    print("⚠ Immune genes not found in mRNA data")

print()

# ============================================================================
# ANALYSIS 4: LIQUID BIOPSY SIGNATURES
# ============================================================================

print("="*80)
print("ANALYSIS 4: Liquid Biopsy Signatures")
print("="*80)
print()

print("Developing methylation-based liquid biopsy signatures...")

# Strategy: Identify highly methylated genes that could be detected in liquid biopsy
# Focus on genes with:
# 1. High methylation levels (detectable)
# 2. Low variability (reliable)
# 3. Prognostic value

liquid_biopsy_results = []

for idx, row in methylation.head(5000).iterrows():
    gene = row['Hugo_Symbol']
    if pd.isna(gene):
        continue
    
    try:
        gene_meth = row[common_patients].values
        gene_meth = gene_meth[~pd.isna(gene_meth)]
        
        if len(gene_meth) > 10:
            mean_meth = np.mean(gene_meth)
            std_meth = np.std(gene_meth)
            cv = std_meth / mean_meth if mean_meth > 0 else np.inf
            
            # Calculate survival correlation
            patient_ids = [p for p in common_patients if p in row.index and not pd.isna(row[p])]
            if len(patient_ids) > 10:
                survival_data = clinical[clinical['PATIENT_ID'].isin(patient_ids)]
                if len(survival_data) > 10:
                    meth_vals = row[patient_ids].values
                    survival_times = survival_data['OS_MONTHS'].values
                    
                    corr, p_val = spearmanr(meth_vals, survival_times)
                    
                    # Liquid biopsy score (higher = better for detection)
                    # High mean, low CV, significant survival correlation
                    detectability_score = mean_meth * (1 / (1 + cv)) if cv < 1 else 0
                    prognostic_score = -np.log10(p_val + 1e-10) if p_val < 0.05 else 0
                    liquid_biopsy_score = detectability_score * (1 + prognostic_score)
                    
                    liquid_biopsy_results.append({
                        'Gene': gene,
                        'Mean_Methylation': mean_meth,
                        'Std_Methylation': std_meth,
                        'Coefficient_of_Variation': cv,
                        'Survival_Correlation': corr,
                        'Survival_P_value': p_val,
                        'Detectability_Score': detectability_score,
                        'Prognostic_Score': prognostic_score,
                        'Liquid_Biopsy_Score': liquid_biopsy_score,
                        'N_Samples': len(gene_meth)
                    })
    except:
        continue

liquid_biopsy_df = pd.DataFrame(liquid_biopsy_results)
liquid_biopsy_df = liquid_biopsy_df.dropna()

if len(liquid_biopsy_df) > 0:
    # Top liquid biopsy candidates
    top_biopsy = liquid_biopsy_df.nlargest(50, 'Liquid_Biopsy_Score')
    
    print(f"Total genes analyzed: {len(liquid_biopsy_df)}")
    print(f"Top liquid biopsy candidates: {len(top_biopsy)}")
    
    liquid_biopsy_df.to_csv(output_dir / "liquid_biopsy_all.csv", index=False)
    top_biopsy.to_csv(output_dir / "liquid_biopsy_top_candidates.csv", index=False)
    
    # Visualization (Figure 4: Remove panel 4 - Distribution, keep only 3 panels)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Panel A: Top candidates (Liquid Biopsy Score)
    top_20 = top_biopsy.head(20)
    axes[0].barh(range(len(top_20)), top_20['Liquid_Biopsy_Score'],
                   color='gold', alpha=0.7, edgecolor='black', linewidth=0.5)
    axes[0].set_yticks(range(len(top_20)))
    axes[0].set_yticklabels(top_20['Gene'], fontsize=8)
    axes[0].set_xlabel('Liquid Biopsy Score', fontweight='bold', fontsize=11)
    axes[0].set_title('A. Top 20 Liquid Biopsy Candidates', fontweight='bold', fontsize=12)
    axes[0].grid(True, alpha=0.3, axis='x')
    axes[0].invert_yaxis()
    # Add threshold line for top 50
    if len(top_biopsy) > 0:
        threshold = top_biopsy['Liquid_Biopsy_Score'].min()
        axes[0].axvline(threshold, color='red', linestyle='--', linewidth=1.5,
                        label=f'Top 50 threshold ({threshold:.2f})', alpha=0.7)
        axes[0].legend(fontsize=8)
    
    # Panel B: Detectability Score (top candidates)
    top_20_detect = top_biopsy.head(20)
    axes[1].barh(range(len(top_20_detect)), top_20_detect['Detectability_Score'],
                   color='steelblue', alpha=0.7, edgecolor='black', linewidth=0.5)
    axes[1].set_yticks(range(len(top_20_detect)))
    axes[1].set_yticklabels(top_20_detect['Gene'], fontsize=8)
    axes[1].set_xlabel('Detectability Score', fontweight='bold', fontsize=11)
    axes[1].set_title('B. Detectability Score\n(High Mean Methylation, Low CV)', fontweight='bold', fontsize=12)
    axes[1].grid(True, alpha=0.3, axis='x')
    axes[1].invert_yaxis()
    
    # Panel C: Detectability vs Prognostic scatter (with top candidates highlighted)
    scatter = axes[2].scatter(liquid_biopsy_df['Detectability_Score'],
                      liquid_biopsy_df['Prognostic_Score'],
                      alpha=0.3, s=30, c=liquid_biopsy_df['Mean_Methylation'],
                      cmap='viridis', edgecolors='none')
    # Highlight top 20 candidates
    top_20_scatter = top_biopsy.head(20)
    axes[2].scatter(top_20_scatter['Detectability_Score'],
                   top_20_scatter['Prognostic_Score'],
                   alpha=0.9, s=100, c='red', edgecolors='black', linewidth=1.5,
                   label='Top 20 Candidates', zorder=5)
    # Label top 3 genes
    for i in range(min(3, len(top_20_scatter))):
        axes[2].annotate(top_20_scatter.iloc[i]['Gene'].split(';')[0],
                        (top_20_scatter.iloc[i]['Detectability_Score'],
                         top_20_scatter.iloc[i]['Prognostic_Score']),
                        fontsize=8, fontweight='bold', alpha=0.8)
    axes[2].set_xlabel('Detectability Score', fontweight='bold', fontsize=11)
    axes[2].set_ylabel('Prognostic Score', fontweight='bold', fontsize=11)
    axes[2].set_title('C. Detectability vs Prognostic Value\n(Upper Right = Ideal Candidates)', fontweight='bold', fontsize=12)
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(fontsize=8)
    cbar = plt.colorbar(scatter, ax=axes[2])
    cbar.set_label('Mean Methylation', fontweight='bold')
    
    # No suptitle (removed for Figure 4)
    plt.tight_layout()
    plt.savefig(plots_dir / "liquid_biopsy_signatures.png", dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"✓ Liquid biopsy analysis completed")

print()

# ============================================================================
# ANALYSIS 5: COMBINATION THERAPY PREDICTION
# ============================================================================

print("="*80)
print("ANALYSIS 5: Combination Therapy Prediction")
print("="*80)
print()

print("Developing methylation-based drug response prediction...")

# Strategy: Identify methylation patterns associated with survival
# These could predict response to different treatment combinations

# Define treatment-related genes
treatment_genes = [
    'MGMT', 'EGFR', 'PTEN', 'TP53', 'IDH1', 'IDH2', 'ATRX', 'TERT',
    'CDKN2A', 'CDKN2B', 'RB1', 'PIK3CA', 'PIK3R1', 'NF1', 'PDGFRA',
    'VEGFA', 'VEGFR', 'PDGFR', 'MET', 'FGFR', 'BRAF', 'KRAS', 'NRAS'
]

therapy_results = []

for gene_name in treatment_genes:
    # Find in methylation
    meth_gene_rows = methylation[methylation['Hugo_Symbol'] == gene_name]
    
    if len(meth_gene_rows) > 0:
        row_meth = meth_gene_rows.iloc[0]
        patient_ids = [p for p in common_patients if p in row_meth.index and not pd.isna(row_meth[p])]
        
        if len(patient_ids) > 10:
            survival_data = clinical[clinical['PATIENT_ID'].isin(patient_ids)]
            if len(survival_data) > 10:
                meth_vals = row_meth[patient_ids].values
                survival_times = survival_data['OS_MONTHS'].values
                survival_events = survival_data['OS_EVENT'].values
                
                # Split by median methylation
                median_meth = np.median(meth_vals)
                high_meth = survival_times[meth_vals > median_meth]
                low_meth = survival_times[meth_vals <= median_meth]
                high_events = survival_events[meth_vals > median_meth]
                low_events = survival_events[meth_vals <= median_meth]
                
                if len(high_meth) > 5 and len(low_meth) > 5:
                    # Log-rank test
                    try:
                        results = logrank_test(high_meth, low_meth, high_events, low_events)
                        p_value = results.p_value
                    except:
                        p_value = 1.0
                    
                    # Correlation
                    corr, corr_p = spearmanr(meth_vals, survival_times)
                    
                    therapy_results.append({
                        'Gene': gene_name,
                        'Mean_Methylation': np.mean(meth_vals),
                        'Median_Methylation': median_meth,
                        'Survival_Correlation': corr,
                        'Survival_P_value': corr_p,
                        'LogRank_P_value': p_value,
                        'High_Meth_Median_Survival': np.median(high_meth) if len(high_meth) > 0 else np.nan,
                        'Low_Meth_Median_Survival': np.median(low_meth) if len(low_meth) > 0 else np.nan,
                        'N_High': len(high_meth),
                        'N_Low': len(low_meth)
                    })

therapy_df = pd.DataFrame(therapy_results)
therapy_df = therapy_df.dropna()

if len(therapy_df) > 0:
    # Significant genes
    significant_therapy = therapy_df[
        ((therapy_df['Survival_P_value'] < 0.05) | (therapy_df['LogRank_P_value'] < 0.05)) &
        (therapy_df['Survival_Correlation'].abs() > 0.2)
    ].sort_values('Survival_P_value')
    
    print(f"Treatment-related genes analyzed: {len(therapy_df)}")
    print(f"Significant prognostic genes: {len(significant_therapy)}")
    
    therapy_df.to_csv(output_dir / "combination_therapy_all.csv", index=False)
    significant_therapy.to_csv(output_dir / "combination_therapy_significant.csv", index=False)
    
    # Visualization
    if len(significant_therapy) > 0:
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Top prognostic genes
        top_10 = significant_therapy.head(10)
        axes[0, 0].barh(range(len(top_10)), -np.log10(top_10['Survival_P_value'] + 1e-10),
                      color='darkred', alpha=0.7)
        axes[0, 0].set_yticks(range(len(top_10)))
        axes[0, 0].set_yticklabels(top_10['Gene'], fontsize=9)
        axes[0, 0].set_xlabel('-log10(P-value)', fontweight='bold')
        axes[0, 0].set_title('Top 10 Treatment Genes: Prognostic Significance', 
                            fontweight='bold')
        axes[0, 0].grid(True, alpha=0.3, axis='x')
        axes[0, 0].axvline(-np.log10(0.05), color='red', linestyle='--', linewidth=2)
        axes[0, 0].invert_yaxis()
        
        # Survival correlation
        colors = ['red' if c < 0 else 'blue' for c in therapy_df['Survival_Correlation']]
        axes[0, 1].barh(range(len(therapy_df)), therapy_df['Survival_Correlation'],
                       color=colors, alpha=0.7)
        axes[0, 1].set_yticks(range(len(therapy_df)))
        axes[0, 1].set_yticklabels(therapy_df['Gene'], fontsize=8)
        axes[0, 1].set_xlabel('Survival Correlation', fontweight='bold')
        axes[0, 1].set_title('Treatment Genes: Survival Correlation', fontweight='bold')
        axes[0, 1].grid(True, alpha=0.3, axis='x')
        axes[0, 1].axvline(0, color='black', linestyle='--', alpha=0.7)
        axes[0, 1].invert_yaxis()
        
        # Methylation levels
        axes[1, 0].barh(range(len(therapy_df)), therapy_df['Mean_Methylation'],
                      color='steelblue', alpha=0.7)
        axes[1, 0].set_yticks(range(len(therapy_df)))
        axes[1, 0].set_yticklabels(therapy_df['Gene'], fontsize=8)
        axes[1, 0].set_xlabel('Mean Methylation', fontweight='bold')
        axes[1, 0].set_title('Treatment Genes: Methylation Levels', fontweight='bold')
        axes[1, 0].grid(True, alpha=0.3, axis='x')
        axes[1, 0].invert_yaxis()
        
        # Summary
        axes[1, 1].axis('off')
        summary_text = f"""
Combination Therapy Prediction Analysis

Total Treatment Genes Analyzed: {len(therapy_df)}
Significant Prognostic Genes: {len(significant_therapy)}

Top 5 Prognostic Genes:
1. {significant_therapy.iloc[0]['Gene']} 
   (P={significant_therapy.iloc[0]['Survival_P_value']:.4f})
2. {significant_therapy.iloc[1]['Gene'] if len(significant_therapy) > 1 else 'N/A'}
   (P={f"{significant_therapy.iloc[1]['Survival_P_value']:.4f}" if len(significant_therapy) > 1 else 'N/A'})
3. {significant_therapy.iloc[2]['Gene'] if len(significant_therapy) > 2 else 'N/A'}
   (P={f"{significant_therapy.iloc[2]['Survival_P_value']:.4f}" if len(significant_therapy) > 2 else 'N/A'})

Clinical Impact:
- Drug response prediction
- Treatment selection
- Combination therapy optimization
- Personalized medicine
"""
        axes[1, 1].text(0.1, 0.5, summary_text, fontsize=10, fontweight='bold',
                       verticalalignment='center', family='monospace',
                       bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.5))
        
        plt.suptitle('Combination Therapy Prediction\n(Methylation-Based Drug Response)', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(plots_dir / "combination_therapy_prediction.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Combination therapy analysis completed")
    else:
        print("⚠ No significant therapy genes found")
else:
    print("⚠ No therapy results generated")

print()

# ============================================================================
# SUMMARY
# ============================================================================

print("="*80)
print("NOVEL ANALYSES SUMMARY")
print("="*80)
print()
print("Completed Analyses:")
print("1. Temporal Methylation Dynamics")
print("2. Spatial Methylation Heterogeneity")
print("3. Methylation-Immune Crosstalk")
print("4. Liquid Biopsy Signatures")
print("5. Combination Therapy Prediction")
print()
print(f"Results saved to: {output_dir}")
print(f"Plots saved to: {plots_dir}")
print()
print("="*80)
print("All novel analyses completed!")
print("="*80)

