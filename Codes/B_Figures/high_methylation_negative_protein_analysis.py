#!/usr/bin/env python3
"""
High Methylation - Negative Protein Impact Analysis
Yüksek metilasyonlu genlerden hangilerinin protein miktarını negatif etkilediğini bulma
Biomarker potansiyeli olan metilasyon profillerini belirleme
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from repo_paths import REPO_ROOT, data_file


# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

# Create output directories
output_dir = Path("results")
plots_dir = Path("plots")
output_dir.mkdir(exist_ok=True)
plots_dir.mkdir(exist_ok=True)

print("="*80)
print("High Methylation - Negative Protein Impact Analysis")
print("="*80)
print()

# ============================================================================
# 1. DATA LOADING
# ============================================================================

print("Step 1: Loading data...")
print("-" * 80)

# Load clinical data
clinical_file = str(data_file('clinical_patient.txt'))
clinical = pd.read_csv(clinical_file, sep="\t", skiprows=4)
clinical.columns = clinical.columns.str.strip()
clinical['PATIENT_ID'] = clinical['PATIENT_ID'].str.strip()
# Fix survival data: use "Deceased" (case-insensitive) and proper time calculation
clinical['VITAL_STATUS'] = clinical['VITAL_STATUS'].str.strip()
clinical['OS_EVENT'] = (clinical['VITAL_STATUS'].str.lower() == 'deceased').astype(int)
# Use time-to-death for deceased, last-contact for living (censored)
death_days = pd.to_numeric(clinical['PATH_DIAG_TO_DEATH_DAYS'], errors='coerce')
contact_days = pd.to_numeric(clinical['PATH_DIAG_TO_LAST_CONTACT_DAYS'], errors='coerce')
clinical['OS_DAYS'] = np.where(clinical['OS_EVENT'] == 1, death_days, contact_days)
clinical = clinical.dropna(subset=['OS_DAYS', 'OS_EVENT'])
clinical = clinical[clinical['OS_DAYS'] > 0]

# Load mRNA data
mrna_file = str(data_file("mrna_seq_fpkm.txt"))
mrna = pd.read_csv(mrna_file, sep="\t")
mrna.columns = mrna.columns.str.strip()

# Load protein data
protein_file = str(data_file('protein_quantification.txt'))
protein = pd.read_csv(protein_file, sep="\t")
protein.columns = protein.columns.str.strip()

# Load methylation data
print("Loading methylation data...")
methylation_file = str(data_file('methylation_epic.txt'))
methylation = pd.read_csv(methylation_file, sep="\t", nrows=20000)  # Increased for better coverage
print(f"Methylation data loaded: {methylation.shape}")

# Get common patients
clinical_patients = set(clinical['PATIENT_ID'].unique())
mrna_patients = set([col for col in mrna.columns if col.startswith('C3')])
protein_patients = set([col for col in protein.columns if col.startswith('C3')])
methylation_patients = set([col for col in methylation.columns if col.startswith('C3')])

common_patients = sorted(list(clinical_patients & mrna_patients & protein_patients & methylation_patients))
print(f"Common patients: {len(common_patients)}")

# Filter datasets
clinical_filtered = clinical[clinical['PATIENT_ID'].isin(common_patients)].copy()
mrna_filtered = mrna[['Hugo_Symbol'] + [p for p in common_patients if p in mrna.columns]].copy()
protein_filtered = protein[[col for col in protein.columns if col == 'Composite.Element.REF' or col in common_patients]].copy()
methylation_filtered = methylation[[col for col in methylation.columns if col == 'Hugo_Symbol' or col in common_patients]].copy()

print("Data loading completed!")
print()

# ============================================================================
# 2. GENE MAPPING AND NORMALIZATION
# ============================================================================

print("Step 2: Gene mapping and normalization...")
print("-" * 80)

# Extract genes
mrna_genes = set(mrna_filtered['Hugo_Symbol'].unique())
protein_genes = set(protein_filtered['Composite.Element.REF'].str.split('|').str[0].unique())
meth_genes = set(methylation_filtered['Hugo_Symbol'].dropna().unique())

common_genes = sorted(list(mrna_genes & protein_genes & meth_genes))
print(f"Common genes across all datasets: {len(common_genes)}")

# Normalize data
mrna_values = mrna_filtered[common_patients].values
mrna_values = np.nan_to_num(mrna_values, nan=0.0, posinf=0.0, neginf=0.0)
mrna_normalized = np.log2(mrna_values + 1)
mrna_df = pd.DataFrame(mrna_normalized, index=mrna_filtered['Hugo_Symbol'], columns=common_patients)
mrna_df = mrna_df.fillna(0)

protein_values = protein_filtered[common_patients].values
protein_values = np.nan_to_num(protein_values, nan=0.0, posinf=0.0, neginf=0.0)
protein_df = pd.DataFrame(protein_values, index=protein_filtered['Composite.Element.REF'].str.split('|').str[0], columns=common_patients)
protein_df = protein_df.fillna(0)

meth_values = methylation_filtered[common_patients].values
meth_values = np.nan_to_num(meth_values, nan=0.5, posinf=0.5, neginf=0.5)
meth_df = pd.DataFrame(meth_values, index=methylation_filtered['Hugo_Symbol'], columns=common_patients)
meth_df = meth_df.fillna(0.5)

print("Normalization completed!")
print()

# ============================================================================
# 3. HIGH METHYLATION - NEGATIVE PROTEIN IMPACT ANALYSIS
# ============================================================================

print("Step 3: Identifying high methylation genes with negative protein impact...")
print("-" * 80)

results = []

for gene in common_genes:
    if gene not in meth_df.index or gene not in protein_df.index:
        continue
    
    try:
        # Get data
        meth_data = meth_df.loc[gene, common_patients].values
        protein_data = protein_df.loc[gene, common_patients].values
        
        # Remove NaN
        valid_idx = ~(np.isnan(meth_data) | np.isnan(protein_data))
        if valid_idx.sum() < 10:  # At least 10 patients
            continue
        
        meth_data = meth_data[valid_idx]
        protein_data = protein_data[valid_idx]
        
        # Calculate statistics
        meth_mean = np.mean(meth_data)
        meth_median = np.median(meth_data)
        meth_std = np.std(meth_data)
        protein_mean = np.mean(protein_data)
        protein_median = np.median(protein_data)
        
        # Calculate correlation
        try:
            corr, p_val = spearmanr(meth_data, protein_data)
        except:
            corr, p_val = np.nan, np.nan
        
        # Define high methylation (beta > 0.7)
        high_meth_threshold = 0.7
        is_high_methylated = meth_mean > high_meth_threshold
        
        # Define negative impact (negative correlation)
        is_negative_impact = corr < 0
        
        # Calculate high methylation patients
        high_meth_patients = (meth_data > high_meth_threshold).sum()
        high_meth_percentage = (high_meth_patients / len(meth_data)) * 100
        
        # Calculate protein difference (high meth vs low meth)
        high_meth_protein = protein_data[meth_data > high_meth_threshold]
        low_meth_protein = protein_data[meth_data <= high_meth_threshold]
        
        if len(high_meth_protein) > 0 and len(low_meth_protein) > 0:
            protein_diff = np.mean(high_meth_protein) - np.mean(low_meth_protein)
        else:
            protein_diff = np.nan
        
        results.append({
            'Gene': gene,
            'Methylation_Mean': meth_mean,
            'Methylation_Median': meth_median,
            'Methylation_STD': meth_std,
            'Protein_Mean': protein_mean,
            'Protein_Median': protein_median,
            'Meth_Protein_Corr': corr,
            'Meth_Protein_P': p_val,
            'Is_High_Methylated': is_high_methylated,
            'Is_Negative_Impact': is_negative_impact,
            'High_Meth_Patients': high_meth_patients,
            'High_Meth_Percentage': high_meth_percentage,
            'Protein_Diff_High_vs_Low_Meth': protein_diff,
            'High_Meth_Protein_Mean': np.mean(high_meth_protein) if len(high_meth_protein) > 0 else np.nan,
            'Low_Meth_Protein_Mean': np.mean(low_meth_protein) if len(low_meth_protein) > 0 else np.nan
        })
    except Exception as e:
        continue

results_df = pd.DataFrame(results)

# Filter for high methylation + negative impact
high_meth_negative = results_df[
    (results_df['Is_High_Methylated'] == True) & 
    (results_df['Is_Negative_Impact'] == True) &
    (results_df['Meth_Protein_P'] < 0.05)  # Significant correlation
].copy()

# Sort by correlation (most negative first)
high_meth_negative = high_meth_negative.sort_values('Meth_Protein_Corr')

print(f"Total genes analyzed: {len(results_df)}")
print(f"High methylation genes (mean > 0.7): {results_df['Is_High_Methylated'].sum()}")
print(f"Negative impact genes (corr < 0): {results_df['Is_Negative_Impact'].sum()}")
print(f"High methylation + Negative impact + Significant: {len(high_meth_negative)}")
print()

# Save all results
results_df.to_csv(output_dir / "methylation_protein_impact_analysis.csv", index=False)
high_meth_negative.to_csv(output_dir / "high_methylation_negative_protein_genes.csv", index=False)

print("Top 20 High Methylation Genes with Negative Protein Impact:")
print("="*80)
print(high_meth_negative[['Gene', 'Methylation_Mean', 'Meth_Protein_Corr', 'Meth_Protein_P', 
                          'Protein_Diff_High_vs_Low_Meth', 'High_Meth_Percentage']].head(20).to_string(index=False))
print()

# ============================================================================
# 4. BIOMARKER POTENTIAL ANALYSIS
# ============================================================================

print("Step 4: Biomarker Potential Analysis...")
print("-" * 80)

# Calculate biomarker potential score
# Factors: 1) Strong negative correlation, 2) High methylation, 3) Significant p-value, 4) Large protein difference

high_meth_negative['Biomarker_Score'] = (
    -high_meth_negative['Meth_Protein_Corr'] * 0.4 +  # Stronger negative correlation = higher score
    high_meth_negative['Methylation_Mean'] * 0.2 +     # Higher methylation = higher score
    (-np.log10(high_meth_negative['Meth_Protein_P'] + 1e-10)) * 0.2 +  # More significant = higher score
    abs(high_meth_negative['Protein_Diff_High_vs_Low_Meth']) * 0.2  # Larger difference = higher score
)

# Sort by biomarker score
high_meth_negative = high_meth_negative.sort_values('Biomarker_Score', ascending=False)

print("Top 20 Biomarker Candidates (High Methylation → Negative Protein Impact):")
print("="*80)
biomarker_candidates = high_meth_negative[['Gene', 'Methylation_Mean', 'Meth_Protein_Corr', 
                                           'Meth_Protein_P', 'Protein_Diff_High_vs_Low_Meth',
                                           'High_Meth_Percentage', 'Biomarker_Score']].head(20)
print(biomarker_candidates.to_string(index=False))
print()

# Save biomarker candidates
biomarker_candidates.to_csv(output_dir / "biomarker_candidates_high_meth_negative_protein.csv", index=False)

# ============================================================================
# 5. SURVIVAL ANALYSIS FOR TOP BIOMARKER CANDIDATES
# ============================================================================

print("Step 5: Survival Analysis for Top Biomarker Candidates...")
print("-" * 80)

# Select top 10 biomarker candidates
top_biomarkers = high_meth_negative.head(10)['Gene'].tolist()

survival_results = []

for gene in top_biomarkers:
    try:
        # Get methylation data
        meth_data = meth_df.loc[gene, common_patients]
        
        # Classify patients: High methylation (beta > 0.7) vs Low methylation
        high_meth_patients = meth_data[meth_data > 0.7].index.tolist()
        low_meth_patients = meth_data[meth_data <= 0.7].index.tolist()
        
        if len(high_meth_patients) < 5 or len(low_meth_patients) < 5:
            continue
        
        # Get survival data
        high_meth_survival = clinical_filtered[clinical_filtered['PATIENT_ID'].isin(high_meth_patients)][['OS_DAYS', 'OS_EVENT']].dropna()
        low_meth_survival = clinical_filtered[clinical_filtered['PATIENT_ID'].isin(low_meth_patients)][['OS_DAYS', 'OS_EVENT']].dropna()
        
        high_meth_survival = high_meth_survival[high_meth_survival['OS_DAYS'] > 0]
        low_meth_survival = low_meth_survival[low_meth_survival['OS_DAYS'] > 0]
        
        if len(high_meth_survival) < 5 or len(low_meth_survival) < 5:
            continue
        
        # Calculate median survival
        high_meth_median = high_meth_survival['OS_DAYS'].median()
        low_meth_median = low_meth_survival['OS_DAYS'].median()
        
        # Log-rank test
        from lifelines.statistics import logrank_test
        try:
            lr_result = logrank_test(high_meth_survival['OS_DAYS'], low_meth_survival['OS_DAYS'],
                                     high_meth_survival['OS_EVENT'], low_meth_survival['OS_EVENT'])
            p_value = lr_result.p_value
        except:
            p_value = np.nan
        
        survival_results.append({
            'Gene': gene,
            'High_Meth_Count': len(high_meth_patients),
            'Low_Meth_Count': len(low_meth_patients),
            'High_Meth_Median_Survival': high_meth_median,
            'Low_Meth_Median_Survival': low_meth_median,
            'Survival_Diff': high_meth_median - low_meth_median,
            'LogRank_P': p_value,
            'Is_Significant': p_value < 0.05 if not np.isnan(p_value) else False
        })
    except Exception as e:
        continue

survival_df = pd.DataFrame(survival_results)
survival_df = survival_df.sort_values('LogRank_P')

print("Survival Analysis Results for Top Biomarker Candidates:")
print("="*80)
print(survival_df.to_string(index=False))
print()

# Save survival results
survival_df.to_csv(output_dir / "biomarker_survival_analysis.csv", index=False)

# ============================================================================
# 6. VISUALIZATION
# ============================================================================

print("Step 6: Creating visualizations...")
print("-" * 80)

# Plot 1: Top biomarker candidates
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: Methylation vs Protein Correlation for top candidates
top_20 = high_meth_negative.head(20)
scatter = axes[0, 0].scatter(top_20['Methylation_Mean'], top_20['Meth_Protein_Corr'], 
                   s=100, alpha=0.7, c=-np.log10(top_20['Meth_Protein_P'] + 1e-10), cmap='viridis', edgecolors='black', linewidth=0.5)
axes[0, 0].set_xlabel('Mean Methylation (Beta)', fontweight='bold')
axes[0, 0].set_ylabel('Methylation-Protein Correlation', fontweight='bold')
axes[0, 0].set_title('A. Top 20 Biomarker Candidates: Methylation vs Correlation', fontweight='bold')
axes[0, 0].axhline(0, color='red', linestyle='--', alpha=0.7, linewidth=1)
axes[0, 0].axvline(0.7, color='red', linestyle='--', alpha=0.7, linewidth=1, label='High Meth Threshold')
# Label top genes
for i, (idx, row) in enumerate(top_20.iterrows()):
    if row['Gene'] in ['EGFR', 'B4GALT1', 'CTSZ', 'VASP', 'PHF11', 'TSPYL5']:
        axes[0, 0].annotate(row['Gene'], (row['Methylation_Mean'], row['Meth_Protein_Corr']),
                           fontsize=8, alpha=0.8, fontweight='bold')
axes[0, 0].legend(fontsize=8)
axes[0, 0].grid(True, alpha=0.3)
cbar = plt.colorbar(scatter, ax=axes[0, 0])
cbar.set_label('-log10(P-value)', fontweight='bold')

# Plot 2: Protein difference (High vs Low Meth) - color by sign
top_20_sorted = top_20.sort_values('Protein_Diff_High_vs_Low_Meth', ascending=True)
colors_diff = ['red' if x < 0 else 'blue' for x in top_20_sorted['Protein_Diff_High_vs_Low_Meth']]
axes[0, 1].barh(range(len(top_20_sorted)), top_20_sorted['Protein_Diff_High_vs_Low_Meth'], color=colors_diff, alpha=0.7)
axes[0, 1].set_yticks(range(len(top_20_sorted)))
axes[0, 1].set_yticklabels(top_20_sorted['Gene'], fontsize=8)
axes[0, 1].set_xlabel('Protein Difference (High Meth - Low Meth)', fontweight='bold')
axes[0, 1].set_title('B. Protein Level Difference: High vs Low Methylation', fontweight='bold')
axes[0, 1].axvline(0, color='black', linestyle='-', alpha=0.5, linewidth=1)
axes[0, 1].grid(True, alpha=0.3, axis='x')
# Add value labels for clarity
for i, (idx, row) in enumerate(top_20_sorted.iterrows()):
    val = row['Protein_Diff_High_vs_Low_Meth']
    axes[0, 1].text(val + (0.05 if val >= 0 else -0.05), i, f'{val:.3f}', 
                    va='center', ha='left' if val >= 0 else 'right', fontsize=7)

# Plot 3: Biomarker score distribution
top_20_score_sorted = top_20.sort_values('Biomarker_Score', ascending=True)
axes[1, 0].barh(range(len(top_20_score_sorted)), top_20_score_sorted['Biomarker_Score'], color='steelblue', alpha=0.7, edgecolor='black', linewidth=0.5)
axes[1, 0].set_yticks(range(len(top_20_score_sorted)))
axes[1, 0].set_yticklabels(top_20_score_sorted['Gene'], fontsize=8)
axes[1, 0].set_xlabel('Biomarker Score', fontweight='bold')
axes[1, 0].set_title('C. Biomarker Potential Score (Top 20)', fontweight='bold')
axes[1, 0].grid(True, alpha=0.3, axis='x')
# Highlight survival-significant genes
if len(survival_df) > 0:
    sig_genes_set = set(survival_df[survival_df['Is_Significant'] == True]['Gene'].tolist())
    for i, (idx, row) in enumerate(top_20_score_sorted.iterrows()):
        if row['Gene'] in sig_genes_set:
            axes[1, 0].barh(i, row['Biomarker_Score'], color='green', alpha=0.5)

# Plot 4: Survival analysis results (only significant genes, sorted by p-value)
if len(survival_df) > 0:
    sig_genes = survival_df[survival_df['Is_Significant'] == True].sort_values('LogRank_P', ascending=True)
    if len(sig_genes) > 0:
        axes[1, 1].barh(range(len(sig_genes)), -np.log10(sig_genes['LogRank_P'] + 1e-10), 
                       color='green', alpha=0.7, edgecolor='black', linewidth=0.5)
        axes[1, 1].set_yticks(range(len(sig_genes)))
        axes[1, 1].set_yticklabels(sig_genes['Gene'], fontsize=8)
        axes[1, 1].set_xlabel('-log10(P-value)', fontweight='bold')
        axes[1, 1].set_title('D. Survival Significance (Log-Rank Test)', fontweight='bold')
        axes[1, 1].axvline(-np.log10(0.05), color='red', linestyle='--', linewidth=1, label='P=0.05')
        axes[1, 1].axvline(-np.log10(0.01), color='orange', linestyle='--', linewidth=1, label='P=0.01')
        axes[1, 1].legend(fontsize=8)
        axes[1, 1].grid(True, alpha=0.3, axis='x')
        # Add p-value labels
        for i, (idx, row) in enumerate(sig_genes.iterrows()):
            pval = row['LogRank_P']
            axes[1, 1].text(-np.log10(pval + 1e-10) + 0.1, i, f'p={pval:.4f}', 
                           va='center', fontsize=7)
    else:
        axes[1, 1].text(0.5, 0.5, 'No significant survival differences found', 
                       ha='center', va='center', transform=axes[1, 1].transAxes, fontsize=11)
        axes[1, 1].set_title('D. Survival Analysis', fontweight='bold')
else:
    axes[1, 1].text(0.5, 0.5, 'No survival data available', 
                   ha='center', va='center', transform=axes[1, 1].transAxes, fontsize=11)
    axes[1, 1].set_title('D. Survival Analysis', fontweight='bold')

plt.suptitle('High-Methylation Negative-Protein Biomarkers with Clinical Impact', 
            fontsize=14, fontweight='bold', y=0.995)
plt.tight_layout()
plt.savefig(plots_dir / "high_methylation_negative_protein_biomarkers.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print(f"Visualization saved to {plots_dir / 'high_methylation_negative_protein_biomarkers.png'}")

# ============================================================================
# SUMMARY
# ============================================================================

print()
print("="*80)
print("Analysis Summary")
print("="*80)
print()
print(f"Total genes analyzed: {len(results_df)}")
print(f"High methylation genes (mean > 0.7): {results_df['Is_High_Methylated'].sum()}")
print(f"High methylation + Negative impact + Significant: {len(high_meth_negative)}")
print()
print("Top 10 Biomarker Candidates:")
print("="*80)
for i, row in high_meth_negative.head(10).iterrows():
    print(f"{i+1}. {row['Gene']}: Meth={row['Methylation_Mean']:.3f}, "
          f"Corr={row['Meth_Protein_Corr']:.3f}, P={row['Meth_Protein_P']:.2e}, "
          f"Score={row['Biomarker_Score']:.3f}")
print()
print("Files saved:")
print(f"- results/methylation_protein_impact_analysis.csv")
print(f"- results/high_methylation_negative_protein_genes.csv")
print(f"- results/biomarker_candidates_high_meth_negative_protein.csv")
print(f"- results/biomarker_survival_analysis.csv")
print(f"- plots/high_methylation_negative_protein_biomarkers.png")
print()
print("="*80)
print("Analysis completed successfully!")
print("="*80)











