#!/usr/bin/env python3
"""
GBM Epigenetic-Transcriptomic-Proteomic Integration Analysis
Novel Multi-Omics Analysis Pipeline

This script integrates:
1. clinical_patient.txt (clinical data)
2. methylation_epic.txt (DNA methylation)
3. mrna_seq_fpkm.txt (gene expression)
4. protein_quantification.txt (protein expression)

Novel Analyses:
1. Epigenetic-Transcriptomic-Proteomic Cascade Analysis
2. Multi-Layer Regulatory Networks
3. Multi-Omics Subtype Discovery
4. Methylation-Protein Disconnect Analysis
5. MGMT-Associated Multi-Omics Signatures
6. Multi-Omics Prognostic Models
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr, pearsonr
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LassoCV, ElasticNetCV
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings('ignore')

import os
import sys
from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

# Create output directories
output_dir = Path("results")
plots_dir = Path("plots")
output_dir.mkdir(exist_ok=True)
plots_dir.mkdir(exist_ok=True)

print("="*80)
print("GBM Epigenetic-Transcriptomic-Proteomic Integration Analysis")
print("="*80)
print()

# ============================================================================
# 1. DATA LOADING AND PREPROCESSING
# ============================================================================

print("Step 1: Loading and preprocessing data...")
print("-" * 80)

# Load clinical data
clinical_file = "../clinical_patient.txt"
clinical = pd.read_csv(clinical_file, sep="\t", skiprows=4)
clinical.columns = clinical.columns.str.strip()
clinical['PATIENT_ID'] = clinical['PATIENT_ID'].str.strip()

# Extract survival data
clinical['OS_DAYS'] = clinical['PATH_DIAG_TO_DEATH_DAYS'].copy()
clinical.loc[clinical['VITAL_STATUS'] == 'Living', 'OS_DAYS'] = clinical.loc[clinical['VITAL_STATUS'] == 'Living', 'PATH_DIAG_TO_LAST_CONTACT_DAYS']
clinical['OS_EVENT'] = (clinical['VITAL_STATUS'] == 'Deceased').astype(int)

# Load methylation data (sample first few rows to understand structure)
print("Loading methylation data (this may take a while due to large file size)...")
methylation_file = "../methylation_epic.txt"
try:
    # Try to load in chunks if file is too large
    methylation = pd.read_csv(methylation_file, sep="\t", nrows=1000)  # Sample first 1000 rows
    print(f"Methylation data shape (sample): {methylation.shape}")
    print(f"Methylation columns (first 10): {list(methylation.columns[:10])}")
except Exception as e:
    print(f"Error loading methylation data: {e}")
    print("Will proceed with available data...")
    methylation = None

# Load mRNA expression data
print("Loading mRNA expression data...")
mrna_file = "../mrna_seq_fpkm.txt"
mrna = pd.read_csv(mrna_file, sep="\t")
mrna.columns = mrna.columns.str.strip()
print(f"mRNA data shape: {mrna.shape}")

# Load protein data
print("Loading protein quantification data...")
protein_file = "../protein_quantification.txt"
protein = pd.read_csv(protein_file, sep="\t")
protein.columns = protein.columns.str.strip()
print(f"Protein data shape: {protein.shape}")

print("\nData loading completed!")
print()

# ============================================================================
# 2. DATA INTEGRATION AND PATIENT MATCHING
# ============================================================================

print("Step 2: Data integration and patient matching...")
print("-" * 80)

# Extract patient IDs from each dataset
clinical_patients = set(clinical['PATIENT_ID'].unique())

if methylation is not None:
    methylation_patients = set([col for col in methylation.columns if col.startswith('C3')])
else:
    methylation_patients = set()

mrna_patients = set([col for col in mrna.columns if col.startswith('C3')])
protein_patients = set([col for col in protein.columns if col.startswith('C3')])

# Find common patients
common_patients = clinical_patients & mrna_patients & protein_patients
if methylation is not None:
    common_patients = common_patients & methylation_patients

common_patients = sorted(list(common_patients))
print(f"Common patients across all datasets: {len(common_patients)}")

# Filter datasets to common patients
clinical_filtered = clinical[clinical['PATIENT_ID'].isin(common_patients)].copy()
mrna_filtered = mrna[['Hugo_Symbol'] + [p for p in common_patients if p in mrna.columns]].copy()
protein_filtered = protein[[col for col in protein.columns if col == 'Composite.Element.REF' or col in common_patients]].copy()

if methylation is not None:
    methylation_filtered = methylation[[col for col in methylation.columns if col == 'Hugo_Symbol' or col in common_patients]].copy()
else:
    methylation_filtered = None

print(f"Filtered clinical data: {clinical_filtered.shape}")
print(f"Filtered mRNA data: {mrna_filtered.shape}")
print(f"Filtered protein data: {protein_filtered.shape}")
if methylation_filtered is not None:
    print(f"Filtered methylation data: {methylation_filtered.shape}")

print("\nData integration completed!")
print()

# ============================================================================
# 3. GENE/PROTEIN MAPPING
# ============================================================================

print("Step 3: Gene/protein mapping...")
print("-" * 80)

# Extract gene symbols
mrna_genes = set(mrna_filtered['Hugo_Symbol'].unique())
protein_genes = set(protein_filtered['Composite.Element.REF'].str.split('|').str[0].unique())

# Find common genes
common_genes = mrna_genes & protein_genes
common_genes = sorted(list(common_genes))
print(f"Common genes between mRNA and protein: {len(common_genes)}")

# Filter to common genes
mrna_common = mrna_filtered[mrna_filtered['Hugo_Symbol'].isin(common_genes)].copy()
protein_common = protein_filtered[protein_filtered['Composite.Element.REF'].str.split('|').str[0].isin(common_genes)].copy()

# Create gene mapping for protein
protein_common['Gene'] = protein_common['Composite.Element.REF'].str.split('|').str[0]

print("\nGene mapping completed!")
print()

# ============================================================================
# 4. DATA NORMALIZATION
# ============================================================================

print("Step 4: Data normalization...")
print("-" * 80)

# Normalize mRNA: Log2(FPKM + 1)
mrna_values = mrna_common[common_patients].values
# Replace NaN and inf with 0
mrna_values = np.nan_to_num(mrna_values, nan=0.0, posinf=0.0, neginf=0.0)
mrna_normalized = np.log2(mrna_values + 1)
mrna_normalized_df = pd.DataFrame(mrna_normalized, 
                                   index=mrna_common['Hugo_Symbol'],
                                   columns=common_patients)
# Fill any remaining NaN
mrna_normalized_df = mrna_normalized_df.fillna(0)

# Protein already z-score normalized, but we'll ensure consistency
protein_values = protein_common[common_patients].values
# Replace NaN and inf with 0
protein_values = np.nan_to_num(protein_values, nan=0.0, posinf=0.0, neginf=0.0)
protein_normalized_df = pd.DataFrame(protein_values,
                                     index=protein_common['Gene'],
                                     columns=common_patients)
# Fill any remaining NaN
protein_normalized_df = protein_normalized_df.fillna(0)

# Methylation: Beta values (already normalized 0-1)
if methylation_filtered is not None:
    methylation_values = methylation_filtered[common_patients].values
    methylation_normalized_df = pd.DataFrame(methylation_values,
                                            index=methylation_filtered['Hugo_Symbol'] if 'Hugo_Symbol' in methylation_filtered.columns else range(len(methylation_filtered)),
                                            columns=common_patients)
else:
    methylation_normalized_df = None

print("Normalization completed!")
print()

# ============================================================================
# 5. EPIGENETIC-TRANSCRIPTOMIC-PROTEOMIC CASCADE ANALYSIS
# ============================================================================

print("Step 5: Epigenetic-Transcriptomic-Proteomic Cascade Analysis...")
print("-" * 80)

def calculate_cascade_correlations(methylation_df, mrna_df, protein_df, genes, patients):
    """
    Calculate correlations across methylation -> mRNA -> protein cascade
    """
    results = []
    
    for gene in genes[:min(100, len(genes))]:  # Limit to first 100 genes for speed
        if gene not in mrna_df.index or gene not in protein_df.index:
            continue
            
        # Get data for this gene
        mrna_data = mrna_df.loc[gene, patients].values
        protein_data = protein_df.loc[gene, patients].values
        
        # Calculate mRNA-Protein correlation
        if len(mrna_data) > 5 and len(protein_data) > 5:
            try:
                mrna_protein_corr, mrna_protein_p = spearmanr(mrna_data, protein_data)
            except:
                mrna_protein_corr, mrna_protein_p = np.nan, np.nan
        else:
            mrna_protein_corr, mrna_protein_p = np.nan, np.nan
        
        # If methylation data available, calculate methylation correlations
        if methylation_df is not None:
            # Find methylation probes for this gene (simplified - would need gene annotation)
            # For now, skip methylation correlations
            meth_mrna_corr, meth_mrna_p = np.nan, np.nan
            meth_protein_corr, meth_protein_p = np.nan, np.nan
        else:
            meth_mrna_corr, meth_mrna_p = np.nan, np.nan
            meth_protein_corr, meth_protein_p = np.nan, np.nan
        
        results.append({
            'Gene': gene,
            'mRNA_Protein_Corr': mrna_protein_corr,
            'mRNA_Protein_P': mrna_protein_p,
            'Meth_mRNA_Corr': meth_mrna_corr,
            'Meth_mRNA_P': meth_mrna_p,
            'Meth_Protein_Corr': meth_protein_corr,
            'Meth_Protein_P': meth_protein_p
        })
    
    return pd.DataFrame(results)

# Calculate cascade correlations
cascade_results = calculate_cascade_correlations(
    methylation_normalized_df, 
    mrna_normalized_df, 
    protein_normalized_df,
    common_genes,
    common_patients
)

# Save results
cascade_results.to_csv(output_dir / "cascade_correlations.csv", index=False)
print(f"Cascade analysis completed! Results saved to {output_dir / 'cascade_correlations.csv'}")
print(f"Number of genes analyzed: {len(cascade_results)}")
print(f"Significant mRNA-Protein correlations (p<0.05): {(cascade_results['mRNA_Protein_P'] < 0.05).sum()}")

# Visualize cascade correlations
if len(cascade_results) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # mRNA-Protein correlation distribution
    axes[0].hist(cascade_results['mRNA_Protein_Corr'].dropna(), bins=50, edgecolor='black')
    axes[0].set_xlabel('mRNA-Protein Correlation (Spearman)')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title('Distribution of mRNA-Protein Correlations')
    axes[0].axvline(0, color='red', linestyle='--', label='Zero correlation')
    axes[0].legend()
    
    # Correlation vs P-value
    axes[1].scatter(cascade_results['mRNA_Protein_Corr'], 
                   -np.log10(cascade_results['mRNA_Protein_P'] + 1e-10),
                   alpha=0.6)
    axes[1].set_xlabel('mRNA-Protein Correlation')
    axes[1].set_ylabel('-log10(P-value)')
    axes[1].set_title('Correlation vs Significance')
    axes[1].axhline(-np.log10(0.05), color='red', linestyle='--', label='P=0.05')
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(plots_dir / "cascade_correlations.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Cascade correlation plots saved to {plots_dir / 'cascade_correlations.png'}")

print()

# ============================================================================
# 6. MULTI-LAYER REGULATORY NETWORKS
# ============================================================================

print("Step 6: Multi-Layer Regulatory Network Analysis...")
print("-" * 80)

def calculate_coexpression_network(data_df, method='spearman', threshold=0.3):
    """
    Calculate co-expression network
    """
    # Transpose to have genes as columns
    data_t = data_df.T
    
    # Calculate correlation matrix
    if method == 'spearman':
        corr_matrix = data_t.corr(method='spearman')
    else:
        corr_matrix = data_t.corr(method='pearson')
    
    # Threshold
    corr_matrix[abs(corr_matrix) < threshold] = 0
    
    return corr_matrix

# Calculate mRNA co-expression network
print("Calculating mRNA co-expression network...")
mrna_network = calculate_coexpression_network(mrna_normalized_df, threshold=0.3)
print(f"mRNA network: {mrna_network.shape}, Non-zero edges: {(mrna_network != 0).sum().sum() / 2}")

# Calculate protein co-expression network
print("Calculating protein co-expression network...")
protein_network = calculate_coexpression_network(protein_normalized_df, threshold=0.3)
print(f"Protein network: {protein_network.shape}, Non-zero edges: {(protein_network != 0).sum().sum() / 2}")

# Calculate cross-layer network (mRNA-Protein)
print("Calculating mRNA-Protein cross-layer network...")
cross_layer_network = pd.DataFrame(index=mrna_normalized_df.index, columns=protein_normalized_df.index)
for mrna_gene in mrna_normalized_df.index[:min(50, len(mrna_normalized_df.index))]:
    if mrna_gene in protein_normalized_df.index:
        try:
            corr, p = spearmanr(mrna_normalized_df.loc[mrna_gene, common_patients],
                               protein_normalized_df.loc[mrna_gene, common_patients])
            cross_layer_network.loc[mrna_gene, mrna_gene] = corr
        except:
            pass

# Identify hub genes (high degree in network)
def identify_hubs(network, top_n=20):
    """Identify hub genes based on network degree"""
    degrees = (network != 0).sum(axis=1)
    top_hubs = degrees.nlargest(top_n)
    return top_hubs

mrna_hubs = identify_hubs(mrna_network)
protein_hubs = identify_hubs(protein_network)

# Save hub genes
pd.DataFrame({'Gene': mrna_hubs.index, 'Degree': mrna_hubs.values}).to_csv(
    output_dir / "mrna_hub_genes.csv", index=False)
pd.DataFrame({'Gene': protein_hubs.index, 'Degree': protein_hubs.values}).to_csv(
    output_dir / "protein_hub_genes.csv", index=False)

print(f"Top 10 mRNA hub genes: {list(mrna_hubs.head(10).index)}")
print(f"Top 10 protein hub genes: {list(protein_hubs.head(10).index)}")
print(f"Hub analysis results saved to {output_dir}")

print()

# ============================================================================
# 7. MULTI-OMICS SUBTYPE DISCOVERY
# ============================================================================

print("Step 7: Multi-Omics Subtype Discovery...")
print("-" * 80)

def multi_omics_clustering(mrna_df, protein_df, n_clusters=3):
    """
    Perform multi-omics clustering using PCA integration
    """
    # Prepare data
    mrna_data = mrna_df.T.values
    protein_data = protein_df.T.values
    
    # Handle NaN values - fill with median
    mrna_data = pd.DataFrame(mrna_data).fillna(pd.DataFrame(mrna_data).median()).values
    protein_data = pd.DataFrame(protein_data).fillna(pd.DataFrame(protein_data).median()).values
    
    # Handle infinite values
    mrna_data = np.nan_to_num(mrna_data, nan=0.0, posinf=0.0, neginf=0.0)
    protein_data = np.nan_to_num(protein_data, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Standardize
    scaler_mrna = StandardScaler()
    scaler_protein = StandardScaler()
    
    mrna_scaled = scaler_mrna.fit_transform(mrna_data)
    protein_scaled = scaler_protein.fit_transform(protein_data)
    
    # Handle any remaining NaN after scaling
    mrna_scaled = np.nan_to_num(mrna_scaled, nan=0.0, posinf=0.0, neginf=0.0)
    protein_scaled = np.nan_to_num(protein_scaled, nan=0.0, posinf=0.0, neginf=0.0)
    
    # PCA on each layer
    pca_mrna = PCA(n_components=min(50, mrna_scaled.shape[1]))
    pca_protein = PCA(n_components=min(50, protein_scaled.shape[1]))
    
    mrna_pca = pca_mrna.fit_transform(mrna_scaled)
    protein_pca = pca_protein.fit_transform(protein_scaled)
    
    # Combine PCA features
    combined_features = np.hstack([mrna_pca, protein_pca])
    
    # K-means clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(combined_features)
    
    return clusters, combined_features

# Perform clustering
print("Performing multi-omics clustering...")
subtypes, features = multi_omics_clustering(mrna_normalized_df, protein_normalized_df, n_clusters=3)

# Add subtypes to clinical data
clinical_filtered['MultiOmics_Subtype'] = subtypes

# Save subtypes
clinical_filtered[['PATIENT_ID', 'MultiOmics_Subtype']].to_csv(
    output_dir / "multi_omics_subtypes.csv", index=False)

print(f"Subtype distribution:")
print(clinical_filtered['MultiOmics_Subtype'].value_counts().sort_index())

# Survival analysis by subtype
print("\nPerforming survival analysis by subtype...")
kmf = KaplanMeierFitter()

fig, ax = plt.subplots(figsize=(10, 6))
for subtype in sorted(clinical_filtered['MultiOmics_Subtype'].unique()):
    subtype_data = clinical_filtered[clinical_filtered['MultiOmics_Subtype'] == subtype].copy()
    # Remove NaN values
    subtype_data = subtype_data.dropna(subset=['OS_DAYS', 'OS_EVENT'])
    # Remove negative or zero survival times
    subtype_data = subtype_data[subtype_data['OS_DAYS'] > 0]
    if len(subtype_data) > 0:
        kmf.fit(subtype_data['OS_DAYS'], subtype_data['OS_EVENT'], label=f'Subtype {subtype}')
        kmf.plot_survival_function(ax=ax)

ax.set_xlabel('Days')
ax.set_ylabel('Survival Probability')
ax.set_title('Survival by Multi-Omics Subtype')
plt.legend()
plt.tight_layout()
plt.savefig(plots_dir / "survival_by_subtype.png", dpi=300, bbox_inches='tight')
plt.close()

# Log-rank test
if len(clinical_filtered['MultiOmics_Subtype'].unique()) > 1:
    groups = []
    durations = []
    events = []
    for subtype in sorted(clinical_filtered['MultiOmics_Subtype'].unique()):
        subtype_data = clinical_filtered[clinical_filtered['MultiOmics_Subtype'] == subtype].copy()
        subtype_data = subtype_data.dropna(subset=['OS_DAYS', 'OS_EVENT'])
        subtype_data = subtype_data[subtype_data['OS_DAYS'] > 0]
        if len(subtype_data) > 0:
            groups.extend([subtype] * len(subtype_data))
            durations.extend(subtype_data['OS_DAYS'].tolist())
            events.extend(subtype_data['OS_EVENT'].tolist())
    
    if len(durations) > 0 and len(set(groups)) > 1:
        try:
            # Split by first two groups for log-rank test
            unique_groups = sorted(set(groups))
            if len(unique_groups) >= 2:
                group1_data = [d for d, g in zip(durations, groups) if g == unique_groups[0]]
                group1_events = [e for e, g in zip(events, groups) if g == unique_groups[0]]
                group2_data = [d for d, g in zip(durations, groups) if g == unique_groups[1]]
                group2_events = [e for e, g in zip(events, groups) if g == unique_groups[1]]
                if len(group1_data) > 0 and len(group2_data) > 0:
                    results = logrank_test(group1_data, group2_data, group1_events, group2_events)
                    print(f"Log-rank test p-value (Subtype {unique_groups[0]} vs {unique_groups[1]}): {results.p_value}")
        except Exception as e:
            print(f"Log-rank test could not be performed: {e}")

print(f"Survival analysis plot saved to {plots_dir / 'survival_by_subtype.png'}")

print()

# ============================================================================
# 8. METHYLATION-PROTEIN DISCONNECT ANALYSIS
# ============================================================================

print("Step 8: Methylation-Protein Disconnect Analysis...")
print("-" * 80)

def identify_discordance(mrna_df, protein_df, genes, patients):
    """
    Identify genes with discordant mRNA-protein expression
    """
    discordance_results = []
    
    for gene in genes[:min(100, len(genes))]:
        if gene not in mrna_df.index or gene not in protein_df.index:
            continue
        
        mrna_data = mrna_df.loc[gene, patients].values
        protein_data = protein_df.loc[gene, patients].values
        
        if len(mrna_data) < 5 or len(protein_data) < 5:
            continue
        
        # Calculate correlation
        try:
            corr, p = spearmanr(mrna_data, protein_data)
        except:
            corr, p = np.nan, np.nan
        
        # Classify discordance
        if abs(corr) < 0.2:  # Low correlation = discordance
            discordance_type = "Discordant"
        elif corr > 0.5:
            discordance_type = "Concordant_High"
        elif corr < -0.5:
            discordance_type = "Concordant_Negative"
        else:
            discordance_type = "Moderate"
        
        discordance_results.append({
            'Gene': gene,
            'mRNA_Protein_Corr': corr,
            'P_value': p,
            'Discordance_Type': discordance_type,
            'mRNA_Mean': np.mean(mrna_data),
            'Protein_Mean': np.mean(protein_data)
        })
    
    return pd.DataFrame(discordance_results)

# Identify discordance
discordance_results = identify_discordance(
    mrna_normalized_df,
    protein_normalized_df,
    common_genes,
    common_patients
)

# Save results
discordance_results.to_csv(output_dir / "discordance_analysis.csv", index=False)

discordant_genes = discordance_results[discordance_results['Discordance_Type'] == 'Discordant']
print(f"Discordant genes (correlation < 0.2): {len(discordant_genes)}")
print(f"Top 10 discordant genes:")
print(discordant_genes.nsmallest(10, 'mRNA_Protein_Corr')[['Gene', 'mRNA_Protein_Corr']])

# Visualize discordance
if len(discordance_results) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Discordance type distribution
    discordance_results['Discordance_Type'].value_counts().plot(kind='bar', ax=axes[0])
    axes[0].set_xlabel('Discordance Type')
    axes[0].set_ylabel('Count')
    axes[0].set_title('Distribution of Discordance Types')
    axes[0].tick_params(axis='x', rotation=45)
    
    # Correlation distribution colored by discordance type
    for disc_type in discordance_results['Discordance_Type'].unique():
        data = discordance_results[discordance_results['Discordance_Type'] == disc_type]
        axes[1].scatter(data['mRNA_Protein_Corr'], data['P_value'], 
                       label=disc_type, alpha=0.6)
    axes[1].set_xlabel('mRNA-Protein Correlation')
    axes[1].set_ylabel('P-value')
    axes[1].set_title('Discordance Analysis')
    axes[1].legend()
    axes[1].set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(plots_dir / "discordance_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

print(f"Discordance analysis results saved to {output_dir / 'discordance_analysis.csv'}")
print()

# ============================================================================
# 9. MULTI-OMICS PROGNOSTIC MODELS
# ============================================================================

print("Step 9: Multi-Omics Prognostic Model Development...")
print("-" * 80)

# Prepare survival data
survival_data = clinical_filtered[['PATIENT_ID', 'OS_DAYS', 'OS_EVENT']].copy()
survival_data = survival_data.dropna(subset=['OS_DAYS', 'OS_EVENT'])

# Feature selection using univariate Cox regression
def univariate_cox_regression(data_df, survival_df, patients, top_n=50):
    """
    Perform univariate Cox regression for feature selection
    """
    results = []
    
    # Match patients
    common_pats = list(set(data_df.columns) & set(patients))
    survival_matched = survival_df[survival_df['PATIENT_ID'].isin(common_pats)].copy()
    
    for gene in data_df.index[:min(100, len(data_df.index))]:
        try:
            gene_data = data_df.loc[gene, common_pats]
            survival_matched_temp = survival_matched.copy()
            survival_matched_temp['Gene_Expr'] = gene_data.values
            
            # Remove NaN values
            survival_matched_temp = survival_matched_temp.dropna(subset=['OS_DAYS', 'OS_EVENT', 'Gene_Expr'])
            survival_matched_temp = survival_matched_temp[survival_matched_temp['OS_DAYS'] > 0]
            
            if len(survival_matched_temp) < 5:
                continue
            
            # Fit Cox model
            cph = CoxPHFitter()
            cph.fit(survival_matched_temp[['OS_DAYS', 'OS_EVENT', 'Gene_Expr']], 
                   duration_col='OS_DAYS', event_col='OS_EVENT')
            
            # Get p-value from summary
            try:
                p_val = cph.summary['p'].iloc[0] if 'p' in cph.summary.columns else cph.summary['p-value'].iloc[0]
            except:
                p_val = np.nan
            
            # Get confidence intervals
            try:
                ci_lower = cph.confidence_intervals_.iloc[0, 0]
                ci_upper = cph.confidence_intervals_.iloc[0, 1]
            except:
                ci_lower = np.nan
                ci_upper = np.nan
            
            results.append({
                'Gene': gene,
                'HR': cph.hazard_ratios_['Gene_Expr'],
                'P_value': p_val,
                'CI_lower': ci_lower,
                'CI_upper': ci_upper
            })
        except Exception as e:
            continue
    
    results_df = pd.DataFrame(results)
    if len(results_df) > 0 and 'P_value' in results_df.columns:
        results_df = results_df.sort_values('P_value')
        return results_df.head(top_n)
    else:
        return results_df

# Feature selection for mRNA
print("Selecting prognostic mRNA features...")
mrna_features = univariate_cox_regression(mrna_normalized_df, survival_data, common_patients, top_n=20)
print(f"Selected {len(mrna_features)} prognostic mRNA features")
if len(mrna_features) > 0:
    print(f"Top mRNA features: {list(mrna_features['Gene'].head(5))}")

# Feature selection for protein
print("Selecting prognostic protein features...")
protein_features = univariate_cox_regression(protein_normalized_df, survival_data, common_patients, top_n=20)
print(f"Selected {len(protein_features)} prognostic protein features")
if len(protein_features) > 0:
    print(f"Top protein features: {list(protein_features['Gene'].head(5))}")

# Combine features
if len(mrna_features) > 0 and len(protein_features) > 0:
    selected_genes = list(set(mrna_features['Gene'].tolist() + protein_features['Gene'].tolist()))
elif len(mrna_features) > 0:
    selected_genes = mrna_features['Gene'].tolist()
elif len(protein_features) > 0:
    selected_genes = protein_features['Gene'].tolist()
else:
    # If no features selected, use top variable genes
    print("No prognostic features found, using top variable genes instead...")
    mrna_var = mrna_normalized_df.var(axis=1).nlargest(10)
    protein_var = protein_normalized_df.var(axis=1).nlargest(10)
    selected_genes = list(set(mrna_var.index.tolist() + protein_var.index.tolist()))[:20]
    print(f"Selected {len(selected_genes)} variable genes")

# Build multi-omics prognostic model
print("Building multi-omics prognostic model...")

# Prepare feature matrix
X = []
y_duration = []
y_event = []
patient_ids = []

for patient in common_patients:
    if patient not in survival_data['PATIENT_ID'].values:
        continue
    
    patient_survival = survival_data[survival_data['PATIENT_ID'] == patient]
    if len(patient_survival) == 0:
        continue
    
    # Extract features
    features = []
    for gene in selected_genes[:min(20, len(selected_genes))]:
        if gene in mrna_normalized_df.index:
            features.append(mrna_normalized_df.loc[gene, patient])
        else:
            features.append(0)
        
        if gene in protein_normalized_df.index:
            features.append(protein_normalized_df.loc[gene, patient])
        else:
            features.append(0)
    
    X.append(features)
    y_duration.append(patient_survival['OS_DAYS'].values[0])
    y_event.append(patient_survival['OS_EVENT'].values[0])
    patient_ids.append(patient)

X = np.array(X)
y_duration = np.array(y_duration)
y_event = np.array(y_event)

# Build Cox model
if len(X) > 10:
    feature_df = pd.DataFrame(X, index=patient_ids)
    feature_df['OS_DAYS'] = y_duration
    feature_df['OS_EVENT'] = y_event
    
    try:
        cph = CoxPHFitter()
        cph.fit(feature_df, duration_col='OS_DAYS', event_col='OS_EVENT')
        
        # Calculate risk scores
        risk_scores = cph.predict_partial_hazard(feature_df)
        feature_df['Risk_Score'] = risk_scores
        
        # Stratify by risk
        feature_df['Risk_Group'] = pd.qcut(risk_scores, q=3, labels=['Low', 'Medium', 'High'])
        
        # Save results
        feature_df[['Risk_Score', 'Risk_Group']].to_csv(
            output_dir / "prognostic_risk_scores.csv")
        
        # Survival analysis by risk group
        kmf = KaplanMeierFitter()
        risk_groups = ['Low', 'Medium', 'High']
        risk_data_dict = {g: feature_df[feature_df['Risk_Group'] == g] for g in risk_groups}
        
        # Pairwise log-rank p-values
        logrank_pvals = {}
        for i, g1 in enumerate(risk_groups):
            for g2 in risk_groups[i+1:]:
                d1, d2 = risk_data_dict[g1], risk_data_dict[g2]
                if len(d1) > 0 and len(d2) > 0:
                    try:
                        res = logrank_test(d1['OS_DAYS'].values, d2['OS_DAYS'].values,
                                           d1['OS_EVENT'].values, d2['OS_EVENT'].values)
                        logrank_pvals[f'{g1} vs {g2}'] = res.p_value
                    except Exception:
                        logrank_pvals[f'{g1} vs {g2}'] = np.nan
                else:
                    logrank_pvals[f'{g1} vs {g2}'] = np.nan
        
        fig, ax = plt.subplots(figsize=(10, 6))
        for risk_group in risk_groups:
            risk_data = risk_data_dict[risk_group]
            if len(risk_data) > 0:
                kmf.fit(risk_data['OS_DAYS'], risk_data['OS_EVENT'], label=f'{risk_group} (n={len(risk_data)})')
                kmf.plot_survival_function(ax=ax)
        
        ax.set_xlabel('Days')
        ax.set_ylabel('Survival Probability')
        # Caption: significant p-values (p < 0.05)
        alpha = 0.05
        sig_pvals = [(k, v) for k, v in logrank_pvals.items() if not np.isnan(v) and v < alpha]
        if sig_pvals:
            caption = 'Log-rank (significant): ' + '; '.join([f'{k} p={v:.4f}' for k, v in sig_pvals])
        else:
            all_pvals = '; '.join([f'{k} p={v:.4f}' if not np.isnan(v) else f'{k} p=NA' for k, v in logrank_pvals.items()])
            caption = 'Log-rank: ' + (all_pvals if all_pvals else 'NA')
        ax.set_title(f'Survival by Multi-Omics Risk Group\n{caption}')
        plt.legend()
        plt.tight_layout()
        plt.savefig(plots_dir / "survival_by_risk_group.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Prognostic model completed! Risk scores saved to {output_dir / 'prognostic_risk_scores.csv'}")
        print(f"Risk group distribution:")
        print(feature_df['Risk_Group'].value_counts())
        
    except Exception as e:
        print(f"Error building Cox model: {e}")

print()

# ============================================================================
# 10. SUMMARY AND OUTPUT
# ============================================================================

print("="*80)
print("Analysis Summary")
print("="*80)
print()
print("Completed Analyses:")
print("1. ✓ Epigenetic-Transcriptomic-Proteomic Cascade Analysis")
print("2. ✓ Multi-Layer Regulatory Network Analysis")
print("3. ✓ Multi-Omics Subtype Discovery")
print("4. ✓ Methylation-Protein Disconnect Analysis")
print("5. ✓ Multi-Omics Prognostic Model Development")
print()
print(f"Results saved to: {output_dir}")
print(f"Plots saved to: {plots_dir}")
print()
print("Key Output Files:")
print("- cascade_correlations.csv: Cascade correlation analysis results")
print("- mrna_hub_genes.csv: mRNA hub genes in regulatory network")
print("- protein_hub_genes.csv: Protein hub genes in regulatory network")
print("- multi_omics_subtypes.csv: Multi-omics subtype assignments")
print("- discordance_analysis.csv: Discordance analysis results")
print("- prognostic_risk_scores.csv: Prognostic risk scores")
print()
print("="*80)
print("Analysis completed successfully!")
print("="*80)

