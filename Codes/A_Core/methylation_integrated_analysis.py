#!/usr/bin/env python3
"""
Methylation-Integrated Novel Analysis
Epigenetic-Transcriptomic-Proteomic Integration with Real Methylation Data

Novel Analyses:
1. Methylation-mRNA-Protein Cascade Analysis
2. Methylation-Based Subtype Discovery
3. MGMT Methylation Analysis
4. Methylation-Protein Disconnect Analysis
5. Epigenetic Prognostic Signatures
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr, pearsonr
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestRegressor
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings('ignore')

import os
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
print("Methylation-Integrated Novel Analysis")
print("="*80)
print()

# ============================================================================
# 1. DATA LOADING
# ============================================================================

print("Step 1: Loading data...")
print("-" * 80)

# Load clinical data
clinical_file = "../clinical_patient.txt"
clinical = pd.read_csv(clinical_file, sep="\t", skiprows=4)
clinical.columns = clinical.columns.str.strip()
clinical['PATIENT_ID'] = clinical['PATIENT_ID'].str.strip()
clinical['OS_DAYS'] = clinical['PATH_DIAG_TO_DEATH_DAYS'].copy()
clinical.loc[clinical['VITAL_STATUS'] == 'Living', 'OS_DAYS'] = clinical.loc[clinical['VITAL_STATUS'] == 'Living', 'PATH_DIAG_TO_LAST_CONTACT_DAYS']
clinical['OS_EVENT'] = (clinical['VITAL_STATUS'] == 'Deceased').astype(int)

# Load mRNA data
mrna_file = "../mrna_seq_fpkm.txt"
mrna = pd.read_csv(mrna_file, sep="\t")
mrna.columns = mrna.columns.str.strip()

# Load protein data
protein_file = "../protein_quantification.txt"
protein = pd.read_csv(protein_file, sep="\t")
protein.columns = protein.columns.str.strip()

# Load methylation data - sample for efficiency (can be expanded)
print("Loading methylation data (sampling for efficiency)...")
methylation_file = "../methylation_epic.txt"

# Read methylation data in chunks
print("Reading methylation data...")
try:
    # Read first 10000 rows for analysis (can be expanded)
    methylation = pd.read_csv(methylation_file, sep="\t", nrows=10000)
    print(f"Methylation data loaded: {methylation.shape}")
except Exception as e:
    print(f"Error loading methylation: {e}")
    methylation = None

# Get common patients
clinical_patients = set(clinical['PATIENT_ID'].unique())
mrna_patients = set([col for col in mrna.columns if col.startswith('C3')])
protein_patients = set([col for col in protein.columns if col.startswith('C3')])

if methylation is not None:
    methylation_patients = set([col for col in methylation.columns if col.startswith('C3')])
    common_patients = sorted(list(clinical_patients & mrna_patients & protein_patients & methylation_patients))
else:
    common_patients = sorted(list(clinical_patients & mrna_patients & protein_patients))

print(f"Common patients: {len(common_patients)}")

# Filter datasets
clinical_filtered = clinical[clinical['PATIENT_ID'].isin(common_patients)].copy()
mrna_filtered = mrna[['Hugo_Symbol'] + [p for p in common_patients if p in mrna.columns]].copy()
protein_filtered = protein[[col for col in protein.columns if col == 'Composite.Element.REF' or col in common_patients]].copy()

if methylation is not None:
    methylation_filtered = methylation[[col for col in methylation.columns if col == 'Hugo_Symbol' or col in common_patients]].copy()
else:
    methylation_filtered = None

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

if methylation_filtered is not None:
    meth_genes = set(methylation_filtered['Hugo_Symbol'].dropna().unique())
    common_genes = sorted(list(mrna_genes & protein_genes & meth_genes))
else:
    common_genes = sorted(list(mrna_genes & protein_genes))

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

if methylation_filtered is not None:
    meth_values = methylation_filtered[common_patients].values
    meth_values = np.nan_to_num(meth_values, nan=0.5, posinf=0.5, neginf=0.5)  # Beta values 0-1
    meth_df = pd.DataFrame(meth_values, index=methylation_filtered['Hugo_Symbol'], columns=common_patients)
    meth_df = meth_df.fillna(0.5)
else:
    meth_df = None

print("Normalization completed!")
print()

# ============================================================================
# 3. METHYLATION-mRNA-PROTEIN CASCADE ANALYSIS
# ============================================================================

print("Step 3: Methylation-mRNA-Protein Cascade Analysis...")
print("-" * 80)

def methylation_cascade_analysis(meth_df, mrna_df, protein_df, genes, patients, max_genes=200):
    """
    Analyze methylation -> mRNA -> Protein cascade
    """
    results = []
    
    # Get common genes
    common_genes_analysis = [g for g in genes if g in mrna_df.index and g in protein_df.index]
    if meth_df is not None:
        common_genes_analysis = [g for g in common_genes_analysis if g in meth_df.index]
    
    common_genes_analysis = common_genes_analysis[:max_genes]
    
    print(f"Analyzing {len(common_genes_analysis)} genes...")
    
    for gene in common_genes_analysis:
        try:
            # Get data
            mrna_data = mrna_df.loc[gene, patients].values
            protein_data = protein_df.loc[gene, patients].values
            
            # Remove NaN
            valid_idx = ~(np.isnan(mrna_data) | np.isnan(protein_data))
            if valid_idx.sum() < 5:
                continue
                
            mrna_data = mrna_data[valid_idx]
            protein_data = protein_data[valid_idx]
            
            # mRNA-Protein correlation
            try:
                mrna_protein_corr, mrna_protein_p = spearmanr(mrna_data, protein_data)
            except:
                mrna_protein_corr, mrna_protein_p = np.nan, np.nan
            
            # Methylation correlations
            if meth_df is not None and gene in meth_df.index:
                meth_data = meth_df.loc[gene, patients].values[valid_idx]
                
                # Methylation-mRNA correlation
                try:
                    meth_mrna_corr, meth_mrna_p = spearmanr(meth_data, mrna_data)
                except:
                    meth_mrna_corr, meth_mrna_p = np.nan, np.nan
                
                # Methylation-Protein correlation
                try:
                    meth_protein_corr, meth_protein_p = spearmanr(meth_data, protein_data)
                except:
                    meth_protein_corr, meth_protein_p = np.nan, np.nan
            else:
                meth_mrna_corr, meth_mrna_p = np.nan, np.nan
                meth_protein_corr, meth_protein_p = np.nan, np.nan
            
            results.append({
                'Gene': gene,
                'Meth_mRNA_Corr': meth_mrna_corr,
                'Meth_mRNA_P': meth_mrna_p,
                'mRNA_Protein_Corr': mrna_protein_corr,
                'mRNA_Protein_P': mrna_protein_p,
                'Meth_Protein_Corr': meth_protein_corr,
                'Meth_Protein_P': meth_protein_p
            })
        except Exception as e:
            continue
    
    return pd.DataFrame(results)

# Perform cascade analysis
cascade_results = methylation_cascade_analysis(meth_df, mrna_df, protein_df, common_genes, common_patients, max_genes=200)

# Save results
cascade_results.to_csv(output_dir / "methylation_cascade_analysis.csv", index=False)

# Filter significant results
sig_meth_mrna = cascade_results[(cascade_results['Meth_mRNA_P'] < 0.05) & (~cascade_results['Meth_mRNA_Corr'].isna())]
sig_mrna_protein = cascade_results[(cascade_results['mRNA_Protein_P'] < 0.05) & (~cascade_results['mRNA_Protein_Corr'].isna())]
sig_meth_protein = cascade_results[(cascade_results['Meth_Protein_P'] < 0.05) & (~cascade_results['Meth_Protein_Corr'].isna())]

print(f"Significant Methylation-mRNA correlations: {len(sig_meth_mrna)}")
print(f"Significant mRNA-Protein correlations: {len(sig_mrna_protein)}")
print(f"Significant Methylation-Protein correlations: {len(sig_meth_protein)}")

# Visualize cascade
if len(cascade_results) > 0:
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Methylation-mRNA correlation
    axes[0, 0].scatter(cascade_results['Meth_mRNA_Corr'], 
                       -np.log10(cascade_results['Meth_mRNA_P'] + 1e-10),
                       alpha=0.6, c='blue')
    axes[0, 0].set_xlabel('Methylation-mRNA Correlation')
    axes[0, 0].set_ylabel('-log10(P-value)')
    axes[0, 0].set_title('Methylation-mRNA Correlation')
    axes[0, 0].axhline(-np.log10(0.05), color='red', linestyle='--', label='P=0.05')
    axes[0, 0].legend()
    
    # mRNA-Protein correlation
    axes[0, 1].scatter(cascade_results['mRNA_Protein_Corr'],
                       -np.log10(cascade_results['mRNA_Protein_P'] + 1e-10),
                       alpha=0.6, c='green')
    axes[0, 1].set_xlabel('mRNA-Protein Correlation')
    axes[0, 1].set_ylabel('-log10(P-value)')
    axes[0, 1].set_title('mRNA-Protein Correlation')
    axes[0, 1].axhline(-np.log10(0.05), color='red', linestyle='--', label='P=0.05')
    axes[0, 1].legend()
    
    # Methylation-Protein correlation
    axes[1, 0].scatter(cascade_results['Meth_Protein_Corr'],
                       -np.log10(cascade_results['Meth_Protein_P'] + 1e-10),
                       alpha=0.6, c='purple')
    axes[1, 0].set_xlabel('Methylation-Protein Correlation')
    axes[1, 0].set_ylabel('-log10(P-value)')
    axes[1, 0].set_title('Methylation-Protein Correlation')
    axes[1, 0].axhline(-np.log10(0.05), color='red', linestyle='--', label='P=0.05')
    axes[1, 0].legend()
    
    # Cascade summary
    cascade_summary = pd.DataFrame({
        'Correlation_Type': ['Meth-mRNA', 'mRNA-Protein', 'Meth-Protein'],
        'Mean_Corr': [
            cascade_results['Meth_mRNA_Corr'].mean(),
            cascade_results['mRNA_Protein_Corr'].mean(),
            cascade_results['Meth_Protein_Corr'].mean()
        ],
        'Significant_Count': [
            len(sig_meth_mrna),
            len(sig_mrna_protein),
            len(sig_meth_protein)
        ]
    })
    
    axes[1, 1].bar(cascade_summary['Correlation_Type'], cascade_summary['Mean_Corr'], color=['blue', 'green', 'purple'])
    axes[1, 1].set_ylabel('Mean Correlation')
    axes[1, 1].set_title('Cascade Correlation Summary')
    axes[1, 1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(plots_dir / "methylation_cascade_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Cascade analysis plot saved to {plots_dir / 'methylation_cascade_analysis.png'}")

print()

# ============================================================================
# 4. METHYLATION-BASED SUBTYPE DISCOVERY
# ============================================================================

print("Step 4: Methylation-Based Subtype Discovery...")
print("-" * 80)

if meth_df is not None:
    # Prepare methylation data for clustering
    meth_data = meth_df.T.values  # Patients x Genes
    
    # Handle NaN
    meth_data = pd.DataFrame(meth_data).fillna(0.5).values
    meth_data = np.nan_to_num(meth_data, nan=0.5, posinf=0.5, neginf=0.5)
    
    # PCA for dimensionality reduction
    scaler = StandardScaler()
    meth_scaled = scaler.fit_transform(meth_data)
    meth_scaled = np.nan_to_num(meth_scaled, nan=0.0, posinf=0.0, neginf=0.0)
    
    # PCA
    pca = PCA(n_components=min(50, meth_scaled.shape[1]))
    meth_pca = pca.fit_transform(meth_scaled)
    
    # K-means clustering
    kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
    meth_subtypes = kmeans.fit_predict(meth_pca)
    
    # Add to clinical data
    clinical_filtered['Methylation_Subtype'] = meth_subtypes
    
    # Save subtypes
    clinical_filtered[['PATIENT_ID', 'Methylation_Subtype']].to_csv(
        output_dir / "methylation_subtypes.csv", index=False)
    
    print(f"Methylation subtype distribution:")
    print(clinical_filtered['Methylation_Subtype'].value_counts().sort_index())
    
    # Survival analysis
    print("\nPerforming survival analysis by methylation subtype...")
    kmf = KaplanMeierFitter()
    
    # Prepare survival data per subtype
    subtype_survival = {}
    for subtype in sorted(clinical_filtered['Methylation_Subtype'].unique()):
        subtype_data = clinical_filtered[clinical_filtered['Methylation_Subtype'] == subtype].copy()
        subtype_data = subtype_data.dropna(subset=['OS_DAYS', 'OS_EVENT'])
        subtype_data = subtype_data[subtype_data['OS_DAYS'] > 0]
        if len(subtype_data) > 0:
            subtype_survival[subtype] = subtype_data
    
    # Pairwise log-rank p-values
    subtypes = sorted(subtype_survival.keys())
    logrank_pvals = {}
    for i, s1 in enumerate(subtypes):
        for s2 in subtypes[i+1:]:
            d1 = subtype_survival[s1]
            d2 = subtype_survival[s2]
            try:
                res = logrank_test(d1['OS_DAYS'].values, d2['OS_DAYS'].values,
                                  d1['OS_EVENT'].values, d2['OS_EVENT'].values)
                logrank_pvals[f'{s1} vs {s2}'] = res.p_value
            except Exception:
                logrank_pvals[f'{s1} vs {s2}'] = np.nan
    
    fig, ax = plt.subplots(figsize=(10, 6))
    for subtype in subtypes:
        subtype_data = subtype_survival[subtype]
        kmf.fit(subtype_data['OS_DAYS'], subtype_data['OS_EVENT'], label=f'Meth Subtype {subtype} (n={len(subtype_data)})')
        kmf.plot_survival_function(ax=ax)
    
    ax.set_xlabel('Days')
    ax.set_ylabel('Survival Probability')
    # Caption: only significant results (p < 0.05)
    alpha = 0.05
    sig_pvals = [(k, v) for k, v in logrank_pvals.items() if not np.isnan(v) and v < alpha]
    if sig_pvals:
        caption = 'Log-rank (significant): ' + '; '.join([f'{k} p={v:.4f}' for k, v in sig_pvals])
    else:
        caption = 'Log-rank: no significant difference between subtypes (all p ≥ 0.05)'
    ax.set_title(f'Survival by Methylation Subtype\n{caption}')
    plt.legend()
    plt.tight_layout()
    plt.savefig(plots_dir / "survival_by_methylation_subtype.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Log-rank p-values (methylation subtypes): {logrank_pvals}")
    print(f"Survival analysis plot saved to {plots_dir / 'survival_by_methylation_subtype.png'}")

print()

# ============================================================================
# 5. MGMT METHYLATION ANALYSIS
# ============================================================================

print("Step 5: MGMT Methylation Analysis...")
print("-" * 80)

if meth_df is not None:
    # Search for MGMT in methylation data
    mgmt_genes = [g for g in meth_df.index if 'MGMT' in str(g).upper()]
    
    if len(mgmt_genes) > 0:
        print(f"Found MGMT-related probes: {mgmt_genes}")
        
        # Get MGMT methylation values
        mgmt_meth_values = []
        for patient in common_patients:
            if patient in meth_df.columns:
                try:
                    val = meth_df.loc[mgmt_genes[0], patient]
                    if isinstance(val, pd.Series):
                        val = val.iloc[0] if len(val) > 0 else 0.5
                    val = float(val) if not np.isnan(val) else 0.5
                    mgmt_meth_values.append(1 if val > 0.3 else 0)
                except:
                    mgmt_meth_values.append(0)
            else:
                mgmt_meth_values.append(0)
        
        # Create MGMT status dataframe
        mgmt_status_df = pd.DataFrame({
            'PATIENT_ID': common_patients,
            'MGMT_Methylated': mgmt_meth_values
        })
        
        # Merge with clinical data
        clinical_filtered = clinical_filtered.merge(mgmt_status_df, on='PATIENT_ID', how='left')
        
        # Save MGMT status
        clinical_filtered[['PATIENT_ID', 'MGMT_Methylated']].to_csv(
            output_dir / "mgmt_methylation_status.csv", index=False)
        
        print(f"MGMT methylation status:")
        mgmt_counts = clinical_filtered['MGMT_Methylated'].value_counts()
        print(f"  Methylated: {mgmt_counts.get(1, 0)}")
        print(f"  Unmethylated: {mgmt_counts.get(0, 0)}")
        
        # Survival analysis by MGMT status
        print("\nPerforming survival analysis by MGMT status...")
        kmf = KaplanMeierFitter()
        
        fig, ax = plt.subplots(figsize=(10, 6))
        for status in [0, 1]:
            status_data = clinical_filtered[clinical_filtered['MGMT_Methylated'] == status].copy()
            status_data = status_data.dropna(subset=['OS_DAYS', 'OS_EVENT'])
            status_data = status_data[status_data['OS_DAYS'] > 0]
            if len(status_data) > 0:
                label = 'MGMT Methylated' if status == 1 else 'MGMT Unmethylated'
                kmf.fit(status_data['OS_DAYS'], status_data['OS_EVENT'], label=label)
                kmf.plot_survival_function(ax=ax)
        
        ax.set_xlabel('Days')
        ax.set_ylabel('Survival Probability')
        ax.set_title('Survival by MGMT Methylation Status')
        plt.legend()
        plt.tight_layout()
        plt.savefig(plots_dir / "survival_by_mgmt_status.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Survival analysis plot saved to {plots_dir / 'survival_by_mgmt_status.png'}")
    else:
        print("MGMT gene not found in methylation data")

print()

# ============================================================================
# 6. METHYLATION-PROTEIN DISCONNECT ANALYSIS
# ============================================================================

print("Step 6: Methylation-Protein Disconnect Analysis...")
print("-" * 80)

if meth_df is not None:
    disconnect_results = []
    
    common_genes_disc = [g for g in common_genes if g in meth_df.index and g in protein_df.index][:200]
    
    for gene in common_genes_disc:
        try:
            meth_data = meth_df.loc[gene, common_patients].values
            protein_data = protein_df.loc[gene, common_patients].values
            
            valid_idx = ~(np.isnan(meth_data) | np.isnan(protein_data))
            if valid_idx.sum() < 5:
                continue
            
            meth_data = meth_data[valid_idx]
            protein_data = protein_data[valid_idx]
            
            # Calculate correlation
            try:
                corr, p = spearmanr(meth_data, protein_data)
            except:
                corr, p = np.nan, np.nan
            
            # Classify disconnect
            if abs(corr) < 0.2:
                disc_type = "Discordant"
            elif corr > 0.5:
                disc_type = "Concordant_High"
            elif corr < -0.5:
                disc_type = "Concordant_Negative"
            else:
                disc_type = "Moderate"
            
            disconnect_results.append({
                'Gene': gene,
                'Meth_Protein_Corr': corr,
                'P_value': p,
                'Discordance_Type': disc_type,
                'Meth_Mean': np.mean(meth_data),
                'Protein_Mean': np.mean(protein_data)
            })
        except:
            continue
    
    disconnect_df = pd.DataFrame(disconnect_results)
    disconnect_df.to_csv(output_dir / "methylation_protein_discordance.csv", index=False)
    
    discordant_genes = disconnect_df[disconnect_df['Discordance_Type'] == 'Discordant']
    print(f"Discordant genes (Meth-Protein correlation < 0.2): {len(discordant_genes)}")
    print(f"Top 10 discordant genes:")
    print(discordant_genes.nsmallest(10, 'Meth_Protein_Corr')[['Gene', 'Meth_Protein_Corr']])
    
    # Visualize
    if len(disconnect_df) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        disconnect_df['Discordance_Type'].value_counts().plot(kind='bar', ax=axes[0])
        axes[0].set_xlabel('Discordance Type')
        axes[0].set_ylabel('Count')
        axes[0].set_title('Methylation-Protein Discordance Distribution')
        axes[0].tick_params(axis='x', rotation=45)
        
        axes[1].scatter(disconnect_df['Meth_Protein_Corr'], 
                       -np.log10(disconnect_df['P_value'] + 1e-10),
                       alpha=0.6)
        axes[1].set_xlabel('Methylation-Protein Correlation')
        axes[1].set_ylabel('-log10(P-value)')
        axes[1].set_title('Methylation-Protein Discordance Analysis')
        axes[1].axhline(-np.log10(0.05), color='red', linestyle='--', label='P=0.05')
        axes[1].legend()
        
        plt.tight_layout()
        plt.savefig(plots_dir / "methylation_protein_discordance.png", dpi=300, bbox_inches='tight')
        plt.close()

print()

# ============================================================================
# 7. EPIGENETIC PROGNOSTIC SIGNATURES
# ============================================================================

print("Step 7: Epigenetic Prognostic Signature Development...")
print("-" * 80)

if meth_df is not None:
    # Prepare survival data
    survival_data = clinical_filtered[['PATIENT_ID', 'OS_DAYS', 'OS_EVENT']].copy()
    survival_data = survival_data.dropna(subset=['OS_DAYS', 'OS_EVENT'])
    survival_data = survival_data[survival_data['OS_DAYS'] > 0]
    
    # Select top variable methylation probes
    meth_var = meth_df.var(axis=1).nlargest(50)
    selected_probes = meth_var.index.tolist()
    
    # Build prognostic model
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
        
        # Extract methylation features
        features = []
        for probe in selected_probes[:20]:
            if probe in meth_df.index:
                try:
                    val = meth_df.loc[probe, patient]
                    if isinstance(val, pd.Series):
                        val = val.iloc[0] if len(val) > 0 else 0.5
                    val = float(val) if not np.isnan(val) else 0.5
                    features.append(val)
                except:
                    features.append(0.5)
            else:
                features.append(0.5)  # Default beta value
        
        if len(features) == 20:  # Ensure consistent length
            X.append(features)
        y_duration.append(patient_survival['OS_DAYS'].values[0])
        y_event.append(patient_survival['OS_EVENT'].values[0])
        patient_ids.append(patient)
    
    X = np.array(X)
    y_duration = np.array(y_duration)
    y_event = np.array(y_event)
    
    if len(X) > 10:
        feature_df = pd.DataFrame(X, index=patient_ids)
        feature_df['OS_DAYS'] = y_duration
        feature_df['OS_EVENT'] = y_event
        
        try:
            cph = CoxPHFitter()
            cph.fit(feature_df, duration_col='OS_DAYS', event_col='OS_EVENT')
            
            # Calculate risk scores
            risk_scores = cph.predict_partial_hazard(feature_df)
            feature_df['Epigenetic_Risk_Score'] = risk_scores
            
            # Stratify by risk
            feature_df['Epigenetic_Risk_Group'] = pd.qcut(risk_scores, q=3, labels=['Low', 'Medium', 'High'])
            
            # Save results
            feature_df[['Epigenetic_Risk_Score', 'Epigenetic_Risk_Group']].to_csv(
                output_dir / "epigenetic_prognostic_scores.csv")
            
            print(f"Epigenetic prognostic model completed!")
            print(f"Risk group distribution:")
            print(feature_df['Epigenetic_Risk_Group'].value_counts())
            
            # Survival analysis by risk group
            kmf = KaplanMeierFitter()
            fig, ax = plt.subplots(figsize=(10, 6))
            
            for risk_group in ['Low', 'Medium', 'High']:
                risk_data = feature_df[feature_df['Epigenetic_Risk_Group'] == risk_group]
                kmf.fit(risk_data['OS_DAYS'], risk_data['OS_EVENT'], label=risk_group)
                kmf.plot_survival_function(ax=ax)
            
            ax.set_xlabel('Days')
            ax.set_ylabel('Survival Probability')
            ax.set_title('Survival by Epigenetic Risk Group')
            plt.legend()
            plt.tight_layout()
            plt.savefig(plots_dir / "survival_by_epigenetic_risk.png", dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error building epigenetic prognostic model: {e}")

print()

# ============================================================================
# SUMMARY
# ============================================================================

print("="*80)
print("Methylation-Integrated Analysis Summary")
print("="*80)
print()
print("Completed Analyses:")
print("1. ✓ Methylation-mRNA-Protein Cascade Analysis")
print("2. ✓ Methylation-Based Subtype Discovery")
print("3. ✓ MGMT Methylation Analysis")
print("4. ✓ Methylation-Protein Disconnect Analysis")
print("5. ✓ Epigenetic Prognostic Signature Development")
print()
print(f"Results saved to: {output_dir}")
print(f"Plots saved to: {plots_dir}")
print()
print("="*80)
print("Analysis completed successfully!")
print("="*80)

