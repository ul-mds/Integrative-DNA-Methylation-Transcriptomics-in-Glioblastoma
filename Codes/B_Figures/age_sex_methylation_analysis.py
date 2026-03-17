#!/usr/bin/env python3
"""
Age and Sex-Based Methylation Analysis
Analyze methylation changes by age groups (Young: 0-40, Mid: 41-55, Old: 56+) 
and sex (Male/Female)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, kruskal, mannwhitneyu
from scipy.stats import linregress
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# Create output directories
output_dir = Path("results")
plots_dir = Path("plots")
output_dir.mkdir(exist_ok=True)
plots_dir.mkdir(exist_ok=True)

print("="*80)
print("Age and Sex-Based Methylation Analysis")
print("="*80)
print()

# ============================================================================
# 1. LOAD CLINICAL DATA
# ============================================================================

print("Step 1: Loading clinical data...")
print("-" * 80)

clinical_file = "../clinical_patient.txt"
clinical = pd.read_csv(clinical_file, sep="\t", skiprows=4)
clinical.columns = clinical.columns.str.strip()
clinical['PATIENT_ID'] = clinical['PATIENT_ID'].str.strip()

print(f"Clinical data loaded: {clinical.shape}")
print(f"Columns: {list(clinical.columns)}")

# Check age and sex columns
print("\nChecking available columns...")
if 'AGE' in clinical.columns:
    age_col = 'AGE'
elif 'AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS' in clinical.columns:
    age_col = 'AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS'
else:
    # Try to find age column
    age_cols = [col for col in clinical.columns if 'AGE' in col.upper()]
    if age_cols:
        age_col = age_cols[0]
        print(f"Using age column: {age_col}")
    else:
        print("ERROR: Age column not found!")
        print("Available columns:", list(clinical.columns))
        exit(1)

if 'SEX' in clinical.columns:
    sex_col = 'SEX'
elif 'GENDER' in clinical.columns:
    sex_col = 'GENDER'
else:
    sex_cols = [col for col in clinical.columns if 'SEX' in col.upper() or 'GENDER' in col.upper()]
    if sex_cols:
        sex_col = sex_cols[0]
        print(f"Using sex column: {sex_col}")
    else:
        print("ERROR: Sex column not found!")
        exit(1)

print(f"Age column: {age_col}")
print(f"Sex column: {sex_col}")

# Create age groups
clinical['Age_Group'] = pd.cut(clinical[age_col], 
                               bins=[0, 40, 55, 200], 
                               labels=['Young', 'Mid', 'Old'],
                               include_lowest=True)

# Clean sex column
clinical[sex_col] = clinical[sex_col].str.strip().str.upper()
clinical = clinical[clinical[sex_col].isin(['MALE', 'FEMALE', 'M', 'F'])]

# Standardize sex values
clinical['Sex'] = clinical[sex_col].map({'MALE': 'Male', 'M': 'Male', 
                                         'FEMALE': 'Female', 'F': 'Female'})

# Remove rows with missing age or sex
clinical_clean = clinical.dropna(subset=[age_col, 'Sex', 'Age_Group'])
clinical_clean = clinical_clean[clinical_clean['PATIENT_ID'].notna()]

print(f"\nClinical data after cleaning: {len(clinical_clean)} patients")
print(f"\nAge group distribution:")
print(clinical_clean['Age_Group'].value_counts().sort_index())
print(f"\nSex distribution:")
print(clinical_clean['Sex'].value_counts())
print(f"\nAge group by Sex:")
print(pd.crosstab(clinical_clean['Age_Group'], clinical_clean['Sex']))
print()

# ============================================================================
# 2. LOAD METHYLATION DATA
# ============================================================================

print("Step 2: Loading methylation data...")
print("-" * 80)

methylation_file = "../methylation_epic.txt"
print("Loading methylation data (this may take a while)...")

# Read methylation data
methylation = pd.read_csv(methylation_file, sep="\t", nrows=50000)  # Adjust as needed
print(f"Methylation data loaded: {methylation.shape}")

# Get patient IDs from methylation data
meth_patients = [col for col in methylation.columns if col.startswith('C3')]
print(f"Patients in methylation data: {len(meth_patients)}")

# Get common patients
common_patients = sorted(list(set(clinical_clean['PATIENT_ID'].unique()) & 
                              set(meth_patients)))
print(f"Common patients: {len(common_patients)}")

# Filter clinical data for common patients
clinical_filtered = clinical_clean[clinical_clean['PATIENT_ID'].isin(common_patients)].copy()
print(f"Clinical data for common patients: {len(clinical_filtered)}")
print()

# ============================================================================
# 3. ANALYZE METHYLATION BY AGE AND SEX
# ============================================================================

print("Step 3: Analyzing methylation by age and sex...")
print("-" * 80)

results = []

# Process genes in batches
batch_size = 1000
total_genes = len(methylation)

print(f"Processing {total_genes} genes in batches of {batch_size}...")

for batch_start in range(0, total_genes, batch_size):
    batch_end = min(batch_start + batch_size, total_genes)
    batch_genes = methylation.iloc[batch_start:batch_end]
    
    for idx, row in batch_genes.iterrows():
        gene = row['Hugo_Symbol']
        
        if pd.isna(gene):
            continue
        
        try:
            # Get gene methylation values for common patients
            gene_meth = row[common_patients].values
            
            # Create DataFrame with methylation, age, sex
            gene_df = pd.DataFrame({
                'PATIENT_ID': common_patients,
                'Methylation': gene_meth
            })
            
            # Merge with clinical data
            gene_df = gene_df.merge(clinical_filtered[['PATIENT_ID', age_col, 'Age_Group', 'Sex']], 
                                   on='PATIENT_ID', how='inner')
            
            # Remove NaN methylation values
            gene_df = gene_df.dropna(subset=['Methylation'])
            
            if len(gene_df) < 10:  # Need at least 10 samples
                continue
            
            # Analyze for each sex separately
            for sex in ['Male', 'Female']:
                sex_data = gene_df[gene_df['Sex'] == sex].copy()
                
                if len(sex_data) < 5:  # Need at least 5 samples per sex
                    continue
                
                # Calculate mean methylation by age group
                age_group_means = sex_data.groupby('Age_Group')['Methylation'].mean()
                
                # Check if we have all three age groups
                if len(age_group_means) < 2:
                    continue
                
                # Calculate trend (increase with age)
                # Use Spearman correlation with age as continuous variable
                try:
                    spearman_corr, spearman_p = spearmanr(sex_data[age_col], sex_data['Methylation'])
                except:
                    spearman_corr, spearman_p = np.nan, np.nan
                
                # Calculate linear regression slope
                try:
                    slope, intercept, r_value, p_value, std_err = linregress(
                        sex_data[age_col], sex_data['Methylation'])
                except:
                    slope, intercept, r_value, p_value, std_err = np.nan, np.nan, np.nan, np.nan, np.nan
                
                # Calculate mean methylation for each age group
                young_mean = age_group_means.get('Young', np.nan)
                mid_mean = age_group_means.get('Mid', np.nan)
                old_mean = age_group_means.get('Old', np.nan)
                
                # Calculate increase from Young to Old
                if not pd.isna(young_mean) and not pd.isna(old_mean):
                    increase_young_to_old = old_mean - young_mean
                else:
                    increase_young_to_old = np.nan
                
                # Calculate increase from Mid to Old
                if not pd.isna(mid_mean) and not pd.isna(old_mean):
                    increase_mid_to_old = old_mean - mid_mean
                else:
                    increase_mid_to_old = np.nan
                
                # Kruskal-Wallis test for differences between age groups
                try:
                    young_vals = sex_data[sex_data['Age_Group'] == 'Young']['Methylation'].values
                    mid_vals = sex_data[sex_data['Age_Group'] == 'Mid']['Methylation'].values
                    old_vals = sex_data[sex_data['Age_Group'] == 'Old']['Methylation'].values
                    
                    # Only test if we have at least 2 groups with data
                    groups_to_test = []
                    group_names = []
                    if len(young_vals) > 0:
                        groups_to_test.append(young_vals)
                        group_names.append('Young')
                    if len(mid_vals) > 0:
                        groups_to_test.append(mid_vals)
                        group_names.append('Mid')
                    if len(old_vals) > 0:
                        groups_to_test.append(old_vals)
                        group_names.append('Old')
                    
                    if len(groups_to_test) >= 2:
                        kruskal_stat, kruskal_p = kruskal(*groups_to_test)
                    else:
                        kruskal_stat, kruskal_p = np.nan, np.nan
                except:
                    kruskal_stat, kruskal_p = np.nan, np.nan
                
                # Store results
                results.append({
                    'Gene': gene,
                    'Sex': sex,
                    'N_Samples': len(sex_data),
                    'Young_Mean': young_mean,
                    'Mid_Mean': mid_mean,
                    'Old_Mean': old_mean,
                    'Increase_Young_to_Old': increase_young_to_old,
                    'Increase_Mid_to_Old': increase_mid_to_old,
                    'Spearman_Correlation': spearman_corr,
                    'Spearman_P': spearman_p,
                    'Linear_Slope': slope,
                    'Linear_R_Squared': r_value**2 if not pd.isna(r_value) else np.nan,
                    'Linear_P': p_value,
                    'Kruskal_Statistic': kruskal_stat,
                    'Kruskal_P': kruskal_p,
                    'Mean_Methylation': sex_data['Methylation'].mean(),
                    'Std_Methylation': sex_data['Methylation'].std()
                })
        
        except Exception as e:
            continue
    
    if (batch_start // batch_size + 1) % 10 == 0:
        print(f"  Processed {batch_end}/{total_genes} genes...")

results_df = pd.DataFrame(results)
print(f"\nTotal gene-sex combinations analyzed: {len(results_df)}")
print()

# ============================================================================
# 4. FILTER FOR GENES WITH AGE-RELATED INCREASE
# ============================================================================

print("Step 4: Identifying genes with age-related increase...")
print("-" * 80)

# Filter for genes with positive correlation and significant p-value
increasing_genes = results_df[
    (results_df['Spearman_Correlation'] > 0) & 
    (results_df['Spearman_P'] < 0.05) &
    (results_df['Increase_Young_to_Old'] > 0)
].copy()

increasing_genes = increasing_genes.sort_values('Spearman_Correlation', ascending=False)

print(f"Genes with significant age-related increase: {len(increasing_genes)}")
print()

# Separate by sex
male_increasing = increasing_genes[increasing_genes['Sex'] == 'Male'].copy()
female_increasing = increasing_genes[increasing_genes['Sex'] == 'Female'].copy()

print(f"Male-specific increasing genes: {len(male_increasing)}")
print(f"Female-specific increasing genes: {len(female_increasing)}")
print()

# Save results
results_df.to_csv(output_dir / "age_sex_methylation_analysis_all.csv", index=False)
increasing_genes.to_csv(output_dir / "age_sex_methylation_increasing_genes.csv", index=False)
male_increasing.to_csv(output_dir / "age_sex_methylation_increasing_male.csv", index=False)
female_increasing.to_csv(output_dir / "age_sex_methylation_increasing_female.csv", index=False)

# ============================================================================
# 5. TOP GENES WITH AGE-RELATED INCREASE
# ============================================================================

print("Step 5: Top genes with age-related increase...")
print("-" * 80)

print("\nTop 20 Genes - Male (Age-Related Increase):")
print("="*80)
if len(male_increasing) > 0:
    top_male = male_increasing.head(20)
    for i, (idx, row) in enumerate(top_male.iterrows(), 1):
        print(f"{i:3d}. {row['Gene']:30s} | Corr: {row['Spearman_Correlation']:7.4f} | "
              f"P: {row['Spearman_P']:.4f} | Increase: {row['Increase_Young_to_Old']:7.4f} | "
              f"Old Mean: {row['Old_Mean']:.4f}")
else:
    print("No significant male-specific increasing genes found.")

print("\nTop 20 Genes - Female (Age-Related Increase):")
print("="*80)
if len(female_increasing) > 0:
    top_female = female_increasing.head(20)
    for i, (idx, row) in enumerate(top_female.iterrows(), 1):
        print(f"{i:3d}. {row['Gene']:30s} | Corr: {row['Spearman_Correlation']:7.4f} | "
              f"P: {row['Spearman_P']:.4f} | Increase: {row['Increase_Young_to_Old']:7.4f} | "
              f"Old Mean: {row['Old_Mean']:.4f}")
else:
    print("No significant female-specific increasing genes found.")

print()

# ============================================================================
# 6. VISUALIZATION
# ============================================================================

print("Step 6: Creating visualizations...")
print("-" * 80)

# Create comprehensive visualization (G and H panels removed for Figure 3)
fig = plt.figure(figsize=(20, 10))
gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

# 1. Top increasing genes - Male (max 2 gene names)
ax1 = fig.add_subplot(gs[0, 0])
if len(male_increasing) > 0:
    top_male_20 = male_increasing.head(20)
    bars = ax1.barh(range(len(top_male_20)), top_male_20['Spearman_Correlation'],
                   color='steelblue', alpha=0.7)
    ax1.set_yticks(range(len(top_male_20)))
    # Show only first 2 gene names, rest empty
    gene_labels = [''] * len(top_male_20)
    if len(top_male_20) >= 1:
        gene_labels[0] = top_male_20.iloc[0]['Gene']
    if len(top_male_20) >= 2:
        gene_labels[1] = top_male_20.iloc[1]['Gene']
    ax1.set_yticklabels(gene_labels, fontsize=8)
    ax1.set_xlabel('Spearman Correlation', fontweight='bold')
    ax1.set_title('A. Top 20 Male Genes\n(Age-Related Increase)', fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='x')
    ax1.invert_yaxis()
else:
    ax1.text(0.5, 0.5, 'No significant\nmale genes', ha='center', va='center',
            transform=ax1.transAxes, fontsize=12)
    ax1.set_title('A. Top 20 Male Genes', fontweight='bold')

# 2. Top increasing genes - Female (max 2 gene names)
ax2 = fig.add_subplot(gs[0, 1])
if len(female_increasing) > 0:
    top_female_20 = female_increasing.head(20)
    bars = ax2.barh(range(len(top_female_20)), top_female_20['Spearman_Correlation'],
                   color='pink', alpha=0.7)
    ax2.set_yticks(range(len(top_female_20)))
    # Show only first 2 gene names, rest empty
    gene_labels = [''] * len(top_female_20)
    if len(top_female_20) >= 1:
        gene_labels[0] = top_female_20.iloc[0]['Gene']
    if len(top_female_20) >= 2:
        gene_labels[1] = top_female_20.iloc[1]['Gene']
    ax2.set_yticklabels(gene_labels, fontsize=8)
    ax2.set_xlabel('Spearman Correlation', fontweight='bold')
    ax2.set_title('B. Top 20 Female Genes\n(Age-Related Increase)', fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='x')
    ax2.invert_yaxis()
else:
    ax2.text(0.5, 0.5, 'No significant\nfemale genes', ha='center', va='center',
            transform=ax2.transAxes, fontsize=12)
    ax2.set_title('B. Top 20 Female Genes', fontweight='bold')

# 3. Increase from Young to Old - Male (max 2 gene names)
ax3 = fig.add_subplot(gs[0, 2])
if len(male_increasing) > 0:
    top_male_20 = male_increasing.head(20)
    bars = ax3.barh(range(len(top_male_20)), top_male_20['Increase_Young_to_Old'],
                   color='steelblue', alpha=0.7)
    ax3.set_yticks(range(len(top_male_20)))
    # Show only first 2 gene names, rest empty
    gene_labels = [''] * len(top_male_20)
    if len(top_male_20) >= 1:
        gene_labels[0] = top_male_20.iloc[0]['Gene']
    if len(top_male_20) >= 2:
        gene_labels[1] = top_male_20.iloc[1]['Gene']
    ax3.set_yticklabels(gene_labels, fontsize=8)
    ax3.set_xlabel('Methylation Increase\n(Old - Young)', fontweight='bold')
    ax3.set_title('C. Methylation Increase\n(Male: Young→Old)', fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='x')
    ax3.invert_yaxis()
else:
    ax3.text(0.5, 0.5, 'No data', ha='center', va='center',
            transform=ax3.transAxes, fontsize=12)
    ax3.set_title('C. Methylation Increase (Male)', fontweight='bold')

# 4. Increase from Young to Old - Female (max 2 gene names)
ax4 = fig.add_subplot(gs[1, 0])
if len(female_increasing) > 0:
    top_female_20 = female_increasing.head(20)
    bars = ax4.barh(range(len(top_female_20)), top_female_20['Increase_Young_to_Old'],
                   color='pink', alpha=0.7)
    ax4.set_yticks(range(len(top_female_20)))
    # Show only first 2 gene names, rest empty
    gene_labels = [''] * len(top_female_20)
    if len(top_female_20) >= 1:
        gene_labels[0] = top_female_20.iloc[0]['Gene']
    if len(top_female_20) >= 2:
        gene_labels[1] = top_female_20.iloc[1]['Gene']
    ax4.set_yticklabels(gene_labels, fontsize=8)
    ax4.set_xlabel('Methylation Increase\n(Old - Young)', fontweight='bold')
    ax4.set_title('D. Methylation Increase\n(Female: Young→Old)', fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='x')
    ax4.invert_yaxis()
else:
    ax4.text(0.5, 0.5, 'No data', ha='center', va='center',
            transform=ax4.transAxes, fontsize=12)
    ax4.set_title('D. Methylation Increase (Female)', fontweight='bold')

# 5. Age group means comparison - Male (max 2 gene names)
ax5 = fig.add_subplot(gs[1, 1])
if len(male_increasing) > 0:
    top_male_10 = male_increasing.head(10)
    x_pos = np.arange(len(top_male_10))
    width = 0.25
    
    young_means = top_male_10['Young_Mean'].fillna(0)
    mid_means = top_male_10['Mid_Mean'].fillna(0)
    old_means = top_male_10['Old_Mean'].fillna(0)
    
    ax5.bar(x_pos - width, young_means, width, label='Young', alpha=0.7, color='lightblue')
    ax5.bar(x_pos, mid_means, width, label='Mid', alpha=0.7, color='steelblue')
    ax5.bar(x_pos + width, old_means, width, label='Old', alpha=0.7, color='darkblue')
    
    ax5.set_xlabel('Genes', fontweight='bold')
    ax5.set_ylabel('Mean Methylation', fontweight='bold')
    ax5.set_title('E. Age Group Means\n(Top 10 Male Genes)', fontweight='bold')
    ax5.set_xticks(x_pos)
    # Show only first 2 gene names, rest empty
    gene_labels = [''] * len(top_male_10)
    if len(top_male_10) >= 1:
        gene_labels[0] = top_male_10.iloc[0]['Gene']
    if len(top_male_10) >= 2:
        gene_labels[1] = top_male_10.iloc[1]['Gene']
    ax5.set_xticklabels(gene_labels, rotation=45, ha='right', fontsize=7)
    ax5.legend()
    ax5.grid(True, alpha=0.3, axis='y')
else:
    ax5.text(0.5, 0.5, 'No data', ha='center', va='center',
            transform=ax5.transAxes, fontsize=12)
    ax5.set_title('E. Age Group Means (Male)', fontweight='bold')

# 6. Age group means comparison - Female (max 2 gene names)
ax6 = fig.add_subplot(gs[1, 2])
if len(female_increasing) > 0:
    top_female_10 = female_increasing.head(10)
    x_pos = np.arange(len(top_female_10))
    width = 0.25
    
    young_means = top_female_10['Young_Mean'].fillna(0)
    mid_means = top_female_10['Mid_Mean'].fillna(0)
    old_means = top_female_10['Old_Mean'].fillna(0)
    
    ax6.bar(x_pos - width, young_means, width, label='Young', alpha=0.7, color='lightpink')
    ax6.bar(x_pos, mid_means, width, label='Mid', alpha=0.7, color='pink')
    ax6.bar(x_pos + width, old_means, width, label='Old', alpha=0.7, color='hotpink')
    
    ax6.set_xlabel('Genes', fontweight='bold')
    ax6.set_ylabel('Mean Methylation', fontweight='bold')
    ax6.set_title('F. Age Group Means\n(Top 10 Female Genes)', fontweight='bold')
    ax6.set_xticks(x_pos)
    # Show only first 2 gene names, rest empty
    gene_labels = [''] * len(top_female_10)
    if len(top_female_10) >= 1:
        gene_labels[0] = top_female_10.iloc[0]['Gene']
    if len(top_female_10) >= 2:
        gene_labels[1] = top_female_10.iloc[1]['Gene']
    ax6.set_xticklabels(gene_labels, rotation=45, ha='right', fontsize=7)
    ax6.legend()
    ax6.grid(True, alpha=0.3, axis='y')
else:
    ax6.text(0.5, 0.5, 'No data', ha='center', va='center',
            transform=ax6.transAxes, fontsize=12)
    ax6.set_title('F. Age Group Means (Female)', fontweight='bold')

# G and H panels removed for Figure 3 (information moved to caption)
plt.savefig(plots_dir / "age_sex_methylation_analysis.png", dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Visualization saved: age_sex_methylation_analysis.png")

# ============================================================================
# SUMMARY
# ============================================================================

print()
print("="*80)
print("Analysis Summary")
print("="*80)
print()
print(f"Total gene-sex combinations analyzed: {len(results_df):,}")
print(f"Genes with significant age-related increase (p<0.05): {len(increasing_genes)}")
print(f"  - Male-specific: {len(male_increasing)}")
print(f"  - Female-specific: {len(female_increasing)}")
print()
print("Files saved:")
print(f"- results/age_sex_methylation_analysis_all.csv")
print(f"- results/age_sex_methylation_increasing_genes.csv")
print(f"- results/age_sex_methylation_increasing_male.csv")
print(f"- results/age_sex_methylation_increasing_female.csv")
print(f"- plots/age_sex_methylation_analysis.png")
print()
print("="*80)
print("Analysis completed successfully!")
print("="*80)











