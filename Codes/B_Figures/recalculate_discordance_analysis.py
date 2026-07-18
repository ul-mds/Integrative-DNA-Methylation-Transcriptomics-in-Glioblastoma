#!/usr/bin/env python3
"""
Recalculate Discordance Analysis with proper sorting
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path

script_dir = Path(__file__).resolve().parent
output_dir = script_dir / "results"
plots_dir = script_dir / "plots"

# Load data
print("Loading data...")
methylation_file = "../methylation_epic.txt"
protein_file = "../protein_quantification.txt"

# Load methylation (sample first 20000 rows for speed)
methylation = pd.read_csv(methylation_file, sep="\t", nrows=20000)
methylation.columns = methylation.columns.str.strip()

# Load protein
protein = pd.read_csv(protein_file, sep="\t")
protein.columns = protein.columns.str.strip()

# Get common patients
meth_patients = [col for col in methylation.columns if col.startswith('C3')]
protein_patients = [col for col in protein.columns if col.startswith('C3')]
common_patients = sorted(list(set(meth_patients) & set(protein_patients)))

print(f"Common patients: {len(common_patients)}")

# Prepare dataframes
meth_df = methylation.set_index('Hugo_Symbol')[common_patients]
protein_df = protein.set_index('Composite.Element.REF').str.split('|').str[0]
protein_df = protein_df.reset_index()
protein_df = protein_df.groupby('Composite.Element.REF')[common_patients].first()
protein_df.index = protein_df.index.str.split('|').str[0]

# Get common genes
common_genes = sorted(list(set(meth_df.index) & set(protein_df.index)))
print(f"Common genes: {len(common_genes)}")

# Calculate discordance
print("Calculating discordance...")
results = []

for gene in common_genes[:500]:  # Analyze first 500 genes
    if gene not in meth_df.index or gene not in protein_df.index:
        continue
    
    try:
        meth_data = meth_df.loc[gene, common_patients].values
        protein_data = protein_df.loc[gene, common_patients].values
        
        valid_idx = ~(np.isnan(meth_data) | np.isnan(protein_data))
        if valid_idx.sum() < 10:
            continue
        
        meth_data = meth_data[valid_idx]
        protein_data = protein_data[valid_idx]
        
        # Calculate correlation
        try:
            corr, p_val = spearmanr(meth_data, protein_data)
        except:
            corr, p_val = np.nan, np.nan
        
        if np.isnan(corr):
            continue
        
        # Classify discordance
        abs_corr = abs(corr)
        if abs_corr < 0.2:
            disc_type = "Discordant"
        elif corr > 0.5:
            disc_type = "Concordant_High"
        elif corr < -0.5:
            disc_type = "Concordant_Negative"
        else:
            disc_type = "Moderate"
        
        results.append({
            'Gene': gene,
            'Meth_Protein_Corr': corr,
            'Abs_Meth_Protein_Corr': abs_corr,
            'P_value': p_val,
            'Discordance_Type': disc_type,
            'Meth_Mean': np.mean(meth_data),
            'Protein_Mean': np.mean(protein_data)
        })
    except Exception as e:
        continue

discordance_df = pd.DataFrame(results)
print(f"Total genes analyzed: {len(discordance_df)}")
print(f"Discordant genes: {len(discordance_df[discordance_df['Discordance_Type'] == 'Discordant'])}")

# Save results
discordance_df.to_csv(output_dir / "methylation_protein_discordance_recalculated.csv", index=False)

# Get top 20 discordant genes (lowest absolute correlation)
discordant_genes = discordance_df[discordance_df['Discordance_Type'] == 'Discordant'].copy()
top_20_discordant = discordant_genes.nsmallest(20, 'Abs_Meth_Protein_Corr')

print("\nTop 20 Discordant Genes:")
print(top_20_discordant[['Gene', 'Meth_Protein_Corr', 'Abs_Meth_Protein_Corr', 'Meth_Mean', 'Protein_Mean']].to_string(index=False))

# Check protein mean range
print(f"\nProtein Mean Range:")
print(f"  Min: {discordance_df['Protein_Mean'].min():.4f}")
print(f"  Max: {discordance_df['Protein_Mean'].max():.4f}")
print(f"  Mean: {discordance_df['Protein_Mean'].mean():.4f}")
print(f"  Std: {discordance_df['Protein_Mean'].std():.4f}")

# Check if all top 20 start with A
top_20_genes = top_20_discordant['Gene'].tolist()
a_starting = [g for g in top_20_genes if g.startswith('A')]
print(f"\nTop 20 genes starting with 'A': {len(a_starting)}/{len(top_20_genes)}")
