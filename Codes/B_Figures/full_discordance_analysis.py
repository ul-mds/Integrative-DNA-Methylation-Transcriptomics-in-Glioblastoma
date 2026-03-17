#!/usr/bin/env python3
"""
Full Discordance Analysis - Analyze ALL genes, not just A-letter genes
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path

script_dir = Path(__file__).resolve().parent
output_dir = script_dir / "results"
plots_dir = script_dir / "plots"

print("="*80)
print("FULL DISCORDANCE ANALYSIS - All Genes")
print("="*80)
print()

# Load data
print("Loading data...")
methylation_file = "../methylation_epic.txt"
protein_file = "../protein_quantification.txt"

# Load methylation (read in chunks to handle large file)
print("Loading methylation data (this may take a while)...")
methylation_chunks = []
chunk_size = 50000
for i, chunk in enumerate(pd.read_csv(methylation_file, sep="\t", chunksize=chunk_size)):
    methylation_chunks.append(chunk)
    if i % 10 == 0:
        print(f"  Loaded {i*chunk_size:,} rows...")
methylation = pd.concat(methylation_chunks, ignore_index=True)
print(f"Methylation data loaded: {methylation.shape}")

# Load protein
print("Loading protein data...")
protein = pd.read_csv(protein_file, sep="\t")
print(f"Protein data loaded: {protein.shape}")

# Clean column names
methylation.columns = methylation.columns.str.strip()
protein.columns = protein.columns.str.strip()

# Get common patients
meth_patients = [col for col in methylation.columns if col.startswith('C3')]
protein_patients = [col for col in protein.columns if col.startswith('C3')]
common_patients = sorted(list(set(meth_patients) & set(protein_patients)))
print(f"Common patients: {len(common_patients)}")

# Prepare methylation dataframe
meth_df = methylation.set_index('Hugo_Symbol')[common_patients]
print(f"Methylation genes: {len(meth_df)}")

# Prepare protein dataframe
protein_df = protein.copy()
if 'Composite.Element.REF' in protein_df.columns:
    protein_df['Gene'] = protein_df['Composite.Element.REF'].str.split('|').str[0]
    protein_df = protein_df.groupby('Gene')[common_patients].mean()  # Average if multiple probes per gene
else:
    protein_df = protein_df.set_index(protein_df.columns[0])
    protein_df.index = protein_df.index.str.split('|').str[0]
    protein_df = protein_df.groupby(protein_df.index)[common_patients].mean()

print(f"Protein genes: {len(protein_df)}")

# Get common genes
common_genes = sorted(list(set(meth_df.index) & set(protein_df.index)))
print(f"Common genes: {len(common_genes)}")
print()

# Calculate discordance for ALL genes
print("Calculating discordance for all genes...")
print("-" * 80)

results = []
processed = 0

for gene in common_genes:
    if gene not in meth_df.index or gene not in protein_df.index:
        continue
    
    try:
        meth_data = meth_df.loc[gene, common_patients].values
        protein_data = protein_df.loc[gene, common_patients].values
        
        valid_idx = ~(np.isnan(meth_data) | np.isnan(protein_data))
        if valid_idx.sum() < 10:  # Need at least 10 patients
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
        
        processed += 1
        if processed % 1000 == 0:
            print(f"  Processed {processed:,} genes...")
            
    except Exception as e:
        continue

discordance_df = pd.DataFrame(results)
print(f"\nTotal genes analyzed: {len(discordance_df):,}")

# Summary statistics
print("\nDiscordance Type Distribution:")
print(discordance_df['Discordance_Type'].value_counts())
print()

discordant_genes = discordance_df[discordance_df['Discordance_Type'] == 'Discordant']
print(f"Discordant genes: {len(discordant_genes):,}")

# Get top 20 discordant genes (lowest absolute correlation)
top_20_discordant = discordant_genes.nsmallest(20, 'Abs_Meth_Protein_Corr')
print("\nTop 20 Discordant Genes (Lowest |Correlation|):")
print("="*80)
print(top_20_discordant[['Gene', 'Meth_Protein_Corr', 'Abs_Meth_Protein_Corr', 'Meth_Mean', 'Protein_Mean']].to_string(index=False))

# Check first letters
first_letters = [g[0] for g in top_20_discordant['Gene'].tolist()]
from collections import Counter
letter_counts = Counter(first_letters)
print(f"\nFirst letter distribution in top 20:")
for letter, count in sorted(letter_counts.items()):
    print(f"  {letter}: {count}")

# Protein statistics
print(f"\nProtein Mean Statistics:")
print(f"  Min: {discordance_df['Protein_Mean'].min():.4f}")
print(f"  Max: {discordance_df['Protein_Mean'].max():.4f}")
print(f"  Mean: {discordance_df['Protein_Mean'].mean():.4f}")
print(f"  Std: {discordance_df['Protein_Mean'].std():.4f}")

# Save results
discordance_df.to_csv(output_dir / "methylation_protein_discordance_full.csv", index=False)
top_20_discordant.to_csv(output_dir / "top_20_discordant_genes_full.csv", index=False)

print(f"\n✓ Results saved:")
print(f"  - {output_dir / 'methylation_protein_discordance_full.csv'}")
print(f"  - {output_dir / 'top_20_discordant_genes_full.csv'}")
