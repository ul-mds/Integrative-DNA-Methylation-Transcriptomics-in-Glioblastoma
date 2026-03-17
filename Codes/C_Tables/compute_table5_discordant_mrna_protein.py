#!/usr/bin/env python3
"""
Compute mRNA-protein correlation for ALL discordant genes (from methylation_protein_discordance_full).
Then select top 10 by mRNA-protein significance - ensuring genes span the alphabet (not just A-genes).
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from pathlib import Path

script_dir = Path(__file__).resolve().parent
base_dir = script_dir.parent
output_dir = script_dir / "results"

# Load discordant genes
meth_full = pd.read_csv(output_dir / "methylation_protein_discordance_full.csv")
discordant = meth_full[meth_full['Discordance_Type'] == 'Discordant'].copy()
discordant_genes = discordant['Gene'].tolist()
print(f"Discordant genes: {len(discordant_genes)}")

# Load mRNA and protein data
print("Loading mRNA data...")
mrna = pd.read_csv(base_dir / "mrna_seq_fpkm.txt", sep="\t")
mrna.columns = mrna.columns.str.strip()
mrna = mrna.set_index('Hugo_Symbol')

print("Loading protein data...")
protein = pd.read_csv(base_dir / "protein_quantification.txt", sep="\t")
protein.columns = protein.columns.str.strip()
protein['Gene'] = protein['Composite.Element.REF'].str.split('|').str[0]
protein = protein.drop(columns=['Composite.Element.REF']).set_index('Gene')

# Get common patients
mrna_patients = [c for c in mrna.columns if c.startswith('C3')]
protein_patients = [c for c in protein.columns if c.startswith('C3')]
common_patients = sorted(list(set(mrna_patients) & set(protein_patients)))
print(f"Common patients: {len(common_patients)}")

# Compute mRNA-protein for each discordant gene
results = []
for gene in discordant_genes:
    if gene not in mrna.index or gene not in protein.index:
        continue
    try:
        mrna_vals = mrna.loc[gene, common_patients].values
        prot_vals = protein.loc[gene, common_patients].values
        if hasattr(mrna_vals, '__len__') and len(mrna_vals.shape) > 1:
            mrna_vals = mrna_vals.flatten()
        valid = ~(np.isnan(mrna_vals.astype(float)) | np.isnan(prot_vals.astype(float)))
        if valid.sum() < 10:
            continue
        mrna_vals = np.array(mrna_vals, dtype=float)[valid]
        prot_vals = np.array(prot_vals, dtype=float)[valid]
        r, p = spearmanr(mrna_vals, prot_vals)
        if np.isnan(p) or p <= 0:
            continue
        row = discordant[discordant['Gene'] == gene].iloc[0]
        results.append({
            'Gene': gene,
            'Meth_Protein_Corr': row['Meth_Protein_Corr'],
            'P_value': row['P_value'],
            'mRNA_Protein_Corr': r,
            'mRNA_Protein_P': p,
            'Meth_Mean': row['Meth_Mean'],
            'Protein_Mean': row['Protein_Mean']
        })
    except Exception as e:
        continue

df = pd.DataFrame(results)
df = df.sort_values('mRNA_Protein_P', ascending=True)
top10 = df.head(10)

print(f"\nGenes with mRNA-protein data: {len(df)}")
print("\nTop 10 by mRNA-protein significance:")
print(top10[['Gene', 'Meth_Protein_Corr', 'mRNA_Protein_Corr', 'mRNA_Protein_P']].to_string())
print("\nFirst letters:", [g[0] for g in top10['Gene'].tolist()])

# Save
top10_out = top10.copy()
top10_out.insert(0, 'Rank', range(1, 11))
top10_out.to_csv(output_dir / "table5_top10_discordant_by_mrna_protein_significance.csv", index=False)
print(f"\nSaved to {output_dir / 'table5_top10_discordant_by_mrna_protein_significance.csv'}")
