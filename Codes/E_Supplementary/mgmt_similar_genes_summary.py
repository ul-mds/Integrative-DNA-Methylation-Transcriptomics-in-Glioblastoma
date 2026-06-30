#!/usr/bin/env python3
"""
Create a clean summary of MGMT-similar genes
"""

import pandas as pd
import numpy as np

# Load results
similarity_all = pd.read_csv("results/mgmt_similarity_analysis_all_genes.csv")

# Filter for genes with at least one valid similarity metric
valid_genes = similarity_all[
    (~similarity_all['Spearman_Correlation'].isna()) | 
    (~similarity_all['Cosine_Similarity'].isna()) |
    (~similarity_all['Pearson_Correlation'].isna())
].copy()

# Sort by composite score (best first)
valid_genes = valid_genes.sort_values('Composite_Similarity_Score', ascending=False)

# Get top 50
top_50 = valid_genes.head(50)

# Create clean summary
summary = top_50[['Gene', 'Pearson_Correlation', 'Spearman_Correlation', 
                  'Cosine_Similarity', 'Composite_Similarity_Score',
                  'Mean_Difference', 'Gene_Methylation_Mean', 'MGMT_Methylation_Mean']].copy()

# Round values
for col in ['Pearson_Correlation', 'Spearman_Correlation', 'Cosine_Similarity', 
           'Composite_Similarity_Score', 'Mean_Difference', 'Gene_Methylation_Mean', 
           'MGMT_Methylation_Mean']:
    summary[col] = summary[col].round(4)

# Replace NaN with "N/A"
summary = summary.fillna('N/A')

print("="*100)
print("MGMT BENZERLİK ANALİZİ - TOP 50 GEN")
print("="*100)
print()
print(f"Toplam analiz edilen gen: {len(similarity_all):,}")
print(f"Geçerli benzerlik metriklerine sahip gen: {len(valid_genes):,}")
print()
print("TOP 50 GEN - MGMT'ye En Benzer Metilasyon Profilleri:")
print("="*100)
print()

for i, (idx, row) in enumerate(top_50.iterrows(), 1):
    pearson = f"{row['Pearson_Correlation']:.4f}" if not pd.isna(row['Pearson_Correlation']) else "N/A"
    spearman = f"{row['Spearman_Correlation']:.4f}" if not pd.isna(row['Spearman_Correlation']) else "N/A"
    cosine = f"{row['Cosine_Similarity']:.4f}" if not pd.isna(row['Cosine_Similarity']) else "N/A"
    score = f"{row['Composite_Similarity_Score']:.4f}" if not pd.isna(row['Composite_Similarity_Score']) else "0.0000"
    mean_diff = f"{row['Mean_Difference']:.4f}" if not pd.isna(row['Mean_Difference']) else "N/A"
    
    print(f"{i:3d}. {row['Gene']:30s} | Spearman: {spearman:8s} | Cosine: {cosine:8s} | Score: {score:8s} | Mean Diff: {mean_diff:8s}")

print()
print("="*100)

# Save clean summary
summary.to_csv("results/mgmt_similar_genes_top50_clean_summary.csv", index=False)
print(f"\nClean summary saved to: results/mgmt_similar_genes_top50_clean_summary.csv")











