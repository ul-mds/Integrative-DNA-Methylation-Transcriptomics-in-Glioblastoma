#!/usr/bin/env python3
"""
MGMT Methylation Profile Similarity Analysis
Find genes with methylation profiles most similar to MGMT
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
from scipy.spatial.distance import correlation
from sklearn.metrics.pairwise import cosine_similarity
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
print("MGMT Methylation Profile Similarity Analysis")
print("="*80)
print()

# ============================================================================
# 1. LOAD DATA
# ============================================================================

print("Step 1: Loading methylation data...")
print("-" * 80)

# Load methylation data
methylation_file = "../methylation_epic.txt"
print("Loading methylation data (this may take a while)...")

# Read methylation data - use larger sample for better coverage
# First check how many genes we have
print("Reading methylation data header...")
methylation_sample = pd.read_csv(methylation_file, sep="\t", nrows=1000)
print(f"Sample shape: {methylation_sample.shape}")

# Read full data but limit to reasonable size for efficiency
methylation = pd.read_csv(methylation_file, sep="\t", nrows=100000)  # Increased for better coverage
print(f"Methylation data loaded: {methylation.shape}")

# Get patient IDs
patients = [col for col in methylation.columns if col.startswith('C3')]
print(f"Number of patients: {len(patients)}")

# Check if MGMT exists
mgmt_genes = methylation[methylation['Hugo_Symbol'].str.contains('MGMT', case=False, na=False)]
print(f"MGMT-related genes found: {len(mgmt_genes)}")
if len(mgmt_genes) > 0:
    print(f"MGMT gene symbols: {mgmt_genes['Hugo_Symbol'].unique()}")
    mgmt_gene = mgmt_genes['Hugo_Symbol'].iloc[0]  # Use first MGMT gene
    print(f"Using MGMT gene: {mgmt_gene}")
else:
    print("ERROR: MGMT gene not found in methylation data!")
    exit(1)

print()

# ============================================================================
# 2. EXTRACT MGMT METHYLATION PROFILE
# ============================================================================

print("Step 2: Extracting MGMT methylation profile...")
print("-" * 80)

# Get MGMT methylation values
mgmt_meth = methylation[methylation['Hugo_Symbol'] == mgmt_gene]
if len(mgmt_meth) == 0:
    # Try case-insensitive search
    mgmt_meth = methylation[methylation['Hugo_Symbol'].str.upper() == 'MGMT']
    
if len(mgmt_meth) == 0:
    print("ERROR: Could not find MGMT methylation data!")
    exit(1)

mgmt_profile_raw = mgmt_meth[patients].iloc[0].values
# Fill NaN with median for MGMT (more robust than mean for beta values)
mgmt_median = np.nanmedian(mgmt_profile_raw)
mgmt_profile = np.nan_to_num(mgmt_profile_raw, nan=mgmt_median, posinf=mgmt_median, neginf=mgmt_median)

print(f"MGMT profile: {len(mgmt_profile)} values, {np.sum(~np.isnan(mgmt_profile_raw))} valid (raw)")
print(f"MGMT methylation range: {mgmt_profile.min():.3f} - {mgmt_profile.max():.3f}")
print(f"MGMT methylation mean: {mgmt_profile.mean():.3f}, median: {mgmt_median:.3f}")
print(f"MGMT methylation profile extracted successfully!")
print()

# ============================================================================
# 3. CALCULATE SIMILARITY TO MGMT
# ============================================================================

print("Step 3: Calculating similarity to MGMT profile...")
print("-" * 80)

similarity_results = []

# Process genes in batches for efficiency
batch_size = 1000
total_genes = len(methylation)

print(f"Processing {total_genes} genes in batches of {batch_size}...")

for batch_start in range(0, total_genes, batch_size):
    batch_end = min(batch_start + batch_size, total_genes)
    batch_genes = methylation.iloc[batch_start:batch_end]
    
    for idx, row in batch_genes.iterrows():
        gene = row['Hugo_Symbol']
        
        # Skip if gene is NaN
        if pd.isna(gene):
            continue
        
        # Skip MGMT itself
        try:
            if str(gene).upper() == 'MGMT' or 'MGMT' in str(gene).upper():
                continue
        except:
            continue
        
        try:
            # Get gene methylation profile
            gene_profile_raw = row[patients].values
            
            # Fill NaN values with gene median (more robust)
            gene_median = np.nanmedian(gene_profile_raw)
            if np.isnan(gene_median):
                gene_median = 0.5  # Default beta value
            gene_profile = np.nan_to_num(gene_profile_raw, nan=gene_median, posinf=gene_median, neginf=gene_median)
            
            # Skip if no variation (all values same)
            if np.std(gene_profile) < 0.01:
                continue
            
            # Use filled profiles for correlation
            mgmt_valid = mgmt_profile
            gene_valid = gene_profile
                
            if np.std(gene_valid) < 0.01 or np.std(mgmt_valid) < 0.01:
                continue
            
            # Calculate similarity metrics using valid values only
            # 1. Pearson correlation
            try:
                pearson_corr, pearson_p = pearsonr(mgmt_valid, gene_valid)
            except:
                pearson_corr, pearson_p = np.nan, np.nan
            
            # 2. Spearman correlation
            try:
                spearman_corr, spearman_p = spearmanr(mgmt_valid, gene_valid)
            except:
                spearman_corr, spearman_p = np.nan, np.nan
            
            # 3. Cosine similarity
            try:
                cosine_sim = cosine_similarity([mgmt_valid], [gene_valid])[0][0]
            except:
                cosine_sim = np.nan
            
            # 4. 1 - Correlation distance (higher = more similar)
            try:
                corr_distance = correlation(mgmt_valid, gene_valid)
                similarity_score = 1 - corr_distance  # Convert distance to similarity
            except:
                corr_distance = np.nan
                similarity_score = np.nan
            
            # 5. Mean absolute difference (lower = more similar)
            mae = np.mean(np.abs(mgmt_valid - gene_valid))
            
            # 6. Root mean square difference (lower = more similar)
            rmse = np.sqrt(np.mean((mgmt_valid - gene_valid)**2))
            
            # Calculate statistics
            gene_mean = np.mean(gene_valid)
            gene_std = np.std(gene_valid)
            mgmt_mean = np.mean(mgmt_valid)
            
            # Mean difference
            mean_diff = abs(gene_mean - mgmt_mean)
            
            similarity_results.append({
                'Gene': gene,
                'Pearson_Correlation': pearson_corr,
                'Pearson_P': pearson_p,
                'Spearman_Correlation': spearman_corr,
                'Spearman_P': spearman_p,
                'Cosine_Similarity': cosine_sim,
                'Correlation_Similarity': similarity_score,
                'MAE': mae,
                'RMSE': rmse,
                'Mean_Difference': mean_diff,
                'Gene_Methylation_Mean': gene_mean,
                'Gene_Methylation_STD': gene_std,
                'MGMT_Methylation_Mean': mgmt_mean
            })
        except Exception as e:
            continue
    
    if (batch_start // batch_size + 1) % 10 == 0:
        print(f"  Processed {batch_end}/{total_genes} genes...")

similarity_df = pd.DataFrame(similarity_results)
print(f"\nTotal genes analyzed: {len(similarity_df)}")
print()

# ============================================================================
# 4. RANK GENES BY SIMILARITY
# ============================================================================

print("Step 4: Ranking genes by similarity to MGMT...")
print("-" * 80)

# Create composite similarity score only if we have results
if len(similarity_df) == 0:
    print("WARNING: No genes with sufficient valid values found!")
    print("Please check the data or lower the threshold.")
    exit(1)

# Weight: Pearson (0.3), Spearman (0.3), Cosine (0.2), Correlation Similarity (0.2)
similarity_df['Composite_Similarity_Score'] = (
    similarity_df['Pearson_Correlation'].fillna(0) * 0.3 +
    similarity_df['Spearman_Correlation'].fillna(0) * 0.3 +
    similarity_df['Cosine_Similarity'].fillna(0) * 0.2 +
    similarity_df['Correlation_Similarity'].fillna(0) * 0.2
)

# Also rank by individual metrics
similarity_df = similarity_df.sort_values('Composite_Similarity_Score', ascending=False)

# Filter for significant correlations
sig_pearson = similarity_df[(similarity_df['Pearson_P'] < 0.05) & 
                            (~similarity_df['Pearson_Correlation'].isna())]
sig_spearman = similarity_df[(similarity_df['Spearman_P'] < 0.05) & 
                             (~similarity_df['Spearman_Correlation'].isna())]

print(f"Genes with significant Pearson correlation (p<0.05): {len(sig_pearson)}")
print(f"Genes with significant Spearman correlation (p<0.05): {len(sig_spearman)}")
print()

# ============================================================================
# 5. TOP SIMILAR GENES
# ============================================================================

print("Step 5: Top genes similar to MGMT...")
print("-" * 80)

# Top 50 by composite score
top_50 = similarity_df.head(50)

# Top by Pearson correlation
top_pearson = similarity_df.nlargest(50, 'Pearson_Correlation')

# Top by Spearman correlation
top_spearman = similarity_df.nlargest(50, 'Spearman_Correlation')

# Top by Cosine similarity
top_cosine = similarity_df.nlargest(50, 'Cosine_Similarity')

print("\nTop 20 Genes Most Similar to MGMT (Composite Score):")
print("="*80)
top_20_display = top_50.head(20)[['Gene', 'Pearson_Correlation', 'Spearman_Correlation', 
                                  'Cosine_Similarity', 'Composite_Similarity_Score', 
                                  'Mean_Difference', 'Pearson_P']].copy()
top_20_display['Pearson_Correlation'] = top_20_display['Pearson_Correlation'].round(4)
top_20_display['Spearman_Correlation'] = top_20_display['Spearman_Correlation'].round(4)
top_20_display['Cosine_Similarity'] = top_20_display['Cosine_Similarity'].round(4)
top_20_display['Composite_Similarity_Score'] = top_20_display['Composite_Similarity_Score'].round(4)
top_20_display['Mean_Difference'] = top_20_display['Mean_Difference'].round(4)
top_20_display['Pearson_P'] = top_20_display['Pearson_P'].apply(lambda x: f"{x:.2e}" if not pd.isna(x) else "N/A")
print(top_20_display.to_string(index=False))
print()

# Save results
similarity_df.to_csv(output_dir / "mgmt_similarity_analysis_all_genes.csv", index=False)
top_50.to_csv(output_dir / "mgmt_similar_genes_top50.csv", index=False)
top_pearson.to_csv(output_dir / "mgmt_similar_genes_by_pearson.csv", index=False)
top_spearman.to_csv(output_dir / "mgmt_similar_genes_by_spearman.csv", index=False)
top_cosine.to_csv(output_dir / "mgmt_similar_genes_by_cosine.csv", index=False)

print(f"Results saved to {output_dir}")
print()

# ============================================================================
# 6. VISUALIZATION
# ============================================================================

print("Step 6: Creating visualizations...")
print("-" * 80)

fig = plt.figure(figsize=(18, 12))
gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

# 1. Top 20 genes by composite score
ax1 = fig.add_subplot(gs[0, 0])
top_20 = top_50.head(20)
bars = ax1.barh(range(len(top_20)), top_20['Composite_Similarity_Score'], 
               color=plt.cm.viridis(np.linspace(0, 1, len(top_20))))
ax1.set_yticks(range(len(top_20)))
ax1.set_yticklabels(top_20['Gene'], fontsize=9)
ax1.set_xlabel('Composite Similarity Score', fontweight='bold')
ax1.set_title('A. Top 20 Genes Similar to MGMT\n(Composite Score)', fontweight='bold')
ax1.grid(True, alpha=0.3, axis='x')
ax1.invert_yaxis()

# 2. Pearson correlation distribution
ax2 = fig.add_subplot(gs[0, 1])
ax2.hist(similarity_df['Pearson_Correlation'].dropna(), bins=50, edgecolor='black', alpha=0.7)
ax2.axvline(similarity_df['Pearson_Correlation'].mean(), color='red', linestyle='--', 
           linewidth=2, label=f'Mean: {similarity_df["Pearson_Correlation"].mean():.3f}')
ax2.set_xlabel('Pearson Correlation with MGMT', fontweight='bold')
ax2.set_ylabel('Frequency', fontweight='bold')
ax2.set_title('B. Pearson Correlation Distribution', fontweight='bold')
ax2.legend()
ax2.grid(True, alpha=0.3)

# 3. Spearman correlation distribution
ax3 = fig.add_subplot(gs[0, 2])
ax3.hist(similarity_df['Spearman_Correlation'].dropna(), bins=50, edgecolor='black', alpha=0.7, color='green')
ax3.axvline(similarity_df['Spearman_Correlation'].mean(), color='red', linestyle='--', 
           linewidth=2, label=f'Mean: {similarity_df["Spearman_Correlation"].mean():.3f}')
ax3.set_xlabel('Spearman Correlation with MGMT', fontweight='bold')
ax3.set_ylabel('Frequency', fontweight='bold')
ax3.set_title('C. Spearman Correlation Distribution', fontweight='bold')
ax3.legend()
ax3.grid(True, alpha=0.3)

# 4. Top genes by Pearson correlation
ax4 = fig.add_subplot(gs[1, 0])
top_pearson_20 = top_pearson.head(20)
bars = ax4.barh(range(len(top_pearson_20)), top_pearson_20['Pearson_Correlation'],
               color='steelblue', alpha=0.7)
ax4.set_yticks(range(len(top_pearson_20)))
ax4.set_yticklabels(top_pearson_20['Gene'], fontsize=9)
ax4.set_xlabel('Pearson Correlation', fontweight='bold')
ax4.set_title('D. Top 20 by Pearson Correlation', fontweight='bold')
ax4.grid(True, alpha=0.3, axis='x')
ax4.invert_yaxis()

# 5. Top genes by Spearman correlation
ax5 = fig.add_subplot(gs[1, 1])
top_spearman_20 = top_spearman.head(20)
bars = ax5.barh(range(len(top_spearman_20)), top_spearman_20['Spearman_Correlation'],
               color='green', alpha=0.7)
ax5.set_yticks(range(len(top_spearman_20)))
ax5.set_yticklabels(top_spearman_20['Gene'], fontsize=9)
ax5.set_xlabel('Spearman Correlation', fontweight='bold')
ax5.set_title('E. Top 20 by Spearman Correlation', fontweight='bold')
ax5.grid(True, alpha=0.3, axis='x')
ax5.invert_yaxis()

# 6. Cosine similarity distribution
ax6 = fig.add_subplot(gs[1, 2])
ax6.hist(similarity_df['Cosine_Similarity'].dropna(), bins=50, edgecolor='black', alpha=0.7, color='purple')
ax6.axvline(similarity_df['Cosine_Similarity'].mean(), color='red', linestyle='--', 
           linewidth=2, label=f'Mean: {similarity_df["Cosine_Similarity"].mean():.3f}')
ax6.set_xlabel('Cosine Similarity with MGMT', fontweight='bold')
ax6.set_ylabel('Frequency', fontweight='bold')
ax6.set_title('F. Cosine Similarity Distribution', fontweight='bold')
ax6.legend()
ax6.grid(True, alpha=0.3)

# 7. Correlation comparison (Pearson vs Spearman)
ax7 = fig.add_subplot(gs[2, 0])
ax7.scatter(similarity_df['Pearson_Correlation'], similarity_df['Spearman_Correlation'],
           alpha=0.5, s=30, edgecolors='black', linewidth=0.3)
ax7.set_xlabel('Pearson Correlation', fontweight='bold')
ax7.set_ylabel('Spearman Correlation', fontweight='bold')
ax7.set_title('G. Pearson vs Spearman Correlation', fontweight='bold')
ax7.plot([-1, 1], [-1, 1], 'r--', alpha=0.5, label='y=x')
ax7.grid(True, alpha=0.3)
ax7.legend()

# Highlight top 10
top_10_genes = top_50.head(10)
for _, row in top_10_genes.iterrows():
    if not pd.isna(row['Pearson_Correlation']) and not pd.isna(row['Spearman_Correlation']):
        ax7.scatter(row['Pearson_Correlation'], row['Spearman_Correlation'],
                   s=200, marker='*', color='yellow', edgecolors='black', linewidth=1, zorder=5)
        ax7.annotate(row['Gene'], 
                    (row['Pearson_Correlation'], row['Spearman_Correlation']),
                    fontsize=7, alpha=0.9, fontweight='bold', zorder=6)

# 8. Mean difference vs Similarity
ax8 = fig.add_subplot(gs[2, 1])
ax8.scatter(top_50['Mean_Difference'], top_50['Composite_Similarity_Score'],
           alpha=0.6, s=60, c=top_50['Pearson_Correlation'], cmap='coolwarm',
           edgecolors='black', linewidth=0.3)
ax8.set_xlabel('Mean Methylation Difference', fontweight='bold')
ax8.set_ylabel('Composite Similarity Score', fontweight='bold')
ax8.set_title('H. Similarity vs Mean Difference\n(Top 50 Genes)', fontweight='bold')
ax8.grid(True, alpha=0.3)
cbar8 = plt.colorbar(ax8.collections[0], ax=ax8)
cbar8.set_label('Pearson Correlation', fontweight='bold')

# 9. Summary statistics
ax9 = fig.add_subplot(gs[2, 2])
ax9.axis('off')

summary_text = f"""
MGMT Similarity Analysis Summary

Total Genes Analyzed: {len(similarity_df):,}
MGMT Gene: {mgmt_gene}

Top Similarity Metrics:
  - Highest Pearson: {top_pearson.iloc[0]['Gene']} (r={top_pearson.iloc[0]['Pearson_Correlation']:.4f})
  - Highest Spearman: {top_spearman.iloc[0]['Gene']} (r={top_spearman.iloc[0]['Spearman_Correlation']:.4f})
  - Highest Cosine: {top_cosine.iloc[0]['Gene']} (sim={top_cosine.iloc[0]['Cosine_Similarity']:.4f})

Significant Correlations:
  - Pearson (p<0.05): {len(sig_pearson)}
  - Spearman (p<0.05): {len(sig_spearman)}

Mean Correlations:
  - Pearson: {similarity_df['Pearson_Correlation'].mean():.4f}
  - Spearman: {similarity_df['Spearman_Correlation'].mean():.4f}
  - Cosine: {similarity_df['Cosine_Similarity'].mean():.4f}
"""

ax9.text(0.1, 0.5, summary_text, fontsize=10, fontweight='bold',
        verticalalignment='center', family='monospace',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.suptitle('MGMT Methylation Profile Similarity Analysis', 
            fontsize=16, fontweight='bold', y=0.995)
plt.savefig(plots_dir / "mgmt_similarity_analysis.png", dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Visualization saved: mgmt_similarity_analysis.png")

# ============================================================================
# 7. DETAILED TOP 50 LIST
# ============================================================================

print("\nStep 7: Generating detailed top 50 list...")
print("-" * 80)

top_50_detailed = top_50[['Gene', 'Pearson_Correlation', 'Pearson_P', 
                          'Spearman_Correlation', 'Spearman_P',
                          'Cosine_Similarity', 'Composite_Similarity_Score',
                          'Mean_Difference', 'MAE', 'RMSE',
                          'Gene_Methylation_Mean', 'MGMT_Methylation_Mean']].copy()

# Round values
for col in ['Pearson_Correlation', 'Spearman_Correlation', 'Cosine_Similarity', 
           'Composite_Similarity_Score', 'Mean_Difference', 'MAE', 'RMSE',
           'Gene_Methylation_Mean', 'MGMT_Methylation_Mean']:
    top_50_detailed[col] = top_50_detailed[col].round(4)

top_50_detailed['Pearson_P'] = top_50_detailed['Pearson_P'].apply(
    lambda x: f"{x:.2e}" if not pd.isna(x) else "N/A")
top_50_detailed['Spearman_P'] = top_50_detailed['Spearman_P'].apply(
    lambda x: f"{x:.2e}" if not pd.isna(x) else "N/A")

top_50_detailed.to_csv(output_dir / "mgmt_similar_genes_top50_detailed.csv", index=False)

print("\nTop 50 Genes Most Similar to MGMT (Detailed):")
print("="*80)
print(top_50_detailed.to_string(index=False))
print()

# ============================================================================
# SUMMARY
# ============================================================================

print("="*80)
print("Analysis Summary")
print("="*80)
print()
print(f"MGMT Gene: {mgmt_gene}")
print(f"Total genes analyzed: {len(similarity_df):,}")
print(f"Genes with significant Pearson correlation: {len(sig_pearson)}")
print(f"Genes with significant Spearman correlation: {len(sig_spearman)}")
print()
print("Top 10 Genes Most Similar to MGMT:")
print("-" * 80)
for i, (idx, row) in enumerate(top_50.head(10).iterrows(), 1):
    print(f"{i:2d}. {row['Gene']:15s} | Pearson: {row['Pearson_Correlation']:7.4f} | "
          f"Spearman: {row['Spearman_Correlation']:7.4f} | "
          f"Cosine: {row['Cosine_Similarity']:7.4f} | "
          f"Score: {row['Composite_Similarity_Score']:7.4f}")
print()
print("Files saved:")
print(f"- results/mgmt_similarity_analysis_all_genes.csv")
print(f"- results/mgmt_similar_genes_top50.csv")
print(f"- results/mgmt_similar_genes_top50_detailed.csv")
print(f"- results/mgmt_similar_genes_by_pearson.csv")
print(f"- results/mgmt_similar_genes_by_spearman.csv")
print(f"- results/mgmt_similar_genes_by_cosine.csv")
print(f"- plots/mgmt_similarity_analysis.png")
print()
print("="*80)
print("Analysis completed successfully!")
print("="*80)

