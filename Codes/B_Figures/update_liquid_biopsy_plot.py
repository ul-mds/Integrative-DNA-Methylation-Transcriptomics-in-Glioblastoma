#!/usr/bin/env python3
"""
Update liquid biopsy plot for Figure 4
- Remove panel 4 (Distribution)
- Keep only 3 panels (A, B, C)
- Remove title
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "stix"
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path

script_dir = Path(__file__).resolve().parent
results_dir = script_dir / "results" / "novel_analyses"
plots_dir = script_dir / "plots" / "novel_analyses"
plots_dir.mkdir(parents=True, exist_ok=True)

# Load data
top_csv = results_dir / "liquid_biopsy_top_candidates.csv"
all_csv = results_dir / "liquid_biopsy_all.csv"

top_biopsy = pd.read_csv(top_csv) if top_csv.exists() else pd.DataFrame()
liquid_biopsy_df = pd.read_csv(all_csv) if all_csv.exists() else pd.DataFrame()

print(f"Loaded {len(top_biopsy)} top candidates")
print(f"Loaded {len(liquid_biopsy_df)} total genes")

if len(top_biopsy) == 0:
    print("Error: No data found!")
    exit(1)

# Create plot (1x3 grid, panel 4 removed)
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Clean specific gene labels for display (drop probe-suffix numerals)
def _clean_labels(genes):
    cleaned = []
    for g in genes:
        g = str(g).replace("CTD-2609K8.3", "CTD-2609K8").replace("AC004893.11", "AC004893")
        cleaned.append(g)
    return cleaned

# Panel A: Top candidates (Liquid Biopsy Score)
top_20 = top_biopsy.head(20)
axes[0].barh(range(len(top_20)), top_20['Liquid_Biopsy_Score'],
               color='gold', alpha=0.7, edgecolor='black', linewidth=0.5)
axes[0].set_yticks(range(len(top_20)))
axes[0].set_yticklabels(_clean_labels(top_20['Gene']), fontsize=10, fontweight='bold')
axes[0].set_xlabel('Liquid Biopsy Score', fontweight='bold', fontsize=11)
axes[0].text(-0.05, 1.0, 'a', transform=axes[0].transAxes, ha='right', va='top', fontweight='bold', fontsize=16)
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
axes[1].set_yticklabels(_clean_labels(top_20_detect['Gene']), fontsize=10, fontweight='bold')
axes[1].set_xlabel('Detectability Score', fontweight='bold', fontsize=11)
axes[1].text(-0.05, 1.0, 'b', transform=axes[1].transAxes, ha='right', va='top', fontweight='bold', fontsize=16)
axes[1].grid(True, alpha=0.3, axis='x')
axes[1].invert_yaxis()

# Panel C: Detectability vs Prognostic scatter (with top candidates highlighted)
if len(liquid_biopsy_df) > 0:
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
        gene_name = top_20_scatter.iloc[i]['Gene'].split(';')[0]
        axes[2].annotate(gene_name,
                        (top_20_scatter.iloc[i]['Detectability_Score'],
                         top_20_scatter.iloc[i]['Prognostic_Score']),
                        fontsize=11, fontweight='bold', alpha=0.8)
    axes[2].set_xlabel('Detectability Score', fontweight='bold', fontsize=11)
    axes[2].set_ylabel('Prognostic Score', fontweight='bold', fontsize=11)
    axes[2].text(-0.05, 1.0, 'c', transform=axes[2].transAxes, ha='right', va='top', fontweight='bold', fontsize=16)
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(fontsize=8)
    cbar = plt.colorbar(scatter, ax=axes[2])
    cbar.set_label('Mean Methylation', fontweight='bold')

# No suptitle (removed for Figure 4)
plt.tight_layout()

# Widen panel C toward the left to close the gap between panels b and c
# (keep the right edge fixed so the colorbar position is unaffected)
pos_b = axes[1].get_position()
pos_c = axes[2].get_position()
new_left = pos_b.x1 + 0.035
if new_left < pos_c.x0:
    axes[2].set_position([new_left, pos_c.y0, pos_c.x1 - new_left, pos_c.height])

plt.savefig(plots_dir / "liquid_biopsy_signatures.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print(f"✓ Updated plot: {plots_dir / 'liquid_biopsy_signatures.png'}")

# Print summary
print("\n" + "="*80)
print("Summary:")
print("="*80)
print(f"Top candidate: {top_biopsy.iloc[0]['Gene']}")
print(f"  - Liquid Biopsy Score: {top_biopsy.iloc[0]['Liquid_Biopsy_Score']:.4f}")
print(f"  - Detectability Score: {top_biopsy.iloc[0]['Detectability_Score']:.4f}")
print(f"  - Mean Methylation: {top_biopsy.iloc[0]['Mean_Methylation']:.4f}")
print(f"  - CV: {top_biopsy.iloc[0]['Coefficient_of_Variation']:.4f}")
print(f"  - Survival Correlation: {top_biopsy.iloc[0]['Survival_Correlation']:.4f}")
print(f"  - Survival P-value: {top_biopsy.iloc[0]['Survival_P_value']:.4f}")
print(f"\nTop 50 threshold: {top_biopsy['Liquid_Biopsy_Score'].min():.4f}")
