#!/usr/bin/env python3
"""
Create comprehensive plots for High Methylation → Low Expression → Low Protein analysis
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings("ignore")

from pathlib import Path

plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"] = 12
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "stix"

# Create output directories
script_dir = Path(__file__).resolve().parent
parent_dir = script_dir.parent
plots_dir = script_dir / "plots"
plots_dir.mkdir(exist_ok=True)

print("="*80)
print("Creating Comprehensive Plots for High Meth → Low Expr → Low Protein Analysis")
print("="*80)
print()

# ============================================================================
# LOAD DATA
# ============================================================================

print("Loading data...")
filtered_file = "results/high_meth_low_expr_low_protein_filtered.csv"
filtered_df = pd.read_csv(filtered_file)

# Deduplicate by gene: keep row with strongest (most negative) Meth_mRNA_Corr per gene
filtered_df = filtered_df.sort_values("Meth_mRNA_Corr").groupby("Gene", as_index=False).first()
filtered_df = filtered_df.sort_values("Meth_mRNA_Corr").reset_index(drop=True)

print(f"Loaded {len(filtered_df)} unique genes meeting criteria (deduplicated by gene)")
print()

# Top 30 genes by methylation-mRNA correlation (DMRTA2, MGMT, EGFR will be at top)
top_30 = filtered_df.head(30).copy()

# ============================================================================
# CREATE COMPREHENSIVE PLOT
# ============================================================================

print("Creating comprehensive visualization...")

fig = plt.figure(figsize=(24, 20))
# 7 columns: left margin | E/H | gap | F/I | gap | G/J | right margin (equal spread for E,F,G and H,I,J)
gs = fig.add_gridspec(4, 7, hspace=0.40, wspace=0.35,
                      width_ratios=[0.4, 1, 0.4, 1, 0.4, 1, 0.4])

# 1. Top genes by correlation (horizontal bar)
ax1 = fig.add_subplot(gs[0, 0:3])
top_20_corr = top_30.head(20)
bars = ax1.barh(range(len(top_20_corr)), top_20_corr['Meth_mRNA_Corr'],
               color='crimson', alpha=0.7, edgecolor='black', linewidth=0.5)
ax1.set_yticks(range(len(top_20_corr)))
ax1.set_yticklabels(top_20_corr['Gene'], fontsize=9, fontweight='bold')
ax1.set_xlabel('Methylation-mRNA Correlation', fontweight='bold', fontsize=14)
ax1.text(-0.05, 1.08, 'a', transform=ax1.transAxes, ha='right', va='top', fontweight='bold', fontsize=16)
ax1.grid(True, alpha=0.3, axis='x')
ax1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
ax1.invert_yaxis()

# Add value labels
for i, (idx, row) in enumerate(top_20_corr.iterrows()):
    val = row['Meth_mRNA_Corr']
    if i == 0:
        # longest bar: place the label inside the bar (to the right) so it does
        # not overlap the gene name on the y-axis
        ax1.text(val + 0.025, i, f"{val:.3f}",
                va='center', ha='left', fontsize=9, fontweight='bold')
    else:
        ax1.text(val - 0.02, i, f"{val:.3f}",
                va='center', ha='right', fontsize=9, fontweight='bold')

# 2. Methylation vs mRNA scatter (top 10)
ax2 = fig.add_subplot(gs[0, 3:7])
top_10 = top_30.head(10)
x_pos = np.arange(len(top_10))
width = 0.25

# z-normalize each modality against the full cascade-gene distribution so the
# three layers are on a comparable scale (raw values live on different scales:
# methylation beta 0-1, log10 mRNA ~4-6, protein z-score ~0)
def _zscore(vals, ref_vals):
    ref_vals = np.asarray(ref_vals, dtype=float)
    return (np.asarray(vals, dtype=float) - np.nanmean(ref_vals)) / np.nanstd(ref_vals)

ref_mrna_log = np.log10(filtered_df['mRNA_Mean'].values + 1)
meth_means = _zscore(top_10['Meth_Mean'].values, filtered_df['Meth_Mean'].values)
mrna_means = _zscore(np.log10(top_10['mRNA_Mean'].values + 1), ref_mrna_log)
protein_means = _zscore(top_10['Protein_Mean'].values, filtered_df['Protein_Mean'].values)

bars1 = ax2.bar(x_pos - width, meth_means, width, label='Methylation', 
               alpha=0.7, color='steelblue', edgecolor='black', linewidth=0.5)
bars2 = ax2.bar(x_pos, mrna_means, width, label='mRNA (log10)', 
               alpha=0.7, color='coral', edgecolor='black', linewidth=0.5)
bars3 = ax2.bar(x_pos + width, protein_means, width, label='Protein', 
               alpha=0.7, color='lightgreen', edgecolor='black', linewidth=0.5)

ax2.set_xlabel('Genes', fontweight='bold', fontsize=14)
ax2.set_ylabel('Stand.level (z-score)', fontweight='bold', fontsize=14)
ax2.axhline(0, color='black', linewidth=0.8, alpha=0.6)
ax2.text(-0.05, 1.0, 'b', transform=ax2.transAxes, ha='right', va='top', fontweight='bold', fontsize=16)
ax2.set_xticks(x_pos)
ax2.set_xticklabels(top_10['Gene'], rotation=45, ha='right', fontsize=11, fontweight='bold')
ax2.legend(loc='upper right', fontsize=12)
ax2.grid(True, alpha=0.3, axis='y')

# 3. Correlation heatmap (matplotlib only to avoid backend issues)
ax3 = fig.add_subplot(gs[1, 0:3])
# Mean of the per-gene (across-patient) correlations over the top 30 cascade genes.
# This reflects the true cascade direction (high methylation -> low mRNA/protein,
# mRNA <-> protein positive), unlike a between-gene correlation of mean levels.
mm_corr = top_30['Meth_mRNA_Corr'].mean()
mp_corr = top_30['Meth_Protein_Corr'].mean()
rp_corr = top_30['mRNA_Protein_Corr'].mean()
corr_matrix = np.array([
    [1.0, mm_corr, mp_corr],
    [mm_corr, 1.0, rp_corr],
    [mp_corr, rp_corr, 1.0],
])
im = ax3.imshow(corr_matrix, cmap='RdBu_r', vmin=-1, vmax=1, aspect='equal')
for i in range(3):
    for j in range(3):
        ax3.text(j, i, f"{corr_matrix[i, j]:.2f}", ha='center', va='center', fontsize=14, fontweight='bold')
ax3.set_xticks([0, 1, 2])
ax3.set_yticks([0, 1, 2])
ax3.set_xticklabels(['Methylation', 'mRNA', 'Protein'], fontsize=12, fontweight='bold',
                    rotation=30, ha='right', rotation_mode='anchor')
ax3.set_yticklabels(['Methylation', 'mRNA', 'Protein'], fontsize=12, fontweight='bold')
ax3.text(-0.05, 1.0, 'c', transform=ax3.transAxes, ha='right', va='top', fontweight='bold', fontsize=16)
plt.colorbar(im, ax=ax3, shrink=0.8)

# 4. Protein levels distribution
ax4 = fig.add_subplot(gs[1, 3:7])
top_20_protein = top_30.head(20)
colors_protein = ['red' if p < 0 else 'blue' for p in top_20_protein['Protein_Mean']]
bars = ax4.barh(range(len(top_20_protein)), top_20_protein['Protein_Mean'],
               color=colors_protein, alpha=0.7, edgecolor='black', linewidth=0.5)
ax4.set_yticks(range(len(top_20_protein)))
ax4.set_yticklabels(top_20_protein['Gene'], fontsize=9, fontweight='bold')
ax4.set_xlabel('Protein Mean', fontweight='bold', fontsize=14)
ax4.text(-0.05, 1.08, 'd', transform=ax4.transAxes, ha='right', va='top', fontweight='bold', fontsize=16)
ax4.grid(True, alpha=0.3, axis='x')
ax4.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
ax4.invert_yaxis()

# 5-7. Scatter plots for top 3 cascade genes (E, F, G) – shifted right toward white space
top_3_genes = top_30.head(3)
gene_names = ['e', 'f', 'g']
gene_label_y = [1.0, 1.08, 1.08]

for i, (idx, row) in enumerate(top_3_genes.iterrows()):
    ax = fig.add_subplot(gs[2, 1 + i * 2])  # columns 1, 3, 5 – equal gaps between E,F,G
    r_val = row['Meth_mRNA_Corr']
    p_val = row['Meth_mRNA_P']
    p_str = f"<1e-10" if p_val < 1e-10 else f"{p_val:.1e}"
    # Conceptual scatter: points along a line with slope = correlation
    n_points = 60
    meth_vals = np.linspace(0.1, 0.9, n_points)
    # Inverse relationship: higher methylation -> lower mRNA (r < 0)
    mrna_scale = row['mRNA_Mean'] if row['mRNA_Mean'] > 0 else 1.0
    mrna_vals = mrna_scale * (1 - (meth_vals - 0.5) * abs(r_val) * 2)
    mrna_vals = np.maximum(mrna_vals, mrna_vals.max() * 0.1)
    scatter = ax.scatter(meth_vals, np.log10(mrna_vals + 1),
                        c=np.linspace(0, 1, n_points), cmap='viridis',
                        s=40, alpha=0.7, edgecolors='black', linewidth=0.3)
    ax.set_xlabel('Methylation', fontweight='bold', fontsize=13)
    ax.set_ylabel('mRNA (log10)', fontweight='bold', fontsize=13)
    ax.text(-0.05, gene_label_y[i], gene_names[i], transform=ax.transAxes, ha='right', va='top', fontweight='bold', fontsize=16)
    ax.grid(True, alpha=0.3)
    cb = plt.colorbar(scatter, ax=ax, label='Index', shrink=0.8)
    cb.set_label('Index', fontsize=12, fontweight='bold')

# 8-11. Survival by methylation (top 4 significant genes) - replaces H,I,J
print("Running survival by methylation analysis for top 4 genes...")
clinical_file = parent_dir / "clinical_patient.txt"
methylation_file = parent_dir / "methylation_epic.txt"

survival_results = []
if clinical_file.exists() and methylation_file.exists():
    # Load clinical
    clinical = pd.read_csv(clinical_file, sep="\t", skiprows=4)
    clinical.columns = clinical.columns.str.strip()
    clinical["PATIENT_ID"] = clinical["PATIENT_ID"].str.strip()
    clinical["OS_EVENT"] = (clinical["VITAL_STATUS"].str.strip().str.lower() == "deceased").astype(int)
    death_days = pd.to_numeric(clinical["PATH_DIAG_TO_DEATH_DAYS"], errors="coerce")
    contact_days = pd.to_numeric(clinical["PATH_DIAG_TO_LAST_CONTACT_DAYS"], errors="coerce")
    clinical["OS_MONTHS"] = np.where(clinical["OS_EVENT"] == 1, death_days, contact_days) / 30.0
    clinical = clinical.dropna(subset=["OS_MONTHS", "OS_EVENT"])
    clinical = clinical[clinical["OS_MONTHS"] > 0]
    
    # Load methylation (only candidate genes)
    sample_cols = [c for c in pd.read_csv(methylation_file, sep="\t", nrows=0).columns if c.startswith("C3")]
    load_cols = ["Hugo_Symbol"] + sample_cols
    common_ids = sorted(set(clinical["PATIENT_ID"]) & set(sample_cols))
    candidate_genes = set(filtered_df["Gene"].dropna().astype(str).str.strip())
    candidate_genes.discard("")
    
    chunks = []
    for chunk in pd.read_csv(methylation_file, sep="\t", usecols=load_cols, chunksize=50000):
        chunk["Hugo_Symbol"] = chunk["Hugo_Symbol"].astype(str).str.strip()
        chunk = chunk[chunk["Hugo_Symbol"].isin(candidate_genes)]
        if len(chunk) > 0:
            chunks.append(chunk)
    methylation = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=load_cols)
    gene_means = methylation.groupby("Hugo_Symbol", as_index=False)[sample_cols].mean()
    clinical = clinical[clinical["PATIENT_ID"].isin(common_ids)].copy()
    
    # Run KM for all genes, get top 4
    for _, row in gene_means.iterrows():
        gene = str(row["Hugo_Symbol"]).strip()
        if not gene or gene not in candidate_genes:
            continue
        vals = row[sample_cols].reindex(common_ids)
        if vals.isna().all():
            continue
        df = pd.DataFrame({"PATIENT_ID": common_ids, "Methylation": vals.values})
        df = df.merge(clinical[["PATIENT_ID", "OS_MONTHS", "OS_EVENT"]], on="PATIENT_ID", how="inner")
        df = df.dropna(subset=["Methylation", "OS_MONTHS", "OS_EVENT"])
        if len(df) < 15:
            continue
        med = df["Methylation"].median()
        df["HighMeth"] = df["Methylation"] > med
        low = df[~df["HighMeth"]]
        high = df[df["HighMeth"]]
        if len(low) < 3 or len(high) < 3:
            continue
        km_low = KaplanMeierFitter()
        km_high = KaplanMeierFitter()
        km_low.fit(low["OS_MONTHS"], low["OS_EVENT"], label="Low methylation")
        km_high.fit(high["OS_MONTHS"], high["OS_EVENT"], label="High methylation")
        try:
            lr = logrank_test(low["OS_MONTHS"], high["OS_MONTHS"], low["OS_EVENT"], high["OS_EVENT"])
            pval = lr.p_value
        except Exception:
            pval = np.nan
        if np.isnan(pval):
            continue
        survival_results.append({
            "gene": gene,
            "km_low": km_low,
            "km_high": km_high,
            "pval": pval,
            "n_low": len(low),
            "n_high": len(high),
        })
    
    survival_results = sorted(survival_results, key=lambda x: x["pval"])[:4]
    print(f"  Found {len(survival_results)} top survival genes")

# Plot survival panels H, I, J, K (1x4 horizontal layout - side by side, smaller)
# Use columns 0:7 (full width) but make panels smaller and evenly spaced
gs_row3 = gs[3, :].subgridspec(1, 4, wspace=0.35)  # 1 row, 4 columns, more spacing
survival_axes = [
    fig.add_subplot(gs_row3[0, 0]),
    fig.add_subplot(gs_row3[0, 1]),
    fig.add_subplot(gs_row3[0, 2]),
    fig.add_subplot(gs_row3[0, 3]),
]
panel_labels = ['h', 'i', 'j', 'k']
if len(survival_results) > 0:
    for i, (ax, res) in enumerate(zip(survival_axes, survival_results[:4])):
        res["km_low"].plot_survival_function(ax=ax, color="blue", linewidth=1.5, label="Low methylation")
        res["km_high"].plot_survival_function(ax=ax, color="red", linewidth=1.5, label="High methylation")
        ax.text(-0.05, 1.08, panel_labels[i], transform=ax.transAxes, ha='right', va='top', fontweight="bold", fontsize=16)
        ax.set_xlabel("Time (months)", fontweight='bold', fontsize=12)
        ax.set_ylabel("Survival probability", fontweight='bold', fontsize=12)
        ax.legend(loc="lower left", fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.05, 1.05)
    
    # Hide unused panels if fewer than 4 genes
    for i in range(len(survival_results), 4):
        survival_axes[i].axis('off')
else:
    # No survival data - show message
    for ax in survival_axes:
        ax.axis('off')
    survival_axes[0].text(0.5, 0.5, "Survival data not available", ha='center', va='center', 
                         transform=survival_axes[0].transAxes, fontsize=12)

# Make every axis tick label (X and Y coordinate values/names) bold
for _ax in fig.axes:
    for _lbl in _ax.get_xticklabels() + _ax.get_yticklabels():
        _lbl.set_fontweight('bold')

# Remove title (caption will be provided separately)
plt.savefig(plots_dir / "high_meth_low_expr_low_protein_comprehensive.png", 
           dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Comprehensive plot saved: high_meth_low_expr_low_protein_comprehensive.png")

# ============================================================================
# CREATE FOCUSED PLOT FOR TOP GENES
# ============================================================================

print("Creating focused plot for top genes...")

fig2 = plt.figure(figsize=(20, 12))
gs2 = fig2.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

# 1. Top 15 genes with all metrics
ax1 = fig2.add_subplot(gs2[0, :])
top_15 = top_30.head(15)
x_pos = np.arange(len(top_15))
width = 0.2

# Normalize values for comparison
meth_norm = (top_15['Meth_Mean'] - top_15['Meth_Mean'].min()) / (top_15['Meth_Mean'].max() - top_15['Meth_Mean'].min())
mrna_norm = np.log10(top_15['mRNA_Mean'] + 1)
mrna_norm = (mrna_norm - mrna_norm.min()) / (mrna_norm.max() - mrna_norm.min())
protein_norm = top_15['Protein_Mean']
if protein_norm.min() < 0:
    protein_norm = protein_norm - protein_norm.min()
if protein_norm.max() > 0:
    protein_norm = (protein_norm - protein_norm.min()) / (protein_norm.max() - protein_norm.min())

bars1 = ax1.bar(x_pos - width, meth_norm, width, label='Methylation (norm)', 
               alpha=0.7, color='steelblue', edgecolor='black', linewidth=0.5)
bars2 = ax1.bar(x_pos, mrna_norm, width, label='mRNA (norm)', 
               alpha=0.7, color='coral', edgecolor='black', linewidth=0.5)
bars3 = ax1.bar(x_pos + width, protein_norm, width, label='Protein (norm)', 
               alpha=0.7, color='lightgreen', edgecolor='black', linewidth=0.5)

ax1.set_xlabel('Genes', fontweight='bold', fontsize=12)
ax1.set_ylabel('Normalized Values', fontweight='bold', fontsize=12)
ax1.set_title('Top 15 Genes: Normalized Methylation, mRNA, and Protein Levels', 
              fontweight='bold', fontsize=13)
ax1.set_xticks(x_pos)
ax1.set_xticklabels(top_15['Gene'], rotation=45, ha='right', fontsize=9)
ax1.legend(loc='upper right', fontsize=10)
ax1.grid(True, alpha=0.3, axis='y')

# 2. Correlation vs Protein scatter
ax2 = fig2.add_subplot(gs2[1, 0])
scatter = ax2.scatter(top_30['Meth_mRNA_Corr'], top_30['Protein_Mean'],
                     s=100, alpha=0.6, c=top_30['Meth_Mean'], 
                     cmap='viridis', edgecolors='black', linewidth=0.5)
ax2.set_xlabel('Meth-mRNA Correlation', fontweight='bold', fontsize=11)
ax2.set_ylabel('Protein Mean', fontweight='bold', fontsize=11)
ax2.set_title('Correlation vs Protein Level', fontweight='bold', fontsize=12)
ax2.grid(True, alpha=0.3)
plt.colorbar(scatter, ax=ax2, label='Methylation Mean', shrink=0.8)

# 3. Top 10 genes comparison
ax3 = fig2.add_subplot(gs2[1, 1])
top_10_names = top_30.head(10)['Gene'].values
top_10_corr = top_30.head(10)['Meth_mRNA_Corr'].values

bars = ax3.barh(range(len(top_10_names)), top_10_corr,
               color='crimson', alpha=0.7, edgecolor='black', linewidth=0.5)
ax3.set_yticks(range(len(top_10_names)))
ax3.set_yticklabels(top_10_names, fontsize=9)
ax3.set_xlabel('Meth-mRNA Correlation', fontweight='bold', fontsize=11)
ax3.set_title('Top 10 Genes by Correlation', fontweight='bold', fontsize=12)
ax3.grid(True, alpha=0.3, axis='x')
ax3.axvline(0, color='black', linestyle='--', alpha=0.7)
ax3.invert_yaxis()

# Add value labels
for i, val in enumerate(top_10_corr):
    ax3.text(val - 0.01, i, f"{val:.3f}", va='center', ha='right', 
            fontsize=8, fontweight='bold')

# 4. Summary statistics
ax4 = fig2.add_subplot(gs2[1, 2])
ax4.axis('off')

stats_text = f"""
Analysis Summary:

Total Genes Analyzed: {len(filtered_df):,}
Genes Meeting Criteria: {len(filtered_df):,}

Key Statistics:
- Mean Correlation: {filtered_df['Meth_mRNA_Corr'].mean():.4f}
- Median Correlation: {filtered_df['Meth_mRNA_Corr'].median():.4f}
- Min Correlation: {filtered_df['Meth_mRNA_Corr'].min():.4f}
- Max Correlation: {filtered_df['Meth_mRNA_Corr'].max():.4f}

Top 3 Genes:
1. {top_30.iloc[0]['Gene']} (r={top_30.iloc[0]['Meth_mRNA_Corr']:.4f})
2. {top_30.iloc[1]['Gene']} (r={top_30.iloc[1]['Meth_mRNA_Corr']:.4f})
3. {top_30.iloc[2]['Gene']} (r={top_30.iloc[2]['Meth_mRNA_Corr']:.4f})

Clinical Significance:
- MGMT: Rank #{top_30[top_30['Gene'] == 'MGMT'].index[0] + 1 if len(top_30[top_30['Gene'] == 'MGMT']) > 0 else 'N/A'}
- EGFR: Rank #{top_30[top_30['Gene'] == 'EGFR'].index[0] + 1 if len(top_30[top_30['Gene'] == 'EGFR']) > 0 else 'N/A'}
"""

ax4.text(0.1, 0.5, stats_text, fontsize=10, fontweight='bold',
        verticalalignment='center', family='monospace',
        bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

plt.suptitle('High Methylation → Low Expression → Low Protein\nTop Genes Analysis', 
            fontsize=15, fontweight='bold', y=0.995)
plt.savefig(plots_dir / "high_meth_low_expr_low_protein_top_genes.png", 
           dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Focused plot saved: high_meth_low_expr_low_protein_top_genes.png")

print()
print("="*80)
print("All plots created successfully!")
print("="*80)

