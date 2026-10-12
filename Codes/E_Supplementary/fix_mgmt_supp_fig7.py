#!/usr/bin/env python3
"""
Fix Supplementary Figure 7 (MGMT methylation similarity).

The stored CSV was corrupt: Pearson was globally NaN and Spearman/Cosine/
Mean_Difference were populated for only 4 genes (so panels A, B, C, D, E, G, H
were empty or showed only ~4 genes). This script recomputes ALL similarity
metrics vs MGMT for every gene (vectorized), rebuilds the composite score with
the original weights, rewrites the CSV, and regenerates the figure.
"""

import numpy as np
import pandas as pd
from pathlib import Path
import sys
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from repo_paths import REPO_ROOT, data_file
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "stix"

BASE = REPO_ROOT
METH = data_file("methylation_epic.txt")
CSV = BASE / "results" / "mgmt_similarity_analysis_all_genes.csv"
NROWS = 100000  # match original analysis window

print("Reading methylation matrix (first 100k rows) ...")
meth = pd.read_csv(METH, sep="\t", nrows=NROWS)
patients = [c for c in meth.columns if c.startswith("C3")]
n = len(patients)
print(f"  {meth.shape[0]} probes x {n} patients")

# ---- MGMT profile (median-filled, matching original) --------------------
mgmt_rows = meth[meth["Hugo_Symbol"].astype(str).str.upper() == "MGMT"]
mgmt_raw = mgmt_rows[patients].iloc[0].values.astype(float)
mgmt_med = np.nanmedian(mgmt_raw)
mgmt = np.nan_to_num(mgmt_raw, nan=mgmt_med, posinf=mgmt_med, neginf=mgmt_med)
mgmt_mean = mgmt.mean()

# ---- gene matrix, NaN filled with per-gene median -----------------------
X = meth[patients].to_numpy(dtype=float)
row_med = np.nanmedian(X, axis=1)
row_med = np.where(np.isnan(row_med), 0.5, row_med)
nan_idx = np.where(np.isnan(X))
X[nan_idx] = np.take(row_med, nan_idx[0])

gene_std = X.std(axis=1)
valid = gene_std >= 0.01  # skip near-constant genes (matches original filter)

def pearson_vs(target):
    tc = target - target.mean()
    t_ss = np.sqrt((tc ** 2).sum())
    Xc = X - X.mean(axis=1, keepdims=True)
    num = Xc @ tc
    den = np.sqrt((Xc ** 2).sum(axis=1)) * t_ss
    with np.errstate(divide="ignore", invalid="ignore"):
        rr = num / den
    rr[~np.isfinite(rr)] = np.nan
    return rr

print("Computing Pearson ...")
pear = pearson_vs(mgmt)

print("Computing Spearman (rank-based, vectorized) ...")
# average-rank not needed for continuous beta; ordinal ranks are fine
rank_X = np.argsort(np.argsort(X, axis=1), axis=1).astype(float)
rank_mgmt = np.argsort(np.argsort(mgmt)).astype(float)
rmc = rank_mgmt - rank_mgmt.mean()
rm_ss = np.sqrt((rmc ** 2).sum())
rXc = rank_X - rank_X.mean(axis=1, keepdims=True)
num = rXc @ rmc
den = np.sqrt((rXc ** 2).sum(axis=1)) * rm_ss
with np.errstate(divide="ignore", invalid="ignore"):
    spear = num / den
spear[~np.isfinite(spear)] = np.nan

print("Computing Cosine, Mean_Difference, MAE, RMSE ...")
cos = (X @ mgmt) / (np.linalg.norm(X, axis=1) * np.linalg.norm(mgmt))
mean_diff = np.abs(X.mean(axis=1) - mgmt_mean)
mae = np.abs(X - mgmt).mean(axis=1)
rmse = np.sqrt(((X - mgmt) ** 2).mean(axis=1))

# two-sided p-values from t-distribution
def pvals(rv):
    with np.errstate(divide="ignore", invalid="ignore"):
        t = rv * np.sqrt((n - 2) / (1 - rv ** 2))
    p = 2 * stats.t.sf(np.abs(t), n - 2)
    return np.where(np.isnan(rv), np.nan, p)

pear_p = pvals(pear)
spear_p = pvals(spear)

# Correlation_Similarity = 1 - correlation distance = Pearson (centered)
corr_sim = pear.copy()

df_all = pd.DataFrame({
    "Gene": meth["Hugo_Symbol"].astype(str).values,
    "Pearson_Correlation": pear,
    "Pearson_P": pear_p,
    "Spearman_Correlation": spear,
    "Spearman_P": spear_p,
    "Cosine_Similarity": cos,
    "Correlation_Similarity": corr_sim,
    "MAE": mae,
    "RMSE": rmse,
    "Mean_Difference": mean_diff,
    "Gene_Methylation_Mean": X.mean(axis=1),
    "Gene_Methylation_STD": gene_std,
    "MGMT_Methylation_Mean": mgmt_mean,
})
# drop near-constant genes and MGMT itself
df_all = df_all[valid]
df_all = df_all[~df_all["Gene"].str.upper().str.contains("MGMT", na=False)]
# first occurrence per gene
df_all = df_all.drop_duplicates(subset="Gene", keep="first").reset_index(drop=True)

# Composite score with original weights
df_all["Composite_Similarity_Score"] = (
    df_all["Pearson_Correlation"].fillna(0) * 0.3
    + df_all["Spearman_Correlation"].fillna(0) * 0.3
    + df_all["Cosine_Similarity"].fillna(0) * 0.2
    + df_all["Correlation_Similarity"].fillna(0) * 0.2
)

sim = df_all
print(f"  Recomputed metrics for {len(sim)} genes")
for c in ["Pearson_Correlation", "Spearman_Correlation", "Cosine_Similarity", "Mean_Difference"]:
    print(f"    {c}: non-null={sim[c].notna().sum()}")

sim.to_csv(CSV, index=False)
print(f"  Updated CSV saved: {CSV}")

pear_vals = sim["Pearson_Correlation"].dropna()

# ---------------------------------------------------------------- figure
plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams["savefig.dpi"] = 300

top_20 = sim.nlargest(20, "Composite_Similarity_Score")
top_pearson_20 = sim.nlargest(20, "Pearson_Correlation")
top_spearman_20 = sim.nlargest(20, "Spearman_Correlation")
top_50 = sim.nlargest(50, "Composite_Similarity_Score")

LABEL_FS = 7  # smaller gene-label font

def short_gene(g):
    """Keep at most 2 gene symbols (semicolon/comma separated)."""
    parts = [p for p in str(g).replace(",", ";").split(";") if p]
    return ";".join(parts[:2])

fig = plt.figure(figsize=(20, 16))
gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.3)

# A. Composite
ax1 = fig.add_subplot(gs[0, 0])
ax1.barh(range(len(top_20)), top_20["Composite_Similarity_Score"], color="#4B0082", alpha=0.85)
ax1.set_yticks(range(len(top_20)))
ax1.set_yticklabels([short_gene(g) for g in top_20["Gene"]], fontsize=LABEL_FS)
ax1.set_xlabel("Composite Similarity Score", fontweight="bold")
ax1.set_title("a. Top 20 Genes Similar to MGMT\n(Composite Score)", fontweight="bold")
ax1.grid(True, alpha=0.3, axis="x")
ax1.invert_yaxis()

# B. Pearson distribution
ax2 = fig.add_subplot(gs[0, 1])
pear_vals = sim["Pearson_Correlation"].dropna()
ax2.hist(pear_vals, bins=50, edgecolor="black", alpha=0.7, color="#1f77b4")
ax2.axvline(pear_vals.mean(), color="red", linestyle="--", linewidth=2,
            label=f"Mean: {pear_vals.mean():.3f}")
ax2.set_xlabel("Pearson Correlation with MGMT", fontweight="bold")
ax2.set_ylabel("Frequency", fontweight="bold")
ax2.set_title("b. Pearson Correlation Distribution", fontweight="bold")
ax2.legend()
ax2.grid(True, alpha=0.3)

# C. Spearman distribution
ax3 = fig.add_subplot(gs[0, 2])
sp_vals = sim["Spearman_Correlation"].dropna()
ax3.hist(sp_vals, bins=50, edgecolor="black", alpha=0.7, color="green")
ax3.axvline(sp_vals.mean(), color="red", linestyle="--", linewidth=2,
            label=f"Mean: {sp_vals.mean():.3f}")
ax3.set_xlabel("Spearman Correlation with MGMT", fontweight="bold")
ax3.set_ylabel("Frequency", fontweight="bold")
ax3.set_title("c. Spearman Correlation Distribution", fontweight="bold")
ax3.legend()
ax3.grid(True, alpha=0.3)

# D. Top 20 by Pearson
ax4 = fig.add_subplot(gs[1, 0])
ax4.barh(range(len(top_pearson_20)), top_pearson_20["Pearson_Correlation"], color="#1f77b4", alpha=0.85)
ax4.set_yticks(range(len(top_pearson_20)))
ax4.set_yticklabels([short_gene(g) for g in top_pearson_20["Gene"]], fontsize=LABEL_FS)
ax4.set_xlabel("Pearson Correlation", fontweight="bold")
ax4.set_title("d. Top 20 by Pearson Correlation", fontweight="bold")
ax4.grid(True, alpha=0.3, axis="x")
ax4.invert_yaxis()

# E. Top 20 by Spearman
ax5 = fig.add_subplot(gs[1, 1])
ax5.barh(range(len(top_spearman_20)), top_spearman_20["Spearman_Correlation"], color="green", alpha=0.85)
ax5.set_yticks(range(len(top_spearman_20)))
ax5.set_yticklabels([short_gene(g) for g in top_spearman_20["Gene"]], fontsize=LABEL_FS)
ax5.set_xlabel("Spearman Correlation", fontweight="bold")
ax5.set_title("e. Top 20 by Spearman Correlation", fontweight="bold")
ax5.grid(True, alpha=0.3, axis="x")
ax5.invert_yaxis()

# F. Cosine distribution
ax6 = fig.add_subplot(gs[1, 2])
cos_vals = sim["Cosine_Similarity"].dropna()
ax6.hist(cos_vals, bins=50, edgecolor="black", alpha=0.7, color="purple")
ax6.axvline(cos_vals.mean(), color="red", linestyle="--", linewidth=2,
            label=f"Mean: {cos_vals.mean():.3f}")
ax6.set_xlabel("Cosine Similarity with MGMT", fontweight="bold")
ax6.set_ylabel("Frequency", fontweight="bold")
ax6.set_title("f. Cosine Similarity Distribution", fontweight="bold")
ax6.legend()
ax6.grid(True, alpha=0.3)

# G. Pearson vs Spearman
ax7 = fig.add_subplot(gs[2, 0])
both = sim[["Pearson_Correlation", "Spearman_Correlation"]].dropna()
ax7.scatter(both["Pearson_Correlation"], both["Spearman_Correlation"],
            alpha=0.3, s=12, c="#4B0082")
ax7.plot([-1, 1], [-1, 1], "r--", alpha=0.6, label="y = x")
for _, row in top_20.iterrows():
    if pd.notna(row["Pearson_Correlation"]) and pd.notna(row["Spearman_Correlation"]):
        ax7.scatter(row["Pearson_Correlation"], row["Spearman_Correlation"],
                    s=40, c="#DD8452", edgecolor="black", zorder=5)
lim_lo = min(both.min().min(), -0.1) - 0.05
lim_hi = max(both.max().max(), 0.3) + 0.05
ax7.set_xlim(lim_lo, lim_hi); ax7.set_ylim(lim_lo, lim_hi)
ax7.set_xlabel("Pearson Correlation", fontweight="bold")
ax7.set_ylabel("Spearman Correlation", fontweight="bold")
ax7.set_title("g. Pearson vs Spearman Correlation", fontweight="bold")
ax7.grid(True, alpha=0.3)
ax7.legend()

# H. Mean diff vs composite (top 50), colored by Pearson
ax8 = fig.add_subplot(gs[2, 1])
sc = ax8.scatter(top_50["Mean_Difference"], top_50["Composite_Similarity_Score"],
                 alpha=0.75, s=60, c=top_50["Pearson_Correlation"], cmap="coolwarm",
                 edgecolor="black", linewidth=0.4)
ax8.set_xlabel("Mean Methylation Difference", fontweight="bold")
ax8.set_ylabel("Composite Similarity Score", fontweight="bold")
ax8.set_title("h. Similarity vs Mean Difference\n(Top 50 Genes)", fontweight="bold")
ax8.grid(True, alpha=0.3)
cbar8 = plt.colorbar(sc, ax=ax8)
cbar8.set_label("Pearson Correlation", fontweight="bold")

# I. Summary
ax9 = fig.add_subplot(gs[2, 2])
ax9.axis("off")
sig_pearson = int(((sim["Pearson_P"] < 0.05) & sim["Pearson_Correlation"].notna()).sum())
sig_spearman = int(((sim["Spearman_P"] < 0.05) & sim["Spearman_Correlation"].notna()).sum())
tp = top_pearson_20.iloc[0]
ts = top_spearman_20.iloc[0]
tc = sim.nlargest(1, "Cosine_Similarity").iloc[0]
summary_text = (
    "MGMT Similarity Analysis Summary\n\n"
    f"Total Genes Analyzed: {len(sim):,}\n"
    "MGMT Gene: MGMT\n\n"
    "Top Similarity Metrics:\n"
    f"  - Highest Pearson: {short_gene(tp['Gene'])} (r={tp['Pearson_Correlation']:.4f})\n"
    f"  - Highest Spearman: {short_gene(ts['Gene'])} (r={ts['Spearman_Correlation']:.4f})\n"
    f"  - Highest Cosine: {short_gene(tc['Gene'])} (sim={tc['Cosine_Similarity']:.4f})\n\n"
    "Significant Correlations:\n"
    f"  - Pearson (p<0.05): {sig_pearson}\n"
    f"  - Spearman (p<0.05): {sig_spearman}\n\n"
    "Mean Correlations:\n"
    f"  - Pearson: {pear_vals.mean():.4f}\n"
    f"  - Spearman: {sp_vals.mean():.4f}\n"
    f"  - Cosine: {cos_vals.mean():.4f}"
)
ax9.text(0.05, 0.5, summary_text, fontsize=11, fontweight="bold", va="center",
         bbox=dict(boxstyle="round,pad=0.6", facecolor="#FBE9C8", edgecolor="#C8A04B"))

fig.suptitle("MGMT Methylation Profile Similarity Analysis", fontsize=16, fontweight="bold", y=0.995)

plots_dir = BASE / "plots"
plots_dir.mkdir(exist_ok=True)
out1 = plots_dir / "mgmt_similarity_analysis.png"
out2 = BASE / "figures" / "Supplementary" / "Supplementary_Figure_7.png"
plt.savefig(out1, dpi=300, bbox_inches="tight", facecolor="white")
plt.savefig(out2, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print("Saved:", out1)
print("Saved:", out2)
