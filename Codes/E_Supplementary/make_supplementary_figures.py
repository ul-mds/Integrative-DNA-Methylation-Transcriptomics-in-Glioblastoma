#!/usr/bin/env python3
"""
Generate distinct supplementary figures into figures/Supplementary/.

- Supplementary Figure 2: comprehensive biomarker analysis with a corrected
  Panel F title ("F. Survival Significance" wrapped, with "(Log-Rank Test)").
- Supplementary Figure 4: detailed liquid-biopsy analysis, intentionally
  different panels from main Figure 5.
- Supplementary Figure 5: detailed methylation-protein discordance analysis,
  intentionally different panels from main Figure 6.

Main manuscript figures are NOT modified by this script.
"""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings

warnings.filterwarnings("ignore")

plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"] = 10
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "stix"

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS = SCRIPT_DIR / "results"
SUPP_DIR = SCRIPT_DIR / "figures" / "Supplementary"
SUPP_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================================
# Supplementary Figure 2 - Biomarker analysis (Panel F title fixed)
# ============================================================================
def make_supp_fig2() -> None:
    biomarker_df = pd.read_csv(
        RESULTS / "biomarker_candidates_high_meth_negative_protein.csv"
    )
    survival_df = pd.read_csv(RESULTS / "biomarker_survival_analysis.csv")

    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.45, wspace=0.32)
    top_20 = biomarker_df.head(20)

    # A. Methylation vs Protein correlation
    ax1 = fig.add_subplot(gs[0, 0])
    sc = ax1.scatter(
        top_20["Methylation_Mean"], top_20["Meth_Protein_Corr"], s=150, alpha=0.7,
        c=-np.log10(top_20["Meth_Protein_P"] + 1e-10), cmap="viridis",
        edgecolors="black", linewidth=0.5,
    )
    ax1.set_xlabel("Mean Methylation (Beta Value)", fontweight="bold")
    ax1.set_ylabel("Methylation-Protein Correlation", fontweight="bold")
    ax1.set_title("a. Methylation vs Protein Correlation\n(Top 20 Biomarker Candidates)", fontweight="bold")
    ax1.axhline(0, color="red", linestyle="--", alpha=0.5, linewidth=1)
    ax1.axvline(0.7, color="red", linestyle="--", alpha=0.5, linewidth=1, label="High Meth Threshold")
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    plt.colorbar(sc, ax=ax1).set_label("-log10(P-value)", fontweight="bold")
    for _, row in top_20.head(5).iterrows():
        ax1.annotate(row["Gene"], (row["Methylation_Mean"], row["Meth_Protein_Corr"]),
                     fontsize=8, alpha=0.8, fontweight="bold")

    # B. Biomarker score ranking
    ax2 = fig.add_subplot(gs[0, 1])
    colors = plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, len(top_20)))
    ax2.barh(range(len(top_20)), top_20["Biomarker_Score"], color=colors)
    ax2.set_yticks(range(len(top_20)))
    ax2.set_yticklabels(top_20["Gene"], fontsize=8)
    ax2.set_xlabel("Biomarker Score", fontweight="bold")
    ax2.set_title("b. Biomarker Potential Score\n(Ranked by Score)", fontweight="bold")
    ax2.grid(True, alpha=0.3, axis="x")
    ax2.invert_yaxis()

    # C. Protein difference
    ax3 = fig.add_subplot(gs[0, 2])
    colors_diff = ["coral" if x < 0 else "steelblue" for x in top_20["Protein_Diff_High_vs_Low_Meth"]]
    ax3.barh(range(len(top_20)), top_20["Protein_Diff_High_vs_Low_Meth"], color=colors_diff)
    ax3.set_yticks(range(len(top_20)))
    ax3.set_yticklabels(top_20["Gene"], fontsize=8)
    ax3.set_xlabel("Protein Difference\n(High Meth - Low Meth)", fontweight="bold")
    ax3.set_title("c. Protein Level Difference\n(High vs Low Methylation)", fontweight="bold")
    ax3.axvline(0, color="black", linestyle="-", alpha=0.5, linewidth=1)
    ax3.grid(True, alpha=0.3, axis="x")
    ax3.invert_yaxis()

    # D. Correlation vs significance (volcano-like)
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.scatter(
        biomarker_df["Meth_Protein_Corr"], -np.log10(biomarker_df["Meth_Protein_P"] + 1e-10),
        alpha=0.6, s=80, c=biomarker_df["Methylation_Mean"], cmap="coolwarm",
        edgecolors="black", linewidth=0.3,
    )
    ax4.set_xlabel("Methylation-Protein Correlation", fontweight="bold")
    ax4.set_ylabel("-log10(P-value)", fontweight="bold")
    ax4.set_title("d. Correlation vs Significance\n(All Biomarker Candidates)", fontweight="bold")
    ax4.axhline(-np.log10(0.05), color="red", linestyle="--", alpha=0.7, label="P=0.05")
    ax4.axvline(0, color="black", linestyle="-", alpha=0.3)
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    plt.colorbar(ax4.collections[0], ax=ax4).set_label("Mean Methylation", fontweight="bold")
    for _, row in top_20.head(5).iterrows():
        ax4.scatter(row["Meth_Protein_Corr"], -np.log10(row["Meth_Protein_P"] + 1e-10),
                    s=200, marker="*", color="yellow", edgecolors="black", linewidth=1, zorder=5)
        ax4.annotate(row["Gene"], (row["Meth_Protein_Corr"], -np.log10(row["Meth_Protein_P"] + 1e-10)),
                     fontsize=7, alpha=0.9, fontweight="bold", zorder=6)

    # E. High methylation prevalence
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.barh(range(len(top_20)), top_20["High_Meth_Percentage"], color="teal", alpha=0.7)
    ax5.set_yticks(range(len(top_20)))
    ax5.set_yticklabels(top_20["Gene"], fontsize=8)
    ax5.set_xlabel("High Methylation Patients (%)", fontweight="bold")
    ax5.set_title("e. High Methylation Prevalence\n(% Patients with Meth > 0.7)", fontweight="bold")
    ax5.set_xlim(0, 100)
    ax5.grid(True, alpha=0.3, axis="x")
    ax5.invert_yaxis()

    # F. Survival significance (FIXED TITLE - wrapped + Log-Rank Test on its own line)
    ax6 = fig.add_subplot(gs[1, 2])
    sig_genes = survival_df[survival_df["Is_Significant"] == True].sort_values("LogRank_P")
    if len(sig_genes) > 0:
        colors_surv = ["green" if p < 0.01 else "orange" for p in sig_genes["LogRank_P"]]
        ax6.barh(range(len(sig_genes)), -np.log10(sig_genes["LogRank_P"] + 1e-10),
                 color=colors_surv, alpha=0.7)
        ax6.set_yticks(range(len(sig_genes)))
        ax6.set_yticklabels(sig_genes["Gene"], fontsize=8)
        ax6.set_xlabel("-log10(P-value)", fontweight="bold")
        ax6.axvline(-np.log10(0.05), color="red", linestyle="--", alpha=0.7, label="P=0.05")
        ax6.grid(True, alpha=0.3, axis="x")
        ax6.legend(fontsize=8, loc="lower right")
        ax6.invert_yaxis()
    else:
        ax6.text(0.5, 0.5, "No significant\nsurvival differences",
                 ha="center", va="center", transform=ax6.transAxes, fontsize=12)
    ax6.set_title("f. Survival Significance\n(Log-Rank Test)", fontweight="bold")

    # G. Top 10 summary table
    ax7 = fig.add_subplot(gs[2, :])
    ax7.axis("off")
    summary = top_20.head(10)[["Gene", "Methylation_Mean", "Meth_Protein_Corr",
                               "Meth_Protein_P", "Biomarker_Score"]].copy()
    summary["Methylation_Mean"] = summary["Methylation_Mean"].round(3)
    summary["Meth_Protein_Corr"] = summary["Meth_Protein_Corr"].round(3)
    summary["Meth_Protein_P"] = summary["Meth_Protein_P"].apply(lambda x: f"{x:.2e}")
    summary["Biomarker_Score"] = summary["Biomarker_Score"].round(3)
    summary.columns = ["Gene", "Methylation\n(Mean)", "Correlation", "P-value", "Biomarker\nScore"]
    table = ax7.table(cellText=summary.values, colLabels=summary.columns,
                      cellLoc="center", loc="center", bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 2)
    for j in range(len(summary.columns)):
        table[(0, j)].set_facecolor("#4CAF50")
        table[(0, j)].set_text_props(weight="bold", color="white")
    for i in range(1, len(summary) + 1):
        for j in range(len(summary.columns)):
            if i % 2 == 0:
                table[(i, j)].set_facecolor("#f0f0f0")
    ax7.set_title("g. Top 10 Biomarker Candidates Summary", fontweight="bold", pad=20)

    plt.suptitle("High Methylation \u2192 Negative Protein Impact: Biomarker Candidates",
                 fontsize=14, fontweight="bold", y=0.995)
    out = SUPP_DIR / "Supplementary_Figure_2.png"
    plt.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved {out}")


# ============================================================================
# Supplementary Figure 4 - Liquid biopsy (distinct from main Figure 5)
# ============================================================================
def make_supp_fig4() -> None:
    all_df = pd.read_csv(RESULTS / "novel_analyses" / "liquid_biopsy_all.csv")
    top_df = pd.read_csv(RESULTS / "novel_analyses" / "liquid_biopsy_top_candidates.csv")
    top20 = top_df.head(20)
    threshold = top_df["Liquid_Biopsy_Score"].min()

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.subplots_adjust(hspace=0.32, wspace=0.28)

    # A. Detectability components: Mean Methylation vs CV (all genes), top highlighted
    ax = axes[0, 0]
    sc = ax.scatter(all_df["Mean_Methylation"], all_df["Coefficient_of_Variation"],
                    c=all_df["Liquid_Biopsy_Score"], cmap="viridis", s=18, alpha=0.35,
                    edgecolors="none")
    ax.scatter(top20["Mean_Methylation"], top20["Coefficient_of_Variation"],
               s=90, c="red", edgecolors="black", linewidth=1.2, zorder=5,
               label="Top 20 candidates")
    for i in range(min(3, len(top20))):
        ax.annotate(top20.iloc[i]["Gene"].split(";")[0],
                    (top20.iloc[i]["Mean_Methylation"], top20.iloc[i]["Coefficient_of_Variation"]),
                    fontsize=9, fontweight="bold")
    ax.set_xlabel("Mean Methylation (Beta)", fontweight="bold")
    ax.set_ylabel("Coefficient of Variation", fontweight="bold")
    ax.set_title("a. Detectability Components\n(High Methylation, Low Variability Preferred)", fontweight="bold")
    ax.set_ylim(0, min(2.0, all_df["Coefficient_of_Variation"].quantile(0.99)))
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    plt.colorbar(sc, ax=ax).set_label("Liquid Biopsy Score", fontweight="bold")

    # B. Distribution of liquid biopsy scores (all 4,111 genes)
    ax = axes[0, 1]
    ax.hist(all_df["Liquid_Biopsy_Score"], bins=60, color="gold", alpha=0.8, edgecolor="black", linewidth=0.3)
    ax.axvline(threshold, color="red", linestyle="--", linewidth=1.5,
               label=f"Top 50 threshold ({threshold:.2f})")
    ax.set_xlabel("Liquid Biopsy Score", fontweight="bold")
    ax.set_ylabel("Number of Genes", fontweight="bold")
    ax.set_title(f"b. Genome-wide Score Distribution\n(n = {len(all_df):,} genes)", fontweight="bold")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3, axis="y")
    ax.legend(fontsize=8)

    # C. Prognostic landscape: survival correlation vs significance (all genes)
    ax = axes[1, 0]
    ax.scatter(all_df["Survival_Correlation"], -np.log10(all_df["Survival_P_value"] + 1e-10),
               s=18, alpha=0.3, c="steelblue", edgecolors="none")
    ax.scatter(top20["Survival_Correlation"], -np.log10(top20["Survival_P_value"] + 1e-10),
               s=90, c="red", edgecolors="black", linewidth=1.2, zorder=5, label="Top 20 candidates")
    ax.axhline(-np.log10(0.05), color="red", linestyle="--", alpha=0.7, label="P = 0.05")
    ax.axvline(0, color="black", linestyle="-", alpha=0.3)
    for i in range(min(3, len(top20))):
        ax.annotate(top20.iloc[i]["Gene"].split(";")[0],
                    (top20.iloc[i]["Survival_Correlation"], -np.log10(top20.iloc[i]["Survival_P_value"] + 1e-10)),
                    fontsize=9, fontweight="bold")
    ax.set_xlabel("Survival Correlation", fontweight="bold")
    ax.set_ylabel("-log10(Survival P-value)", fontweight="bold")
    ax.set_title("c. Prognostic Landscape\n(Survival Correlation vs Significance)", fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    # D. Top 20 candidates - methylation level with variability annotation
    ax = axes[1, 1]
    y = range(len(top20))
    ax.barh(list(y), top20["Mean_Methylation"], color="mediumseagreen", alpha=0.8,
            edgecolor="black", linewidth=0.4)
    ax.set_yticks(list(y))
    ax.set_yticklabels(top20["Gene"], fontsize=9, fontweight="bold")
    ax.set_xlabel("Mean Methylation (Beta)", fontweight="bold")
    ax.set_title("d. Top 20 Candidates: Methylation Level\n(annotated with CV)", fontweight="bold")
    ax.set_xlim(0, 1.05)
    for yi, (_, row) in zip(y, top20.iterrows()):
        ax.text(row["Mean_Methylation"] + 0.01, yi, f"CV={row['Coefficient_of_Variation']:.2f}",
                va="center", fontsize=7)
    ax.grid(True, alpha=0.3, axis="x")
    ax.invert_yaxis()

    out = SUPP_DIR / "Supplementary_Figure_4.png"
    plt.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved {out}")


# ============================================================================
# Supplementary Figure 5 - Discordance (distinct from main Figure 6)
# ============================================================================
def make_supp_fig5() -> None:
    disc_file = RESULTS / "methylation_protein_discordance_full.csv"
    if not disc_file.exists():
        disc_file = RESULTS / "methylation_protein_discordance.csv"
    df = pd.read_csv(disc_file)
    df["Abs_Corr"] = df["Meth_Protein_Corr"].abs()

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.subplots_adjust(hspace=0.33, wspace=0.28)

    # A. Distribution of methylation-protein correlations (all genes)
    ax = axes[0, 0]
    ax.hist(df["Meth_Protein_Corr"], bins=40, color="slateblue", alpha=0.8,
            edgecolor="black", linewidth=0.3)
    ax.axvspan(-0.2, 0.2, color="coral", alpha=0.2, label="Discordant band (|r| < 0.2)")
    ax.axvline(0, color="black", linestyle="-", alpha=0.5)
    ax.set_xlabel("Methylation-Protein Correlation (r)", fontweight="bold")
    ax.set_ylabel("Number of Genes", fontweight="bold")
    ax.set_title(f"a. Correlation Distribution\n(n = {len(df)} genes)", fontweight="bold")
    ax.grid(True, alpha=0.3, axis="y")
    ax.legend(fontsize=8)

    # B. Discordance category composition (pie)
    ax = axes[0, 1]
    counts = df["Discordance_Type"].value_counts()
    palette = {"Discordant": "coral", "Moderate": "gray",
               "Concordant_High": "green", "Concordant_Negative": "red"}
    pie_colors = [palette.get(k, "steelblue") for k in counts.index]
    ax.pie(counts.values, labels=counts.index, colors=pie_colors, autopct="%1.1f%%",
           startangle=90, wedgeprops={"edgecolor": "white", "linewidth": 1.5},
           textprops={"fontsize": 10, "fontweight": "bold"})
    ax.set_title("b. Discordance Category Composition", fontweight="bold")

    # C. Absolute correlation vs significance (weak-correlation predominance)
    ax = axes[1, 0]
    for dt in df["Discordance_Type"].unique():
        sub = df[df["Discordance_Type"] == dt]
        ax.scatter(sub["Abs_Corr"], -np.log10(sub["P_value"] + 1e-10),
                   s=55, alpha=0.6, label=dt, color=palette.get(dt, "steelblue"),
                   edgecolors="black", linewidth=0.3)
    ax.axvline(0.2, color="orange", linestyle="--", alpha=0.7, label="Discordant threshold")
    ax.axhline(-np.log10(0.05), color="red", linestyle="--", alpha=0.7, label="P = 0.05")
    ax.set_xlabel("|Methylation-Protein Correlation|", fontweight="bold")
    ax.set_ylabel("-log10(P-value)", fontweight="bold")
    ax.set_title("c. Weak-Correlation Predominance\n(Absolute Correlation vs Significance)", fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    # D. Mean methylation distribution by discordance type (boxplot)
    ax = axes[1, 1]
    order = [t for t in ["Discordant", "Moderate", "Concordant_High", "Concordant_Negative"]
             if t in df["Discordance_Type"].unique()]
    box_data = [df.loc[df["Discordance_Type"] == t, "Meth_Mean"].values for t in order]
    bp = ax.boxplot(box_data, patch_artist=True, labels=order, showmeans=True)
    for patch, t in zip(bp["boxes"], order):
        patch.set_facecolor(palette.get(t, "steelblue"))
        patch.set_alpha(0.6)
    ax.set_ylabel("Mean Methylation (Beta)", fontweight="bold")
    ax.set_title("d. Methylation Level by Discordance Type", fontweight="bold")
    ax.tick_params(axis="x", rotation=20)
    ax.grid(True, alpha=0.3, axis="y")

    out = SUPP_DIR / "Supplementary_Figure_5.png"
    plt.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved {out}")


# ============================================================================
# Supplementary Figure 3 - Age/sex methylation (distinct from main Figure 4)
# ============================================================================
def make_supp_fig3() -> None:
    male = pd.read_csv(RESULTS / "age_sex_methylation_increasing_male.csv")
    female = pd.read_csv(RESULTS / "age_sex_methylation_increasing_female.csv")
    male_c, female_c = "#2c7fb8", "#dd3497"

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.subplots_adjust(hspace=0.32, wspace=0.26)

    # A. Distribution of age-methylation correlations by sex
    ax = axes[0, 0]
    bins = np.linspace(0, max(male["Spearman_Correlation"].max(),
                              female["Spearman_Correlation"].max()), 35)
    ax.hist(male["Spearman_Correlation"], bins=bins, alpha=0.6, color=male_c,
            label=f"Male (n = {len(male)})", edgecolor="white", linewidth=0.3)
    ax.hist(female["Spearman_Correlation"], bins=bins, alpha=0.6, color=female_c,
            label=f"Female (n = {len(female)})", edgecolor="white", linewidth=0.3)
    ax.set_xlabel("Age-Methylation Spearman Correlation", fontweight="bold")
    ax.set_ylabel("Number of Genes", fontweight="bold")
    ax.set_title("a. Age-Correlation Distribution by Sex\n(Genes with Age-Related Increase)", fontweight="bold")
    ax.grid(True, alpha=0.3, axis="y")
    ax.legend(fontsize=9)

    # B. Volcano: correlation vs significance, by sex (strongest gene labeled)
    ax = axes[0, 1]
    ax.scatter(male["Spearman_Correlation"], -np.log10(male["Spearman_P"] + 1e-10),
               s=22, alpha=0.4, color=male_c, edgecolors="none", label="Male")
    ax.scatter(female["Spearman_Correlation"], -np.log10(female["Spearman_P"] + 1e-10),
               s=22, alpha=0.4, color=female_c, edgecolors="none", label="Female")
    ax.axhline(-np.log10(0.05), color="red", linestyle="--", alpha=0.7, label="P = 0.05")
    for df_s, col in [(male, male_c), (female, female_c)]:
        top = df_s.sort_values("Spearman_Correlation", ascending=False).iloc[0]
        ax.scatter(top["Spearman_Correlation"], -np.log10(top["Spearman_P"] + 1e-10),
                   s=140, marker="*", color="yellow", edgecolors=col, linewidth=1.5, zorder=5)
        ax.annotate(str(top["Gene"]).split(";")[0],
                    (top["Spearman_Correlation"], -np.log10(top["Spearman_P"] + 1e-10)),
                    fontsize=9, fontweight="bold")
    ax.set_xlabel("Age-Methylation Spearman Correlation", fontweight="bold")
    ax.set_ylabel("-log10(Spearman P-value)", fontweight="bold")
    ax.set_title("b. Significance Landscape by Sex\n(Volcano of Age-Methylation Association)", fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    # C. Age trajectories (Young -> Mid -> Old) for top 5 genes per sex
    ax = axes[1, 0]
    stages = ["Young", "Mid", "Old"]
    xpos = [0, 1, 2]
    for _, row in male.sort_values("Spearman_Correlation", ascending=False).head(5).iterrows():
        ax.plot(xpos, [row["Young_Mean"], row["Mid_Mean"], row["Old_Mean"]],
                marker="o", color=male_c, alpha=0.8, linewidth=1.8)
        ax.annotate(str(row["Gene"]).split(";")[0], (2, row["Old_Mean"]),
                    fontsize=7, color=male_c, fontweight="bold", xytext=(4, 0),
                    textcoords="offset points", va="center")
    for _, row in female.sort_values("Spearman_Correlation", ascending=False).head(5).iterrows():
        ax.plot(xpos, [row["Young_Mean"], row["Mid_Mean"], row["Old_Mean"]],
                marker="s", color=female_c, alpha=0.8, linewidth=1.8)
        ax.annotate(str(row["Gene"]).split(";")[0], (2, row["Old_Mean"]),
                    fontsize=7, color=female_c, fontweight="bold", xytext=(4, 0),
                    textcoords="offset points", va="center")
    ax.set_xticks(xpos)
    ax.set_xticklabels(["Young\n(0-40)", "Mid\n(41-55)", "Old\n(56+)"])
    ax.set_ylabel("Mean Methylation (Beta)", fontweight="bold")
    ax.set_title("c. Age Trajectories of Top Genes\n(Top 5 Male = circles, Top 5 Female = squares)", fontweight="bold")
    ax.set_xlim(-0.2, 2.6)
    ax.grid(True, alpha=0.3)
    handles = [plt.Line2D([0], [0], color=male_c, marker="o", label="Male"),
               plt.Line2D([0], [0], color=female_c, marker="s", label="Female")]
    ax.legend(handles=handles, fontsize=9)

    # D. Magnitude of age-related increase (Old - Young) by sex
    ax = axes[1, 1]
    box_data = [male["Increase_Young_to_Old"].values, female["Increase_Young_to_Old"].values]
    bp = ax.boxplot(box_data, patch_artist=True, labels=["Male", "Female"], showmeans=True)
    for patch, col in zip(bp["boxes"], [male_c, female_c]):
        patch.set_facecolor(col)
        patch.set_alpha(0.6)
    ax.set_ylabel("Methylation Increase (Old - Young)", fontweight="bold")
    ax.set_title("d. Magnitude of Age-Related Increase by Sex", fontweight="bold")
    ax.grid(True, alpha=0.3, axis="y")

    out = SUPP_DIR / "Supplementary_Figure_3.png"
    plt.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved {out}")


if __name__ == "__main__":
    make_supp_fig2()
    make_supp_fig3()
    make_supp_fig4()
    make_supp_fig5()
    print("Done.")
