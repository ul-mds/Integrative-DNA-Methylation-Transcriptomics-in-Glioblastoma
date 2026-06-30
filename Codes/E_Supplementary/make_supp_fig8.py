#!/usr/bin/env python3
"""
Generate Supplementary Figure 8 (EGFR methylation in the context of
amplification and mutation) into figures/Supplementary/.

Every panel carries a lowercase letter label (a-f) with a short title.
Only the figure styling (panel letters, short titles, font) differs from the
original; the underlying data and plot content are unchanged.
"""
import pandas as pd, numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "stix"
from pathlib import Path

np.random.seed(0)
df = pd.read_csv('results/egfr_genetic_analysis/egfr_patient_level_with_cna.csv')
out = Path('figures/Supplementary/Supplementary_Figure_8.png')

plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12, 'savefig.dpi': 300})
fig = plt.figure(figsize=(15, 9))
gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.32)

C_NONAMP, C_AMP = '#4C72B0', '#C44E52'
gistic_colors = ['#55A868', '#DD8452', '#C44E52']


def annotate_p(ax, p):
    txt = 'p < 0.001' if p < 1e-3 else f'p = {p:.3g}'
    ax.text(0.5, 0.97, txt, transform=ax.transAxes, ha='center', va='top', fontsize=10,
            bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='gray', alpha=0.8))


# Row a-c: by GISTIC copy-number category (0,1,2)
labels_g = ['Diploid\n(0)', 'Gain\n(1)', 'Amp\n(2)']
panelABC = [
    ('methylation', 'EGFR methylation (β)', 'a. Methylation by CNA'),
    ('mrna', 'EGFR mRNA (FPKM)', 'b. mRNA by CNA'),
    ('protein', 'EGFR protein (ratio)', 'c. Protein by CNA'),
]
for j, (col, ylab, title) in enumerate(panelABC):
    ax = fig.add_subplot(gs[0, j])
    data = [df.loc[df['gistic'] == g, col].dropna() for g in [0, 1, 2]]
    bp = ax.boxplot(data, patch_artist=True, widths=0.6, tick_labels=labels_g, showfliers=False)
    for patch, c in zip(bp['boxes'], gistic_colors):
        patch.set_facecolor(c); patch.set_alpha(0.6)
    for med in bp['medians']:
        med.set_color('black')
    for i, d in enumerate(data):
        x = np.random.normal(i + 1, 0.06, len(d))
        ax.scatter(x, d, s=10, c='black', alpha=0.35, zorder=3)
    if col == 'mrna':
        ax.set_yscale('log')
    ax.set_ylabel(ylab)
    ax.set_xlabel('GISTIC copy-number')
    valid = [d for d in data if len(d) > 0]
    if len(valid) >= 2:
        _, p = stats.kruskal(*valid)
        annotate_p(ax, p)
    ax.set_title(title, loc='left', fontweight='bold')

# Row d: methylation-protein scatter amp vs non-amp
axB = fig.add_subplot(gs[1, 0])
for val, name, c in [(False, 'Non-amplified', C_NONAMP), (True, 'Amplified', C_AMP)]:
    s = df[df['EGFR_amp'] == val]
    axB.scatter(s['methylation'], s['protein'], label=name, alpha=0.7, c=c,
                edgecolor='white', linewidth=0.5, s=45)
    sub = s[['methylation', 'protein']].dropna()
    if len(sub) >= 8:
        z = np.polyfit(sub['methylation'], sub['protein'], 1)
        xs = np.linspace(sub['methylation'].min(), sub['methylation'].max(), 50)
        axB.plot(xs, np.polyval(z, xs), c=c, lw=2, ls='--')
axB.set_xlabel('EGFR methylation (β)'); axB.set_ylabel('EGFR protein (ratio)')
axB.legend(fontsize=9, loc='upper right')
axB.set_title('d. Methylation–protein', loc='left', fontweight='bold')
ra, _ = stats.spearmanr(*df[df['EGFR_amp']][['methylation', 'protein']].dropna().values.T)
rn, _ = stats.spearmanr(*df[~df['EGFR_amp']][['methylation', 'protein']].dropna().values.T)
axB.text(0.03, 0.04, f'Amp: ρ={ra:.2f}\nNon-amp: ρ={rn:.2f}', transform=axB.transAxes, fontsize=9,
         va='bottom', bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='gray', alpha=0.85))

# Row e-f: mutant vs WT for protein and mRNA
panelEF = [
    ('protein', 'EGFR protein (ratio)', gs[1, 1], 'e. Protein by mutation'),
    ('mrna', 'EGFR mRNA (FPKM)', gs[1, 2], 'f. mRNA by mutation'),
]
for col, ylab, pos, title in panelEF:
    ax = fig.add_subplot(pos)
    data = [df.loc[~df['EGFR_mut'], col].dropna(), df.loc[df['EGFR_mut'], col].dropna()]
    bp = ax.boxplot(data, patch_artist=True, widths=0.55, tick_labels=['WT', 'Mutant'], showfliers=False)
    for patch, c in zip(bp['boxes'], [C_NONAMP, C_AMP]):
        patch.set_facecolor(c); patch.set_alpha(0.6)
    for med in bp['medians']:
        med.set_color('black')
    for i, d in enumerate(data):
        x = np.random.normal(i + 1, 0.06, len(d))
        ax.scatter(x, d, s=12, c='black', alpha=0.35, zorder=3)
    if col == 'mrna':
        ax.set_yscale('log')
    ax.set_ylabel(ylab); ax.set_xlabel('EGFR mutation status')
    _, p = stats.mannwhitneyu(data[1], data[0], alternative='two-sided')
    annotate_p(ax, p)
    ax.set_title(title, loc='left', fontweight='bold')

fig.suptitle('Supplementary Figure 8. EGFR amplification and mutation (CPTAC GBM, n=97)',
             fontsize=13, fontweight='bold', y=0.98)
plt.savefig(out, bbox_inches='tight', facecolor='white')
print('Saved', out)
from PIL import Image
print('Size', Image.open(out).size)
