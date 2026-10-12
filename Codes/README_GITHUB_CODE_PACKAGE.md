# Minimum Publishable Code Package (BJC — GBM Multi-Omics)

This folder is the **reviewer-facing code package** for the British Journal of Cancer manuscript on integrative DNA methylation, transcriptomics, and proteomics in glioblastoma.

**Public repository:** [ul-mds/Integrative-DNA-Methylation-Transcriptomics-in-Glioblastoma](https://github.com/ul-mds/Integrative-DNA-Methylation-Transcriptomics-in-Glioblastoma) (`Codes/` directory)

Shared path helpers: `repo_paths.py` (`REPO_ROOT`, `data_file()`).

## Data source (required)

Primary analyses use **TCGA-GBM / CPTAC GBM (2021)** via cBioPortal: [`gbm_cptac_2021`](https://www.cbioportal.org/study/summary?id=gbm_cptac_2021).

Place downloaded exports at the **repository root** (not inside `Codes/`):

- `methylation_epic.txt`
- `protein_quantification.txt`
- `data_mrna_seq_v2_rsem.txt`
- `clinical_patient.txt`

---

## A) Core multi-omics analysis

- `A_Core/gbm_epigenetic_transcriptomic_proteomic_analysis.py`
- `A_Core/methylation_integrated_analysis.py`
- `A_Core/novel_analyses_pipeline.py` *(optional orchestration)*

```bash
python Codes/A_Core/methylation_integrated_analysis.py
```

---

## B) Main figures (BJC Figures 1–3)

| BJC figure | Script |
|------------|--------|
| Figure 1 — cascade + survival | `B_Figures/create_high_meth_low_expr_low_protein_plots.py` |
| Figure 2 — biomarkers | `B_Figures/high_methylation_negative_protein_analysis.py` |
| Figure 3 — discordance | `B_Figures/full_discordance_analysis.py` |

**Moved to supplementary (same scripts):**

- Age/sex — `B_Figures/age_sex_methylation_analysis.py` → SF5, ST22
- Liquid biopsy — `B_Figures/update_liquid_biopsy_plot.py` → SF6, ST23

```bash
python Codes/B_Figures/create_high_meth_low_expr_low_protein_plots.py
python Codes/B_Figures/high_methylation_negative_protein_analysis.py
python Codes/B_Figures/full_discordance_analysis.py
```

---

## C) Main tables (BJC Tables 1–3)

| BJC table | Script |
|-----------|--------|
| Table 1 — cascade genes | `C_Tables/sync_manuscript_table1_from_csv.py` |
| Table 2 — biomarkers | `C_Tables/final_table2_update.py` |
| Table 3 — discordant genes | `C_Tables/compute_table5_discordant_mrna_protein.py` |

```bash
python Codes/C_Tables/sync_manuscript_table1_from_csv.py
python Codes/C_Tables/final_table2_update.py
python Codes/C_Tables/compute_table5_discordant_mrna_protein.py
```

---

## D) Validation (internal + external)

Under `D_Validation/validation/` (also linked as `validation/` at repo root):

- `scripts/run_tcga_internal_validation.py`
- `scripts/run_external_dataset_validation.py`
- `README.md`, `VALIDATION_REPORT.md`

```bash
python Codes/D_Validation/validation/scripts/run_tcga_internal_validation.py
python Codes/D_Validation/validation/scripts/run_external_dataset_validation.py
```

---

## E) Supplementary materials (SF1–SF11, ST1–ST23)

- `E_Supplementary/make_supplementary_figures.py`
- `E_Supplementary/fix_mgmt_supp_fig7.py` (SF4 MGMT in BJC numbering)
- `E_Supplementary/make_supp_fig8.py` (SF11 EGFR)
- `E_Supplementary/mgmt_similarity_analysis.py`
- `E_Supplementary/mgmt_similar_genes_summary.py`
- `E_Supplementary/update_supplementary_folder.py` — builds `Supplementary/` at repo root
- `E_Supplementary/sync_bjc_supplementary_parallel.py` — local BJC submission folder sync

```bash
python Codes/E_Supplementary/update_supplementary_folder.py
```

---

## Excluded from GitHub (local manuscript tooling)

Word automation, reference renumbering, debug/inspect scripts, and full `results/` / `plots/` outputs remain outside this package; regenerate outputs with the commands above.

---

## Code availability statement (manuscript)

Use in the BJC manuscript:

> Analysis code is available at https://github.com/ul-mds/Integrative-DNA-Methylation-Transcriptomics-in-Glioblastoma/tree/main/Codes
