# Minimum Publishable Code Package (GBM Multi-Omics Manuscript)

This folder documents the **minimal, editor-/reviewer-facing code package** required to reproduce the key analyses, figures, tables, and validation reported in the manuscript. The goal is to provide a **concise, scientifically relevant** set of scripts (without repository noise such as Word editing utilities, formatting helpers, or one-off debug scripts).

## Data source (required)

All primary analyses are based on the **TCGA-GBM / CPTAC GBM (2021) multi-omics study** as distributed via cBioPortal:

- Study page: [`gbm_cptac_2021`](https://www.cbioportal.org/study/summary?id=gbm_cptac_2021)

Download the study data (methylation, mRNA, protein, and clinical/survival) and place the resulting files into your local project according to the paths expected by the scripts (see “How to run” below). If you choose a different folder structure, update the input paths in the scripts accordingly.

---

## A) Core multi-omics analysis (primary findings)

Include these scripts because they implement the **main analysis logic** (multi-omics integration and cascade discovery) used throughout the Results section:

- `gbm_epigenetic_transcriptomic_proteomic_analysis.py`
- `methylation_integrated_analysis.py`

*(Optional, only if it is the top-level driver you actually use to run multiple modules end-to-end):*

- `novel_analyses_pipeline.py`

### How to run (A)

From the repository root:

```bash
python gbm_epigenetic_transcriptomic_proteomic_analysis.py
python methylation_integrated_analysis.py
```

---

## B) Main figure generation (Figures 1–5)

Include these scripts because they generate the **primary manuscript figures** from the processed/filtered data:

- **Figure 1 (cascade discovery/visualization)**: `create_high_meth_low_expr_low_protein_plots.py`
- **Figure 2 (biomarker analysis/plots)**: `high_methylation_negative_protein_analysis.py`
- **Figure 3 (age/sex methylation analysis/plots)**: `age_sex_methylation_analysis.py`
- **Figure 4 (liquid biopsy prioritization/plots)**: `update_liquid_biopsy_plot.py`
- **Figure 5 (discordance analysis/plots)**: `full_discordance_analysis.py`

### How to run (B)

```bash
python create_high_meth_low_expr_low_protein_plots.py
python high_methylation_negative_protein_analysis.py
python age_sex_methylation_analysis.py
python update_liquid_biopsy_plot.py
python full_discordance_analysis.py
```

---

## C) Main table generation (Tables 1–5)

Include these scripts because they build or synchronize the **core tables** reported in the manuscript:

- **Table 1 (top cascade genes)**: `sync_manuscript_table1_from_csv.py`
- **Table 2 (top biomarkers)**: `final_table2_update.py` *(use the single “final” Table 2 script only)*
- **Table 5 (discordant genes)**: `compute_table5_discordant_mrna_protein.py`

### How to run (C)

```bash
python sync_manuscript_table1_from_csv.py
python final_table2_update.py
python compute_table5_discordant_mrna_protein.py
```

---

## D) Validation framework (internal + external)

Include these files because they demonstrate **robustness and generalizability**, which is typically one of the most persuasive components for editors/reviewers:

- `validation/scripts/run_tcga_internal_validation.py` *(internal validation; train/test split)*
- `validation/scripts/run_external_dataset_validation.py` *(external validation; GEO/other cohorts)*
- `validation/README.md` *(how to run validation)*
- `validation/VALIDATION_REPORT.md` *(summary of validation outcomes and outputs)*

*(Optional, only if it is the canonical one-command entrypoint you use):*

- `validation/scripts/run_validation_pipeline.py`

### How to run (D)

```bash
python validation/scripts/run_tcga_internal_validation.py
python validation/scripts/run_external_dataset_validation.py
```

---

## E) Supplementary materials generation (Figures SF1–SF8, Tables, Methods)

Include these scripts because they generate the **supplementary figures** and package
the submission-ready `Supplementary/` folder (Supplementary Figures, Supplementary
Tables, and Supplementary Methods) reported alongside the manuscript:

- **Supplementary Figures 2, 4, 5 (detailed biomarker / liquid-biopsy / discordance panels)**: `make_supplementary_figures.py`
- **Supplementary Figure 7 (MGMT methylation similarity)**: `fix_mgmt_supp_fig7.py`
- **Supplementary Figure 8 (EGFR methylation vs. amplification and mutation)**: `make_supp_fig8.py`
- **MGMT-similarity analysis (basis of SF7)**: `mgmt_similarity_analysis.py`
- **MGMT-similar gene summary table**: `mgmt_similar_genes_summary.py`
- **Supplementary folder packaging (Tables S1–S20, Figures SF1–SF8, Methods)**: `update_supplementary_folder.py`

### How to run (E)

```bash
# 1) Regenerate supplementary figures
python E_Supplementary/make_supplementary_figures.py
python E_Supplementary/fix_mgmt_supp_fig7.py
python E_Supplementary/make_supp_fig8.py

# 2) (Optional) Recompute the MGMT-similarity analysis behind SF7
python E_Supplementary/mgmt_similarity_analysis.py
python E_Supplementary/mgmt_similar_genes_summary.py

# 3) Package the submission-ready Supplementary/ folder
#    (Supplementary_Tables.xlsx, Supplementary_Figures.docx, Supplementary_Methods.docx)
python E_Supplementary/update_supplementary_folder.py
```

---

## Explicitly excluded (to avoid unnecessary repository clutter)

Do **not** include the following categories in the GitHub “minimum publishable package”, unless a journal specifically requests them:

- **Word/Docx automation and highlighting utilities**: scripts containing `docx`, Word sync, yellow highlighting, etc.
- **Reference renumbering / formatting utilities**: `renumber_*`, `fix_reference_*`, `optimize_reference_order.py`, and similar.
- **One-off inspection/debug scripts**: `inspect_*`, `diagnose_*`, `debug_*`, `minimal_test.py`, `test_write.py`.
- **Text humanization / phrasing rewrite scripts**: `apply_humanized_*`, `apply_ai_*`.
- **Large outputs**: do not commit full `results/` or `plots/` unless required; instead, provide commands to regenerate them and optionally a small representative example.

---

## Practical GitHub recommendation (what editors expect)

- **Reproducibility**: add a single environment file (`requirements.txt` or `environment.yml`) and a short top-level README describing:
  - data source(s) (e.g., TCGA/CPTAC identifiers),
  - the exact command(s) to run A–D,
  - expected output filenames for Figures 1–5 and Tables 1–5.
- **Clarity**: keep this “minimum package” small and focused; it should be obvious how a reviewer can reproduce the headline results.

## Notes for editors/reviewers (practical expectations)

- These scripts are intended to be run **after** downloading the `gbm_cptac_2021` study data from cBioPortal and placing the data into the expected input locations.
- Reproducibility is best supported by providing:
  - a pinned environment file (`requirements.txt` or `environment.yml`),
  - a short top-level `README.md` describing input files, run commands, and expected outputs.

