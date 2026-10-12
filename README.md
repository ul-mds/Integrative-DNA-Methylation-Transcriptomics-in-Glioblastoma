# Integrative DNA Methylation–Transcriptomics in Glioblastoma

Integrative analysis of DNA methylation, mRNA expression, and protein abundance in **TCGA-CPTAC glioblastoma** ([`gbm_cptac_2021`](https://www.cbioportal.org/study/summary?id=gbm_cptac_2021)).

This repository supports the **British Journal of Cancer (BJC)** submission: three main figures, three main tables, and supplementary materials SF1–SF11 / ST1–ST23.

**Code package:** [Codes/](Codes/) on GitHub  
**Repository:** [ul-mds/Integrative-DNA-Methylation-Transcriptomics-in-Glioblastoma](https://github.com/ul-mds/Integrative-DNA-Methylation-Transcriptomics-in-Glioblastoma)

## Repository layout

| Path | Description |
|------|-------------|
| [`Codes/A_Core/`](Codes/A_Core/) | Core multi-omics integration and cascade discovery |
| [`Codes/B_Figures/`](Codes/B_Figures/) | Main figure generation (BJC Figures 1–3) |
| [`Codes/C_Tables/`](Codes/C_Tables/) | Main table generation (BJC Tables 1–3) |
| [`Codes/D_Validation/`](Codes/D_Validation/) | Internal (TCGA train/test) and external (GEO) validation |
| [`Codes/E_Supplementary/`](Codes/E_Supplementary/) | Supplementary figures (SF1–SF11) and `Supplementary/` packaging |
| [`Supplementary/`](Supplementary/) | Submission supplementary tables, figures, and extended methods |
| [`Codes/D_Validation/validation/`](Codes/D_Validation/validation/) | Validation scripts, CSV summaries, and plots (canonical path for GitHub) |

## Input data (repository root)

Download the `gbm_cptac_2021` study from cBioPortal and place these files in the **repository root** (see `Codes/README_GITHUB_CODE_PACKAGE.md` for details):

- `methylation_epic.txt`
- `protein_quantification.txt`
- `data_mrna_seq_v2_rsem.txt` (or the mRNA file used in your export)
- `clinical_patient.txt`

## Quick start

```bash
pip install -r requirements.txt

# Core analyses (run from repository root)
python Codes/A_Core/methylation_integrated_analysis.py

# Validation
python Codes/D_Validation/validation/scripts/run_tcga_internal_validation.py
python Codes/D_Validation/validation/scripts/run_external_dataset_validation.py

# Regenerate supplementary materials
python Codes/E_Supplementary/update_supplementary_folder.py
```

## BJC main-text mapping

| BJC item | Script / output |
|----------|-----------------|
| Figure 1 | Cascade + survival — `Codes/B_Figures/create_high_meth_low_expr_low_protein_plots.py` |
| Figure 2 | Biomarkers — `Codes/B_Figures/high_methylation_negative_protein_analysis.py` |
| Figure 3 | Discordance — `Codes/B_Figures/full_discordance_analysis.py` |
| Table 1 | `Codes/C_Tables/sync_manuscript_table1_from_csv.py` |
| Table 2 | `Codes/C_Tables/final_table2_update.py` |
| Table 3 | Discordant genes — `Codes/C_Tables/compute_table5_discordant_mrna_protein.py` |

Age/sex and liquid-biopsy analyses (former main Figures 4–5 / Tables 3–4) are in supplementary materials ST22–ST23 and SF5–SF6; see `Codes/README_GITHUB_CODE_PACKAGE.md`.

## Citation

If you use this code or supplementary data, please cite the associated BJC manuscript (in submission).

## License

See repository settings on GitHub for license terms.
