# Methylation in Glioblastoma — Multi-Omics Integration

Integrative analysis of DNA methylation, mRNA expression, and protein abundance in **TCGA-CPTAC glioblastoma** ([`gbm_cptac_2021`](https://www.cbioportal.org/study/summary?id=gbm_cptac_2021)).

## Repository layout

| Path | Description |
|------|-------------|
| [`Codes/A_Core/`](Codes/A_Core/) | Core multi-omics integration and cascade discovery |
| [`Codes/B_Figures/`](Codes/B_Figures/) | Main figure generation (Figures 1–5) |
| [`Codes/C_Tables/`](Codes/C_Tables/) | Main table generation (Tables 1–5) |
| [`Codes/D_Validation/`](Codes/D_Validation/) | Internal (TCGA train/test) and external (GSE7696, GSE230770) validation |
| [`Codes/E_Supplementary/`](Codes/E_Supplementary/) | Supplementary figure generation (SF1–SF8) and `Supplementary/` packaging |
| [`Supplementary/`](Supplementary/) | Supplementary Tables S1–S20, Supplementary Figures SF1–SF8, extended methods |
| [`validation/`](validation/) | Validation inputs, scripts, and result outputs |

## Quick start

```bash
pip install -r requirements.txt
# Download gbm_cptac_2021 data from cBioPortal and place input files as described in Codes/README_GITHUB_CODE_PACKAGE.md

# Core analyses
python Codes/A_Core/methylation_integrated_analysis.py

# Validation
python validation/scripts/run_tcga_internal_validation.py
python validation/scripts/run_external_dataset_validation.py

# Regenerate supplementary materials (Tables, Figures, Methods)
python Codes/E_Supplementary/update_supplementary_folder.py
```

## Supplementary materials

- **Supplementary_Tables.xlsx** — Tables S1–S20 (cohort summary, cascade genes, biomarkers, Cox regression, validation)
- **Supplementary_Figures.docx** — Figures SF1–SF7 with panel-level captions
- **Supplementary_Methods.docx** — Extended methods (archive; main Methods are in the manuscript)

Regenerate after updating analysis outputs:

```bash
python Codes/E_Supplementary/update_supplementary_folder.py
```

## Citation

If you use this code or supplementary data, please cite the associated manuscript (in preparation).

## License

See repository settings on GitHub for license terms.
