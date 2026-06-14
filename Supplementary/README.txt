Supplementary Materials — Multi-Omics GBM Study
============================================

Files:
  1. Supplementary_Tables.xlsx     Tables S1–S20
  2. Supplementary_Figures.docx    Figures SF1–SF7 (captions match manuscript)
  3. Supplementary_Methods.docx    Extended methods (archive; not cited in main text)

npj Precision Oncology submission:
  - Main Methods are in the manuscript body (no Supplementary Methods file in submission PDF).
  - Merge Supplementary Tables + Supplementary Figures into one Supplementary Information PDF.

Regenerate all files:
  python3 update_supplementary_folder.py

Key updates (latest run):
  - S1: TCGA-CPTAC data source
  - S6: All 6 biomarkers + Cox regression (age/sex adjusted)
  - S14: Full discordant gene list (n=232)
  - S18–S20: Validation results from validation/results/*.csv
  - SF1–SF7: Manuscript-aligned captions and plot images

Note: Manuscript list text for Table 14 says "100 genes"; full analysis identifies 232 discordant genes (S14).
