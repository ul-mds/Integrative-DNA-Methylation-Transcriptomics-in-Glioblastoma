#!/usr/bin/env python3
"""
Add validation results to Supplementary_Tables.xlsx.
- S18: TCGA internal validation (train-test split)
- S19: External validation results (GSE7696, GSE230770)
- S20: External validation summary by biomarker
Run from Methylation_Transcriptomics/.
"""
import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parent
SUPP_TABLES = BASE / "Supplementary" / "Supplementary_Tables.xlsx"
VALIDATION_DIR = BASE / "validation" / "results"


def add_validation_tables():
    # Load existing Excel
    with pd.ExcelFile(SUPP_TABLES) as xl:
        existing_sheets = xl.sheet_names

    # Check if validation tables already exist
    if any(s.startswith("S18") or s.startswith("S19") or s.startswith("S20") for s in existing_sheets):
        print("Validation tables may already exist. Overwriting...")

    # Load validation CSVs
    tcga_internal = pd.read_csv(VALIDATION_DIR / "tcga_internal_validation.csv")
    external_results = pd.read_csv(VALIDATION_DIR / "external_validation_results.csv")
    external_summary = pd.read_csv(VALIDATION_DIR / "external_validation_summary.csv")

    # Format TCGA internal validation for S18
    tcga_display = tcga_internal.copy()
    tcga_display = tcga_display.rename(columns={
        "Train_P": "Train p-value",
        "Train_Sig": "Train Significant",
        "Train_Diff": "Train Survival Diff (days)",
        "Test_P": "Test p-value",
        "Test_Sig": "Test Significant",
        "Test_Diff": "Test Survival Diff (days)",
        "Concordant": "Directionally Concordant",
    })
    tcga_display["Train Significant"] = tcga_display["Train Significant"].map({True: "Yes", False: "No"})
    tcga_display["Test Significant"] = tcga_display["Test Significant"].map({True: "Yes", False: "No"})
    tcga_display["Directionally Concordant"] = tcga_display["Directionally Concordant"].map({True: "Yes", False: "No"})
    tcga_display = tcga_display.round(4)

    # Format external results for S19
    ext_display = external_results.copy()
    ext_display["Significant"] = ext_display["Significant"].map({True: "Yes", False: "No"})
    ext_display["Concordant_TCGA"] = ext_display["Concordant_TCGA"].map({True: "Yes", False: "No"})
    ext_display = ext_display.rename(columns={
        "Data_Type": "Data Type",
        "N_High": "N High",
        "N_Low": "N Low",
        "Median_High": "Median Survival High (days)",
        "Median_Low": "Median Survival Low (days)",
        "Survival_Diff": "Survival Diff (days)",
        "LogRank_P": "Log-rank p-value",
        "Concordant_TCGA": "Concordant with TCGA",
        "TCGA_Diff": "TCGA Diff (days)",
        "TCGA_P": "TCGA p-value",
    })
    ext_display = ext_display.round(4)

    # Format external summary for S20
    ext_sum_display = external_summary.copy()
    ext_sum_display = ext_sum_display.rename(columns={
        "Gene": "Gene"
        , "TCGA_Diff": "TCGA Survival Diff (days)"
        , "N_Validated_Datasets": "N Datasets Validated"
        , "N_Significant": "N Significant"
        , "N_Concordant": "N Concordant"
        , "Validation_Rate": "Validation Rate"
    })

    # Write to Excel - append new sheets
    with pd.ExcelWriter(SUPP_TABLES, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        # S18: TCGA Internal Validation
        tcga_display.to_excel(writer, sheet_name="S18_TCGA_Internal", index=False, startrow=1)
        ws18 = writer.sheets["S18_TCGA_Internal"]
        ws18.cell(1, 1, "Supplementary Table 18: TCGA internal validation (70/30 train-test split). Biomarker survival associations in held-out test set.")
        ws18.merge_cells("A1:H1")

        # S19: External Validation Results
        ext_display.to_excel(writer, sheet_name="S19_External_Validation", index=False, startrow=1)
        ws19 = writer.sheets["S19_External_Validation"]
        ws19.cell(1, 1, "Supplementary Table 19: External validation results. GSE7696 (expression, n=80); GSE230770 (methylation, n=51). Survival Diff = High group median - Low group median.")
        ws19.merge_cells("A1:M1")

        # S20: External Validation Summary
        ext_sum_display.to_excel(writer, sheet_name="S20_External_Summary", index=False, startrow=1)
        ws20 = writer.sheets["S20_External_Summary"]
        ws20.cell(1, 1, "Supplementary Table 20: External validation summary by biomarker. Validation Rate = concordant datasets / total datasets validated.")
        ws20.merge_cells("A1:F1")

    print("Added validation tables to Supplementary_Tables.xlsx:")
    print("  - S18_TCGA_Internal: TCGA train-test split results")
    print("  - S19_External_Validation: GSE7696 + GSE230770 detailed results")
    print("  - S20_External_Summary: Biomarker-level validation summary")
    return 0


if __name__ == "__main__":
    exit(add_validation_tables())
