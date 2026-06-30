# Codes (GitHub-ready layout)

This folder contains the **minimal publishable code package**, organized into five sections matching the manuscript’s reproducibility needs.

## A) Core multi-omics analysis

See: `A_Core/`

## B) Main figure generation (Figures 1–5)

See: `B_Figures/`

## C) Main table generation (Tables 1–5)

See: `C_Tables/`

## D) Validation (internal + external)

See: `D_Validation/`

## E) Supplementary materials generation (Figures SF1–SF8, Tables, Methods)

See: `E_Supplementary/`

These scripts regenerate the submission-ready supplementary figures and package
the `../Supplementary/` folder (Supplementary Figures, Tables, and Methods).

For detailed rationale and run commands, see `README_GITHUB_CODE_PACKAGE.md`.

## Supplementary output folder (repository root)

Submission-ready supplementary files live in `../Supplementary/` and are regenerated with:

```bash
python Codes/E_Supplementary/update_supplementary_folder.py
```

See the top-level `README.md` for the full repository layout.
