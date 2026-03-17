#!/usr/bin/env python3
"""
Update Table 1 in MANUSCRIPT_NATURE_COMMUNICATIONS_FINAL.docx with top 20 cascade
genes from results/high_meth_low_expr_low_protein_filtered.csv (DMRTA2, MGMT, EGFR, etc.).
Run from Methylation_Transcriptomics/.
"""
from docx import Document
from docx.shared import Pt
import pandas as pd
from pathlib import Path

doc_path = Path("MANUSCRIPT_NATURE_COMMUNICATIONS_FINAL.docx")
csv_path = Path("results/high_meth_low_expr_low_protein_filtered.csv")

if not doc_path.exists():
    print(f"Not found: {doc_path}")
    exit(1)
if not csv_path.exists():
    print(f"Not found: {csv_path}")
    exit(1)

# Load CSV: dedupe by gene, top 20 by Meth_mRNA_Corr
df = pd.read_csv(csv_path)
df = df.sort_values("Meth_mRNA_Corr").groupby("Gene", as_index=False).first()
df = df.sort_values("Meth_mRNA_Corr").reset_index(drop=True).head(20)

def p_fmt(p):
    try:
        x = float(p)
        if x < 1e-10:
            return "<1e-10"
        if x < 1e-5:
            return "<1e-5"
        if x < 1e-4:
            return "<1e-4"
        if x < 1e-3:
            return "<1e-3"
        return f"{x:.2e}"
    except Exception:
        return str(p)

# Build table rows: Rank, Gene, Meth-mRNA Correlation, Meth-mRNA P-value, Meth-Protein Correlation, Meth-Protein P-value, mRNA-Protein Correlation, mRNA-Protein P-value
headers = ["Rank", "Gene", "Meth-mRNA Correlation", "Meth-mRNA P-value", "Meth-Protein Correlation", "Meth-Protein P-value", "mRNA-Protein Correlation", "mRNA-Protein P-value"]
rows = []
for i, row in df.iterrows():
    rows.append([
        str(i + 1),
        str(row["Gene"]),
        f"{row['Meth_mRNA_Corr']:.4f}",
        p_fmt(row["Meth_mRNA_P"]),
        f"{row['Meth_Protein_Corr']:.4f}",
        p_fmt(row["Meth_Protein_P"]),
        f"{row['mRNA_Protein_Corr']:.4f}",
        p_fmt(row["mRNA_Protein_P"]),
    ])

doc = Document(str(doc_path))
# Find first table that looks like Table 1 (has "Gene" in first row and several columns)
target_table = None
for t in doc.tables:
    if len(t.rows) < 2 or len(t.columns) < 4:
        continue
    first_row_text = " ".join(c.text.strip() for c in t.rows[0].cells).lower()
    if "gene" in first_row_text and ("meth" in first_row_text or "correlation" in first_row_text):
        target_table = t
        break

if target_table is None:
    print("No candidate Table 1 found in document (table with 'Gene' and correlation headers).")
    exit(1)

# Ensure table has enough rows (1 header + 20 data)
nrows_need = 1 + len(rows)
while len(target_table.rows) < nrows_need:
    target_table.add_row()
# Set header
for j, h in enumerate(headers):
    if j < len(target_table.rows[0].cells):
        target_table.rows[0].cells[j].text = h
# Set data rows
for i, row in enumerate(rows):
    r = i + 1
    if r >= len(target_table.rows):
        break
    for j, val in enumerate(row):
        if j < len(target_table.rows[r].cells):
            target_table.rows[r].cells[j].text = val

doc.save(str(doc_path))
print(f"Updated Table 1 in {doc_path} with top 20 cascade genes (DMRTA2, MGMT, EGFR, ...).")
