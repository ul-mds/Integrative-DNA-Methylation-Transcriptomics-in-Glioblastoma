"""Resolve repository root and input data paths for the publishable Codes package."""
from __future__ import annotations

from pathlib import Path

CODES_DIR = Path(__file__).resolve().parent
REPO_ROOT = CODES_DIR.parent

GITHUB_CODES_URL = (
    "https://github.com/ul-mds/"
    "Integrative-DNA-Methylation-Transcriptomics-in-Glioblastoma/tree/main/Codes"
)


def data_file(name: str) -> Path:
    """Path to a cBioPortal export file placed at the repository root."""
    return REPO_ROOT / name
