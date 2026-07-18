#!/usr/bin/env python3
"""
Download GEO datasets for validation (GSE119774, GSE4271).
Uses direct FTP/HTTPS - no GEOparse required.
"""

import urllib.request
import gzip
import os
from pathlib import Path

EXTERNAL_DIR = Path(__file__).resolve().parent.parent / "external_cohort"

GEO_DATASETS = {
    "GSE119774": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119774/matrix/GSE119774_series_matrix.txt.gz",
    "GSE119774_beta": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119774/suppl/GSE119774_GBM_GSC_i850K_beta_formatted.txt.gz",
    "GSE4271": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE4nnn/GSE4271/matrix/GSE4271_series_matrix.txt.gz",
}


def download_file(url, dest_path):
    """Download file with basic error handling."""
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read()
        with open(dest_path, 'wb') as f:
            f.write(data)
        return True
    except Exception as e:
        print(f"Download failed: {e}")
        return False


def main():
    for name, url in GEO_DATASETS.items():
        base_id = name.replace("_beta", "")
        dest = EXTERNAL_DIR / base_id
        dest.mkdir(parents=True, exist_ok=True)
        fname = url.split("/")[-1]
        out_path = dest / fname
        if out_path.exists():
            print(f"{name}: Already exists")
            continue
        print(f"Downloading {name}...")
        if download_file(url, out_path):
            print(f"  -> {out_path}")
        else:
            print(f"  -> Failed (try manual download from GEO)")


if __name__ == "__main__":
    main()
