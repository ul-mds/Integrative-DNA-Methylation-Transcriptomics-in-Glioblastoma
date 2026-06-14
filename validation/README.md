# Validation Framework for Nature Communications Submission

Bu klasör, **MANUSCRIPT_NATURE_COMMUNICATIONS_FINAL.docx** içindeki sonuçların Nature Communications gereksinimlerine uygun şekilde doğrulanması için kurulmuştur.

## Nature Communications Validation Gereksinimleri

- **Data availability**: Verilerin minimum dataset olarak erişilebilir olması
- **Code availability**: Analiz kodlarının Zenodo/Code Ocean üzerinden paylaşılması
- **Source data**: Figürlerin altındaki verilerin sağlanması
- **Reproducibility**: Sonuçların tekrarlanabilirliği
- **External validation**: Bağımsız kohortlarda doğrulama (önerilir)

## Klasör Yapısı

```
validation/
├── README.md                 # Bu dosya
├── VALIDATION_CHECKLIST.md   # Nature Communications doğrulama kontrol listesi
├── external_cohort/          # Bağımsız kohort validasyonu (GEO, CGGA, vb.)
├── reproducibility/          # Kod tekrarlanabilirlik testleri
├── source_data/              # Figürler için kaynak veri (Nature zorunlu)
├── statistical/              # İstatistiksel validasyon (bootstrap, permutation)
├── bootstrap/                # Bootstrap analizleri
├── cross_validation/        # K-fold cross-validation
└── scripts/                  # Validasyon scriptleri
```

## Manuel Sonuçları Doğrulanacak Ana Bulgular

| Bölüm | Ana Bulgular | Validasyon Türü |
|-------|--------------|-----------------|
| **Figure 1** | 1,678 cascade geni; RCN1, ST8SIA4, MYO1G, GNG4 survival | External cohort, reproducibility |
| **Figure 2** | 6 biomarker: B4GALT1, EGFR, VASP, GFPT2, PHF11, TSPYL5 | External cohort, survival validation |
| **Figure 3** | 842 erkek + 1,456 kadın age-dependent gen | External cohort |
| **Figure 4** | 50 liquid biopsy adayı; SEPT8;CCNI2 top | External cohort, detectability |
| **Figure 5** | 85% discordance (232/273 gen) | Reproducibility, statistical |

## Önerilen Dış Veri Kaynakları (External Validation)

- **GSE7696**: Expression (80 GBM), survival mevcut – **validation tamamlandı**
- **GSE230770**: Methylation EPIC (51 IDH-wt GBM), survival mevcut – **validation tamamlandı**
- **GSE119774**: Methylation 850K (50 GBM) – `scripts/download_geo_data.py` ile indir
- **GSE4271**: Gene expression
- **CGGA** (Chinese Glioma Genome Atlas): mRNAseq-325 + Methylation, survival mevcut
- **TCGA**: Internal validation (train/test split) – `scripts/run_tcga_internal_validation.py`

## Validation Scriptleri (Güncel)

| Script | Açıklama |
|--------|----------|
| `scripts/run_tcga_internal_validation.py` | TCGA %70/%30 split, 6 biomarker survival |
| `scripts/run_external_dataset_validation.py` | GSE7696, GSE230770 external validation |
| `scripts/run_validation.py` | GSE119774, GSE4271, CGGA external validation |
| `scripts/download_geo_data.py` | GEO veri indirme |

## Sonuçlar ve Plotlar

- `results/tcga_internal_validation.csv` – TCGA internal validation
- `results/external_validation_results.csv` – GSE7696, GSE230770 external validation
- `results/external_validation_summary.csv` – Özet (validation rate)
- `plots/tcga_internal_validation.png` – Train vs Test survival difference
- `plots/external_validation_comparison.png` – External validation karşılaştırma
- `VALIDATION_REPORT.md` – Detaylı rapor

## Hızlı Başlangıç

1. `python validation/scripts/run_tcga_internal_validation.py` – TCGA internal validation
2. `python validation/scripts/run_external_dataset_validation.py` – GSE7696, GSE230770 external validation
3. `python validation/scripts/download_geo_data.py` – GEO verilerini indir
4. CGGA: cgga.org.cn'dan indir, `external_cohort/CGGA/` klasörüne koy
5. `python validation/scripts/run_validation.py` – External validation
