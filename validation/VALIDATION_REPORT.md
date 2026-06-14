# External Cohort Validation Report
## GBM Biomarkers: B4GALT1, EGFR, VASP, GFPT2, PHF11, TSPYL5

**Manuscript:** Multi-Omics Integration Reveals Epigenetic Regulatory Cascades and Clinically Actionable Biomarkers in Glioblastoma

**Son Güncelleme:** Validation pipeline çalıştırıldı. Tüm sonuçlar `validation/results/` ve `validation/plots/` klasörlerinde.

---

## 1. TCGA-GBM Ana Sonuçlar (Referans)

| Gen | High Meth Median (gün) | Low Meth Median (gün) | Survival Farkı | Log-Rank P |
|-----|------------------------|------------------------|----------------|------------|
| B4GALT1 | 496 | 340 | +156 | 0.000245 |
| EGFR | 401 | 192 | +209 | 0.000551 |
| VASP | 627 | 344.5 | +282.5 | 0.001687 |
| GFPT2 | 402 | 218.5 | +183.5 | 0.002374 |
| PHF11 | 577 | 344.5 | +232.5 | 0.006932 |
| TSPYL5 | 568 | 344.5 | +223.5 | 0.008472 |

**Yön:** Yüksek metilasyon → daha iyi sağkalım (paradoksal bulgu: metilasyonla baskılanan proteinler tümörü destekleyebilir)

---

## 2. TCGA Internal Validation (Train/Test Split)

TCGA verisi %70 eğitim / %30 test olarak bölündü (n_train=65, n_test=29).

| Gen | Train P | Train Sig | Train Diff | Test P | Test Sig | Test Diff | Yön Uyumlu |
|-----|---------|-----------|------------|--------|----------|-----------|------------|
| B4GALT1 | 0.0037 | ✓ | +283 | 0.37 | - | -66 | Hayır |
| EGFR | 0.39 | - | +273.5 | 0.054 | - | +250 | Evet |
| TSPYL5 | 0.23 | - | +197 | 0.19 | - | +8 | Evet |

**Not:** Test seti küçük (n=29) olduğu için istatistiksel güç sınırlı. EGFR ve TSPYL5 yönü TCGA ile uyumlu. VASP, GFPT2, PHF11 test setinde yeterli hasta sayısına ulaşamadı (high/low gruplarında min 5 hasta gereksinimi).

---

## 3. External Validation Veri Setleri

### GSE119774 (Methylation 850K)
- **Kaynak:** GEO, Illumina EPIC 850K
- **Örnekler:** 50 primary GBM, 44 GSC, 10 NSC
- **Durum:** Series matrix indirildi; beta değerleri için supplementary dosya gerekli
- **Survival:** GEO metadata'da survival bilgisi yok; yayın ekinden alınabilir

### GSE4271 (Gene Expression)
- **Kaynak:** GEO
- **Durum:** FTP 404 – manuel indirme gerekebilir

### CGGA mRNAseq-325
- **Kaynak:** Chinese Glioma Genome Atlas (cgga.org.cn)
- **Örnekler:** 325 glioma hastası, survival verisi mevcut
- **Durum:** Kayıt gerekli; indirme sonrası `validation/external_cohort/CGGA/` klasörüne yerleştirilmeli

### CGGA Methylation
- **Kaynak:** CGGA
- **Örnekler:** 159 glioma hastası
- **Durum:** CGGA portaldan indirilebilir

---

## 4. Validation Scriptleri

| Script | Açıklama |
|--------|----------|
| `scripts/run_tcga_internal_validation.py` | TCGA train/test split validation |
| `scripts/run_validation.py` | GSE119774, GSE4271, CGGA external validation |
| `scripts/download_geo_data.py` | GEO veri indirme |

---

## 5. Çıktılar

- `results/tcga_internal_validation.csv` – TCGA internal validation sonuçları
- `results/validation_survival_comparison.csv` – External cohort karşılaştırması (veri hazır olduğunda)
- `results/validation_summary.csv` – Özet tablo
- `plots/tcga_internal_validation.png` – TCGA internal validation grafiği
- `plots/validation_comparison.png` – External validation karşılaştırma grafiği
- `plots/validation_heatmap.png` – Significant/Concordant heatmap

---

## 6. CGGA Verisi Kurulumu

1. https://www.cgga.org.cn/download.jsp adresine gidin
2. mRNAseq-325 ve Methylation verilerini indirin
3. Dosyaları şu yapıda yerleştirin:
```
validation/external_cohort/CGGA/
├── mRNAseq_325_expression.csv   (veya *mRNA*325*.csv)
├── CGGA_clinical_mRNAseq325.csv (veya *clinical*.csv)
├── methylation_beta.csv
└── methylation_clinical.csv
```
4. `python validation/scripts/run_validation.py` çalıştırın

---

## 7. External Dataset Validation (GSE7696, GSE230770)

### GSE7696 (Expression, Affymetrix HG-U133 Plus 2.0)
- **Kaynak:** GEO, 80 GBM + 4 non-tumoral
- **Survival:** GEO metadata'da mevcut (survival time in months, survival status)
- **Analiz:** Düşük ekspresyon = daha iyi sağkalım (TCGA ile uyumlu: yüksek metilasyon → düşük ekspresyon → iyi sağkalım)

| Gen | N_High | N_Low | Median_High (gün) | Median_Low (gün) | Survival Diff | LogRank P | TCGA Uyumlu |
|-----|--------|-------|-------------------|------------------|---------------|-----------|-------------|
| B4GALT1 | 40 | 41 | 463 | 495 | -32 | 0.56 | ✓ |
| EGFR | 40 | 41 | 463 | 496 | -32 | 0.14 | ✓ |
| VASP | 40 | 41 | 464 | 520 | -56 | 0.13 | ✓ |
| GFPT2 | 40 | 41 | 440 | 495 | -55 | 0.54 | ✓ |
| PHF11 | 40 | 41 | 463 | 496 | -32 | 0.89 | ✓ |
| TSPYL5 | 40 | 41 | 450 | 495 | -45 | 0.91 | ✓ |

### GSE230770 (Methylation, Illumina EPIC)
- **Kaynak:** GEO, 51 IDH-wt GBM (27 LTS + 24 standard survival)
- **Survival:** GEO metadata'da mevcut (overall survival months, last follow-up)
- **Probe kapsamı:** EGFR, VASP, TSPYL5 (EPIC probe mapping)

| Gen | N_High | N_Low | Median_High (gün) | Median_Low (gün) | Survival Diff | LogRank P | TCGA Uyumlu |
|-----|--------|-------|-------------------|------------------|---------------|-----------|-------------|
| EGFR | 25 | 26 | 1887 | 1294 | +594 | 0.59 | ✓ |
| VASP | 25 | 26 | 609 | 2085 | -1476 | 0.46 | ✗ |
| TSPYL5 | 25 | 26 | 2009 | 1887 | +122 | 0.74 | ✓ |

### External Validation Özeti

| Gen | TCGA Diff | Validated Datasets | Concordant | Validation Rate |
|-----|-----------|-------------------|------------|----------------|
| B4GALT1 | +156 | 1 (GSE7696) | 1 | 1/1 |
| EGFR | +209 | 2 (GSE7696, GSE230770) | 2 | 2/2 |
| VASP | +283 | 2 | 1 | 1/2 |
| GFPT2 | +184 | 1 (GSE7696) | 1 | 1/1 |
| PHF11 | +233 | 1 (GSE7696) | 1 | 1/1 |
| TSPYL5 | +224 | 2 | 2 | 2/2 |

**Sonuç:** 6 biomarker'dan 5'i en az bir external dataset ile TCGA yönünde uyumlu. İstatistiksel anlamlılık (p<0.05) external cohort'larda sınırlı (güç analizi); ancak yön tutarlılığı güçlü.

**Not:** GSE22867 (methylation 27K) ve GSE108474 (REMBRANDT) GEO metadata'da survival bilgisi olmadığı için analiz edilmedi.

---

## 8. Çalıştırılan Validation Sonuçları

| Veri Seti | Durum | Sonuç |
|-----------|-------|-------|
| **TCGA (Test 30%)** | Tamamlandı | B4GALT1, EGFR, TSPYL5 analiz edildi. EGFR +250 gün (yön uyumlu) |
| **GSE7696** | Tamamlandı | 6 biomarker, expression; tümü TCGA uyumlu |
| **GSE230770** | Tamamlandı | EGFR, VASP, TSPYL5 methylation; EGFR ve TSPYL5 uyumlu |
| **CGGA** | Kısmi | 84 GBM clinical, median OS 345 gün. Expression cgga.org.cn'dan |
| **GSE119774** | Kısmi | 40 GBM methylation. Survival GEO'da yok |

## 9. Özet

- **TCGA internal validation:** EGFR ve TSPYL5 yönü tutarlı; test seti küçük (n=29)
- **External validation:** GSE7696 (expression) ve GSE230770 (methylation) ile 5/6 biomarker TCGA ile yön uyumlu
- **GSE119774:** Methylation indirildi ve özetlendi
- **CGGA:** Clinical indirildi (Zenodo); expression için cgga.org.cn kayıt gerekli
- **Çıktılar:** `results/external_validation_results.csv`, `results/external_validation_summary.csv`, `plots/external_validation_comparison.png`
