# Nature Communications Validation Checklist

Manuscript: **Multi-Omics Integration Reveals Epigenetic Regulatory Cascades and Clinically Actionable Biomarkers in Glioblastoma**

## 1. Data Availability

- [ ] Data availability statement manuscript içinde yazıldı
- [ ] TCGA-GBM verisi erişim bilgisi (GDC Data Portal) belirtildi
- [ ] Minimum dataset tanımı: methylation, mRNA, protein, clinical
- [ ] Türetilmiş veriler (results/*.csv) depo edildi (Zenodo/Figshare)

## 2. Code Availability

- [ ] Code availability statement manuscript içinde
- [ ] Analiz kodları Zenodo veya Code Ocean'a yüklendi
- [ ] requirements.txt / environment.yml mevcut
- [ ] Kod DOI ile referans listesinde

## 3. Source Data (Figürler için)

- [ ] **Figure 1**: Cascade genleri, correlation değerleri, survival KM verileri
- [ ] **Figure 2**: Biomarker skorları, protein farkları, survival p-values
- [ ] **Figure 3**: Age-sex methylation correlation tabloları
- [ ] **Figure 4**: Liquid biopsy skorları, detectability vs prognostic
- [ ] **Figure 5**: Discordance kategorileri, correlation dağılımı
- [ ] Tüm tablolar (Table 1–5) için source data

## 4. Reproducibility

- [ ] `gbm_epigenetic_transcriptomic_proteomic_analysis.py` tekrar çalıştırıldı
- [ ] Novel analyses scriptleri (liquid biopsy, biomarker) doğrulandı
- [ ] Seed değerleri sabitlendi (np.random.seed)
- [ ] Python versiyonu ve paket versiyonları belgelendi

## 5. External Validation (Önerilen)

- [ ] Bağımsız GBM kohortu belirlendi (GEO/CGGA)
- [ ] 6 biomarker (B4GALT1, EGFR, VASP, GFPT2, PHF11, TSPYL5) survival doğrulandı
- [ ] Liquid biopsy signature external cohort'ta test edildi
- [ ] Cascade pattern external cohort'ta kontrol edildi

## 6. Statistical Validation

- [ ] Bootstrap 95% CI hesaplandı (kritik bulgular için)
- [ ] Multiple testing correction (FDR) kontrol edildi
- [ ] Cox model assumptions kontrol edildi

## 7. Reporting Summary

- [ ] Nature Communications reporting summary (PDF) dolduruldu
- [ ] Life sciences reporting summary template kullanıldı

## 8. Supplementary Materials

- [ ] Supplementary Tables 1–17 source data ile uyumlu
- [ ] Supplementary Methods detaylı
- [ ] Supplementary Figures source data

---

**Son Güncelleme:** Validation klasörü kuruldu.  
**Durum:** Başlangıç aşaması – maddeler doldurulacak.
