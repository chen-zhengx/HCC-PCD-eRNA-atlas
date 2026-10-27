# An enhancer RNA atlas of programmed cell death in hepatocellular carcinoma

Analysis code and processed data for the manuscript "An enhancer RNA atlas of programmed cell death in hepatocellular carcinoma" (TGH-2026-0192). 

## Layout

```
data/                         shared processed data (used by several modules)
00_preprocessing/             read trimming, alignment, eRNA quantification, RPM; gene TPM from GEO
01_samples/                   sample identifiers (Table S1)
02_differential_expression/   DE eRNAs and genes in six datasets (Figure 2)
03_pcd_enrichment/            PCD gene sets and hypergeometric screening (Tables S2, S3)
04_erna_gene_correlation/     eRNA–PCD gene correlations (Table 1, Table S6), tissue-status sensitivity analysis
05_prognostic_model/          univariate Cox + LASSO-Cox models, KM, time-dependent ROC, internal validation
06_clinical/                  clinical characteristics and Cox regression
07_molecular_subtypes/        consensus clustering and correspondence with risk groups
08_risk_group_function/       DE genes between risk groups, GO/KEGG enrichment
09_immune/                    ssGSEA immune cell scores, immune-rank subtypes
10_immunotherapy/             tigeR immunotherapy-response signatures
11_drug/                      oncoPredict drug-sensitivity prediction and paired tests
12_network/                   multilayer co-expression association networks
13_eRNA_QTL/                  eRNA-QTL mapping (MatrixEQTL)
```

## Data

Six paired HCC RNA-seq datasets from GEO (458 samples, 229 tumor/non-tumor pairs):

| Dataset | BioProject | Pairs |
|---|---|---|
| GSE148355 | PRJNA624047 | 18 |
| GSE77509 | PRJNA310654 | 20 |
| GSE94660 | PRJNA371753 | 21 |
| GSE124535 | PRJNA512451 | 35 |
| GSE214846 | PRJNA887382 | 65 |
| GSE144269 | PRJNA603277 | 70 |


## Software

R 4.5, Python 3.11, and for read processing fastp 1.4.0, HISAT2 2.2.3, SAMtools 1.24, bedtools 2.31.1 and BCFtools 1.24.
