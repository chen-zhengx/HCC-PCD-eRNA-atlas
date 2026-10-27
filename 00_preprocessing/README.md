# 00_preprocessing: read processing and expression quantification

This module contains the code that converts the raw RNA-seq reads into the expression matrices in
`../data/expression/`: read trimming with fastp, alignment with HISAT2 and SAMtools, counting of reads in eRNA regions
with bedtools multicov (`batch_scripts/`), normalization to reads per million mapped reads (`code/01_erna_rpm.R`),
and the extraction of protein-coding gene TPM values from the NCBI-generated RNA-seq quantification distributed by GEO
(`code/02_gene_tpm_from_geo.R`).
