# =========================================================
# 对应结果：编码基因 TPM 矩阵（../data/expression/gene_TPM/<GSE>_gene_TPM.csv.gz）
# 输入：GEO 下载（需联网）、../data/annotation/protein_coding_genes_19416.csv、../data/samples/Table_S1_sample_accessions.csv
# 输出：output/<GSE>_gene_TPM.csv
# 运行：在模块目录下 Rscript code/02_gene_tpm_from_geo.R
# =========================================================
datasets <- c("GSE148355", "GSE77509", "GSE94660", "GSE124535", "GSE214846", "GSE144269")
dir.create("output", showWarnings = FALSE)
genes <- read.csv("../data/annotation/protein_coding_genes_19416.csv")
s1 <- read.csv("../data/samples/Table_S1_sample_accessions.csv", stringsAsFactors = FALSE)

for (dataset in datasets) {
  f <- file.path("output", paste0(dataset, "_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"))
  if (!file.exists(f)) {
    url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=", dataset,
                  "&format=file&file=", dataset, "_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
    download.file(url, f, mode = "wb")
  }
  geo <- read.delim(gzfile(f), row.names = 1, check.names = FALSE)
  s <- s1[s1$dataset == dataset, ]
  s <- s[order(s$column_order), ]
  tpm <- geo[as.character(genes$GeneID), s$GSM]
  dimnames(tpm) <- list(genes$Symbol, s$matrix_column)
  write.csv(data.frame(GeneSymbol = rownames(tpm), tpm, check.names = FALSE),
            file.path("output", paste0(dataset, "_gene_TPM.csv")), row.names = FALSE)
  cat(dataset, ": ", nrow(tpm), " genes x ", ncol(tpm), " samples\n", sep = "")
}
