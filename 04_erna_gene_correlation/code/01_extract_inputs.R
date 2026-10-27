# =========================================================
# 对应结果：相关性分析的输入矩阵（559 个共有差异 eRNA 的 RPM；PCD 基因集与 649 个共有差异基因交集的 TPM）
# 输入：../data/expression/eRNA_RPM/<GSE>_eRNA_RPM.csv.gz、../data/expression/gene_TPM/<GSE>_gene_TPM.csv.gz
#       ../data/annotation/protein_coding_genes_19416.csv
#       ../02_differential_expression/input/common_DE_eRNAs_559.csv
#       ../03_pcd_enrichment/input/pcd_genesets_TableS2.csv、../03_pcd_enrichment/input/common_DE_genes_649.csv
# 输出：output/correlation/<PCD>/<GSE>/input/<GSE>_DE_eRNA0.5_matrix.csv、<GSE>_DE_<PCD>_genes.csv
# 运行：在模块目录下 Rscript code/01_extract_inputs.R
# =========================================================
datasets <- c("GSE148355", "GSE77509", "GSE94660", "GSE124535", "GSE214846", "GSE144269")
pcd_types <- c(Anoikis = "Anoikis", Apoptosis = "Apoptosis", Ferroptosis = "Ferroptosis")   # Table S2 中的类别名

de_ernas <- read.csv("../02_differential_expression/input/common_DE_eRNAs_559.csv")$eRNA
de_genes <- read.csv("../03_pcd_enrichment/input/common_DE_genes_649.csv")$gene
genesets <- read.csv("../03_pcd_enrichment/input/pcd_genesets_TableS2.csv", stringsAsFactors = FALSE)
gene_ids <- read.csv("../data/annotation/protein_coding_genes_19416.csv")

for (dataset in datasets) {
  rpm <- read.csv(sprintf("../data/expression/eRNA_RPM/%s_eRNA_RPM.csv.gz", dataset), row.names = 1, check.names = FALSE)
  tpm <- read.csv(sprintf("../data/expression/gene_TPM/%s_gene_TPM.csv.gz", dataset), row.names = 1, check.names = FALSE)
  # GSE214846 的 TPM 矩阵行名为 NCBI GeneID，换成基因符号
  if (dataset == "GSE214846") rownames(tpm) <- gene_ids$Symbol[match(rownames(tpm), gene_ids$GeneID)]
  samples <- read.csv(sprintf("../data/samples/%s_samples.csv", dataset))$name
  samples <- colnames(rpm)[colnames(rpm) %in% samples]        # 保持表达矩阵的列顺序
  erna_mat <- rpm[rownames(rpm) %in% de_ernas, samples]
  for (pcd in names(pcd_types)) {
    pcd_genes <- intersect(genesets$Symbol[genesets$Type == pcd_types[[pcd]]], de_genes)
    gene_mat <- tpm[rownames(tpm) %in% pcd_genes, samples]
    out_dir <- file.path("output", "correlation", pcd, dataset, "input")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(erna_mat, file.path(out_dir, paste0(dataset, "_DE_eRNA0.5_matrix.csv")))
    write.csv(data.frame(Symbol = rownames(gene_mat), gene_mat, check.names = FALSE),
              file.path(out_dir, paste0(dataset, "_DE_", pcd, "_genes.csv")), row.names = FALSE)
    cat(dataset, pcd, ": eRNA", nrow(erna_mat), "x", ncol(erna_mat), "| genes", nrow(gene_mat), "\n")
  }
}
