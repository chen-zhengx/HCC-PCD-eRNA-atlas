# =========================================================
# 对应结果：Table 1（显著 eRNA–PCD 基因对 707 / 1,083 / 1,541；PCD 相关 eRNA 196 / 240 / 276；PCD 相关基因 20 / 40 / 50）
# 输入：output/correlation/<PCD>/<GSE>/output/<GSE>_All_Results.csv.gz（02_spearman_correlation.R 的输出）
# 输出：output/Table1_shared_pairs_summary.csv、output/shared_pairs_<PCD>.csv、output/PCD_related_eRNAs_<PCD>.csv
# 运行：在模块目录下 Rscript code/04_shared_pairs_table1.R
# =========================================================
library(data.table)

pcd_list <- c("Anoikis", "Apoptosis", "Ferroptosis")
selection_datasets <- c("GSE148355", "GSE77509", "GSE94660", "GSE214846", "GSE144269")
fdr_cut <- 1e-5

table1 <- data.frame()
for (pcd in pcd_list) {
  shared <- NULL
  for (dataset in selection_datasets) {
    res <- fread(file.path("output", "correlation", pcd, dataset, "output", paste0(dataset, "_All_Results.csv.gz")))
    sig <- res[FDR < fdr_cut, .(eRNA, Gene)]
    shared <- if (is.null(shared)) sig else merge(shared, sig, by = c("eRNA", "Gene"))
  }
  # 附上 6 个数据集各自的相关系数与 FDR
  for (dataset in c(selection_datasets, "GSE124535")) {
    res <- fread(file.path("output", "correlation", pcd, dataset, "output", paste0(dataset, "_All_Results.csv.gz")))
    setnames(res, c("Correlation", "FDR"), paste0(c("rho_", "FDR_"), dataset))
    shared <- merge(shared, res[, c("eRNA", "Gene", paste0(c("rho_", "FDR_"), dataset)), with = FALSE],
                    by = c("eRNA", "Gene"), all.x = TRUE)
  }
  fwrite(shared, file.path("output", paste0("shared_pairs_", pcd, ".csv")))
  fwrite(data.table(eRNA = unique(shared$eRNA)), file.path("output", paste0("PCD_related_eRNAs_", pcd, ".csv")))
  table1 <- rbind(table1, data.frame(PCD = pcd, Significant_pairs = nrow(shared),
                                     PCD_related_eRNAs = uniqueN(shared$eRNA), PCD_related_genes = uniqueN(shared$Gene)))
}
fwrite(table1, file.path("output", "Table1_shared_pairs_summary.csv"))
print(table1)
