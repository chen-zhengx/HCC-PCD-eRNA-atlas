# =========================================================
# 对应结果：Figure 2A（6 个数据集共有 559 个差异 eRNA）、Figure 2B（共有 649 个差异基因）
# 输入：output/<eRNA|gene>/<GSE>/Condition2_差异分析完整结果.csv、mcnemar_sig_results_FDR0.01_only.csv（02、03 的输出）
# 输出：output/common_DE_eRNAs.csv、output/common_DE_genes.csv、output/common_DE_summary.csv（每个数据集的差异数与累计交集数）
# 运行：在模块目录下 Rscript code/04_common_DE.R
# =========================================================
datasets <- c("GSE148355", "GSE77509", "GSE94660", "GSE124535", "GSE214846", "GSE144269")
logfc_cut <- c(eRNA = 0.5, gene = 1)

summary_all <- data.frame()
for (feature in c("eRNA", "gene")) {
  common <- NULL
  for (dataset in datasets) {
    dir_ds <- file.path("output", feature, dataset)
    # limma 差异（条件 2）
    limma_res <- read.csv(file.path(dir_ds, "Condition2_差异分析完整结果.csv"), row.names = 1, check.names = FALSE)
    limma_sig <- rownames(limma_res)[limma_res$adj.P.Val < 0.01 & abs(limma_res$logFC) > logfc_cut[[feature]]]
    # McNemar 差异（条件 3）
    mc_file <- file.path(dir_ds, "mcnemar_sig_results_FDR0.01_only.csv")
    mc_sig <- if (file.exists(mc_file)) read.csv(mc_file, check.names = FALSE)$Gene else character(0)
    # 合并起来的列表
    de_list <- union(limma_sig, mc_sig)
    write.csv(data.frame(ID = de_list), file.path(dir_ds, paste0("DE_", feature, "_merged_list.csv")), row.names = FALSE)
    # 依次取交集
    common <- if (is.null(common)) de_list else intersect(common, de_list)
    summary_all <- rbind(summary_all, data.frame(Feature = feature, Dataset = dataset,
                                                 limma = length(limma_sig), McNemar = length(mc_sig),
                                                 DE_total = length(de_list), Cumulative_intersection = length(common)))
    cat(feature, dataset, "差异数:", length(de_list), "| 累计交集:", length(common), "\n")
  }
  write.csv(data.frame(ID = common), file.path("output", paste0("common_DE_", feature, "s.csv")), row.names = FALSE)
}
write.csv(summary_all, file.path("output", "common_DE_summary.csv"), row.names = FALSE)
print(summary_all)
