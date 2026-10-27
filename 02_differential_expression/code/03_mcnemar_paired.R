# =========================================================
# 对应结果：差异分析第 3 步——条件 3（a+b+c > n/2 且 a ≤ n/2）的 eRNA/基因按表达状态（≥1 / <1）做配对 McNemar 检验
#           BH 校正，FDR<0.01 为差异（不限 logFC）
# 输入：output/<eRNA|gene>/<GSE>/Condition3_abc_gt_half_表达矩阵.csv（01 的输出）、../data/samples/<GSE>_samples.csv
# 输出：output/<eRNA|gene>/<GSE>/mcnemar_all_results.csv、mcnemar_sig_results_FDR0.01_only.csv、mcnemar_sig_results0.5.csv
# 运行：在模块目录下 Rscript code/03_mcnemar_paired.R
# =========================================================
library(readr)
library(dplyr)
library(tibble)

datasets <- c("GSE148355", "GSE77509", "GSE94660", "GSE124535", "GSE214846", "GSE144269")
for (feature in c("eRNA", "gene")) {
for (dataset in datasets) {
cat("\n==========", feature, dataset, "==========\n")
# 设置文件路径
input_dir <- file.path("output", feature, dataset)
output_dir <- file.path("output", feature, dataset)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# 读取表达矩阵（安全方式）
expr_matrix <- read_csv(
  file.path(input_dir, "Condition3_abc_gt_half_表达矩阵.csv"),
  show_col_types = FALSE
) %>%
  # 确保首列转为行名
  column_to_rownames(colnames(.)[1]) %>%
  as.data.frame()
# 验证数据结构
cat("\n表达矩阵维度:\n")
cat("基因数(行):", nrow(expr_matrix), "\n")
cat("样本数(列):", ncol(expr_matrix), "\n")
# 读取样本信息（统一使用read_csv）
sample_info <- read_csv(
  sprintf("../data/samples/%s_samples.csv", dataset),
  show_col_types = FALSE
)
# 自动对齐表达矩阵列顺序以匹配样本信息
expr_matrix <- expr_matrix[, sample_info$name]  # 修正：使用已读取的expr_matrix对齐列
# 添加分组和配对信息
n_pairs <- nrow(sample_info) / 2
sample_info$PatientID <- factor(rep(1:n_pairs, times = 2))  # 每个病人一对样本
sample_info$Group <- factor(sample_info$Treatment,
                            levels = c("nontumor", "tumor"),
                            labels = c("Normal", "Tumor"))
# 提取正常组和肿瘤组的列索引（新增关键代码）
normal_samples <- which(sample_info$Group == "Normal")  # 正常组样本列索引
tumor_samples <- which(sample_info$Group == "Tumor")    # 肿瘤组样本列索引
# 二值化表达矩阵 (1为阈值)
expr_binary <- ifelse(as.matrix(expr_matrix) < 1, "Low", "High") %>%
  `rownames<-`(rownames(expr_matrix))
# 执行McNemar检验
results <- lapply(rownames(expr_binary), function(gene) {
  tab <- table(
    factor(expr_binary[gene, normal_samples], levels = c("Low", "High")),
    factor(expr_binary[gene, tumor_samples], levels = c("Low", "High"))
  )
  if (all(dim(tab) == c(2,2))) {  # 确保列联表为2x2
    test <- mcnemar.test(tab)
    data.frame(
      Gene = gene,
      McNemar_ChiSq = test$statistic,
      p_value = test$p.value
    )
  } else {
    NULL  # 非2x2表格不检验
  }
}) %>%
  bind_rows() %>%  # 合并所有基因的检验结果
  mutate(
    adj_p_value = p.adjust(p_value, method = "BH"),  # BH法校正p值
    logFC = log2((rowMeans(expr_matrix[, tumor_samples]) + 1) /  # 计算logFC（加1避免除0）
                   (rowMeans(expr_matrix[, normal_samples]) + 1)),
    abs_logFC = abs(logFC)  # 绝对logFC值
  )
# 保存完整检验结果
write_csv(results, file.path(output_dir, "mcnemar_all_results.csv"))
# 筛选仅需FDR<0.01的基因（不限制logFC）
sig_genes_FDR_only <- results %>%
  filter(adj_p_value < 0.01)
write_csv(sig_genes_FDR_only, file.path(output_dir, "mcnemar_sig_results_FDR0.01_only.csv"))
cat("显著差异基因(FDR<0.01):", nrow(sig_genes_FDR_only), "\n")
# 不同logFC阈值的筛选（示例仅保留一个，可根据需要扩展）
sig_genes0.5 <- results %>%
  filter(adj_p_value < 0.01 & abs_logFC > 0.5)
write_csv(sig_genes0.5, file.path(output_dir, "mcnemar_sig_results0.5.csv"))
# 结果统计
cat("\n===== 分析结果 =====\n")
cat("总基因数:", nrow(expr_matrix), "\n")
cat("有效检验基因数:", nrow(results), "\n")
cat("显著差异基因(FDR<0.01 & |logFC|>0.5):", nrow(sig_genes0.5), "\n")
}
}
