# =========================================================
# 对应结果：差异分析第 2 步——条件 2（a > n/2）的 eRNA/基因用 limma 配对设计做差异分析
#           差异 eRNA：FDR<0.01 且 |logFC|>0.5；差异基因：FDR<0.01 且 |logFC|>1（结果按 0.5–1.0 各阈值分表保存）
# 输入：output/<eRNA|gene>/<GSE>/Condition2_a_gt_half_表达矩阵.csv（01 的输出）、../data/samples/<GSE>_samples.csv
# 输出：output/<eRNA|gene>/<GSE>/Condition2_差异分析完整结果.csv、Condition2_差异分析分阈值结果.xlsx、QC_plots/
# 运行：在模块目录下 Rscript code/02_limma_paired.R
# =========================================================
# ===================== Step 0: 加载包与路径设定 =====================
library(limma)
library(openxlsx)
library(ggplot2)

datasets <- c("GSE148355", "GSE77509", "GSE94660", "GSE124535", "GSE214846", "GSE144269")
for (feature in c("eRNA", "gene")) {
for (dataset in datasets) {
cat("\n==========", feature, dataset, "==========\n")
# 输入输出路径
input_dir <- file.path("output", feature, dataset)
output_dir <- file.path("output", feature, dataset)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# ===================== Step 1: 数据加载与对齐 =====================
# 读取表达矩阵
expr_cond2 <- read.csv(file.path(input_dir, "Condition2_a_gt_half_表达矩阵.csv"), row.names = 1, check.names = FALSE)
expr_cond2 <- as.matrix(expr_cond2)
expr_cond2_log <- log2(expr_cond2 + 1)
# 读取样本信息
sample_info <- read.csv(sprintf("../data/samples/%s_samples.csv", dataset), stringsAsFactors = FALSE)
# 自动对齐表达矩阵列顺序以匹配样本信息
expr_cond2_log <- expr_cond2_log[, sample_info$name]
# 添加分组和配对信息
n_pairs <- nrow(sample_info) / 2
sample_info$PatientID <- factor(rep(1:n_pairs, times = 2))  # 每个病人一对样本
sample_info$Group <- factor(sample_info$Treatment,
                            levels = c("nontumor", "tumor"),
                            labels = c("Normal", "Tumor"))
# ===================== Step 2: QC 直方图绘制 =====================
plot_dir <- file.path(output_dir, "QC_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir)
hist_plot <- ggplot(data.frame(Value = as.vector(expr_cond2_log)), aes(x = Value)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "对数转换后表达值整体分布", x = "log2(RPM+1)", y = "频数") +
  theme_bw()
ggsave(file.path(plot_dir, "对数转换直方图.png"), hist_plot, width = 8, height = 6)
# ===================== Step 3: 差异分析（配对设计） =====================
design <- model.matrix(~ PatientID + Group, data = sample_info)
fit <- lmFit(expr_cond2_log, design)
fit <- eBayes(fit)
degs_cond2 <- topTable(fit, coef = "GroupTumor", number = Inf, adjust.method = "BH")
# 保存完整差异分析结果
write.csv(degs_cond2, file.path(output_dir, "Condition2_差异分析完整结果.csv"))
# ===================== Step 4: 全阈值统计汇总 =====================
fc_thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
count_sig_genes <- function(fc) {
  sig <- degs_cond2[degs_cond2$adj.P.Val < 0.01 & abs(degs_cond2$logFC) > fc, ]
  list(
    Total = nrow(sig),
    Up = sum(sig$logFC > 0),
    Down = sum(sig$logFC < 0)
  )
}
stats_list <- lapply(fc_thresholds, count_sig_genes)
stats_df <- data.frame(
  Threshold = paste0("|logFC|>", fc_thresholds),
  Total = sapply(stats_list, function(x) x$Total),
  Up = sapply(stats_list, function(x) x$Up),
  Down = sapply(stats_list, function(x) x$Down)
)
stats_df <- rbind(
  data.frame(
    Threshold = "FDR<0.01 (任意FC)",
    Total = nrow(degs_cond2[degs_cond2$adj.P.Val < 0.01, ]),
    Up = NA,
    Down = NA
  ),
  stats_df
)
# 控制台输出
cat("\n========== Condition2 差异分析报告 ==========\n")
cat("总检测基因数:", nrow(degs_cond2), "\n")
cat("----------------------------------------\n")
print(stats_df)
cat("\n* 注：上调基因 = 肿瘤表达 > 正常；下调基因 = 肿瘤表达 < 正常\n")
# ===================== Step 5: 保存分层结果到 Excel =====================
result_list <- list()
result_list[["汇总统计"]] <- stats_df
for (fc in fc_thresholds) {
  sig_genes <- degs_cond2[degs_cond2$adj.P.Val < 0.01 & abs(degs_cond2$logFC) > fc, ]
  sheet_name <- paste0("FC>", fc, "_eRNA")
  result_list[[sheet_name]] <- sig_genes
}
write.xlsx(
  result_list,
  file = file.path(output_dir, "Condition2_差异分析分阈值结果.xlsx"),
  rowNames = TRUE
)
cat("\n✅ 全部结果保存完成：", output_dir, "\n")
}
}
