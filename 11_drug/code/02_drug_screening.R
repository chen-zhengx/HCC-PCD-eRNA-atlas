# =========================================================
# 对应结果：19 种肿瘤选择性候选药物（Table S11 的药物清单；正文 "19 of 198 drugs"）：对 70 对样本的预测 IC50 做单侧配对
#           Wilcoxon 符号秩检验（肿瘤 < 癌旁），198 种药物 BH 校正，候选 = FDR < 0.01 且平均差值（癌旁 − 肿瘤）> 0
# 输入：input/drug_predictions_GDSC2_140samples.csv
# 输出：output/screening/AllDrugs_StatisticalScreening.csv（198 种药物）、SensitiveDrugs_StatisticalScreening.csv（19 种）
#       DetailedReport_SensitiveDrugs.csv、Sensitivity_VolcanoPlot.png、Boxplot_<药物>.png（后两类论文未使用）
# 运行：在模块目录下 Rscript code/02_drug_screening.R
# =========================================================

# 设置输出目录
output_dir <- "output/screening"

# 创建输出目录（如果不存在）
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 读取数据
data <- read.csv("input/drug_predictions_GDSC2_140samples.csv", header = TRUE, stringsAsFactors = FALSE)

# 分离肿瘤样本和正常样本
tumor_samples <- data[data$Treatment == "tumor", ]
normal_samples <- data[data$Treatment == "nontumor", ]

# 按照PatientID匹配配对样本
matched_data <- merge(tumor_samples, normal_samples, by = "PatientID", suffixes = c("_tumor", "_normal"))

# 获取药物列名（从第4列开始）
drug_columns <- colnames(data)[4:ncol(data)]

# 初始化结果数据框
results_df <- data.frame(
  Drug = drug_columns,
  Mean_Tumor = numeric(length(drug_columns)),
  Mean_Normal = numeric(length(drug_columns)),
  Median_Tumor = numeric(length(drug_columns)),
  Median_Normal = numeric(length(drug_columns)),
  Mean_Difference = numeric(length(drug_columns)),
  Wilcoxon_Pvalue = numeric(length(drug_columns)),
  Wilcoxon_FDR = numeric(length(drug_columns)),
  SensitiveDrug = logical(length(drug_columns)),
  stringsAsFactors = FALSE
)

# 对每种药物进行分析
for (i in 1:length(drug_columns)) {
  drug <- drug_columns[i]
  
  # 构建配对样本的药物列名
  tumor_col <- paste0(drug, "_tumor")
  normal_col <- paste0(drug, "_normal")
  
  # 提取肿瘤和正常样本的药物预测分数
  tumor_scores <- matched_data[[tumor_col]]
  normal_scores <- matched_data[[normal_col]]
  
  # 计算描述性统计量
  mean_tumor <- mean(tumor_scores, na.rm = TRUE)
  mean_normal <- mean(normal_scores, na.rm = TRUE)
  median_tumor <- median(tumor_scores, na.rm = TRUE)
  median_normal <- median(normal_scores, na.rm = TRUE)
  mean_difference <- mean_normal - mean_tumor
  
  # 执行Wilcoxon符号秩检验
  wilcoxon_test <- wilcox.test(tumor_scores, normal_scores, 
                               paired = TRUE, alternative = "less")
  wilcoxon_pvalue <- wilcoxon_test$p.value
  
  # 存储结果
  results_df$Mean_Tumor[i] <- round(mean_tumor, 4)
  results_df$Mean_Normal[i] <- round(mean_normal, 4)
  results_df$Median_Tumor[i] <- round(median_tumor, 4)
  results_df$Median_Normal[i] <- round(median_normal, 4)
  results_df$Mean_Difference[i] <- round(mean_difference, 4)
  results_df$Wilcoxon_Pvalue[i] <- wilcoxon_pvalue
}

# 计算FDR校正
results_df$Wilcoxon_FDR <- p.adjust(results_df$Wilcoxon_Pvalue, method = "fdr")

# 标记敏感性药物（Wilcoxon检验FDR < 0.05且平均差异为正）
results_df$SensitiveDrug <- results_df$Wilcoxon_FDR < 0.01 & results_df$Mean_Difference > 0

# 按FDR值从低到高排序
results_df <- results_df[order(results_df$Wilcoxon_FDR, -results_df$Mean_Difference), ]

# 输出结果
print("统计检验筛选分析结果：")
print(results_df)

# 筛选出敏感性药物
sensitive_drugs <- results_df[results_df$SensitiveDrug, ]

print("敏感性药物列表（通过统计检验筛选）：")
print(sensitive_drugs)

# 将结果保存到指定路径
write.csv(results_df, 
          file.path(output_dir, "AllDrugs_StatisticalScreening.csv"), 
          row.names = FALSE)

write.csv(sensitive_drugs, 
          file.path(output_dir, "SensitiveDrugs_StatisticalScreening.csv"), 
          row.names = FALSE)

# 为敏感性药物生成详细报告
if (nrow(sensitive_drugs) > 0) {
  detailed_report <- data.frame()
  
  for (drug in sensitive_drugs$Drug) {
    tumor_col <- paste0(drug, "_tumor")
    normal_col <- paste0(drug, "_normal")
    
    # 提取数据
    tumor_scores <- matched_data[[tumor_col]]
    normal_scores <- matched_data[[normal_col]]
    
    # 计算更多统计量
    drug_stats <- data.frame(
      Drug = drug,
      Mean_Tumor = mean(tumor_scores, na.rm = TRUE),
      Mean_Normal = mean(normal_scores, na.rm = TRUE),
      SD_Tumor = sd(tumor_scores, na.rm = TRUE),
      SD_Normal = sd(normal_scores, na.rm = TRUE),
      Median_Tumor = median(tumor_scores, na.rm = TRUE),
      Median_Normal = median(normal_scores, na.rm = TRUE),
      Mean_Difference = mean(normal_scores - tumor_scores, na.rm = TRUE),
      Wilcoxon_Pvalue = results_df$Wilcoxon_Pvalue[results_df$Drug == drug],
      Wilcoxon_FDR = results_df$Wilcoxon_FDR[results_df$Drug == drug],
      EffectSize = mean(normal_scores - tumor_scores, na.rm = TRUE) / 
        sqrt((sd(tumor_scores, na.rm = TRUE)^2 + sd(normal_scores, na.rm = TRUE)^2)/2)
    )
    
    detailed_report <- rbind(detailed_report, drug_stats)
  }
  
  write.csv(detailed_report, 
            file.path(output_dir, "DetailedReport_SensitiveDrugs.csv"), 
            row.names = FALSE)
}

# 输出完成信息
cat("分析完成！结果已保存到：", output_dir, "\n")
cat("1. AllDrugs_StatisticalScreening.csv - 所有药物的统计检验结果\n")
cat("2. SensitiveDrugs_StatisticalScreening.csv - 敏感性药物列表\n")
cat("3. DetailedReport_SensitiveDrugs.csv - 敏感性药物详细报告\n")

# 绘制敏感性药物分布图（可选）
if (nrow(sensitive_drugs) > 0) {
  library(ggplot2)
  
  # 创建FDR值分布图
  p <- ggplot(results_df, aes(x = -log10(Wilcoxon_FDR), y = Mean_Difference)) +
    geom_point(aes(color = SensitiveDrug), alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    scale_color_manual(values = c("grey", "red")) +
    labs(title = "药物敏感性分析 - 统计检验筛选",
         x = "-log10(FDR)",
         y = "平均差异 (正常 - 肿瘤)",
         color = "敏感性药物") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "Sensitivity_VolcanoPlot.png"), p, width = 8, height = 6)
  
  # 为每个敏感性药物创建箱线图
  for (drug in sensitive_drugs$Drug) {
    tumor_col <- paste0(drug, "_tumor")
    normal_col <- paste0(drug, "_normal")
    
    # 准备数据
    plot_data <- data.frame(
      Group = rep(c("Tumor", "Normal"), each = nrow(matched_data)),
      Score = c(matched_data[[tumor_col]], matched_data[[normal_col]])
    )
    
    p_box <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
      geom_boxplot() +
      labs(title = paste("药物敏感性:", drug),
           subtitle = paste("FDR =", format(results_df$Wilcoxon_FDR[results_df$Drug == drug], scientific = TRUE, digits = 3)),
           y = "预测敏感性分数") +
      theme_minimal() +
      scale_fill_manual(values = c("Tumor" = "lightcoral", "Normal" = "lightblue"))
    
    ggsave(file.path(output_dir, paste0("Boxplot_", drug, ".png")), p_box, width = 6, height = 5)
  }
}

