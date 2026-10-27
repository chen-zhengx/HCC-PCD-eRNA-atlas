# =========================================================
# 对应结果：Figure 8、Figure S8、Figure S9 与 Supplementary Table 12 的多层共表达关联网络——风险 eRNA 与各层节点之间的边：
#           风险 eRNA 与 ±1 Mb 候选靶基因、PCD 相关基因、免疫细胞（ssGSEA 得分）、19 种候选药物（oncoPredict 预测 IC50）
#           免疫治疗相关基因的 Spearman 相关（70 个肿瘤样本），每个 eRNA × 每层分别做 BH 校正，FDR < 0.01 为显著边
# 输入：output/inputs/（02_prepare_network_inputs.R 的输出）、input/immune_cell_ssGSEA_70tumors.csv
#       input/drug_IC50_19drugs_70tumors.csv、input/risk_eRNAs_by_PCD.csv
# 输出：output/correlations/<PCD>/<eRNA>/eRNA_vs_target_genes|eRNA_vs_PCD_genes|eRNA_vs_immune_cells|eRNA_vs_drugs|
#       eRNA_vs_immunotherapy_genes/：*_Spearman_Full_Results.csv 与 FDR_0.05 至 FDR_0.00001 各阈值的显著结果
#       （网络使用 FDR_0.01；eRNA 与候选靶基因、免疫细胞、药物、免疫治疗基因的相关与 PCD 类型无关，三个网络中结果相同）
# 运行：在模块目录下 Rscript code/03_eRNA_correlations.R
# =========================================================

risk <- read.csv("input/risk_eRNAs_by_PCD.csv", stringsAsFactors = FALSE)
for (r in seq_len(nrow(risk))) {
  pcd     <- risk$PCD[r]            # anoikis / apoptosis / ferroptosis
  pcd_cn  <- risk$PCD_cn[r]         # 文件名中使用的 PCD 名称：失巢死亡 / 凋亡 / 铁死亡
  erna    <- risk$eRNA[r]
  erna_id <- gsub(":", "", erna)    # 文件名中去掉冒号
  run_dir <- file.path("output", "correlations", pcd, erna_id)
  pcd_matrix_file <- list.files(file.path("output", "inputs", "PCD_genes"),
                                pattern = paste0("^", pcd_cn, "[0-9]+个高度相关基因表达矩阵\\.csv$"), full.names = TRUE)
  cat("\n==========", pcd, erna, "==========\n")

  # ---------------------------------------------------------
  # 1. eRNA 与 ±1 Mb 候选靶基因
  # ---------------------------------------------------------
  # 1. 加载必要的包
  library(data.table)
  library(dplyr)

  # 2. 设置输入和输出路径
  input_dir <- "output/inputs"
  output_dir <- file.path(run_dir, "eRNA_vs_target_genes")

  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 3. 读取数据
  # 读取eRNA表达矩阵 (只有一个eRNA)
  erna_file <- file.path(input_dir, "eRNA_expression", paste0(erna_id, "表达.csv"))
  erna_data <- fread(erna_file, header = TRUE)
  # 转换为向量格式，方便后续计算
  erna_vector <- as.numeric(erna_data[1, -1]) # 第一行，排除第一列（名称列）
  erna_name <- erna_data[[1,1]] # 获取eRNA的名字

  # 读取靶基因表达矩阵
  gene_file <- file.path(input_dir, "target_genes", paste0(erna_id, "_1MB靶基因矩阵.csv"))
  gene_data <- fread(gene_file, header = TRUE)

  # 4. 数据预处理与检查
  # 提取基因名
  gene_names <- gene_data[[1]]
  # 提取表达矩阵（转换为数值矩阵）
  gene_matrix <- as.matrix(gene_data[, -1])

  # 检查样本是否匹配（非常重要！）
  if (ncol(gene_matrix) != length(erna_vector)) {
    stop("错误：eRNA和基因矩阵的样本数量不匹配！")
  } else {
    cat("数据检查通过：样本数量一致，共", ncol(gene_matrix), "个样本。\n")
  }

  # 5. 计算斯皮尔曼相关性
  # 初始化结果数据框
  results <- data.frame(
    eRNA = character(),
    Gene = character(),
    Correlation = numeric(),
    P.Value = numeric(),
    stringsAsFactors = FALSE
  )

  cat("正在计算 eRNA:", erna_name, "与", length(gene_names), "个靶基因的斯皮尔曼相关性...\n")

  # 循环计算每个基因与eRNA的相关性
  for (i in 1:length(gene_names)) {
    gene_exp <- as.numeric(gene_matrix[i, ]) # 获取第i个基因的表达向量

    # 进行斯皮尔曼相关性检验
    cor_test_result <- cor.test(erna_vector, gene_exp, method = "spearman", exact = FALSE)

    # 将结果添加到数据框
    results[i, ] <- list(
      eRNA = erna_name,
      Gene = gene_names[i],
      Correlation = cor_test_result$estimate,
      P.Value = cor_test_result$p.value
    )
  }

  # 6. 多重检验校正 (FDR)
  results$FDR <- p.adjust(results$P.Value, method = "BH")

  # 7. 按相关性系数绝对值排序
  results_sorted <- results[order(-abs(results$Correlation)), ]

  # 8. 定义筛选条件（基于FDR）
  filtered_results <- list(
    FDR_0.05 = results_sorted %>% filter(FDR < 0.05) %>% arrange(desc(abs(Correlation))),
    FDR_0.01 = results_sorted %>% filter(FDR < 0.01) %>% arrange(desc(abs(Correlation))),
    FDR_0.001 = results_sorted %>% filter(FDR < 0.001) %>% arrange(desc(abs(Correlation))),
    FDR_0.0001 = results_sorted %>% filter(FDR < 0.0001) %>% arrange(desc(abs(Correlation))),
    FDR_0.00001 = results_sorted %>% filter(FDR < 0.00001) %>% arrange(desc(abs(Correlation)))
  )

  # 9. 保存结果
  # 保存完整结果
  output_file_full <- file.path(output_dir, paste0(erna_id, "_vs_1MB_靶基因_Spearman_Full_Results.csv"))
  fwrite(results_sorted, output_file_full)
  cat("完整结果已保存至:", output_file_full, "\n")

  # 保存不同FDR阈值的结果
  for(threshold in names(filtered_results)) {
    output_file <- file.path(output_dir, paste0(erna_id, "_vs_1MB_靶基因_Spearman_", threshold, ".csv"))
    fwrite(filtered_results[[threshold]], output_file)
    cat("已保存", threshold, "结果:", output_file, "\n")
  }

  # 10. 创建统计摘要并输出
  summary_stats <- data.frame(
    FDR_Threshold = names(filtered_results),
    Significant_Genes = sapply(filtered_results, nrow),
    stringsAsFactors = FALSE
  )

  # 打印到控制台
  cat("\n===== 分析完成 =====\n")
  cat("eRNA:", erna_name, "\n")
  cat("分析的靶基因数量:", nrow(results), "\n\n")

  cat("不同FDR阈值下的显著基因数量:\n")
  for(i in 1:nrow(summary_stats)) {
    cat("  ", summary_stats$FDR_Threshold[i], ":", summary_stats$Significant_Genes[i], "\n")
  }

  # 找最强正相关和负相关基因
  if (nrow(results_sorted) > 0) {
    cat("\n最强相关性基因:\n")
    cat("  最强正相关基因:", results_sorted$Gene[1], 
        "(r =", round(results_sorted$Correlation[1], 4), 
        ", FDR =", format(results_sorted$FDR[1], scientific = TRUE, digits = 3), ")\n")

    min_cor_index <- which.min(results_sorted$Correlation)
    cat("  最强负相关基因:", results_sorted$Gene[min_cor_index], 
        "(r =", round(results_sorted$Correlation[min_cor_index], 4), 
        ", FDR =", format(results_sorted$FDR[min_cor_index], scientific = TRUE, digits = 3), ")\n")
  }

  # 11. 保存统计结果到文件
  stats_file <- file.path(output_dir, "correlation_statistics.txt")
  sink(stats_file)
  cat("eRNA与靶基因相关性分析统计摘要\n")
  cat("================================\n")
  cat("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("eRNA名称:", erna_name, "\n")
  cat("分析的靶基因总数:", nrow(results), "\n\n")

  cat("不同FDR阈值下的显著基因数量:\n")
  for(i in 1:nrow(summary_stats)) {
    cat(summary_stats$FDR_Threshold[i], ":", summary_stats$Significant_Genes[i], "\n")
  }

  cat("\n最强相关性基因:\n")
  if (nrow(results_sorted) > 0) {
    cat("最强正相关基因:", results_sorted$Gene[1], 
        "(r =", round(results_sorted$Correlation[1], 4), 
        ", FDR =", format(results_sorted$FDR[1], scientific = TRUE, digits = 3), ")\n")

    min_cor_index <- which.min(results_sorted$Correlation)
    cat("最强负相关基因:", results_sorted$Gene[min_cor_index], 
        "(r =", round(results_sorted$Correlation[min_cor_index], 4), 
        ", FDR =", format(results_sorted$FDR[min_cor_index], scientific = TRUE, digits = 3), ")\n")
  }

  # 添加数据质量检查信息
  cat("\n数据质量检查:\n")
  cat("样本数量:", ncol(gene_matrix), "\n")
  cat("基因数量:", nrow(results), "\n")
  cat("有效相关性计算:", sum(!is.na(results$Correlation)), "/", nrow(results), "\n")

  sink()

  cat("\n统计摘要已保存至:", stats_file, "\n")

  # ---------------------------------------------------------
  # 2. eRNA 与 PCD 相关基因
  # 凋亡、铁死亡版本只有输入文件名和输出文件名不同
  # ---------------------------------------------------------
  # 1. 加载必要的包
  library(data.table)
  library(dplyr)

  # 2. 设置输入和输出路径

  # 2. 设置输入和输出路径
  input_dir <- "output/inputs"
  output_dir <- file.path(run_dir, "eRNA_vs_PCD_genes")

  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("创建输出目录:", output_dir, "\n")
  }

  # 3. 读取数据
  cat("正在读取数据...\n")

  # 读取eRNA表达矩阵
  erna_file <- file.path(input_dir, "eRNA_expression", paste0(erna_id, "表达.csv"))
  if (!file.exists(erna_file)) {
    stop("错误: 找不到eRNA表达文件: ", erna_file)
  }
  erna_data <- fread(erna_file, header = TRUE)

  # 转换为向量格式
  erna_vector <- as.numeric(unlist(erna_data[1, -1]))  # 第一行，排除第一列
  erna_name <- erna_data[[1,1]]  # 获取eRNA的名字

  # 读取失巢死亡基因表达矩阵
  gene_file <- pcd_matrix_file  # 02 脚本提取的 PCD 相关基因（失巢凋亡 20 个、凋亡 40 个、铁死亡 50 个）
  if (!file.exists(gene_file)) {
    stop("错误: 找不到失巢死亡基因表达文件: ", gene_file)
  }
  gene_data <- fread(gene_file, header = TRUE)

  # 4. 数据预处理与检查
  # 提取基因名
  gene_names <- gene_data[[1]]

  # 提取表达矩阵并转换为数值矩阵
  gene_matrix <- as.matrix(gene_data[, -1])
  gene_matrix_numeric <- apply(gene_matrix, 2, as.numeric)

  # 检查样本是否匹配
  if (ncol(gene_matrix_numeric) != length(erna_vector)) {
    stop(sprintf("错误：eRNA和基因矩阵的样本数量不匹配！eRNA: %d个样本, 基因: %d个样本", 
                 length(erna_vector), ncol(gene_matrix_numeric)))
  } else {
    cat(sprintf("数据检查通过：样本数量一致，共 %d 个样本。\n", ncol(gene_matrix_numeric)))
  }

  # 5. 计算斯皮尔曼相关性
  cat(sprintf("正在计算 eRNA %s 与 %d 个失巢死亡基因的斯皮尔曼相关性...\n", 
              erna_name, length(gene_names)))

  # 初始化结果数据框
  results <- data.frame(
    eRNA = character(),
    Gene = character(),
    Correlation = numeric(),
    P.Value = numeric(),
    stringsAsFactors = FALSE
  )

  # 循环计算每个基因与eRNA的相关性
  for (i in 1:length(gene_names)) {
    gene_exp <- gene_matrix_numeric[i, ]  # 获取第i个基因的表达向量

    # 检查并处理缺失值
    valid_data <- !(is.na(erna_vector) | is.na(gene_exp))
    n_valid <- sum(valid_data)

    if (n_valid < 3) {
      # 样本太少，无法计算相关性
      results[i, ] <- list(
        eRNA = erna_name,
        Gene = gene_names[i],
        Correlation = NA,
        P.Value = NA
      )
      next
    }

    # 进行斯皮尔曼相关性检验
    cor_test_result <- tryCatch({
      cor.test(erna_vector[valid_data], gene_exp[valid_data], 
               method = "spearman", exact = FALSE)
    }, error = function(e) {
      return(NULL)
    })

    if (is.null(cor_test_result)) {
      # 计算失败
      results[i, ] <- list(
        eRNA = erna_name,
        Gene = gene_names[i],
        Correlation = NA,
        P.Value = NA
      )
    } else {
      # 计算成功
      results[i, ] <- list(
        eRNA = erna_name,
        Gene = gene_names[i],
        Correlation = cor_test_result$estimate,
        P.Value = cor_test_result$p.value
      )
    }
  }

  # 移除无法计算相关性的基因
  results_valid <- results[!is.na(results$Correlation), ]
  cat(sprintf("成功计算 %d/%d 个基因的相关性\n", 
              nrow(results_valid), nrow(results)))

  # 6. 多重检验校正 (FDR)
  results_valid$FDR <- p.adjust(results_valid$P.Value, method = "BH")

  # 7. 按相关性系数绝对值排序
  results_sorted <- results_valid[order(-abs(results_valid$Correlation)), ]

  # 8. 定义筛选条件（基于FDR）
  filtered_results <- list(
    FDR_0.05 = results_sorted %>% filter(FDR < 0.05) %>% arrange(desc(abs(Correlation))),
    FDR_0.01 = results_sorted %>% filter(FDR < 0.01) %>% arrange(desc(abs(Correlation))),
    FDR_0.001 = results_sorted %>% filter(FDR < 0.001) %>% arrange(desc(abs(Correlation))),
    FDR_0.0001 = results_sorted %>% filter(FDR < 0.0001) %>% arrange(desc(abs(Correlation))),
    FDR_0.00001 = results_sorted %>% filter(FDR < 0.00001) %>% arrange(desc(abs(Correlation)))
  )

  # 9. 保存结果
  # 保存完整结果
  output_file_full <- file.path(output_dir, paste0(erna_id, "_vs_", pcd_cn, "基因_Spearman_Full_Results.csv"))
  fwrite(results_sorted, output_file_full)
  cat("完整结果已保存至:", output_file_full, "\n")

  # 保存不同FDR阈值的结果
  for(threshold in names(filtered_results)) {
    if (nrow(filtered_results[[threshold]]) > 0) {
      output_file <- file.path(output_dir, 
                               paste0(erna_id, "_vs_", pcd_cn, "基因_Spearman_", threshold, ".csv"))
      fwrite(filtered_results[[threshold]], output_file)
      cat("已保存", threshold, "结果 (", nrow(filtered_results[[threshold]]), "个基因):", output_file, "\n")
    } else {
      cat(threshold, "阈值下无显著基因\n")
    }
  }

  # 10. 创建统计摘要
  summary_stats <- data.frame(
    FDR_Threshold = names(filtered_results),
    Significant_Genes = sapply(filtered_results, nrow),
    stringsAsFactors = FALSE
  )

  # 打印到控制台
  cat("\n===== 分析完成 =====\n")
  cat("eRNA:", erna_name, "\n")
  cat("分析的失巢死亡基因总数:", nrow(results), "\n")
  cat("成功计算相关性的基因数:", nrow(results_valid), "\n\n")

  cat("不同FDR阈值下的显著基因数量:\n")
  for(i in 1:nrow(summary_stats)) {
    cat("  ", summary_stats$FDR_Threshold[i], ":", summary_stats$Significant_Genes[i], "\n")
  }

  # 找最强正相关和负相关基因
  if (nrow(results_sorted) > 0) {
    cat("\n最强相关性基因:\n")

    # 正相关
    pos_idx <- which.max(results_sorted$Correlation)
    cat("  最强正相关基因:", results_sorted$Gene[pos_idx], 
        "(r =", round(results_sorted$Correlation[pos_idx], 4), 
        ", FDR =", format(results_sorted$FDR[pos_idx], scientific = TRUE, digits = 3), ")\n")

    # 负相关
    neg_idx <- which.min(results_sorted$Correlation)
    cat("  最强负相关基因:", results_sorted$Gene[neg_idx], 
        "(r =", round(results_sorted$Correlation[neg_idx], 4), 
        ", FDR =", format(results_sorted$FDR[neg_idx], scientific = TRUE, digits = 3), ")\n")
  }

  # 11. 保存统计结果到文件
  stats_file <- file.path(output_dir, "correlation_statistics.txt")
  sink(stats_file)
  cat("eRNA与失巢死亡基因相关性分析统计摘要\n")
  cat("====================================\n")
  cat("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("eRNA名称:", erna_name, "\n")
  cat("分析的失巢死亡基因总数:", nrow(results), "\n")
  cat("成功计算相关性的基因数:", nrow(results_valid), "\n")
  cat("失败计算的基因数:", nrow(results) - nrow(results_valid), "\n\n")

  cat("不同FDR阈值下的显著基因数量:\n")
  for(i in 1:nrow(summary_stats)) {
    cat(summary_stats$FDR_Threshold[i], ":", summary_stats$Significant_Genes[i], "\n")
  }

  if (nrow(results_sorted) > 0) {
    cat("\n最强相关性基因:\n")
    pos_idx <- which.max(results_sorted$Correlation)
    cat("最强正相关基因:", results_sorted$Gene[pos_idx], 
        "(r =", round(results_sorted$Correlation[pos_idx], 4), 
        ", FDR =", format(results_sorted$FDR[pos_idx], scientific = TRUE, digits = 3), ")\n")

    neg_idx <- which.min(results_sorted$Correlation)
    cat("最强负相关基因:", results_sorted$Gene[neg_idx], 
        "(r =", round(results_sorted$Correlation[neg_idx], 4), 
        ", FDR =", format(results_sorted$FDR[neg_idx], scientific = TRUE, digits = 3), ")\n")
  }

  # 数据质量检查信息
  cat("\n数据质量检查:\n")
  cat("样本数量:", ncol(gene_matrix_numeric), "\n")
  cat("失巢死亡基因数量:", nrow(results), "\n")
  cat("成功计算相关性的基因数:", nrow(results_valid), "\n")
  cat("失败率:", round((nrow(results) - nrow(results_valid))/nrow(results)*100, 2), "%\n")

  # 相关性分布统计
  if (nrow(results_valid) > 0) {
    cat("\n相关性分布统计:\n")
    cat("平均相关性:", round(mean(results_valid$Correlation, na.rm = TRUE), 4), "\n")
    cat("相关性中位数:", round(median(results_valid$Correlation, na.rm = TRUE), 4), "\n")
    cat("正相关基因数:", sum(results_valid$Correlation > 0, na.rm = TRUE), "\n")
    cat("负相关基因数:", sum(results_valid$Correlation < 0, na.rm = TRUE), "\n")
  }
  sink()

  cat("\n统计摘要已保存至:", stats_file, "\n")

  # ---------------------------------------------------------
  # 3. eRNA 与免疫细胞（ssGSEA 得分）
  # ---------------------------------------------------------
  # 1. 加载必要的包
  library(data.table)
  library(dplyr)

  # 2. 设置输入和输出路径
  input_dir <- "output/inputs"
  output_dir <- file.path(run_dir, "eRNA_vs_immune_cells")

  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("创建输出目录:", output_dir, "\n")
  }

  # 3. 读取数据（保留样本名信息）
  cat("正在读取数据...\n")

  # 读取eRNA表达矩阵
  erna_file <- file.path(input_dir, "eRNA_expression", paste0(erna_id, "表达.csv"))
  if (!file.exists(erna_file)) {
    stop("错误：找不到eRNA表达文件：", erna_file)
  }
  erna_data <- fread(erna_file, header = TRUE)

  # 检查eRNA数据格式
  if (nrow(erna_data) == 0) {
    stop("错误：eRNA表达矩阵为空！")
  }
  if (ncol(erna_data) < 2) {
    stop("错误：eRNA表达矩阵格式不正确，至少需要2列（基因名+表达值）")
  }

  # 提取eRNA样本名和表达值
  erna_samples <- colnames(erna_data)[-1]  # 获取样本名（排除第一列基因名）
  erna_vector <- as.numeric(erna_data[1, -1])  # 获取表达值向量
  erna_name <- erna_data[[1,1]]  # 获取eRNA的名称

  cat("✓ eRNA表达数据：", erna_name, "，共", length(erna_samples), "个样本\n")

  # 读取免疫细胞表达矩阵
  immune_file <- file.path("input", "immune_cell_ssGSEA_70tumors.csv")
  if (!file.exists(immune_file)) {
    stop("错误：找不到免疫细胞表达文件：", immune_file)
  }
  immune_data <- fread(immune_file, header = TRUE)

  # 检查免疫细胞数据格式
  if (nrow(immune_data) == 0) {
    stop("错误：免疫细胞表达矩阵为空！")
  }
  if (ncol(immune_data) < 2) {
    stop("错误：免疫细胞表达矩阵格式不正确，至少需要2列（免疫细胞名+表达值）")
  }

  # 提取免疫细胞样本名、免疫细胞名和表达矩阵
  immune_samples <- colnames(immune_data)[-1]  # 获取样本名（排除第一列免疫细胞名）
  immune_names <- immune_data[[1]]  # 获取免疫细胞名
  immune_matrix <- as.matrix(immune_data[, -1])  # 获取表达矩阵
  colnames(immune_matrix) <- immune_samples  # 为矩阵添加样本名

  cat("✓ 免疫细胞表达数据：", length(immune_names), "个免疫细胞，", length(immune_samples), "个样本\n")

  # 4. 样本匹配和验证
  cat("正在检查样本匹配情况...\n")

  # 检查是否有共同的样本
  common_samples <- intersect(erna_samples, immune_samples)

  if (length(common_samples) == 0) {
    stop("错误：eRNA和免疫细胞矩阵没有共同的样本！请检查样本名是否一致。")
  }

  cat("✓ 找到", length(common_samples), "个共同样本\n")

  # 检查是否有样本不匹配
  if (length(common_samples) < length(erna_samples)) {
    missing_in_immune <- setdiff(erna_samples, immune_samples)
    cat("警告：有", length(missing_in_immune), "个eRNA样本在免疫细胞矩阵中不存在：", 
        paste(missing_in_immune, collapse=", "), "\n")
  }

  if (length(common_samples) < length(immune_samples)) {
    missing_in_erna <- setdiff(immune_samples, erna_samples)
    cat("警告：有", length(missing_in_erna), "个免疫细胞样本在eRNA矩阵中不存在：", 
        paste(missing_in_erna, collapse=", "), "\n")
  }

  # 5. 按共同样本重排数据
  cat("正在按共同样本重排数据...\n")

  # 重排eRNA向量（按共同样本顺序）
  erna_vector_aligned <- erna_vector[match(common_samples, erna_samples)]

  # 重排免疫细胞矩阵（按共同样本顺序）
  immune_matrix_aligned <- immune_matrix[, match(common_samples, immune_samples), drop = FALSE]

  # 验证重排结果
  if (!all(colnames(immune_matrix_aligned) == common_samples)) {
    stop("错误：免疫细胞矩阵样本重排失败！")
  }

  if (length(erna_vector_aligned) != length(common_samples)) {
    stop("错误：eRNA向量重排后长度不匹配！")
  }

  cat("✓ 数据重排完成，使用", length(common_samples), "个共同样本进行分析\n")

  # 6. 计算斯皮尔曼相关性
  cat(sprintf("正在计算 eRNA %s 与 %d 个免疫细胞的斯皮尔曼相关性...\n", 
              erna_name, length(immune_names)))

  # 初始化结果数据框
  results <- data.frame(
    eRNA = character(),
    immune_cell = character(),
    Correlation = numeric(),
    P.Value = numeric(),
    stringsAsFactors = FALSE
  )

  # 循环计算每个免疫细胞与eRNA的相关性
  for (i in 1:length(immune_names)) {
    immune_exp <- as.numeric(immune_matrix_aligned[i, ])  # 获取第i个免疫细胞的表达向量

    # 检查并处理缺失值
    valid_data <- !(is.na(erna_vector_aligned) | is.na(immune_exp))
    n_valid <- sum(valid_data)

    if (n_valid < 3) {
      # 样本太少，无法计算相关性
      results[i, ] <- list(
        eRNA = erna_name,
        immune_cell = immune_names[i],
        Correlation = NA,
        P.Value = NA
      )
      next
    }

    # 进行斯皮尔曼相关性检验
    cor_test_result <- tryCatch({
      cor.test(erna_vector_aligned[valid_data], immune_exp[valid_data], 
               method = "spearman", exact = FALSE)
    }, error = function(e) {
      return(NULL)
    })

    if (is.null(cor_test_result)) {
      # 计算失败
      results[i, ] <- list(
        eRNA = erna_name,
        immune_cell = immune_names[i],
        Correlation = NA,
        P.Value = NA
      )
    } else {
      # 计算成功
      results[i, ] <- list(
        eRNA = erna_name,
        immune_cell = immune_names[i],
        Correlation = cor_test_result$estimate,
        P.Value = cor_test_result$p.value
      )
    }
  }

  # 移除无法计算相关性的免疫细胞
  results_valid <- results[!is.na(results$Correlation), ]
  cat(sprintf("✓ 成功计算 %d/%d 个免疫细胞的相关性\n", 
              nrow(results_valid), nrow(results)))

  # 7. 多重检验校正 (FDR)
  cat("正在进行多重检验校正 (FDR)...\n")
  results_valid$FDR <- p.adjust(results_valid$P.Value, method = "BH")

  # 8. 按相关性系数绝对值排序
  results_sorted <- results_valid[order(-abs(results_valid$Correlation)), ]

  # 9. 定义筛选条件（基于FDR）
  filtered_results <- list(
    FDR_0.05 = results_sorted %>% filter(FDR < 0.05) %>% arrange(desc(abs(Correlation))),
    FDR_0.01 = results_sorted %>% filter(FDR < 0.01) %>% arrange(desc(abs(Correlation))),
    FDR_0.001 = results_sorted %>% filter(FDR < 0.001) %>% arrange(desc(abs(Correlation))),
    FDR_0.0001 = results_sorted %>% filter(FDR < 0.0001) %>% arrange(desc(abs(Correlation))),
    FDR_0.00001 = results_sorted %>% filter(FDR < 0.00001) %>% arrange(desc(abs(Correlation)))
  )

  # 10. 保存结果
  # 保存完整结果
  output_file_full <- file.path(output_dir, paste0(erna_id, "表达_vs_免疫细胞_Spearman_Full_Results.csv"))
  fwrite(results_sorted, output_file_full)
  cat("完整结果已保存至:", output_file_full, "\n")

  # 保存不同FDR阈值的结果
  cat("保存不同FDR阈值的结果...\n")
  for(threshold in names(filtered_results)) {
    if (nrow(filtered_results[[threshold]]) > 0) {
      output_file <- file.path(output_dir, 
                               paste0(erna_id, "表达_vs_免疫细胞_Spearman_", threshold, ".csv"))
      fwrite(filtered_results[[threshold]], output_file)
      cat("  已保存", threshold, "结果 (", nrow(filtered_results[[threshold]]), "个免疫细胞):", 
          basename(output_file), "\n")
    } else {
      cat("  ", threshold, "阈值下无显著免疫细胞\n")
    }
  }

  # 11. 创建统计摘要
  summary_stats <- data.frame(
    FDR_Threshold = names(filtered_results),
    Significant_Immune_Cells = sapply(filtered_results, nrow),
    stringsAsFactors = FALSE
  )

  # 打印到控制台
  cat("\n===== 分析完成 =====\n")
  cat("eRNA:", erna_name, "\n")
  cat("分析的免疫细胞总数:", nrow(results), "\n")
  cat("成功计算相关性的免疫细胞数:", nrow(results_valid), "\n")
  cat("使用的共同样本数:", length(common_samples), "\n\n")

  cat("不同FDR阈值下的显著免疫细胞数量:\n")
  for(i in 1:nrow(summary_stats)) {
    cat("  ", summary_stats$FDR_Threshold[i], ":", summary_stats$Significant_Immune_Cells[i], "\n")
  }

  # 找最强正相关和负相关免疫细胞
  if (nrow(results_sorted) > 0) {
    cat("\n最强相关性免疫细胞:\n")

    # 正相关
    pos_idx <- which.max(results_sorted$Correlation)
    cat("  最强正相关免疫细胞:", results_sorted$immune_cell[pos_idx], 
        "(r =", round(results_sorted$Correlation[pos_idx], 4), 
        ", FDR =", format(results_sorted$FDR[pos_idx], scientific = TRUE, digits = 3), ")\n")

    # 负相关
    neg_idx <- which.min(results_sorted$Correlation)
    cat("  最强负相关免疫细胞:", results_sorted$immune_cell[neg_idx], 
        "(r =", round(results_sorted$Correlation[neg_idx], 4), 
        ", FDR =", format(results_sorted$FDR[neg_idx], scientific = TRUE, digits = 3), ")\n")
  }

  # 12. 保存统计结果到文件
  stats_file <- file.path(output_dir, "correlation_statistics.txt")
  sink(stats_file)
  cat("eRNA与免疫细胞相关性分析统计摘要\n")
  cat("====================================\n")
  cat("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("eRNA名称:", erna_name, "\n\n")

  cat("样本匹配信息:\n")
  cat("  eRNA样本数:", length(erna_samples), "\n")
  cat("  免疫细胞样本数:", length(immune_samples), "\n")
  cat("  共同样本数:", length(common_samples), "\n")
  if (length(erna_samples) > length(common_samples)) {
    cat("  缺失的eRNA样本数:", length(erna_samples) - length(common_samples), "\n")
  }
  if (length(immune_samples) > length(common_samples)) {
    cat("  缺失的免疫细胞样本数:", length(immune_samples) - length(common_samples), "\n")
  }

  cat("\n分析统计:\n")
  cat("  分析的免疫细胞总数:", nrow(results), "\n")
  cat("  成功计算相关性的免疫细胞数:", nrow(results_valid), "\n")
  cat("  失败计算的免疫细胞数:", nrow(results) - nrow(results_valid), "\n")
  cat("  计算成功率:", round(nrow(results_valid)/nrow(results)*100, 2), "%\n\n")

  cat("不同FDR阈值下的显著免疫细胞数量:\n")
  for(i in 1:nrow(summary_stats)) {
    cat(summary_stats$FDR_Threshold[i], ":", summary_stats$Significant_Immune_Cells[i], "\n")
  }

  if (nrow(results_sorted) > 0) {
    cat("\n最强相关性免疫细胞:\n")
    pos_idx <- which.max(results_sorted$Correlation)
    cat("最强正相关免疫细胞:", results_sorted$immune_cell[pos_idx], 
        "(r =", round(results_sorted$Correlation[pos_idx], 4), 
        ", FDR =", format(results_sorted$FDR[pos_idx], scientific = TRUE, digits = 3), ")\n")

    neg_idx <- which.min(results_sorted$Correlation)
    cat("最强负相关免疫细胞:", results_sorted$immune_cell[neg_idx], 
        "(r =", round(results_sorted$Correlation[neg_idx], 4), 
        ", FDR =", format(results_sorted$FDR[neg_idx], scientific = TRUE, digits = 3), ")\n")
  }

  # 数据质量检查信息
  cat("\n数据质量检查:\n")
  cat("  最终使用样本数量:", length(common_samples), "\n")
  cat("  免疫细胞数量:", nrow(results), "\n")
  cat("  成功计算相关性的免疫细胞数:", nrow(results_valid), "\n")
  cat("  失败率:", round((nrow(results) - nrow(results_valid))/nrow(results)*100, 2), "%\n")

  # 相关性分布统计
  if (nrow(results_valid) > 0) {
    cat("\n相关性分布统计:\n")
    cat("  平均相关性:", round(mean(results_valid$Correlation, na.rm = TRUE), 4), "\n")
    cat("  相关性中位数:", round(median(results_valid$Correlation, na.rm = TRUE), 4), "\n")
    cat("  正相关免疫细胞数:", sum(results_valid$Correlation > 0, na.rm = TRUE), "\n")
    cat("  负相关免疫细胞数:", sum(results_valid$Correlation < 0, na.rm = TRUE), "\n")

    # 添加相关性强度分布
    cat("\n相关性强度分布:\n")
    cat("  |r| > 0.8:", sum(abs(results_valid$Correlation) > 0.8), "个免疫细胞\n")
    cat("  0.6 < |r| ≤ 0.8:", sum(abs(results_valid$Correlation) > 0.6 & abs(results_valid$Correlation) <= 0.8), "个免疫细胞\n")
    cat("  0.4 < |r| ≤ 0.6:", sum(abs(results_valid$Correlation) > 0.4 & abs(results_valid$Correlation) <= 0.6), "个免疫细胞\n")
    cat("  0.2 < |r| ≤ 0.4:", sum(abs(results_valid$Correlation) > 0.2 & abs(results_valid$Correlation) <= 0.4), "个免疫细胞\n")
    cat("  |r| ≤ 0.2:", sum(abs(results_valid$Correlation) <= 0.2), "个免疫细胞\n")
  }
  sink()

  cat("\n统计摘要已保存至:", stats_file, "\n")

  # 13. 输出样本匹配详情（可选）
  sample_match_file <- file.path(output_dir, "样本匹配详情.txt")
  writeLines(
    c(
      "样本匹配详情",
      "=============",
      paste("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste("eRNA名称:", erna_name),
      "",
      "eRNA样本列表:",
      paste(erna_samples, collapse=", "),
      "",
      "免疫细胞样本列表:",
      paste(immune_samples, collapse=", "),
      "",
      "共同样本列表（实际用于分析）:",
      paste(common_samples, collapse=", "),
      "",
      if (length(setdiff(erna_samples, immune_samples)) > 0) {
        c("仅在eRNA中存在的样本:", paste(setdiff(erna_samples, immune_samples), collapse=", "))
      } else {
        "所有eRNA样本都在免疫细胞矩阵中"
      },
      "",
      if (length(setdiff(immune_samples, erna_samples)) > 0) {
        c("仅在免疫细胞中存在的样本:", paste(setdiff(immune_samples, erna_samples), collapse=", "))
      } else {
        "所有免疫细胞样本都在eRNA矩阵中"
      }
    ),
    sample_match_file
  )

  cat("样本匹配详情已保存至:", sample_match_file, "\n")
  cat("\n分析完成！所有结果已保存到:", output_dir, "\n")

  # ---------------------------------------------------------
  # 4. eRNA 与 19 种候选药物（预测 IC50）
  # ---------------------------------------------------------
  # 1. 加载必要的包
  library(data.table)
  library(dplyr)

  # 2. 设置输入和输出路径
  input_dir <- "output/inputs"
  output_dir <- file.path(run_dir, "eRNA_vs_drugs")

  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("创建输出目录:", output_dir, "\n")
  }

  # 3. 读取数据（保留样本名信息）
  cat("正在读取数据...\n")

  # 读取eRNA表达矩阵
  erna_file <- file.path(input_dir, "eRNA_expression", paste0(erna_id, "表达.csv"))
  if (!file.exists(erna_file)) {
    stop("错误：找不到eRNA表达文件：", erna_file)
  }
  erna_data <- fread(erna_file, header = TRUE)

  # 检查eRNA数据格式
  if (nrow(erna_data) == 0) {
    stop("错误：eRNA表达矩阵为空！")
  }
  if (ncol(erna_data) < 2) {
    stop("错误：eRNA表达矩阵格式不正确，至少需要2列（基因名+表达值）")
  }

  # 提取eRNA样本名和表达值
  erna_samples <- colnames(erna_data)[-1]  # 获取样本名（排除第一列基因名）
  erna_vector <- as.numeric(erna_data[1, -1])  # 获取表达值向量
  erna_name <- erna_data[[1,1]]  # 获取eRNA的名称

  cat("✓ eRNA表达数据：", erna_name, "，共", length(erna_samples), "个样本\n")

  # 读取敏感药物表达矩阵
  drug_file <- file.path("input", "drug_IC50_19drugs_70tumors.csv")
  if (!file.exists(drug_file)) {
    stop("错误：找不到敏感药物表达文件：", drug_file)
  }
  drug_data <- fread(drug_file, header = TRUE)

  # 检查敏感药物数据格式
  if (nrow(drug_data) == 0) {
    stop("错误：敏感药物表达矩阵为空！")
  }
  if (ncol(drug_data) < 2) {
    stop("错误：敏感药物表达矩阵格式不正确，至少需要2列（药物名+表达值）")
  }

  # 提取敏感药物样本名、药物名和表达矩阵
  drug_samples <- colnames(drug_data)[-1]  # 获取样本名（排除第一列药物名）
  drug_names <- drug_data[[1]]  # 获取药物名
  drug_matrix <- as.matrix(drug_data[, -1])  # 获取表达矩阵
  colnames(drug_matrix) <- drug_samples  # 为矩阵添加样本名

  cat("✓ 敏感药物表达数据：", length(drug_names), "种药物，", length(drug_samples), "个样本\n")

  # 4. 样本匹配和验证
  cat("正在检查样本匹配情况...\n")

  # 检查是否有共同的样本
  common_samples <- intersect(erna_samples, drug_samples)

  if (length(common_samples) == 0) {
    stop("错误：eRNA和敏感药物矩阵没有共同的样本！请检查样本名是否一致。")
  }

  cat("✓ 找到", length(common_samples), "个共同样本\n")

  # 检查是否有样本不匹配
  if (length(common_samples) < length(erna_samples)) {
    missing_in_drug <- setdiff(erna_samples, drug_samples)
    cat("警告：有", length(missing_in_drug), "个eRNA样本在敏感药物矩阵中不存在：", 
        paste(missing_in_drug, collapse=", "), "\n")
  }

  if (length(common_samples) < length(drug_samples)) {
    missing_in_erna <- setdiff(drug_samples, erna_samples)
    cat("警告：有", length(missing_in_erna), "个敏感药物样本在eRNA矩阵中不存在：", 
        paste(missing_in_erna, collapse=", "), "\n")
  }

  # 5. 按共同样本重排数据
  cat("正在按共同样本重排数据...\n")

  # 重排eRNA向量（按共同样本顺序）
  erna_vector_aligned <- erna_vector[match(common_samples, erna_samples)]

  # 重排敏感药物矩阵（按共同样本顺序）
  drug_matrix_aligned <- drug_matrix[, match(common_samples, drug_samples), drop = FALSE]

  # 验证重排结果
  if (!all(colnames(drug_matrix_aligned) == common_samples)) {
    stop("错误：敏感药物矩阵样本重排失败！")
  }

  if (length(erna_vector_aligned) != length(common_samples)) {
    stop("错误：eRNA向量重排后长度不匹配！")
  }

  cat("✓ 数据重排完成，使用", length(common_samples), "个共同样本进行分析\n")

  # 6. 计算斯皮尔曼相关性
  cat(sprintf("正在计算 eRNA %s 与 %d 种敏感药物的斯皮尔曼相关性...\n", 
              erna_name, length(drug_names)))

  # 初始化结果数据框
  results <- data.frame(
    eRNA = character(),
    drug = character(),
    Correlation = numeric(),
    P.Value = numeric(),
    stringsAsFactors = FALSE
  )

  # 循环计算每种药物与eRNA的相关性
  for (i in 1:length(drug_names)) {
    drug_exp <- as.numeric(drug_matrix_aligned[i, ])  # 获取第i种药物的表达向量

    # 检查并处理缺失值
    valid_data <- !(is.na(erna_vector_aligned) | is.na(drug_exp))
    n_valid <- sum(valid_data)

    if (n_valid < 3) {
      # 样本太少，无法计算相关性
      results[i, ] <- list(
        eRNA = erna_name,
        drug = drug_names[i],
        Correlation = NA,
        P.Value = NA
      )
      next
    }

    # 进行斯皮尔曼相关性检验
    cor_test_result <- tryCatch({
      cor.test(erna_vector_aligned[valid_data], drug_exp[valid_data], 
               method = "spearman", exact = FALSE)
    }, error = function(e) {
      return(NULL)
    })

    if (is.null(cor_test_result)) {
      # 计算失败
      results[i, ] <- list(
        eRNA = erna_name,
        drug = drug_names[i],
        Correlation = NA,
        P.Value = NA
      )
    } else {
      # 计算成功
      results[i, ] <- list(
        eRNA = erna_name,
        drug = drug_names[i],
        Correlation = cor_test_result$estimate,
        P.Value = cor_test_result$p.value
      )
    }
  }

  # 移除无法计算相关性的药物
  results_valid <- results[!is.na(results$Correlation), ]
  cat(sprintf("✓ 成功计算 %d/%d 种药物的相关性\n", 
              nrow(results_valid), nrow(results)))

  # 7. 多重检验校正 (FDR)
  cat("正在进行多重检验校正 (FDR)...\n")
  results_valid$FDR <- p.adjust(results_valid$P.Value, method = "BH")

  # 8. 按相关性系数绝对值排序
  results_sorted <- results_valid[order(-abs(results_valid$Correlation)), ]

  # 9. 定义筛选条件（基于FDR）
  filtered_results <- list(
    FDR_0.05 = results_sorted %>% filter(FDR < 0.05) %>% arrange(desc(abs(Correlation))),
    FDR_0.01 = results_sorted %>% filter(FDR < 0.01) %>% arrange(desc(abs(Correlation))),
    FDR_0.001 = results_sorted %>% filter(FDR < 0.001) %>% arrange(desc(abs(Correlation))),
    FDR_0.0001 = results_sorted %>% filter(FDR < 0.0001) %>% arrange(desc(abs(Correlation))),
    FDR_0.00001 = results_sorted %>% filter(FDR < 0.00001) %>% arrange(desc(abs(Correlation)))
  )

  # 10. 保存结果
  # 保存完整结果
  output_file_full <- file.path(output_dir, paste0(erna_id, "表达_vs_敏感药物_Spearman_Full_Results.csv"))
  fwrite(results_sorted, output_file_full)
  cat("完整结果已保存至:", output_file_full, "\n")

  # 保存不同FDR阈值的结果
  cat("保存不同FDR阈值的结果...\n")
  for(threshold in names(filtered_results)) {
    if (nrow(filtered_results[[threshold]]) > 0) {
      output_file <- file.path(output_dir, 
                               paste0(erna_id, "表达_vs_敏感药物_Spearman_", threshold, ".csv"))
      fwrite(filtered_results[[threshold]], output_file)
      cat("  已保存", threshold, "结果 (", nrow(filtered_results[[threshold]]), "种药物):", 
          basename(output_file), "\n")
    } else {
      cat("  ", threshold, "阈值下无显著药物\n")
    }
  }

  # 11. 创建统计摘要
  summary_stats <- data.frame(
    FDR_Threshold = names(filtered_results),
    Significant_Drugs = sapply(filtered_results, nrow),
    stringsAsFactors = FALSE
  )

  # 打印到控制台
  cat("\n===== 分析完成 =====\n")
  cat("eRNA:", erna_name, "\n")
  cat("分析的敏感药物总数:", nrow(results), "\n")
  cat("成功计算相关性的药物数:", nrow(results_valid), "\n")
  cat("使用的共同样本数:", length(common_samples), "\n\n")

  cat("不同FDR阈值下的显著药物数量:\n")
  for(i in 1:nrow(summary_stats)) {
    cat("  ", summary_stats$FDR_Threshold[i], ":", summary_stats$Significant_Drugs[i], "\n")
  }

  # 找最强正相关和负相关药物
  if (nrow(results_sorted) > 0) {
    cat("\n最强相关性药物:\n")

    # 正相关
    pos_idx <- which.max(results_sorted$Correlation)
    cat("  最强正相关药物:", results_sorted$drug[pos_idx], 
        "(r =", round(results_sorted$Correlation[pos_idx], 4), 
        ", FDR =", format(results_sorted$FDR[pos_idx], scientific = TRUE, digits = 3), ")\n")

    # 负相关
    neg_idx <- which.min(results_sorted$Correlation)
    cat("  最强负相关药物:", results_sorted$drug[neg_idx], 
        "(r =", round(results_sorted$Correlation[neg_idx], 4), 
        ", FDR =", format(results_sorted$FDR[neg_idx], scientific = TRUE, digits = 3), ")\n")
  }

  # 12. 保存统计结果到文件
  stats_file <- file.path(output_dir, "correlation_statistics.txt")
  sink(stats_file)
  cat("eRNA与敏感药物相关性分析统计摘要\n")
  cat("====================================\n")
  cat("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("eRNA名称:", erna_name, "\n\n")

  cat("样本匹配信息:\n")
  cat("  eRNA样本数:", length(erna_samples), "\n")
  cat("  敏感药物样本数:", length(drug_samples), "\n")
  cat("  共同样本数:", length(common_samples), "\n")
  if (length(erna_samples) > length(common_samples)) {
    cat("  缺失的eRNA样本数:", length(erna_samples) - length(common_samples), "\n")
  }
  if (length(drug_samples) > length(common_samples)) {
    cat("  缺失的敏感药物样本数:", length(drug_samples) - length(common_samples), "\n")
  }

  cat("\n分析统计:\n")
  cat("  分析的敏感药物总数:", nrow(results), "\n")
  cat("  成功计算相关性的药物数:", nrow(results_valid), "\n")
  cat("  失败计算的药物数:", nrow(results) - nrow(results_valid), "\n")
  cat("  计算成功率:", round(nrow(results_valid)/nrow(results)*100, 2), "%\n\n")

  cat("不同FDR阈值下的显著药物数量:\n")
  for(i in 1:nrow(summary_stats)) {
    cat(summary_stats$FDR_Threshold[i], ":", summary_stats$Significant_Drugs[i], "\n")
  }

  if (nrow(results_sorted) > 0) {
    cat("\n最强相关性药物:\n")
    pos_idx <- which.max(results_sorted$Correlation)
    cat("最强正相关药物:", results_sorted$drug[pos_idx], 
        "(r =", round(results_sorted$Correlation[pos_idx], 4), 
        ", FDR =", format(results_sorted$FDR[pos_idx], scientific = TRUE, digits = 3), ")\n")

    neg_idx <- which.min(results_sorted$Correlation)
    cat("最强负相关药物:", results_sorted$drug[neg_idx], 
        "(r =", round(results_sorted$Correlation[neg_idx], 4), 
        ", FDR =", format(results_sorted$FDR[neg_idx], scientific = TRUE, digits = 3), ")\n")
  }

  # 数据质量检查信息
  cat("\n数据质量检查:\n")
  cat("  最终使用样本数量:", length(common_samples), "\n")
  cat("  敏感药物数量:", nrow(results), "\n")
  cat("  成功计算相关性的药物数:", nrow(results_valid), "\n")
  cat("  失败率:", round((nrow(results) - nrow(results_valid))/nrow(results)*100, 2), "%\n")

  # 相关性分布统计
  if (nrow(results_valid) > 0) {
    cat("\n相关性分布统计:\n")
    cat("  平均相关性:", round(mean(results_valid$Correlation, na.rm = TRUE), 4), "\n")
    cat("  相关性中位数:", round(median(results_valid$Correlation, na.rm = TRUE), 4), "\n")
    cat("  正相关药物数:", sum(results_valid$Correlation > 0, na.rm = TRUE), "\n")
    cat("  负相关药物数:", sum(results_valid$Correlation < 0, na.rm = TRUE), "\n")

    # 添加相关性强度分布
    cat("\n相关性强度分布:\n")
    cat("  |r| > 0.8:", sum(abs(results_valid$Correlation) > 0.8), "种药物\n")
    cat("  0.6 < |r| ≤ 0.8:", sum(abs(results_valid$Correlation) > 0.6 & abs(results_valid$Correlation) <= 0.8), "种药物\n")
    cat("  0.4 < |r| ≤ 0.6:", sum(abs(results_valid$Correlation) > 0.4 & abs(results_valid$Correlation) <= 0.6), "种药物\n")
    cat("  0.2 < |r| ≤ 0.4:", sum(abs(results_valid$Correlation) > 0.2 & abs(results_valid$Correlation) <= 0.4), "种药物\n")
    cat("  |r| ≤ 0.2:", sum(abs(results_valid$Correlation) <= 0.2), "种药物\n")
  }
  sink()

  cat("\n统计摘要已保存至:", stats_file, "\n")

  # 13. 输出样本匹配详情（可选）
  sample_match_file <- file.path(output_dir, "样本匹配详情.txt")
  writeLines(
    c(
      "样本匹配详情",
      "=============",
      paste("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste("eRNA名称:", erna_name),
      "",
      "eRNA样本列表:",
      paste(erna_samples, collapse=", "),
      "",
      "敏感药物样本列表:",
      paste(drug_samples, collapse=", "),
      "",
      "共同样本列表（实际用于分析）:",
      paste(common_samples, collapse=", "),
      "",
      if (length(setdiff(erna_samples, drug_samples)) > 0) {
        c("仅在eRNA中存在的样本:", paste(setdiff(erna_samples, drug_samples), collapse=", "))
      } else {
        "所有eRNA样本都在敏感药物矩阵中"
      },
      "",
      if (length(setdiff(drug_samples, erna_samples)) > 0) {
        c("仅在敏感药物中存在的样本:", paste(setdiff(drug_samples, erna_samples), collapse=", "))
      } else {
        "所有敏感药物样本都在eRNA矩阵中"
      }
    ),
    sample_match_file
  )

  cat("样本匹配详情已保存至:", sample_match_file, "\n")
  cat("\n分析完成！所有结果已保存到:", output_dir, "\n")

  # ---------------------------------------------------------
  # 5. eRNA 与免疫治疗相关基因
  # ---------------------------------------------------------
  # 1. 加载必要的包
  library(data.table)
  library(dplyr)

  # 2. 设置输入和输出路径
  input_dir <- "output/inputs"
  output_dir <- file.path(run_dir, "eRNA_vs_immunotherapy_genes")

  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 3. 读取数据
  # 读取eRNA表达矩阵
  erna_file <- file.path(input_dir, "eRNA_expression", paste0(erna_id, "表达.csv"))
  erna_data <- fread(erna_file, header = TRUE)

  # 读取免疫治疗基因表达矩阵
  gene_file <- file.path(input_dir, "免疫治疗基因表达矩阵.csv")
  gene_data <- fread(gene_file, header = TRUE)

  # 显示数据基本信息
  cat("eRNA表达矩阵维度:", dim(erna_data), "\n")
  cat("免疫治疗基因表达矩阵维度:", dim(gene_data), "\n")

  # 4. 自动样本对齐函数
  align_samples_auto <- function(erna_data, gene_data) {
    # 提取样本名称
    erna_samples <- colnames(erna_data)[-1]  # 排除第一列（eRNA名）
    gene_samples <- colnames(gene_data)[-1]  # 排除第一列（基因名）

    cat("eRNA样本数:", length(erna_samples), "\n")
    cat("基因矩阵样本数:", length(gene_samples), "\n")

    # 找共同样本
    common_samples <- intersect(erna_samples, gene_samples)

    if (length(common_samples) == 0) {
      stop("错误：没有共同的样本！无法进行分析。")
    }

    cat("找到共同样本数:", length(common_samples), "\n")
    cat("eRNA特有样本数:", length(setdiff(erna_samples, gene_samples)), "\n")
    cat("基因矩阵特有样本数:", length(setdiff(gene_samples, erna_samples)), "\n")

    if (length(common_samples) < max(length(erna_samples), length(gene_samples))) {
      cat("警告：部分样本不匹配，将只使用共同样本进行分析\n")
    }

    # 按共同样本重新排列
    erna_index <- match(common_samples, erna_samples)
    gene_index <- match(common_samples, gene_samples)

    # 提取对齐后的数据（保留第一列名称）
    aligned_erna <- erna_data[, c(1, erna_index + 1), with = FALSE]
    aligned_gene <- gene_data[, c(1, gene_index + 1), with = FALSE]

    # 更新列名
    colnames(aligned_erna) <- c(colnames(erna_data)[1], common_samples)
    colnames(aligned_gene) <- c(colnames(gene_data)[1], common_samples)

    return(list(
      erna_data = aligned_erna,
      gene_data = aligned_gene,
      common_samples = common_samples
    ))
  }

  # 5. 执行自动样本对齐
  cat("\n=== 执行自动样本对齐 ===\n")
  aligned_data <- align_samples_auto(erna_data, gene_data)

  erna_data <- aligned_data$erna_data
  gene_data <- aligned_data$gene_data
  common_samples <- aligned_data$common_samples

  # 保存对齐信息
  alignment_info_file <- file.path(output_dir, "样本对齐信息.txt")
  sink(alignment_info_file)
  cat("样本对齐信息\n")
  cat("============\n")
  cat("对齐日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("共同样本数:", length(common_samples), "\n")
  cat("共同样本列表:\n")
  cat(paste(common_samples, collapse = "\n"), "\n")
  sink()
  cat("样本对齐信息已保存至:", alignment_info_file, "\n")

  # 6. 数据预处理（使用对齐后的数据）
  # 转换为向量格式，方便后续计算
  erna_vector <- as.numeric(erna_data[1, -1]) # 第一行，排除第一列（名称列）
  erna_name <- erna_data[[1,1]] # 获取eRNA的名字

  # 提取基因名
  gene_names <- gene_data[[1]]
  # 提取表达矩阵（转换为数值矩阵）
  gene_matrix <- as.matrix(gene_data[, -1])

  # 检查样本是否匹配（对齐后应该完全一致）
  if (ncol(gene_matrix) != length(erna_vector)) {
    stop("错误：样本对齐后数量仍然不匹配！")
  } else {
    cat("数据检查通过：样本数量一致，共", ncol(gene_matrix), "个样本。\n")
  }

  # 检查样本名称是否一致
  if (!identical(colnames(erna_data)[-1], colnames(gene_data)[-1])) {
    warning("警告：样本名称不完全一致，但数量相同")
  } else {
    cat("样本名称完全一致！\n")
  }

  # 7. 计算斯皮尔曼相关性
  # 初始化结果数据框
  results <- data.frame(
    eRNA = character(),
    Gene = character(),
    Correlation = numeric(),
    P.Value = numeric(),
    stringsAsFactors = FALSE
  )

  cat("正在计算 eRNA:", erna_name, "与", length(gene_names), "个免疫治疗基因的斯皮尔曼相关性...\n")

  # 循环计算每个基因与eRNA的相关性
  for (i in 1:length(gene_names)) {
    gene_exp <- as.numeric(gene_matrix[i, ]) # 获取第i个基因的表达向量

    # 进行斯皮尔曼相关性检验
    cor_test_result <- cor.test(erna_vector, gene_exp, method = "spearman", exact = FALSE)

    # 将结果添加到数据框
    results[i, ] <- list(
      eRNA = erna_name,
      Gene = gene_names[i],
      Correlation = cor_test_result$estimate,
      P.Value = cor_test_result$p.value
    )
  }

  # 8. 多重检验校正 (FDR)
  results$FDR <- p.adjust(results$P.Value, method = "BH")

  # 9. 按相关性系数绝对值排序
  results_sorted <- results[order(-abs(results$Correlation)), ]

  # 10. 定义筛选条件（基于FDR）
  filtered_results <- list(
    FDR_0.05 = results_sorted %>% filter(FDR < 0.05) %>% arrange(desc(abs(Correlation))),
    FDR_0.01 = results_sorted %>% filter(FDR < 0.01) %>% arrange(desc(abs(Correlation))),
    FDR_0.001 = results_sorted %>% filter(FDR < 0.001) %>% arrange(desc(abs(Correlation))),
    FDR_0.0001 = results_sorted %>% filter(FDR < 0.0001) %>% arrange(desc(abs(Correlation))),
    FDR_0.00001 = results_sorted %>% filter(FDR < 0.00001) %>% arrange(desc(abs(Correlation)))
  )

  # 11. 保存结果
  # 保存完整结果
  output_file_full <- file.path(output_dir, paste0(erna_id, "_vs_免疫基因_Spearman_Full_Results.csv"))
  fwrite(results_sorted, output_file_full)
  cat("完整结果已保存至:", output_file_full, "\n")

  # 保存不同FDR阈值的结果
  for(threshold in names(filtered_results)) {
    output_file <- file.path(output_dir, paste0(erna_id, "_vs_免疫基因_Spearman_", threshold, ".csv"))
    fwrite(filtered_results[[threshold]], output_file)
    cat("已保存", threshold, "结果:", output_file, "\n")
  }

  # 12. 打印简要摘要并写入统计文件
  cat("\n===== 分析完成 =====\n")
  cat("eRNA:", erna_name, "\n")
  cat("分析的免疫治疗基因数量:", nrow(results), "\n")
  cat("使用的样本数量:", length(common_samples), "\n")

  # 创建统计摘要
  summary_stats <- data.frame(
    FDR_Threshold = names(filtered_results),
    Significant_Genes = sapply(filtered_results, nrow),
    stringsAsFactors = FALSE
  )

  # 打印到控制台
  for(i in 1:nrow(summary_stats)) {
    cat("显著相关的免疫治疗基因数量 (", summary_stats$FDR_Threshold[i], "):", summary_stats$Significant_Genes[i], "\n")
  }

  # 保存统计结果到文件
  stats_file <- file.path(output_dir, "correlation_statistics.txt")
  sink(stats_file)
  cat("eRNA与免疫治疗基因相关性分析统计摘要\n")
  cat("====================================\n")
  cat("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("eRNA名称:", erna_name, "\n")
  cat("分析的免疫治疗基因总数:", nrow(results), "\n")
  cat("使用的样本数量:", length(common_samples), "\n\n")

  cat("样本对齐信息:\n")
  cat("共同样本数:", length(common_samples), "\n")
  if (length(common_samples) < 10) {
    cat("样本列表:", paste(common_samples, collapse = ", "), "\n")
  } else {
    cat("前10个样本:", paste(common_samples[1:10], collapse = ", "), "...\n")
  }

  cat("\n不同FDR阈值下的显著基因数量:\n")
  for(i in 1:nrow(summary_stats)) {
    cat(summary_stats$FDR_Threshold[i], ":", summary_stats$Significant_Genes[i], "\n")
  }

  cat("\n最强相关性基因:\n")
  if (nrow(results_sorted) > 0) {
    cat("最强正相关基因:", results_sorted$Gene[1], "(r =", round(results_sorted$Correlation[1], 4), 
        ", FDR =", format(results_sorted$FDR[1], scientific = TRUE, digits = 3), ")\n")
    min_cor_index <- which.min(results_sorted$Correlation)
    cat("最强负相关基因:", results_sorted$Gene[min_cor_index], "(r =", round(results_sorted$Correlation[min_cor_index], 4),
        ", FDR =", format(results_sorted$FDR[min_cor_index], scientific = TRUE, digits = 3), ")\n")
  }
  sink()

  cat("统计摘要已保存至:", stats_file, "\n")
}
