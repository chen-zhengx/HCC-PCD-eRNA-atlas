# =========================================================
# 对应结果：Figure 8、Figure S8、Figure S9 与 Supplementary Table 12 的多层共表达关联网络——PCD 相关基因与其他各层之间的边：
#           PCD 相关基因与候选靶基因、免疫细胞、候选药物、免疫治疗相关基因的 Spearman 相关（70 个肿瘤样本）
#           每个 eRNA × 每层分别做 BH 校正，FDR < 0.01 为显著边
# 输入：output/correlations/<PCD>/<eRNA>/extracted/（04 脚本）、output/inputs/（02 脚本）
#       input/immune_cell_ssGSEA_70tumors.csv、input/drug_IC50_19drugs_70tumors.csv、input/risk_eRNAs_by_PCD.csv
# 输出：output/correlations/<PCD>/<eRNA>/target_genes_vs_PCD_genes|immune_cells_vs_PCD_genes|drugs_vs_PCD_genes|
#       immunotherapy_genes_vs_PCD_genes/：完整结果与各 FDR 阈值的显著结果（网络使用 FDR_0.01）
# 运行：在模块目录下 Rscript code/05_PCD_gene_correlations.R
# =========================================================

fwrite_pcd <- function(x, file, ...) data.table::fwrite(setNames(x, sub("失巢死亡", pcd_cn, names(x))), file, ...)

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

  # 第二层相关的输入：与该 eRNA 显著相关（FDR < 0.01）的节点（04 脚本的提取结果）。
  #        若该 eRNA 在某一层没有显著节点，则用该层的全部节点作为输入（例如失巢凋亡网络中 ENSR00000154915 的免疫细胞层
  #        使用了全部 28 种免疫细胞，ENSR00000250837 的药物层使用了全部 19 种药物，铁死亡网络中 ENSR00000084974 使用了全部
  #        50 个铁死亡相关基因）
  ext_dir <- file.path(run_dir, "extracted")
  use_file <- function(extracted, all_nodes) if (file.exists(extracted)) extracted else all_nodes
  pcd_input    <- use_file(file.path(ext_dir, paste0("提取的", pcd_cn, "基因表达矩阵.csv")), pcd_matrix_file)
  target_input <- use_file(file.path(ext_dir, "提取的靶基因表达矩阵.csv"),
                           file.path("output", "inputs", "target_genes", paste0(erna_id, "_1MB靶基因矩阵.csv")))
  immune_input <- use_file(file.path(ext_dir, "提取的免疫细胞表达矩阵.csv"), file.path("input", "immune_cell_ssGSEA_70tumors.csv"))
  drug_input   <- use_file(file.path(ext_dir, "提取的药物表达矩阵.csv"), file.path("input", "drug_IC50_19drugs_70tumors.csv"))
  immuno_input <- use_file(file.path(ext_dir, "提取的免疫治疗基因表达矩阵.csv"), file.path("output", "inputs", "免疫治疗基因表达矩阵.csv"))
  cat("第二层输入:", basename(c(pcd_input, target_input, immune_input, drug_input, immuno_input)), sep = "\n  ")

  # ---------------------------------------------------------
  # 1. ±1 Mb 候选靶基因与 PCD 相关基因
  # ---------------------------------------------------------
  # 1. 加载必要的包
  library(data.table)
  library(dplyr)

  # 2. 设置输入和输出路径
  input_dir <- file.path(run_dir, "extracted")
  output_dir <- file.path(run_dir, "target_genes_vs_PCD_genes")

  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("创建输出目录:", output_dir, "\n")
  }

  # 3. 读取数据
  cat("正在读取数据...\n")

  # 读取提取靶基因表达矩阵 
  target_gene_file <- target_input
  if (!file.exists(target_gene_file)) {
    stop("错误：找不到靶基因表达矩阵文件 - ", target_gene_file)
  }
  target_gene_data <- fread(target_gene_file, header = TRUE)

  # 检查靶基因表达矩阵是否有数据
  if (nrow(target_gene_data) == 0) {
    stop("错误：靶基因表达矩阵为空！")
  }

  cat("✓ 靶基因表达矩阵: ", nrow(target_gene_data), "个基因, ", 
      ncol(target_gene_data) - 1, "个样本\n")

  # 读取失巢死亡基因表达矩阵
  death_gene_file <- pcd_input
  if (!file.exists(death_gene_file)) {
    stop("错误：找不到失巢死亡基因表达矩阵文件 - ", death_gene_file)
  }
  death_gene_data <- fread(death_gene_file, header = TRUE)

  cat("✓ 失巢死亡基因表达矩阵: ", nrow(death_gene_data), "个基因, ", 
      ncol(death_gene_data) - 1, "个样本\n")

  # 4. 数据预处理与检查
  # 提取基因名
  target_gene_names <- target_gene_data[[1]]
  death_gene_names <- death_gene_data[[1]]

  # 提取表达矩阵（转换为数值矩阵）
  target_gene_matrix <- as.matrix(target_gene_data[, -1])
  death_gene_matrix <- as.matrix(death_gene_data[, -1])

  # 检查样本是否匹配
  if (ncol(target_gene_matrix) != ncol(death_gene_matrix)) {
    stop("错误：靶基因和失巢死亡基因矩阵的样本数量不匹配！")
  }

  cat("✓ 数据检查通过：样本数量一致，共", ncol(target_gene_matrix), "个样本\n")

  # 5. 计算斯皮尔曼相关性（所有基因对）
  total_pairs <- length(target_gene_names) * length(death_gene_names)
  cat("\n正在计算所有靶基因与失巢死亡基因的斯皮尔曼相关性...\n")
  cat("  靶基因数量:", length(target_gene_names), "\n")
  cat("  失巢死亡基因数量:", length(death_gene_names), "\n")
  cat("  总共需要计算的相关性数量:", total_pairs, "\n")

  # 预分配结果矩阵（提高性能）
  results <- data.frame(
    靶基因 = character(total_pairs),
    失巢死亡基因 = character(total_pairs),
    相关性 = numeric(total_pairs),
    P值 = numeric(total_pairs),
    stringsAsFactors = FALSE
  )

  # 初始化进度条
  pb <- txtProgressBar(min = 0, max = total_pairs, style = 3)
  counter <- 0

  # 循环计算每个靶基因与每个失巢死亡基因的相关性
  for (i in 1:length(target_gene_names)) {
    target_gene_vector <- as.numeric(target_gene_matrix[i, ])

    for (j in 1:length(death_gene_names)) {
      death_gene_exp <- as.numeric(death_gene_matrix[j, ])

      # 检查并处理缺失值
      valid_data <- !(is.na(target_gene_vector) | is.na(death_gene_exp))
      n_valid <- sum(valid_data)

      if (n_valid < 3) {
        # 样本太少，无法计算相关性
        results[counter + 1, ] <- list(
          靶基因 = target_gene_names[i],
          失巢死亡基因 = death_gene_names[j],
          相关性 = NA,
          P值 = NA
        )
      } else {
        # 进行斯皮尔曼相关性检验
        cor_test_result <- tryCatch({
          cor.test(target_gene_vector[valid_data], 
                   death_gene_exp[valid_data], 
                   method = "spearman", exact = FALSE)
        }, error = function(e) {
          return(NULL)
        })

        if (is.null(cor_test_result)) {
          # 计算失败
          results[counter + 1, ] <- list(
            靶基因 = target_gene_names[i],
            失巢死亡基因 = death_gene_names[j],
            相关性 = NA,
            P值 = NA
          )
        } else {
          # 计算成功
          results[counter + 1, ] <- list(
            靶基因 = target_gene_names[i],
            失巢死亡基因 = death_gene_names[j],
            相关性 = cor_test_result$estimate,
            P值 = cor_test_result$p.value
          )
        }
      }

      # 更新进度条
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
    }
  }

  close(pb)
  cat("\n✓ 相关性计算完成！\n")

  # 移除无法计算相关性的行
  results_valid <- results[!is.na(results$相关性), ]
  cat("  成功计算的相关性数量:", nrow(results_valid), "/", total_pairs, "\n")

  # 6. 多重检验校正 (FDR)
  cat("正在进行多重检验校正 (FDR)...\n")
  results_valid$FDR <- p.adjust(results_valid$P值, method = "BH")

  # 检查FDR计算
  na_fdr_count <- sum(is.na(results_valid$FDR))
  if (na_fdr_count > 0) {
    cat("警告: FDR计算有", na_fdr_count, "个NA值\n")
  } else {
    cat("✓ FDR计算完成\n")
  }

  # 7. 按相关性系数绝对值排序
  results_sorted <- results_valid[order(-abs(results_valid$相关性)), ]

  # 显示前几行结果
  cat("\n前5个最强相关性结果:\n")
  print(head(results_sorted[, c("靶基因", "失巢死亡基因", "相关性", "P值", "FDR")], 5))

  # 8. 定义筛选条件（基于FDR）
  cat("\n设置FDR阈值进行筛选...\n")

  threshold_names <- c("FDR_0.05", "FDR_0.01", "FDR_0.001", "FDR_0.0001", "FDR_0.00001")
  filtered_results <- list()

  for (threshold in threshold_names) {
    fdr_value <- as.numeric(sub("FDR_", "", threshold))
    filtered_results[[threshold]] <- results_sorted %>% 
      filter(FDR < fdr_value) %>% 
      arrange(desc(abs(相关性)))
  }

  # 9. 保存结果
  # 保存完整结果
  output_file_full <- file.path(output_dir, paste0("所有靶基因_vs_", pcd_cn, "基因_Spearman_Full_Results.csv"))
  fwrite_pcd(results_sorted, output_file_full)  # 列名中的 PCD 名称随循环改变
  cat("\n✓ 完整结果已保存至:", output_file_full, "\n")

  # 保存不同FDR阈值的结果
  cat("\n保存不同FDR阈值的结果...\n")
  for(threshold in threshold_names) {
    if (nrow(filtered_results[[threshold]]) > 0) {
      output_file <- file.path(output_dir, 
                               paste0("所有靶基因_vs_", pcd_cn, "基因_Spearman_", threshold, ".csv"))
      fwrite_pcd(filtered_results[[threshold]], output_file)
      cat("  ✓ ", threshold, ": ", nrow(filtered_results[[threshold]]), "个基因对 -> ", 
          basename(output_file), "\n")
    } else {
      cat("  ⚠ ", threshold, ": 无显著基因对\n")
    }
  }

  # 10. 生成统计摘要
  cat("\n===== 分析完成 =====\n")
  cat("靶基因数量:", length(target_gene_names), "\n")
  cat("失巢死亡基因数量:", length(death_gene_names), "\n")
  cat("分析的总基因对数量:", total_pairs, "\n")
  cat("成功计算的相关性数量:", nrow(results_valid), "\n")

  # 创建统计摘要
  summary_stats <- data.frame(
    FDR阈值 = threshold_names,
    基因对数量 = sapply(filtered_results, nrow),
    百分比 = sapply(filtered_results, function(x) {
      if (nrow(results_valid) > 0) {
        round(nrow(x)/nrow(results_valid)*100, 2)
      } else {
        0
      }
    }),
    stringsAsFactors = FALSE
  )

  # 打印到控制台
  cat("\n不同FDR阈值下的显著基因对数量:\n")
  print(summary_stats)

  # 找最强正相关和负相关基因对
  if (nrow(results_sorted) > 0) {
    cat("\n最强相关性基因对:\n")

    # 正相关
    pos_idx <- which.max(results_sorted$相关性)
    cat("  最强正相关: ", results_sorted$靶基因[pos_idx], "-", 
        results_sorted$失巢死亡基因[pos_idx], 
        " (r = ", round(results_sorted$相关性[pos_idx], 4), 
        ", FDR = ", format(results_sorted$FDR[pos_idx], scientific = TRUE, digits = 2), ")\n", sep="")

    # 负相关
    neg_idx <- which.min(results_sorted$相关性)
    cat("  最强负相关: ", results_sorted$靶基因[neg_idx], "-", 
        results_sorted$失巢死亡基因[neg_idx], 
        " (r = ", round(results_sorted$相关性[neg_idx], 4), 
        ", FDR = ", format(results_sorted$FDR[neg_idx], scientific = TRUE, digits = 2), ")\n", sep="")
  }

  # 11. 保存统计结果到文件
  stats_file <- file.path(output_dir, "所有靶基因相关性分析统计摘要.txt")
  sink(stats_file)
  cat("所有靶基因与失巢死亡基因相关性分析统计摘要\n")
  cat("==========================================\n")
  cat("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("靶基因数量:", length(target_gene_names), "\n")
  cat("失巢死亡基因数量:", length(death_gene_names), "\n")
  cat("总基因对数量:", total_pairs, "\n")
  cat("成功计算的相关性数量:", nrow(results_valid), "\n")
  cat("计算成功率:", round(nrow(results_valid)/total_pairs*100, 2), "%\n\n")

  cat("不同FDR阈值下的显著基因对数量:\n")
  cat("================================\n")
  for(i in 1:nrow(summary_stats)) {
    cat(summary_stats$FDR阈值[i], ": ", 
        summary_stats$基因对数量[i], " (", summary_stats$百分比[i], "%)\n", sep="")
  }

  if (nrow(results_sorted) > 0) {
    cat("\n最强相关性基因对:\n")
    cat("最强正相关基因对: ", results_sorted$靶基因[1], "-", results_sorted$失巢死亡基因[1], 
        "\n  相关性系数 = ", round(results_sorted$相关性[1], 4), 
        "\n  P值 = ", format(results_sorted$P值[1], scientific = TRUE, digits = 3),
        "\n  FDR = ", format(results_sorted$FDR[1], scientific = TRUE, digits = 3), "\n", sep="")

    min_cor_idx <- which.min(results_sorted$相关性)
    cat("\n最强负相关基因对: ", results_sorted$靶基因[min_cor_idx], "-", 
        results_sorted$失巢死亡基因[min_cor_idx], 
        "\n  相关性系数 = ", round(results_sorted$相关性[min_cor_idx], 4), 
        "\n  P值 = ", format(results_sorted$P值[min_cor_idx], scientific = TRUE, digits = 3),
        "\n  FDR = ", format(results_sorted$FDR[min_cor_idx], scientific = TRUE, digits = 3), "\n", sep="")
  }

  cat("\n数据质量信息:\n")
  cat("样本数量:", ncol(target_gene_matrix), "\n")
  cat("总基因对:", total_pairs, "\n")
  cat("有效相关性计算:", nrow(results_valid), "\n")
  cat("无效计算:", total_pairs - nrow(results_valid), "\n")
  sink()

  cat("\n✓ 统计摘要已保存至:", stats_file, "\n")

  # 12. 最终总结
  cat("\n=== 最终输出文件列表 ===\n")
  cat("1. 所有靶基因_vs_失巢死亡基因_Spearman_Full_Results.csv\n")
  for(threshold in threshold_names) {
    if (nrow(filtered_results[[threshold]]) > 0) {
      cat("2. 所有靶基因_vs_失巢死亡基因_Spearman_", threshold, ".csv\n", sep="")
    }
  }
  cat("3. 所有靶基因相关性分析统计摘要.txt\n")
  cat("\n输出目录:", output_dir, "\n")
  cat("分析完成！\n")

  # ---------------------------------------------------------
  # 2. 免疫细胞与 PCD 相关基因
  # ---------------------------------------------------------
  # 1. 加载必要的包
  library(data.table)
  library(dplyr)

  # 2. 设置输入和输出路径
  input_dir <- file.path(run_dir, "extracted")
  output_dir <- file.path(run_dir, "immune_cells_vs_PCD_genes")

  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("创建输出目录:", output_dir, "\n")
  }

  # 3. 读取数据
  cat("正在读取数据...\n")

  # 读取提取的免疫细胞表达矩阵 
  immune_cell_file <- immune_input
  if (!file.exists(immune_cell_file)) {
    stop("错误：找不到免疫细胞表达矩阵文件 - ", immune_cell_file)
  }
  immune_cell_data <- fread(immune_cell_file, header = TRUE)

  # 检查免疫细胞表达矩阵是否有数据
  if (nrow(immune_cell_data) == 0) {
    stop("错误：免疫细胞表达矩阵为空！")
  }

  # 提取免疫细胞样本名
  immune_cell_samples <- colnames(immune_cell_data)[-1]  # 排除第一列（免疫细胞名）
  cat("✓ 免疫细胞表达矩阵: ", nrow(immune_cell_data), "个免疫细胞, ", 
      length(immune_cell_samples), "个样本\n")

  # 读取失巢死亡基因表达矩阵
  anoikis_gene_file <- pcd_input
  if (!file.exists(anoikis_gene_file)) {
    stop("错误：找不到失巢死亡基因表达矩阵文件 - ", anoikis_gene_file)
  }
  anoikis_gene_data <- fread(anoikis_gene_file, header = TRUE)

  # 提取失巢死亡基因样本名
  anoikis_gene_samples <- colnames(anoikis_gene_data)[-1]  # 排除第一列（基因名）
  cat("✓ 失巢死亡基因表达矩阵: ", nrow(anoikis_gene_data), "个基因, ", 
      length(anoikis_gene_samples), "个样本\n")

  # 4. 样本匹配与排序
  cat("\n正在检查样本匹配情况...\n")

  # 检查样本名是否一致
  if (!identical(immune_cell_samples, anoikis_gene_samples)) {
    cat("注意：样本名不完全一致，正在检查...\n")

    # 找出共同样本
    common_samples <- intersect(immune_cell_samples, anoikis_gene_samples)

    if (length(common_samples) == 0) {
      stop("错误：免疫细胞和失巢死亡基因矩阵没有共同的样本！请检查样本名是否一致。")
    }

    cat("✓ 找到", length(common_samples), "个共同样本\n")

    # 检查是否有样本不匹配
    if (length(common_samples) < length(immune_cell_samples)) {
      missing_in_anoikis <- setdiff(immune_cell_samples, anoikis_gene_samples)
      cat("警告：有", length(missing_in_anoikis), "个免疫细胞样本在失巢死亡基因矩阵中不存在：", 
          paste(missing_in_anoikis, collapse=", "), "\n")
    }

    if (length(common_samples) < length(anoikis_gene_samples)) {
      missing_in_immune <- setdiff(anoikis_gene_samples, immune_cell_samples)
      cat("警告：有", length(missing_in_immune), "个失巢死亡基因样本在免疫细胞矩阵中不存在：", 
          paste(missing_in_immune, collapse=", "), "\n")
    }

    # 按共同样本重排数据
    cat("正在按共同样本重排数据...\n")

    # 提取基因名/免疫细胞名
    immune_cell_names <- immune_cell_data[[1]]
    anoikis_gene_names <- anoikis_gene_data[[1]]

    # 重排免疫细胞矩阵（按共同样本顺序）
    immune_cell_matrix <- as.matrix(immune_cell_data[, -1])
    immune_cell_matrix <- immune_cell_matrix[, match(common_samples, immune_cell_samples), drop = FALSE]
    colnames(immune_cell_matrix) <- common_samples

    # 重排失巢死亡基因矩阵（按共同样本顺序）
    anoikis_gene_matrix <- as.matrix(anoikis_gene_data[, -1])
    anoikis_gene_matrix <- anoikis_gene_matrix[, match(common_samples, anoikis_gene_samples), drop = FALSE]
    colnames(anoikis_gene_matrix) <- common_samples

    cat("✓ 数据重排完成，使用", length(common_samples), "个共同样本进行分析\n")

    # 保存样本匹配信息
    sample_info_file <- file.path(output_dir, "样本匹配信息.txt")
    sink(sample_info_file)
    cat("样本匹配信息\n")
    cat("=============\n")
    cat("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("免疫细胞样本数:", length(immune_cell_samples), "\n")
    cat("失巢死亡基因样本数:", length(anoikis_gene_samples), "\n")
    cat("共同样本数:", length(common_samples), "\n")
    cat("\n免疫细胞样本列表:\n")
    cat(paste(immune_cell_samples, collapse=", "), "\n")
    cat("\n失巢死亡基因样本列表:\n")
    cat(paste(anoikis_gene_samples, collapse=", "), "\n")
    cat("\n共同样本列表（实际用于分析）:\n")
    cat(paste(common_samples, collapse=", "), "\n")

    if (length(setdiff(immune_cell_samples, anoikis_gene_samples)) > 0) {
      cat("\n仅在免疫细胞中存在的样本:\n")
      cat(paste(setdiff(immune_cell_samples, anoikis_gene_samples), collapse=", "), "\n")
    }

    if (length(setdiff(anoikis_gene_samples, immune_cell_samples)) > 0) {
      cat("\n仅在失巢死亡基因中存在的样本:\n")
      cat(paste(setdiff(anoikis_gene_samples, immune_cell_samples), collapse=", "), "\n")
    }
    sink()
    cat("✓ 样本匹配信息已保存至:", sample_info_file, "\n")

  } else {
    # 样本完全一致
    cat("✓ 样本完全一致，共", length(immune_cell_samples), "个样本\n")

    # 提取基因名/免疫细胞名
    immune_cell_names <- immune_cell_data[[1]]
    anoikis_gene_names <- anoikis_gene_data[[1]]

    # 提取表达矩阵
    immune_cell_matrix <- as.matrix(immune_cell_data[, -1])
    anoikis_gene_matrix <- as.matrix(anoikis_gene_data[, -1])
  }

  # 5. 数据预处理与检查
  # 检查样本数量是否匹配（经过重排后应该一致）
  if (ncol(immune_cell_matrix) != ncol(anoikis_gene_matrix)) {
    stop("错误：免疫细胞和失巢死亡基因矩阵的样本数量不匹配！")
  }

  cat("✓ 数据检查通过：样本数量一致，共", ncol(immune_cell_matrix), "个样本\n")

  # 6. 计算斯皮尔曼相关性（所有配对）
  total_pairs <- length(immune_cell_names) * length(anoikis_gene_names)
  cat("\n正在计算所有免疫细胞与失巢死亡基因的斯皮尔曼相关性...\n")
  cat("  免疫细胞数量:", length(immune_cell_names), "\n")
  cat("  失巢死亡基因数量:", length(anoikis_gene_names), "\n")
  cat("  总共需要计算的相关性数量:", total_pairs, "\n")

  # 预分配结果矩阵（提高性能）
  results <- data.frame(
    免疫细胞 = character(total_pairs),
    失巢死亡基因 = character(total_pairs),
    相关性 = numeric(total_pairs),
    P值 = numeric(total_pairs),
    stringsAsFactors = FALSE
  )

  # 初始化进度条
  pb <- txtProgressBar(min = 0, max = total_pairs, style = 3)
  counter <- 0

  # 循环计算每个免疫细胞与每个失巢死亡基因的相关性
  for (i in 1:length(immune_cell_names)) {
    immune_cell_vector <- as.numeric(immune_cell_matrix[i, ])

    for (j in 1:length(anoikis_gene_names)) {
      anoikis_gene_exp <- as.numeric(anoikis_gene_matrix[j, ])

      # 检查并处理缺失值
      valid_data <- !(is.na(immune_cell_vector) | is.na(anoikis_gene_exp))
      n_valid <- sum(valid_data)

      if (n_valid < 3) {
        # 样本太少，无法计算相关性
        results[counter + 1, ] <- list(
          免疫细胞 = immune_cell_names[i],
          失巢死亡基因 = anoikis_gene_names[j],
          相关性 = NA,
          P值 = NA
        )
      } else {
        # 进行斯皮尔曼相关性检验
        cor_test_result <- tryCatch({
          cor.test(immune_cell_vector[valid_data], 
                   anoikis_gene_exp[valid_data], 
                   method = "spearman", exact = FALSE)
        }, error = function(e) {
          return(NULL)
        })

        if (is.null(cor_test_result)) {
          # 计算失败
          results[counter + 1, ] <- list(
            免疫细胞 = immune_cell_names[i],
            失巢死亡基因 = anoikis_gene_names[j],
            相关性 = NA,
            P值 = NA
          )
        } else {
          # 计算成功
          results[counter + 1, ] <- list(
            免疫细胞 = immune_cell_names[i],
            失巢死亡基因 = anoikis_gene_names[j],
            相关性 = cor_test_result$estimate,
            P值 = cor_test_result$p.value
          )
        }
      }

      # 更新进度条
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
    }
  }

  close(pb)
  cat("\n✓ 相关性计算完成！\n")

  # 移除无法计算相关性的行
  results_valid <- results[!is.na(results$相关性), ]
  cat("  成功计算的相关性数量:", nrow(results_valid), "/", total_pairs, "\n")

  # 7. 多重检验校正 (FDR)
  cat("正在进行多重检验校正 (FDR)...\n")
  results_valid$FDR <- p.adjust(results_valid$P值, method = "BH")

  # 检查FDR计算
  na_fdr_count <- sum(is.na(results_valid$FDR))
  if (na_fdr_count > 0) {
    cat("警告: FDR计算有", na_fdr_count, "个NA值\n")
  } else {
    cat("✓ FDR计算完成\n")
  }

  # 8. 按相关性系数绝对值排序
  results_sorted <- results_valid[order(-abs(results_valid$相关性)), ]

  # 显示前几行结果
  cat("\n前5个最强相关性结果:\n")
  print(head(results_sorted[, c("免疫细胞", "失巢死亡基因", "相关性", "P值", "FDR")], 5))

  # 9. 定义筛选条件（基于FDR）
  cat("\n设置FDR阈值进行筛选...\n")

  threshold_names <- c("FDR_0.05", "FDR_0.01", "FDR_0.001", "FDR_0.0001", "FDR_0.00001")
  filtered_results <- list()

  for (threshold in threshold_names) {
    fdr_value <- as.numeric(sub("FDR_", "", threshold))
    filtered_results[[threshold]] <- results_sorted %>% 
      filter(FDR < fdr_value) %>% 
      arrange(desc(abs(相关性)))
  }

  # 10. 保存结果
  # 保存完整结果
  output_file_full <- file.path(output_dir, paste0("所有免疫细胞_vs_", pcd_cn, "基因_Spearman_Full_Results.csv"))
  fwrite_pcd(results_sorted, output_file_full)  # 列名中的 PCD 名称随循环改变
  cat("\n✓ 完整结果已保存至:", output_file_full, "\n")

  # 保存不同FDR阈值的结果
  cat("\n保存不同FDR阈值的结果...\n")
  for(threshold in threshold_names) {
    if (nrow(filtered_results[[threshold]]) > 0) {
      output_file <- file.path(output_dir, 
                               paste0("所有免疫细胞_vs_", pcd_cn, "基因_Spearman_", threshold, ".csv"))
      fwrite_pcd(filtered_results[[threshold]], output_file)
      cat("  ✓ ", threshold, ": ", nrow(filtered_results[[threshold]]), "个配对 -> ", 
          basename(output_file), "\n")
    } else {
      cat("  ⚠ ", threshold, ": 无显著配对\n")
    }
  }

  # 11. 生成统计摘要
  cat("\n===== 分析完成 =====\n")
  cat("免疫细胞数量:", length(immune_cell_names), "\n")
  cat("失巢死亡基因数量:", length(anoikis_gene_names), "\n")
  cat("分析的总配对数量:", total_pairs, "\n")
  cat("成功计算的相关性数量:", nrow(results_valid), "\n")

  # 创建统计摘要
  summary_stats <- data.frame(
    FDR阈值 = threshold_names,
    配对数量 = sapply(filtered_results, nrow),
    百分比 = sapply(filtered_results, function(x) {
      if (nrow(results_valid) > 0) {
        round(nrow(x)/nrow(results_valid)*100, 2)
      } else {
        0
      }
    }),
    stringsAsFactors = FALSE
  )

  # 打印到控制台
  cat("\n不同FDR阈值下的显著配对数量:\n")
  print(summary_stats)

  # 找最强正相关和负相关配对
  if (nrow(results_sorted) > 0) {
    cat("\n最强相关性配对:\n")

    # 正相关
    pos_idx <- which.max(results_sorted$相关性)
    cat("  最强正相关: ", results_sorted$免疫细胞[pos_idx], "-", 
        results_sorted$失巢死亡基因[pos_idx], 
        " (r = ", round(results_sorted$相关性[pos_idx], 4), 
        ", FDR = ", format(results_sorted$FDR[pos_idx], scientific = TRUE, digits = 2), ")\n", sep="")

    # 负相关
    neg_idx <- which.min(results_sorted$相关性)
    cat("  最强负相关: ", results_sorted$免疫细胞[neg_idx], "-", 
        results_sorted$失巢死亡基因[neg_idx], 
        " (r = ", round(results_sorted$相关性[neg_idx], 4), 
        ", FDR = ", format(results_sorted$FDR[neg_idx], scientific = TRUE, digits = 2), ")\n", sep="")
  }

  # 12. 保存统计结果到文件
  stats_file <- file.path(output_dir, "所有免疫细胞相关性分析统计摘要.txt")
  sink(stats_file)
  cat("所有免疫细胞与失巢死亡基因相关性分析统计摘要\n")
  cat("==========================================\n")
  cat("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("免疫细胞数量:", length(immune_cell_names), "\n")
  cat("失巢死亡基因数量:", length(anoikis_gene_names), "\n")
  cat("总配对数量:", total_pairs, "\n")
  cat("成功计算的相关性数量:", nrow(results_valid), "\n")
  cat("计算成功率:", round(nrow(results_valid)/total_pairs*100, 2), "%\n\n")

  cat("样本匹配信息:\n")
  cat("  免疫细胞样本数:", length(immune_cell_samples), "\n")
  cat("  失巢死亡基因样本数:", length(anoikis_gene_samples), "\n")
  if (exists("common_samples")) {
    cat("  共同样本数:", length(common_samples), "\n")
  } else {
    cat("  共同样本数:", length(immune_cell_samples), "（完全匹配）\n")
  }
  cat("\n不同FDR阈值下的显著配对数量:\n")
  cat("================================\n")
  for(i in 1:nrow(summary_stats)) {
    cat(summary_stats$FDR阈值[i], ": ", 
        summary_stats$配对数量[i], " (", summary_stats$百分比[i], "%)\n", sep="")
  }

  if (nrow(results_sorted) > 0) {
    cat("\n最强相关性配对:\n")
    cat("最强正相关配对: ", results_sorted$免疫细胞[1], "-", results_sorted$失巢死亡基因[1], 
        "\n  相关性系数 = ", round(results_sorted$相关性[1], 4), 
        "\n  P值 = ", format(results_sorted$P值[1], scientific = TRUE, digits = 3),
        "\n  FDR = ", format(results_sorted$FDR[1], scientific = TRUE, digits = 3), "\n", sep="")

    min_cor_idx <- which.min(results_sorted$相关性)
    cat("\n最强负相关配对: ", results_sorted$免疫细胞[min_cor_idx], "-", 
        results_sorted$失巢死亡基因[min_cor_idx], 
        "\n  相关性系数 = ", round(results_sorted$相关性[min_cor_idx], 4), 
        "\n  P值 = ", format(results_sorted$P值[min_cor_idx], scientific = TRUE, digits = 3),
        "\n  FDR = ", format(results_sorted$FDR[min_cor_idx], scientific = TRUE, digits = 3), "\n", sep="")
  }

  cat("\n数据质量信息:\n")
  cat("  最终使用样本数量:", ncol(immune_cell_matrix), "\n")
  cat("  总配对:", total_pairs, "\n")
  cat("  有效相关性计算:", nrow(results_valid), "\n")
  cat("  无效计算:", total_pairs - nrow(results_valid), "\n")
  sink()

  cat("\n✓ 统计摘要已保存至:", stats_file, "\n")

  # 13. 最终总结
  cat("\n=== 最终输出文件列表 ===\n")
  cat("1. 样本匹配信息.txt\n")
  cat("2. 所有免疫细胞_vs_失巢死亡基因_Spearman_Full_Results.csv\n")
  for(threshold in threshold_names) {
    if (nrow(filtered_results[[threshold]]) > 0) {
      cat("3. 所有免疫细胞_vs_失巢死亡基因_Spearman_", threshold, ".csv\n", sep="")
    }
  }
  cat("4. 所有免疫细胞相关性分析统计摘要.txt\n")
  cat("\n输出目录:", output_dir, "\n")
  cat("分析完成！\n")

  # ---------------------------------------------------------
  # 3. 候选药物与 PCD 相关基因
  # ---------------------------------------------------------
  # 1. 加载必要的包
  library(data.table)
  library(dplyr)

  # 2. 设置输入和输出路径
  input_dir <- file.path(run_dir, "extracted")
  output_dir <- file.path(run_dir, "drugs_vs_PCD_genes")

  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("创建输出目录:", output_dir, "\n")
  }

  # 3. 读取数据
  cat("正在读取数据...\n")

  # 读取提取的药物表达矩阵 
  drug_file <- drug_input
  if (!file.exists(drug_file)) {
    stop("错误：找不到药物表达矩阵文件 - ", drug_file)
  }
  drug_data <- fread(drug_file, header = TRUE)

  # 检查药物表达矩阵是否有数据
  if (nrow(drug_data) == 0) {
    stop("错误：药物表达矩阵为空！")
  }

  # 提取药物样本名
  drug_samples <- colnames(drug_data)[-1]  # 排除第一列（药物名）
  cat("✓ 药物表达矩阵: ", nrow(drug_data), "种药物, ", 
      length(drug_samples), "个样本\n")

  # 读取失巢死亡基因表达矩阵
  anoikis_gene_file <- pcd_input
  if (!file.exists(anoikis_gene_file)) {
    stop("错误：找不到失巢死亡基因表达矩阵文件 - ", anoikis_gene_file)
  }
  anoikis_gene_data <- fread(anoikis_gene_file, header = TRUE)

  # 提取失巢死亡基因样本名
  anoikis_gene_samples <- colnames(anoikis_gene_data)[-1]  # 排除第一列（基因名）
  cat("✓ 失巢死亡基因表达矩阵: ", nrow(anoikis_gene_data), "个基因, ", 
      length(anoikis_gene_samples), "个样本\n")

  # 4. 样本匹配与排序
  cat("\n正在检查样本匹配情况...\n")

  # 检查样本名是否一致
  if (!identical(drug_samples, anoikis_gene_samples)) {
    cat("注意：样本名不完全一致，正在检查...\n")

    # 找出共同样本
    common_samples <- intersect(drug_samples, anoikis_gene_samples)

    if (length(common_samples) == 0) {
      stop("错误：药物和失巢死亡基因矩阵没有共同的样本！请检查样本名是否一致。")
    }

    cat("✓ 找到", length(common_samples), "个共同样本\n")

    # 检查是否有样本不匹配
    if (length(common_samples) < length(drug_samples)) {
      missing_in_anoikis <- setdiff(drug_samples, anoikis_gene_samples)
      cat("警告：有", length(missing_in_anoikis), "个药物样本在失巢死亡基因矩阵中不存在：", 
          paste(missing_in_anoikis, collapse=", "), "\n")
    }

    if (length(common_samples) < length(anoikis_gene_samples)) {
      missing_in_drug <- setdiff(anoikis_gene_samples, drug_samples)
      cat("警告：有", length(missing_in_drug), "个失巢死亡基因样本在药物矩阵中不存在：", 
          paste(missing_in_drug, collapse=", "), "\n")
    }

    # 按共同样本重排数据
    cat("正在按共同样本重排数据...\n")

    # 提取基因名/药物名
    drug_names <- drug_data[[1]]
    anoikis_gene_names <- anoikis_gene_data[[1]]

    # 重排药物矩阵（按共同样本顺序）
    drug_matrix <- as.matrix(drug_data[, -1])
    drug_matrix <- drug_matrix[, match(common_samples, drug_samples), drop = FALSE]
    colnames(drug_matrix) <- common_samples

    # 重排失巢死亡基因矩阵（按共同样本顺序）
    anoikis_gene_matrix <- as.matrix(anoikis_gene_data[, -1])
    anoikis_gene_matrix <- anoikis_gene_matrix[, match(common_samples, anoikis_gene_samples), drop = FALSE]
    colnames(anoikis_gene_matrix) <- common_samples

    cat("✓ 数据重排完成，使用", length(common_samples), "个共同样本进行分析\n")

    # 保存样本匹配信息
    sample_info_file <- file.path(output_dir, "样本匹配信息.txt")
    sink(sample_info_file)
    cat("样本匹配信息\n")
    cat("=============\n")
    cat("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("药物样本数:", length(drug_samples), "\n")
    cat("失巢死亡基因样本数:", length(anoikis_gene_samples), "\n")
    cat("共同样本数:", length(common_samples), "\n")
    cat("\n药物样本列表:\n")
    cat(paste(drug_samples, collapse=", "), "\n")
    cat("\n失巢死亡基因样本列表:\n")
    cat(paste(anoikis_gene_samples, collapse=", "), "\n")
    cat("\n共同样本列表（实际用于分析）:\n")
    cat(paste(common_samples, collapse=", "), "\n")

    if (length(setdiff(drug_samples, anoikis_gene_samples)) > 0) {
      cat("\n仅在药物中存在的样本:\n")
      cat(paste(setdiff(drug_samples, anoikis_gene_samples), collapse=", "), "\n")
    }

    if (length(setdiff(anoikis_gene_samples, drug_samples)) > 0) {
      cat("\n仅在失巢死亡基因中存在的样本:\n")
      cat(paste(setdiff(anoikis_gene_samples, drug_samples), collapse=", "), "\n")
    }
    sink()
    cat("✓ 样本匹配信息已保存至:", sample_info_file, "\n")

  } else {
    # 样本完全一致
    cat("✓ 样本完全一致，共", length(drug_samples), "个样本\n")

    # 提取基因名/药物名
    drug_names <- drug_data[[1]]
    anoikis_gene_names <- anoikis_gene_data[[1]]

    # 提取表达矩阵
    drug_matrix <- as.matrix(drug_data[, -1])
    anoikis_gene_matrix <- as.matrix(anoikis_gene_data[, -1])
  }

  # 5. 数据预处理与检查
  # 检查样本数量是否匹配（经过重排后应该一致）
  if (ncol(drug_matrix) != ncol(anoikis_gene_matrix)) {
    stop("错误：药物和失巢死亡基因矩阵的样本数量不匹配！")
  }

  cat("✓ 数据检查通过：样本数量一致，共", ncol(drug_matrix), "个样本\n")

  # 6. 计算斯皮尔曼相关性（所有配对）
  total_pairs <- length(drug_names) * length(anoikis_gene_names)
  cat("\n正在计算所有药物与失巢死亡基因的斯皮尔曼相关性...\n")
  cat("  药物数量:", length(drug_names), "\n")
  cat("  失巢死亡基因数量:", length(anoikis_gene_names), "\n")
  cat("  总共需要计算的相关性数量:", total_pairs, "\n")

  # 预分配结果矩阵（提高性能）
  results <- data.frame(
    药物 = character(total_pairs),
    失巢死亡基因 = character(total_pairs),
    相关性 = numeric(total_pairs),
    P值 = numeric(total_pairs),
    stringsAsFactors = FALSE
  )

  # 初始化进度条
  pb <- txtProgressBar(min = 0, max = total_pairs, style = 3)
  counter <- 0

  # 循环计算每个药物与每个失巢死亡基因的相关性
  for (i in 1:length(drug_names)) {
    drug_vector <- as.numeric(drug_matrix[i, ])

    for (j in 1:length(anoikis_gene_names)) {
      anoikis_gene_exp <- as.numeric(anoikis_gene_matrix[j, ])

      # 检查并处理缺失值
      valid_data <- !(is.na(drug_vector) | is.na(anoikis_gene_exp))
      n_valid <- sum(valid_data)

      if (n_valid < 3) {
        # 样本太少，无法计算相关性
        results[counter + 1, ] <- list(
          药物 = drug_names[i],
          失巢死亡基因 = anoikis_gene_names[j],
          相关性 = NA,
          P值 = NA
        )
      } else {
        # 进行斯皮尔曼相关性检验
        cor_test_result <- tryCatch({
          cor.test(drug_vector[valid_data], 
                   anoikis_gene_exp[valid_data], 
                   method = "spearman", exact = FALSE)
        }, error = function(e) {
          return(NULL)
        })

        if (is.null(cor_test_result)) {
          # 计算失败
          results[counter + 1, ] <- list(
            药物 = drug_names[i],
            失巢死亡基因 = anoikis_gene_names[j],
            相关性 = NA,
            P值 = NA
          )
        } else {
          # 计算成功
          results[counter + 1, ] <- list(
            药物 = drug_names[i],
            失巢死亡基因 = anoikis_gene_names[j],
            相关性 = cor_test_result$estimate,
            P值 = cor_test_result$p.value
          )
        }
      }

      # 更新进度条
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
    }
  }

  close(pb)
  cat("\n✓ 相关性计算完成！\n")

  # 移除无法计算相关性的行
  results_valid <- results[!is.na(results$相关性), ]
  cat("  成功计算的相关性数量:", nrow(results_valid), "/", total_pairs, "\n")

  # 7. 多重检验校正 (FDR)
  cat("正在进行多重检验校正 (FDR)...\n")
  results_valid$FDR <- p.adjust(results_valid$P值, method = "BH")

  # 检查FDR计算
  na_fdr_count <- sum(is.na(results_valid$FDR))
  if (na_fdr_count > 0) {
    cat("警告: FDR计算有", na_fdr_count, "个NA值\n")
  } else {
    cat("✓ FDR计算完成\n")
  }

  # 8. 按相关性系数绝对值排序
  results_sorted <- results_valid[order(-abs(results_valid$相关性)), ]

  # 显示前几行结果
  cat("\n前5个最强相关性结果:\n")
  print(head(results_sorted[, c("药物", "失巢死亡基因", "相关性", "P值", "FDR")], 5))

  # 9. 定义筛选条件（基于FDR）
  cat("\n设置FDR阈值进行筛选...\n")

  threshold_names <- c("FDR_0.05", "FDR_0.01", "FDR_0.001", "FDR_0.0001", "FDR_0.00001")
  filtered_results <- list()

  for (threshold in threshold_names) {
    fdr_value <- as.numeric(sub("FDR_", "", threshold))
    filtered_results[[threshold]] <- results_sorted %>% 
      filter(FDR < fdr_value) %>% 
      arrange(desc(abs(相关性)))
  }

  # 10. 保存结果
  # 保存完整结果
  output_file_full <- file.path(output_dir, paste0("所有药物_vs_", pcd_cn, "基因_Spearman_Full_Results.csv"))
  fwrite_pcd(results_sorted, output_file_full)  # 列名中的 PCD 名称随循环改变
  cat("\n✓ 完整结果已保存至:", output_file_full, "\n")

  # 保存不同FDR阈值的结果
  cat("\n保存不同FDR阈值的结果...\n")
  for(threshold in threshold_names) {
    if (nrow(filtered_results[[threshold]]) > 0) {
      output_file <- file.path(output_dir, 
                               paste0("所有药物_vs_", pcd_cn, "基因_Spearman_", threshold, ".csv"))
      fwrite_pcd(filtered_results[[threshold]], output_file)
      cat("  ✓ ", threshold, ": ", nrow(filtered_results[[threshold]]), "个配对 -> ", 
          basename(output_file), "\n")
    } else {
      cat("  ⚠ ", threshold, ": 无显著配对\n")
    }
  }

  # 11. 生成统计摘要
  cat("\n===== 分析完成 =====\n")
  cat("药物数量:", length(drug_names), "\n")
  cat("失巢死亡基因数量:", length(anoikis_gene_names), "\n")
  cat("分析的总配对数量:", total_pairs, "\n")
  cat("成功计算的相关性数量:", nrow(results_valid), "\n")

  # 创建统计摘要
  summary_stats <- data.frame(
    FDR阈值 = threshold_names,
    配对数量 = sapply(filtered_results, nrow),
    百分比 = sapply(filtered_results, function(x) {
      if (nrow(results_valid) > 0) {
        round(nrow(x)/nrow(results_valid)*100, 2)
      } else {
        0
      }
    }),
    stringsAsFactors = FALSE
  )

  # 打印到控制台
  cat("\n不同FDR阈值下的显著配对数量:\n")
  print(summary_stats)

  # 找最强正相关和负相关配对
  if (nrow(results_sorted) > 0) {
    cat("\n最强相关性配对:\n")

    # 正相关
    pos_idx <- which.max(results_sorted$相关性)
    cat("  最强正相关: ", results_sorted$药物[pos_idx], "-", 
        results_sorted$失巢死亡基因[pos_idx], 
        " (r = ", round(results_sorted$相关性[pos_idx], 4), 
        ", FDR = ", format(results_sorted$FDR[pos_idx], scientific = TRUE, digits = 2), ")\n", sep="")

    # 负相关
    neg_idx <- which.min(results_sorted$相关性)
    cat("  最强负相关: ", results_sorted$药物[neg_idx], "-", 
        results_sorted$失巢死亡基因[neg_idx], 
        " (r = ", round(results_sorted$相关性[neg_idx], 4), 
        ", FDR = ", format(results_sorted$FDR[neg_idx], scientific = TRUE, digits = 2), ")\n", sep="")
  }

  # 12. 保存统计结果到文件
  stats_file <- file.path(output_dir, "所有药物相关性分析统计摘要.txt")
  sink(stats_file)
  cat("所有药物与失巢死亡基因相关性分析统计摘要\n")
  cat("==========================================\n")
  cat("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("药物数量:", length(drug_names), "\n")
  cat("失巢死亡基因数量:", length(anoikis_gene_names), "\n")
  cat("总配对数量:", total_pairs, "\n")
  cat("成功计算的相关性数量:", nrow(results_valid), "\n")
  cat("计算成功率:", round(nrow(results_valid)/total_pairs*100, 2), "%\n\n")

  cat("样本匹配信息:\n")
  cat("  药物样本数:", length(drug_samples), "\n")
  cat("  失巢死亡基因样本数:", length(anoikis_gene_samples), "\n")
  if (exists("common_samples")) {
    cat("  共同样本数:", length(common_samples), "\n")
  } else {
    cat("  共同样本数:", length(drug_samples), "（完全匹配）\n")
  }
  cat("\n不同FDR阈值下的显著配对数量:\n")
  cat("================================\n")
  for(i in 1:nrow(summary_stats)) {
    cat(summary_stats$FDR阈值[i], ": ", 
        summary_stats$配对数量[i], " (", summary_stats$百分比[i], "%)\n", sep="")
  }

  if (nrow(results_sorted) > 0) {
    cat("\n最强相关性配对:\n")
    cat("最强正相关配对: ", results_sorted$药物[1], "-", results_sorted$失巢死亡基因[1], 
        "\n  相关性系数 = ", round(results_sorted$相关性[1], 4), 
        "\n  P值 = ", format(results_sorted$P值[1], scientific = TRUE, digits = 3),
        "\n  FDR = ", format(results_sorted$FDR[1], scientific = TRUE, digits = 3), "\n", sep="")

    min_cor_idx <- which.min(results_sorted$相关性)
    cat("\n最强负相关配对: ", results_sorted$药物[min_cor_idx], "-", 
        results_sorted$失巢死亡基因[min_cor_idx], 
        "\n  相关性系数 = ", round(results_sorted$相关性[min_cor_idx], 4), 
        "\n  P值 = ", format(results_sorted$P值[min_cor_idx], scientific = TRUE, digits = 3),
        "\n  FDR = ", format(results_sorted$FDR[min_cor_idx], scientific = TRUE, digits = 3), "\n", sep="")
  }

  cat("\n数据质量信息:\n")
  cat("  最终使用样本数量:", ncol(drug_matrix), "\n")
  cat("  总配对:", total_pairs, "\n")
  cat("  有效相关性计算:", nrow(results_valid), "\n")
  cat("  无效计算:", total_pairs - nrow(results_valid), "\n")
  sink()

  cat("\n✓ 统计摘要已保存至:", stats_file, "\n")

  # 13. 最终总结
  cat("\n=== 最终输出文件列表 ===\n")
  cat("1. 样本匹配信息.txt\n")
  cat("2. 所有药物_vs_失巢死亡基因_Spearman_Full_Results.csv\n")
  for(threshold in threshold_names) {
    if (nrow(filtered_results[[threshold]]) > 0) {
      cat("3. 所有药物_vs_失巢死亡基因_Spearman_", threshold, ".csv\n", sep="")
    }
  }
  cat("4. 所有药物相关性分析统计摘要.txt\n")
  cat("\n输出目录:", output_dir, "\n")
  cat("分析完成！\n")

  # ---------------------------------------------------------
  # 4. 免疫治疗相关基因与 PCD 相关基因
  # ---------------------------------------------------------
  # 1. 加载必要的包
  library(data.table)
  library(dplyr)

  # 2. 设置输入和输出路径
  input_dir <- file.path(run_dir, "extracted")
  output_dir <- file.path(run_dir, "immunotherapy_genes_vs_PCD_genes")

  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 3. 读取数据
  # 读取提取免疫治疗基因表达矩阵 
  immune_gene_file <- immuno_input
  immune_gene_data <- fread(immune_gene_file, header = TRUE)

  # 检查免疫治疗基因表达矩阵是否有数据
  if (nrow(immune_gene_data) == 0) {
    stop("错误：免疫治疗基因表达矩阵为空！")
  }

  cat("免疫治疗基因表达矩阵基因数:", nrow(immune_gene_data), "\n")
  cat("免疫治疗基因表达矩阵样本数:", ncol(immune_gene_data) - 1, "\n")

  # 读取失巢死亡基因表达矩阵
  death_gene_file <- pcd_input
  death_gene_data <- fread(death_gene_file, header = TRUE)

  cat("失巢死亡基因表达矩阵基因数:", nrow(death_gene_data), "\n")
  cat("失巢死亡基因表达矩阵样本数:", ncol(death_gene_data) - 1, "\n")

  # 4. 自动样本对齐函数
  align_samples_auto <- function(gene_data1, gene_data2) {
    # 提取样本名称
    samples1 <- colnames(gene_data1)[-1]  # 排除第一列（基因名）
    samples2 <- colnames(gene_data2)[-1]  # 排除第一列（基因名）

    cat("数据集1样本数:", length(samples1), "\n")
    cat("数据集2样本数:", length(samples2), "\n")

    # 找共同样本
    common_samples <- intersect(samples1, samples2)

    if (length(common_samples) == 0) {
      stop("错误：没有共同的样本！无法进行分析。")
    }

    cat("找到共同样本数:", length(common_samples), "\n")
    cat("数据集1特有样本数:", length(setdiff(samples1, samples2)), "\n")
    cat("数据集2特有样本数:", length(setdiff(samples2, samples1)), "\n")

    if (length(common_samples) < max(length(samples1), length(samples2))) {
      cat("警告：部分样本不匹配，将只使用共同样本进行分析\n")
    }

    # 按共同样本重新排列
    gene1_index <- match(common_samples, samples1)
    gene2_index <- match(common_samples, samples2)

    # 提取对齐后的数据（保留第一列基因名）
    aligned_gene1 <- gene_data1[, c(1, gene1_index + 1), with = FALSE]
    aligned_gene2 <- gene_data2[, c(1, gene2_index + 1), with = FALSE]

    # 更新列名
    colnames(aligned_gene1) <- c(colnames(gene_data1)[1], common_samples)
    colnames(aligned_gene2) <- c(colnames(gene_data2)[1], common_samples)

    return(list(
      gene1 = aligned_gene1,
      gene2 = aligned_gene2,
      common_samples = common_samples
    ))
  }

  # 5. 执行自动样本对齐
  cat("\n=== 执行自动样本对齐 ===\n")
  aligned_data <- align_samples_auto(immune_gene_data, death_gene_data)

  immune_gene_data <- aligned_data$gene1
  death_gene_data <- aligned_data$gene2
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

  # 6. 数据预处理与检查（使用对齐后的数据）
  # 提取基因名
  immune_gene_names <- immune_gene_data[[1]]
  death_gene_names <- death_gene_data[[1]]

  # 提取表达矩阵（转换为数值矩阵）
  immune_gene_matrix <- as.matrix(immune_gene_data[, -1])
  death_gene_matrix <- as.matrix(death_gene_data[, -1])

  # 检查样本是否匹配（对齐后应该完全一致）
  if (ncol(immune_gene_matrix) != ncol(death_gene_matrix)) {
    stop("错误：样本对齐后数量仍然不匹配！")
  } else {
    cat("数据检查通过：样本数量一致，共", ncol(immune_gene_matrix), "个样本。\n")
  }

  # 检查样本名称是否一致
  if (!identical(colnames(immune_gene_data)[-1], colnames(death_gene_data)[-1])) {
    warning("警告：样本名称不完全一致，但数量相同")
  } else {
    cat("样本名称完全一致！\n")
  }

  # 7. 计算斯皮尔曼相关性（所有基因对）
  # 初始化结果数据框
  results <- data.frame(
    免疫治疗基因 = character(),
    失巢死亡基因 = character(),
    相关性 = numeric(),
    P值 = numeric(),
    stringsAsFactors = FALSE
  )

  cat("正在计算所有免疫治疗基因与失巢死亡基因的斯皮尔曼相关性...\n")
  cat("免疫治疗基因数量:", length(immune_gene_names), "\n")
  cat("失巢死亡基因数量:", length(death_gene_names), "\n")
  cat("总共需要计算的相关性数量:", length(immune_gene_names) * length(death_gene_names), "\n")
  cat("使用的样本数量:", ncol(immune_gene_matrix), "\n")

  # 初始化进度条
  pb <- txtProgressBar(min = 0, max = length(immune_gene_names) * length(death_gene_names), style = 3)

  counter <- 0

  # 循环计算每个免疫治疗基因与每个失巢死亡基因的相关性
  for (i in 1:length(immune_gene_names)) {
    immune_gene_vector <- as.numeric(immune_gene_matrix[i, ]) # 获取第i个免疫治疗基因的表达向量

    for (j in 1:length(death_gene_names)) {
      death_gene_exp <- as.numeric(death_gene_matrix[j, ]) # 获取第j个失巢死亡基因的表达向量

      # 进行斯皮尔曼相关性检验
      cor_test_result <- cor.test(immune_gene_vector, death_gene_exp, method = "spearman", exact = FALSE)

      # 将结果添加到数据框
      results <- rbind(results, data.frame(
        免疫治疗基因 = immune_gene_names[i],
        失巢死亡基因 = death_gene_names[j],
        相关性 = cor_test_result$estimate,
        P值 = cor_test_result$p.value,
        stringsAsFactors = FALSE
      ))

      # 更新进度条
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
    }
  }

  close(pb)
  cat("相关性计算完成！\n")

  # 8. 多重检验校正 (FDR) - 对所有P值一起进行校正
  cat("正在进行多重检验校正 (FDR)...\n")
  results$FDR <- p.adjust(results$P值, method = "BH")

  # 检查FDR计算是否成功
  if (all(is.na(results$FDR))) {
    cat("警告: FDR计算返回NA值，可能是由于P值为NA\n")
    cat("P值为NA的数量:", sum(is.na(results$P值)), "\n")
  } else {
    cat("FDR计算完成，FDR范围:", range(results$FDR, na.rm = TRUE), "\n")
  }

  # 9. 按相关性系数绝对值排序
  results_sorted <- results[order(-abs(results$相关性)), ]

  # 显示前几行结果（包含FDR列）
  cat("\n前5个最强相关性结果（包含FDR列）:\n")
  print(head(results_sorted[, c("免疫治疗基因", "失巢死亡基因", "相关性", "P值", "FDR")], 5))

  # 10. 定义筛选条件（基于FDR）- 只保留指定的阈值
  cat("\n设置FDR阈值进行筛选...\n")

  # 只保留指定的阈值
  filtered_results <- list(
    FDR_0.05 = results_sorted %>% filter(FDR < 0.05) %>% arrange(desc(abs(相关性))),
    FDR_0.01 = results_sorted %>% filter(FDR < 0.01) %>% arrange(desc(abs(相关性))),
    FDR_0.001 = results_sorted %>% filter(FDR < 0.001) %>% arrange(desc(abs(相关性))),
    FDR_0.0001 = results_sorted %>% filter(FDR < 0.0001) %>% arrange(desc(abs(相关性))),
    FDR_0.00001 = results_sorted %>% filter(FDR < 0.00001) %>% arrange(desc(abs(相关性)))
  )

  # 11. 保存结果
  # 保存完整结果
  output_file_full <- file.path(output_dir, paste0("免疫治疗基因_vs_", pcd_cn, "基因_Spearman_Full_Results.csv"))
  fwrite_pcd(results_sorted, output_file_full)  # 列名中的 PCD 名称随循环改变
  cat("\n完整结果已保存至:", output_file_full, "\n")
  cat("文件包含以下列:", paste(colnames(results_sorted), collapse = ", "), "\n")

  # 保存不同FDR阈值的结果
  cat("\n保存不同FDR阈值的结果...\n")
  threshold_names <- c("FDR_0.05", "FDR_0.01", "FDR_0.001", "FDR_0.0001", "FDR_0.00001")
  for(threshold in threshold_names) {
    output_file <- file.path(output_dir, paste0("免疫治疗基因_vs_", pcd_cn, "基因_Spearman_", threshold, ".csv"))
    fwrite_pcd(filtered_results[[threshold]], output_file)
    cat("已保存", threshold, "结果:", output_file, 
        " (", nrow(filtered_results[[threshold]]), "个基因对)\n")
  }

  # 12. 生成统计摘要并保存到文件
  cat("\n===== 分析完成 =====\n")
  cat("免疫治疗基因数量:", length(immune_gene_names), "\n")
  cat("失巢死亡基因数量:", length(death_gene_names), "\n")
  cat("分析的总基因对数量:", nrow(results), "\n")
  cat("使用的样本数量:", length(common_samples), "\n")

  # 创建统计摘要
  summary_stats <- data.frame(
    FDR阈值 = threshold_names,
    基因对数量 = sapply(filtered_results, nrow),
    百分比 = sapply(filtered_results, function(x) round(nrow(x)/nrow(results)*100, 2)),
    stringsAsFactors = FALSE
  )

  # 打印到控制台
  cat("\n不同FDR阈值下的显著基因对数量:\n")
  print(summary_stats)

  # 找最强正相关和负相关基因对
  if (nrow(results_sorted) > 0) {
    cat("\n最强正相关基因对:", 
        results_sorted$免疫治疗基因[1], "-", results_sorted$失巢死亡基因[1], 
        "(r =", round(results_sorted$相关性[1], 4), 
        ", FDR =", format(results_sorted$FDR[1], scientific = TRUE, digits = 2), ")\n")

    min_cor_index <- which.min(results_sorted$相关性)
    cat("最强负相关基因对:", 
        results_sorted$免疫治疗基因[min_cor_index], "-", results_sorted$失巢死亡基因[min_cor_index], 
        "(r =", round(results_sorted$相关性[min_cor_index], 4), 
        ", FDR =", format(results_sorted$FDR[min_cor_index], scientific = TRUE, digits = 2), ")\n")
  }

  # 保存统计结果到文件
  stats_file <- file.path(output_dir, "免疫治疗基因相关性分析统计摘要.txt")
  sink(stats_file)
  cat("免疫治疗基因与失巢死亡基因相关性分析统计摘要\n")
  cat("============================================\n")
  cat("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("免疫治疗基因数量:", length(immune_gene_names), "\n")
  cat("失巢死亡基因数量:", length(death_gene_names), "\n")
  cat("使用的样本数量:", length(common_samples), "\n")
  cat("分析的总基因对数量:", nrow(results), "\n\n")

  cat("样本对齐信息:\n")
  cat("共同样本数:", length(common_samples), "\n")
  if (length(common_samples) < 10) {
    cat("样本列表:", paste(common_samples, collapse = ", "), "\n")
  } else {
    cat("前10个样本:", paste(common_samples[1:10], collapse = ", "), "...\n")
  }

  cat("\n不同FDR阈值下的显著基因对数量:\n")
  cat("================================\n")
  for(i in 1:nrow(summary_stats)) {
    cat(summary_stats$FDR阈值[i], ":", 
        summary_stats$基因对数量[i], "(", summary_stats$百分比[i], "%)\n")
  }

  cat("\n最强相关性基因对:\n")
  if (nrow(results_sorted) > 0) {
    cat("最强正相关基因对:", results_sorted$免疫治疗基因[1], "-", results_sorted$失巢死亡基因[1], 
        "\n  相关性系数 =", round(results_sorted$相关性[1], 4), 
        "\n  P值 =", format(results_sorted$P值[1], scientific = TRUE, digits = 3),
        "\n  FDR =", format(results_sorted$FDR[1], scientific = TRUE, digits = 3), "\n")

    min_cor_index <- which.min(results_sorted$相关性)
    cat("\n最强负相关基因对:", results_sorted$免疫治疗基因[min_cor_index], "-", results_sorted$失巢死亡基因[min_cor_index], 
        "\n  相关性系数 =", round(results_sorted$相关性[min_cor_index], 4), 
        "\n  P值 =", format(results_sorted$P值[min_cor_index], scientific = TRUE, digits = 3),
        "\n  FDR =", format(results_sorted$FDR[min_cor_index], scientific = TRUE, digits = 3), "\n")
  }

  cat("\n前10个最强相关性基因对:\n")
  print(head(results_sorted[, c("免疫治疗基因", "失巢死亡基因", "相关性", "P值", "FDR")], 10))

  # 添加FDR分布信息
  cat("\nFDR分布:\n")
  fdr_summary <- summary(results_sorted$FDR)
  cat("最小值:", fdr_summary[1], "\n")
  cat("第一四分位数:", fdr_summary[2], "\n")
  cat("中位数:", fdr_summary[3], "\n")
  cat("平均值:", fdr_summary[4], "\n")
  cat("第三四分位数:", fdr_summary[5], "\n")
  cat("最大值:", fdr_summary[6], "\n")

  cat("\nFDR < 0.05的基因对比例:", mean(results_sorted$FDR < 0.05, na.rm = TRUE) * 100, "%\n")
  cat("FDR < 0.01的基因对比例:", mean(results_sorted$FDR < 0.01, na.rm = TRUE) * 100, "%\n")
  cat("FDR < 0.001的基因对比例:", mean(results_sorted$FDR < 0.001, na.rm = TRUE) * 100, "%\n")
  cat("FDR < 0.0001的基因对比例:", mean(results_sorted$FDR < 0.0001, na.rm = TRUE) * 100, "%\n")
  cat("FDR < 0.00001的基因对比例:", mean(results_sorted$FDR < 0.00001, na.rm = TRUE) * 100, "%\n")

  sink()

  cat("\n统计摘要已保存至:", stats_file, "\n")

  # 13. 检查输出文件中的FDR列
  cat("\n=== 检查输出文件 ===")
  cat("\n打开完整结果文件查看列名和FDR列...\n")
  if (file.exists(output_file_full)) {
    # 读取保存的文件前几行检查
    check_file <- fread(output_file_full, nrows = 3)
    cat("输出文件的前3行:\n")
    print(check_file)
    cat("\n输出文件的列名:", paste(colnames(check_file), collapse = ", "), "\n")

    if ("FDR" %in% colnames(check_file)) {
      cat("✓ FDR列存在于输出文件中\n")
    } else {
      cat("✗ 警告: FDR列不存在于输出文件中\n")
    }
  }

  cat("\n=== 分析完成 ===\n")
  cat("共分析了", nrow(results), "个基因对\n")
  cat("使用的共同样本数:", length(common_samples), "\n")
  cat("生成了以下文件:\n")
  cat("1. 样本对齐信息.txt\n")
  cat("2. 免疫治疗基因_vs_失巢死亡基因_Spearman_Full_Results.csv (完整结果)\n")
  for(threshold in threshold_names) {
    cat(paste0("3. 免疫治疗基因_vs_失巢死亡基因_Spearman_", threshold, ".csv"), "\n")
  }
  cat("4. 免疫治疗基因相关性分析统计摘要.txt\n")
  cat("输出目录:", output_dir, "\n")
}
