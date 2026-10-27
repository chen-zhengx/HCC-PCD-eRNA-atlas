# =========================================================
# 对应结果：Figure 8、Figure S8、Figure S9 与 Supplementary Table 12 的多层共表达关联网络——第二层相关的输入：
#           从各层中提取与该风险 eRNA 显著相关（FDR < 0.01）的节点（候选靶基因、PCD 相关基因、免疫细胞、药物、免疫治疗基因）
# 输入：output/correlations/<PCD>/<eRNA>/eRNA_vs_*/ 中 FDR_0.01 的结果（03 脚本）、output/inputs/（02 脚本）
#       input/immune_cell_ssGSEA_70tumors.csv、input/drug_IC50_19drugs_70tumors.csv、input/risk_eRNAs_by_PCD.csv
# 输出：output/correlations/<PCD>/<eRNA>/extracted/提取的靶基因|<PCD>基因|免疫细胞|药物|免疫治疗基因表达矩阵.csv 及匹配列表
# 运行：在模块目录下 Rscript code/04_extract_correlated_nodes.R
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
  # 1. 与该 eRNA 显著相关的 ±1 Mb 候选靶基因
  # ---------------------------------------------------------
  corr_file <- file.path(run_dir, "eRNA_vs_target_genes", paste0(erna_id, "_vs_1MB_靶基因_Spearman_FDR_0.01.csv"))
  if (file.exists(corr_file)) {
    # 简化版本 - 输入输出同路径
    work_dir <- file.path(run_dir, "extracted")
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

    # 检查目录是否存在
    if (!dir.exists(work_dir)) {
      stop("错误：工作目录不存在 - ", work_dir)
    }

    cat("工作目录:", work_dir, "\n")
    cat("正在读取文件...\n")

    # 读取文件
    exp_matrix <- read.csv(file.path("output", "inputs", "target_genes", paste0(erna_id, "_1MB靶基因矩阵.csv")),
                           stringsAsFactors = FALSE, check.names = FALSE)

    # 检查相关性表列名，可能是中文列名
    corr_table <- read.csv(corr_file,  # 相关性结果直接读 03 脚本的输出
                           stringsAsFactors = FALSE, check.names = FALSE)

    # 提取唯一的Gene列表 - 检查列名可能是中文
    if ("Gene" %in% colnames(corr_table)) {
      gene_list <- unique(corr_table$Gene)
    } else if ("靶基因" %in% colnames(corr_table)) {
      gene_list <- unique(corr_table$靶基因)
      # 重命名列为Gene以便后续使用
      corr_table$Gene <- corr_table$靶基因
    } else {
      # 如果两种列名都没有，尝试第一列
      gene_list <- unique(corr_table[[1]])
      warning("相关性表列名不是'Gene'或'靶基因'，使用第一列作为基因列")
    }

    # 检查列名
    if (!"Symbol" %in% colnames(exp_matrix)) {
      colnames(exp_matrix)[1] <- "Symbol"
    }

    # 提取匹配的基因
    matched_genes <- gene_list[gene_list %in% exp_matrix$Symbol]
    not_matched <- gene_list[!gene_list %in% exp_matrix$Symbol]

    if (length(matched_genes) == 0) {
      warning("警告：没有匹配到任何基因！")
      extracted_matrix <- data.frame()
    } else {
      extracted_matrix <- exp_matrix[exp_matrix$Symbol %in% matched_genes, ]
    }

    # 保存结果
    if (length(matched_genes) > 0) {
      write.csv(extracted_matrix, 
                file.path(work_dir, "提取的靶基因表达矩阵.csv"),
                row.names = FALSE, quote = FALSE)
    }

    cat("\n处理完成!\n")
    cat("工作目录:", work_dir, "\n")
    cat("原始表达矩阵基因数:", nrow(exp_matrix), "\n")
    cat("相关性表中Gene数:", length(gene_list), "\n")
    cat("匹配到的Gene数:", length(matched_genes), "\n")
    if (length(matched_genes) > 0) {
      cat("提取的矩阵维度:", dim(extracted_matrix), "\n")
    }

    if (length(not_matched) > 0) {
      cat("\n未匹配的基因（", length(not_matched), "个）:\n")
      cat(paste(not_matched, collapse = ", "), "\n")
    }
  } else {
    cat("本层没有 FDR < 0.01 的节点，跳过提取:", corr_file, "\n")
  }

  # ---------------------------------------------------------
  # 2. 与该 eRNA 显著相关的 PCD 相关基因
  # ---------------------------------------------------------
  corr_file <- file.path(run_dir, "eRNA_vs_PCD_genes", paste0(erna_id, "_vs_", pcd_cn, "基因_Spearman_FDR_0.01.csv"))
  if (file.exists(corr_file)) {
    # 简化版本 - 输入输出同路径
    work_dir <- file.path(run_dir, "extracted")
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

    # 检查目录是否存在
    if (!dir.exists(work_dir)) {
      stop("错误：工作目录不存在 - ", work_dir)
    }

    cat("工作目录:", work_dir, "\n")
    cat("正在读取文件...\n")

    # 读取文件
    exp_matrix <- read.csv(pcd_matrix_file,
                           stringsAsFactors = FALSE, check.names = FALSE)

    # 检查相关性表列名，可能是中文列名
    corr_table <- read.csv(corr_file,  # 相关性结果直接读 03 脚本的输出
                           stringsAsFactors = FALSE, check.names = FALSE)

    # 提取唯一的Gene列表 - 检查列名可能是中文
    if ("Gene" %in% colnames(corr_table)) {
      gene_list <- unique(corr_table$Gene)
    } else if ("靶基因" %in% colnames(corr_table)) {
      gene_list <- unique(corr_table$靶基因)
      # 重命名列为Gene以便后续使用
      corr_table$Gene <- corr_table$靶基因
    } else {
      # 如果两种列名都没有，尝试第一列
      gene_list <- unique(corr_table[[1]])
      warning("相关性表列名不是'Gene'或'靶基因'，使用第一列作为基因列")
    }

    # 检查列名
    if (!"Symbol" %in% colnames(exp_matrix)) {
      colnames(exp_matrix)[1] <- "Symbol"
    }

    # 提取匹配的基因
    matched_genes <- gene_list[gene_list %in% exp_matrix$Symbol]
    not_matched <- gene_list[!gene_list %in% exp_matrix$Symbol]

    if (length(matched_genes) == 0) {
      warning("警告：没有匹配到任何基因！")
      extracted_matrix <- data.frame()
    } else {
      extracted_matrix <- exp_matrix[exp_matrix$Symbol %in% matched_genes, ]
    }

    # 保存结果
    if (length(matched_genes) > 0) {
      write.csv(extracted_matrix, 
                file.path(work_dir, paste0("提取的", pcd_cn, "基因表达矩阵.csv")),
                row.names = FALSE, quote = FALSE)
    }

    cat("\n处理完成!\n")
    cat("工作目录:", work_dir, "\n")
    cat("原始表达矩阵基因数:", nrow(exp_matrix), "\n")
    cat("相关性表中Gene数:", length(gene_list), "\n")
    cat("匹配到的Gene数:", length(matched_genes), "\n")
    if (length(matched_genes) > 0) {
      cat("提取的矩阵维度:", dim(extracted_matrix), "\n")
    }

    if (length(not_matched) > 0) {
      cat("\n未匹配的基因（", length(not_matched), "个）:\n")
      cat(paste(not_matched, collapse = ", "), "\n")
    }
  } else {
    cat("本层没有 FDR < 0.01 的节点，跳过提取:", corr_file, "\n")
  }

  # ---------------------------------------------------------
  # 3. 与该 eRNA 显著相关的免疫细胞
  # ---------------------------------------------------------
  corr_file <- file.path(run_dir, "eRNA_vs_immune_cells", paste0(erna_id, "表达_vs_免疫细胞_Spearman_FDR_0.01.csv"))
  if (file.exists(corr_file)) {
    # 简化版本 - 输入输出同路径
    work_dir <- file.path(run_dir, "extracted")
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

    # 检查目录是否存在
    if (!dir.exists(work_dir)) {
      stop("错误：工作目录不存在 - ", work_dir)
    }

    cat("工作目录:", work_dir, "\n")
    cat("正在读取文件...\n")

    # 读取文件
    exp_matrix <- read.csv(file.path("input", "immune_cell_ssGSEA_70tumors.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE)

    # 检查相关性表列名，可能是中文列名
    corr_table <- read.csv(corr_file,  # 相关性结果直接读 03 脚本的输出
                           stringsAsFactors = FALSE, check.names = FALSE)

    # 提取唯一的免疫细胞列表 - 检查列名可能是中文
    if ("immune_cell" %in% colnames(corr_table)) {
      immune_cell_list <- unique(corr_table$immune_cell)
    } else if ("免疫细胞" %in% colnames(corr_table)) {
      immune_cell_list <- unique(corr_table$免疫细胞)
      # 重命名为immune_cell以便后续使用
      corr_table$immune_cell <- corr_table$免疫细胞
    } else if ("Gene" %in% colnames(corr_table)) {
      # 兼容旧版本
      immune_cell_list <- unique(corr_table$Gene)
      corr_table$immune_cell <- corr_table$Gene
    } else {
      # 如果以上列名都没有，尝试第一列
      immune_cell_list <- unique(corr_table[[1]])
      warning("相关性表列名不是'immune_cell'、'免疫细胞'或'Gene'，使用第一列作为免疫细胞列")
    }

    # 检查表达矩阵列名
    # 表达矩阵的第一列应该是免疫细胞名称
    if (!"免疫细胞" %in% colnames(exp_matrix) && !"immune_cell" %in% colnames(exp_matrix)) {
      # 将第一列重命名为"免疫细胞"以便后续匹配
      colnames(exp_matrix)[1] <- "免疫细胞"
    }

    # 提取匹配的免疫细胞
    if ("免疫细胞" %in% colnames(exp_matrix)) {
      matched_cells <- immune_cell_list[immune_cell_list %in% exp_matrix$免疫细胞]
      not_matched <- immune_cell_list[!immune_cell_list %in% exp_matrix$免疫细胞]

      if (length(matched_cells) == 0) {
        warning("警告：没有匹配到任何免疫细胞！")
        extracted_matrix <- data.frame()
      } else {
        extracted_matrix <- exp_matrix[exp_matrix$免疫细胞 %in% matched_cells, ]
        # 保持原始列名，只重命名第一列为"immune_cell"以便识别
        colnames(extracted_matrix)[1] <- "immune_cell"
      }
    } else if ("immune_cell" %in% colnames(exp_matrix)) {
      matched_cells <- immune_cell_list[immune_cell_list %in% exp_matrix$immune_cell]
      not_matched <- immune_cell_list[!immune_cell_list %in% exp_matrix$immune_cell]

      if (length(matched_cells) == 0) {
        warning("警告：没有匹配到任何免疫细胞！")
        extracted_matrix <- data.frame()
      } else {
        extracted_matrix <- exp_matrix[exp_matrix$immune_cell %in% matched_cells, ]
      }
    } else {
      # 如果都没有找到合适的列名，使用第一列
      colnames(exp_matrix)[1] <- "immune_cell"
      matched_cells <- immune_cell_list[immune_cell_list %in% exp_matrix$immune_cell]
      not_matched <- immune_cell_list[!immune_cell_list %in% exp_matrix$immune_cell]

      if (length(matched_cells) == 0) {
        warning("警告：没有匹配到任何免疫细胞！")
        extracted_matrix <- data.frame()
      } else {
        extracted_matrix <- exp_matrix[exp_matrix$immune_cell %in% matched_cells, ]
      }
    }

    # 保存结果
    if (length(matched_cells) > 0) {
      output_file <- file.path(work_dir, "提取的免疫细胞表达矩阵.csv")
      write.csv(extracted_matrix, output_file, row.names = FALSE, quote = FALSE)
      cat("提取的免疫细胞表达矩阵已保存至:", output_file, "\n")
    }

    cat("\n处理完成!\n")
    cat("工作目录:", work_dir, "\n")
    cat("原始表达矩阵免疫细胞数:", nrow(exp_matrix), "\n")
    cat("相关性表中免疫细胞数:", length(immune_cell_list), "\n")
    cat("匹配到的免疫细胞数:", length(matched_cells), "\n")
    if (length(matched_cells) > 0) {
      cat("提取的矩阵维度:", dim(extracted_matrix), "\n")
    }

    if (length(not_matched) > 0) {
      cat("\n未匹配的免疫细胞（", length(not_matched), "个）:\n")
      cat(paste(not_matched, collapse = ", "), "\n")

      # 保存未匹配列表到文件
      not_matched_file <- file.path(work_dir, "未匹配的免疫细胞列表.txt")
      writeLines(paste("未匹配的免疫细胞列表（共", length(not_matched), "个）：\n",
                       paste(not_matched, collapse = "\n")), 
                 not_matched_file)
      cat("未匹配列表已保存至:", not_matched_file, "\n")
    }

    # 显示匹配的免疫细胞
    if (length(matched_cells) > 0) {
      cat("\n匹配到的免疫细胞（", length(matched_cells), "个）:\n")
      if (length(matched_cells) <= 20) {
        cat(paste(matched_cells, collapse = ", "), "\n")
      } else {
        cat("前20个:", paste(matched_cells[1:20], collapse = ", "), "\n")
        cat("... 以及", length(matched_cells) - 20, "个更多\n")
      }

      # 保存匹配列表到文件
      matched_file <- file.path(work_dir, "匹配的免疫细胞列表.txt")
      writeLines(paste("匹配的免疫细胞列表（共", length(matched_cells), "个）：\n",
                       paste(matched_cells, collapse = "\n")), 
                 matched_file)
      cat("匹配列表已保存至:", matched_file, "\n")
    }
  } else {
    cat("本层没有 FDR < 0.01 的节点，跳过提取:", corr_file, "\n")
  }

  # ---------------------------------------------------------
  # 4. 与该 eRNA 显著相关的药物
  # ---------------------------------------------------------
  corr_file <- file.path(run_dir, "eRNA_vs_drugs", paste0(erna_id, "表达_vs_敏感药物_Spearman_FDR_0.01.csv"))
  if (file.exists(corr_file)) {
    # 简化版本 - 输入输出同路径
    work_dir <- file.path(run_dir, "extracted")
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

    # 检查目录是否存在
    if (!dir.exists(work_dir)) {
      stop("错误：工作目录不存在 - ", work_dir)
    }

    cat("工作目录:", work_dir, "\n")
    cat("正在读取文件...\n")

    # 读取文件
    exp_matrix <- read.csv(file.path("input", "drug_IC50_19drugs_70tumors.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE)

    corr_table <- read.csv(corr_file,  # 相关性结果直接读 03 脚本的输出
                           stringsAsFactors = FALSE, check.names = FALSE)

    # 提取药物列表 - 从相关性表第二列
    if (ncol(corr_table) >= 2) {
      # 直接使用第二列作为药物列
      drug_list <- unique(corr_table[[2]])

      # 给第二列命名为"drug"以便后续使用
      colnames(corr_table)[2] <- "drug"

      cat("从相关性表第二列提取药物列表，共", length(drug_list), "种药物\n")
    } else {
      stop("相关性表列数不足，至少需要2列")
    }

    # 检查表达矩阵药物列
    if ("药物" %in% colnames(exp_matrix)) {
      # 重命名为"drug"以统一
      colnames(exp_matrix)[which(colnames(exp_matrix) == "药物")] <- "drug"
      cat("识别到表达矩阵中名为'药物'的列，已重命名为'drug'\n")
    } else if (!"drug" %in% colnames(exp_matrix)) {
      # 如果既没有"药物"也没有"drug"，尝试第一列
      colnames(exp_matrix)[1] <- "drug"
      warning("表达矩阵列名不是'药物'或'drug'，使用第一列作为药物列")
    }

    # 提取匹配的药物
    matched_drugs <- drug_list[drug_list %in% exp_matrix$drug]
    not_matched <- drug_list[!drug_list %in% exp_matrix$drug]

    if (length(matched_drugs) == 0) {
      warning("警告：没有匹配到任何药物！")
      extracted_matrix <- data.frame()
    } else {
      extracted_matrix <- exp_matrix[exp_matrix$drug %in% matched_drugs, ]

      # 保存结果
      output_file <- file.path(work_dir, "提取的药物表达矩阵.csv")
      write.csv(extracted_matrix, output_file, row.names = FALSE, quote = FALSE)
      cat("提取的药物表达矩阵已保存至:", output_file, "\n")
    }

    # 输出结果统计
    cat("\n处理完成!\n")
    cat("工作目录:", work_dir, "\n")
    cat("原始表达矩阵药物数:", nrow(exp_matrix), "\n")
    cat("相关性表中药物数:", length(drug_list), "\n")
    cat("匹配到的药物数:", length(matched_drugs), "\n")
    if (length(matched_drugs) > 0) {
      cat("提取的矩阵维度:", dim(extracted_matrix), "\n")
    }

    # 保存匹配和未匹配列表
    if (length(not_matched) > 0) {
      not_matched_file <- file.path(work_dir, "未匹配的药物列表.txt")
      writeLines(paste("未匹配的药物列表（共", length(not_matched), "种）：\n",
                       paste(not_matched, collapse = "\n")), 
                 not_matched_file)
      cat("未匹配列表已保存至:", not_matched_file, "\n")
    }

    if (length(matched_drugs) > 0) {
      matched_file <- file.path(work_dir, "匹配的药物列表.txt")
      writeLines(paste("匹配的药物列表（共", length(matched_drugs), "种）：\n",
                       paste(matched_drugs, collapse = "\n")), 
                 matched_file)
      cat("匹配列表已保存至:", matched_file, "\n")

      # 显示前20个匹配的药物
      cat("\n匹配到的药物（前20种）:\n")
      if (length(matched_drugs) <= 20) {
        cat(paste(matched_drugs, collapse = ", "), "\n")
      } else {
        cat(paste(matched_drugs[1:20], collapse = ", "), "\n")
        cat("... 以及", length(matched_drugs) - 20, "种更多\n")
      }
    }
  } else {
    cat("本层没有 FDR < 0.01 的节点，跳过提取:", corr_file, "\n")
  }

  # ---------------------------------------------------------
  # 5. 与该 eRNA 显著相关的免疫治疗相关基因
  # ---------------------------------------------------------
  corr_file <- file.path(run_dir, "eRNA_vs_immunotherapy_genes", paste0(erna_id, "_vs_免疫基因_Spearman_FDR_0.01.csv"))
  if (file.exists(corr_file)) {
    # 简化版本 - 输入输出同路径
    work_dir <- file.path(run_dir, "extracted")
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

    # 检查目录是否存在
    if (!dir.exists(work_dir)) {
      stop("错误：工作目录不存在 - ", work_dir)
    }

    cat("工作目录:", work_dir, "\n")
    cat("正在读取文件...\n")

    # 读取文件
    exp_matrix <- read.csv(file.path("output", "inputs", "免疫治疗基因表达矩阵.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE)

    corr_table <- read.csv(corr_file,  # 相关性结果直接读 03 脚本的输出
                           stringsAsFactors = FALSE, check.names = FALSE)

    # 显示相关性表的列名
    cat("相关性表的列名:", colnames(corr_table), "\n")

    # 检查是否有Gene列（第二列）
    if ("Gene" %in% colnames(corr_table)) {
      # 使用Gene列提取基因列表
      gene_list <- unique(corr_table[["Gene"]])
      cat("从相关性表Gene列提取基因列表，共", length(gene_list), "个基因\n")
    } else if (ncol(corr_table) >= 2) {
      # 如果第二列没有列名，使用第二列
      gene_list <- unique(corr_table[[2]])
      colnames(corr_table)[2] <- "Gene"  # 给第二列命名为Gene
      cat("从相关性表第二列提取基因列表，共", length(gene_list), "个基因\n")
    } else {
      stop("相关性表没有Gene列且列数不足2列")
    }

    # 检查表达矩阵是否有Gene列
    if ("Gene" %in% colnames(exp_matrix)) {
      cat("表达矩阵有名为'Gene'的列\n")
    } else if ("gene" %in% colnames(exp_matrix)) {
      # 如果有小写的gene列，重命名为Gene以统一
      colnames(exp_matrix)[which(colnames(exp_matrix) == "gene")] <- "Gene"
      cat("表达矩阵有名为'gene'的列，已重命名为'Gene'\n")
    } else if ("基因" %in% colnames(exp_matrix)) {
      # 如果有中文列名，重命名为Gene
      colnames(exp_matrix)[which(colnames(exp_matrix) == "基因")] <- "Gene"
      cat("表达矩阵有名为'基因'的列，已重命名为'Gene'\n")
    } else {
      stop("表达矩阵没有找到Gene/gene/基因列")
    }

    # 提取匹配的基因
    matched_genes <- gene_list[gene_list %in% exp_matrix$Gene]
    not_matched <- gene_list[!gene_list %in% exp_matrix$Gene]

    if (length(matched_genes) == 0) {
      warning("警告：没有匹配到任何基因！")
      extracted_matrix <- data.frame()
    } else {
      extracted_matrix <- exp_matrix[exp_matrix$Gene %in% matched_genes, ]

      # 保存结果
      output_file <- file.path(work_dir, "提取的免疫治疗基因表达矩阵.csv")
      write.csv(extracted_matrix, output_file, row.names = FALSE, quote = FALSE)
      cat("提取的免疫治疗基因表达矩阵已保存至:", output_file, "\n")
    }

    # 输出结果统计
    cat("\n处理完成!\n")
    cat("工作目录:", work_dir, "\n")
    cat("原始表达矩阵基因数:", nrow(exp_matrix), "\n")
    cat("相关性表中基因数:", length(gene_list), "\n")
    cat("匹配到的基因数:", length(matched_genes), "\n")
    if (length(matched_genes) > 0) {
      cat("提取的矩阵维度:", dim(extracted_matrix), "\n")
    }

    # 保存匹配和未匹配列表
    if (length(not_matched) > 0) {
      not_matched_file <- file.path(work_dir, "未匹配的基因列表.txt")
      writeLines(paste("未匹配的基因列表（共", length(not_matched), "个）：\n",
                       paste(not_matched, collapse = "\n")), 
                 not_matched_file)
      cat("未匹配列表已保存至:", not_matched_file, "\n")
    }

    if (length(matched_genes) > 0) {
      matched_file <- file.path(work_dir, "匹配的基因列表.txt")
      writeLines(paste("匹配的基因列表（共", length(matched_genes), "个）：\n",
                       paste(matched_genes, collapse = "\n")), 
                 matched_file)
      cat("匹配列表已保存至:", matched_file, "\n")

      # 显示前20个匹配的基因
      cat("\n匹配到的基因（前20个）:\n")
      if (length(matched_genes) <= 20) {
        cat(paste(matched_genes, collapse = ", "), "\n")
      } else {
        cat(paste(matched_genes[1:20], collapse = ", "), "\n")
        cat("... 以及", length(matched_genes) - 20, "个更多\n")
      }
    }
  } else {
    cat("本层没有 FDR < 0.01 的节点，跳过提取:", corr_file, "\n")
  }

}
