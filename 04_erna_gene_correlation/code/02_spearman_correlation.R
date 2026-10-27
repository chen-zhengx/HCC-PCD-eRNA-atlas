# =========================================================
# 对应结果：每个数据集内差异 eRNA 与差异 PCD 相关基因的 Spearman 相关（合并肿瘤与配对非肿瘤样本计算，数据集内 BH 校正）
#           Table 1、Table S6 的上游
# 输入：output/correlation/<PCD>/<GSE>/input/（01_extract_inputs.R 的输出）
# 输出：output/correlation/<PCD>/<GSE>/output/<GSE>_All_Results.csv.gz、<GSE>_FDR_*.csv、Analysis_Report.csv
#       output/correlation/<PCD>/Summary_Results_All_Datasets.csv
# 运行：在模块目录下 Rscript code/02_spearman_correlation.R
# =========================================================
# 加载必要的包
library(dplyr)
library(data.table)
library(foreach)
library(doParallel)
library(matrixStats)  # 新增用于高效计算标准差

pcd_list <- c("Anoikis", "Apoptosis", "Ferroptosis")

# 定义要处理的数据集
datasets <- c(
  "GSE148355",
  "GSE77509",
  "GSE94660",
  "GSE124535",
  "GSE214846",
  "GSE144269"
)

for (pcd in pcd_list) {
# 设置主目录
base_dir <- file.path("output", "correlation", pcd)

# 创建总结果汇总表
summary_results <- data.frame()

# 记录开始时间
overall_start <- Sys.time()

# 遍历每个数据集
for(dataset in datasets) {
  start_time <- Sys.time()
  
  # 1. 设置路径 ------------------------------------------------------------
  input_dir <- file.path(base_dir, dataset, "input")
  output_dir <- file.path(base_dir, dataset, "output")
  
  # 创建输出目录（如果不存在）
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, mode = "0755")
    cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 创建输出目录:", output_dir, "\n")
  }
  
  # 2. 读取数据 ------------------------------------------------------------
  cat("\n[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 开始处理数据集:", dataset, "\n")
  cat("输入目录:", input_dir, "\n")
  
  # 获取输入文件列表（忽略大小写）
  input_files <- list.files(input_dir, pattern = "\\.csv$", full.names = FALSE, ignore.case = TRUE)
  cat("找到的输入文件:", paste(input_files, collapse = ", "), "\n")
  
  # 识别eRNA和基因文件（使用更健壮的正则表达式，忽略大小写）
  erna_file <- grep("_DE_eRNA0\\.5_matrix\\.csv$", input_files, value = TRUE, ignore.case = TRUE)
  gene_file <- grep(paste0("_DE_", pcd, "_genes\\.csv$"), input_files, value = TRUE, ignore.case = TRUE)
  
  # 检查文件是否找到
  if(length(erna_file) == 0) {
    warning("未找到eRNA表达矩阵文件，跳过该数据集")
    next
  }
  if(length(gene_file) == 0) {
    warning("未找到基因表达矩阵文件，跳过该数据集")
    next
  }
  
  # 使用第一个匹配的文件（如果有多个）
  erna_file <- file.path(input_dir, erna_file[1])
  gene_file <- file.path(input_dir, gene_file[1])
  
  cat("使用eRNA文件:", erna_file, "\n")
  cat("使用基因文件:", gene_file, "\n")
  
  # 3. 读取数据文件 ------------------------------------------------------
  cat("正在读取数据文件...\n")
  
  # 使用fread快速读取大数据文件
  tryCatch({
    # 读取eRNA表达矩阵
    erna_matrix <- fread(erna_file, header = TRUE, check.names = FALSE) %>% 
      as.data.frame()
    rownames(erna_matrix) <- erna_matrix[,1]
    erna_values <- erna_matrix[,-1, drop = FALSE] %>% 
      mutate_all(as.numeric) %>% 
      as.matrix()
    
    # 读取基因表达矩阵
    gene_matrix <- fread(gene_file, header = TRUE, check.names = FALSE) %>% 
      as.data.frame()
    rownames(gene_matrix) <- gene_matrix[,1]
    gene_values <- gene_matrix[,-1, drop = FALSE] %>% 
      mutate_all(as.numeric) %>% 
      as.matrix()
  }, error = function(e) {
    cat("读取文件出错:", e$message, "\n")
    next
  })
  
  # 4. 数据预处理 --------------------------------------------------------
  cat("正在进行数据预处理...\n")
  
  # 去除标准差为0的行（使用更高效的方法）
  erna_sd <- matrixStats::rowSds(erna_values, na.rm = TRUE)
  gene_sd <- matrixStats::rowSds(gene_values, na.rm = TRUE)
  
  erna_values <- erna_values[erna_sd > 0 & !is.na(erna_sd), , drop = FALSE]
  gene_values <- gene_values[gene_sd > 0 & !is.na(gene_sd), , drop = FALSE]
  
  # 检查样本顺序一致性
  if(!identical(colnames(erna_values), colnames(gene_values))) {
    # 尝试按列名排序
    common_samples <- intersect(colnames(erna_values), colnames(gene_values))
    if(length(common_samples) > 0) {
      erna_values <- erna_values[, common_samples, drop = FALSE]
      gene_values <- gene_values[, common_samples, drop = FALSE]
      cat("样本顺序不一致，已调整为共同样本:", length(common_samples), "个\n")
    } else {
      cat("错误：eRNA和基因表达矩阵没有共同的样本列！跳过该数据集\n")
      next
    }
  }
  
  cat("处理后数据维度:\n")
  cat("eRNA矩阵:", nrow(erna_values), "行(特征)", ncol(erna_values), "列(样本)\n")
  cat("基因矩阵:", nrow(gene_values), "行(特征)", ncol(gene_values), "列(样本)\n")
  
  # 5. 并行计算相关性 ----------------------------------------------------
  cat("开始计算相关性（Spearman）...\n")
  
  # 初始化并行计算（使用服务器可用核心的75%）
  total_cores <- detectCores()
  n_cores <- max(1, floor(total_cores * 0.75))
  registerDoParallel(n_cores)
  cat("使用", n_cores, "个CPU核心进行并行计算\n")
  
  # 生成所有组合的索引
  combinations <- expand.grid(
    erna_idx = 1:nrow(erna_values),
    gene_idx = 1:nrow(gene_values)
  )
  
  # 计算预计内存需求（MB）
  estimated_mem <- (object.size(erna_values) + object.size(gene_values)) / 1024^2
  cat("预计内存需求:", round(estimated_mem, 2), "MB\n")
  
  # 自动调整分块大小（基于内存考虑）
  chunk_size <- min(1000, max(100, floor(1e6 / (nrow(erna_values) * nrow(gene_values)))))
  n_chunks <- ceiling(nrow(combinations) / chunk_size)
  
  cat("使用分块处理，每块大小:", chunk_size, "，总块数:", n_chunks, "\n")
  
  # 并行计算Spearman相关性
  results <- foreach(chunk = 1:n_chunks, .combine = rbind, .packages = c("dplyr")) %dopar% {
    start <- (chunk - 1) * chunk_size + 1
    end <- min(chunk * chunk_size, nrow(combinations))
    chunk_data <- combinations[start:end, ]
    
    chunk_results <- apply(chunk_data, 1, function(row) {
      i <- row["erna_idx"]
      j <- row["gene_idx"]
      cor_test <- tryCatch({
        cor.test(
          erna_values[i, ], 
          gene_values[j, ],
          method = "spearman",
          use = "pairwise.complete.obs",
          exact = FALSE  # 对于大样本使用近似计算
        )
      }, error = function(e) {
        return(list(estimate = NA, p.value = NA))
      })
      
      data.frame(
        eRNA = rownames(erna_values)[i],
        Gene = rownames(gene_values)[j],
        Correlation = cor_test$estimate,
        P.Value = cor_test$p.value,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, chunk_results)
  }
  
  # 停止并行
  stopImplicitCluster()
  
  # 移除NA值（计算失败的对）
  results <- na.omit(results)
  
  # 6. FDR校正和结果筛选 ------------------------------------------------
  cat("进行FDR校正和结果筛选...\n")
  
  # 计算FDR校正后的p值（使用更快的p.adjust方法）
  results$FDR <- p.adjust(results$P.Value, method = "BH")
  
  # 定义筛选条件（基于FDR）
  filtered_results <- list(
    FDR_0.01 = results %>% filter(FDR < 0.01) %>% arrange(desc(abs(Correlation))),
    FDR_0.001 = results %>% filter(FDR < 0.001) %>% arrange(desc(abs(Correlation))),
    FDR_0.0001 = results %>% filter(FDR < 0.0001) %>% arrange(desc(abs(Correlation))),
    FDR_0.00001 = results %>% filter(FDR < 0.00001) %>% arrange(desc(abs(Correlation)))
  )
  
  # 7. 保存结果 --------------------------------------------------------
  cat("保存结果文件...\n")
  
  # 输出文件名前缀（移除中文"对"字）
  output_prefix <- gsub("对", "", dataset)
  
  # 保存不同FDR阈值的文件
  for(threshold in names(filtered_results)) {
    output_file <- file.path(output_dir, paste0(output_prefix, "_", threshold, ".csv"))
    fwrite(filtered_results[[threshold]], output_file)
    cat("已保存:", output_file, "\n")
  }
  
  # 保存完整结果（使用压缩格式节省空间）
  fwrite(results, file.path(output_dir, paste0(output_prefix, "_All_Results.csv.gz")), compress = "gzip")
  
  # 8. 生成分析报告 ----------------------------------------------------
  report_data <- data.frame(
    Dataset = dataset,
    Timestamp = Sys.time(),
    eRNA_Features = nrow(erna_values),
    Gene_Features = nrow(gene_values),
    Samples = ncol(erna_values),
    Total_Correlations = nrow(results),
    FDR_0.01 = nrow(filtered_results$FDR_0.01),
    FDR_0.001 = nrow(filtered_results$FDR_0.001),
    FDR_0.0001 = nrow(filtered_results$FDR_0.0001),
    FDR_0.00001 = nrow(filtered_results$FDR_0.00001),
    stringsAsFactors = FALSE
  )
  
  # 添加到总结果汇总
  summary_results <- bind_rows(summary_results, report_data)
  
  # 保存单个数据集的报告
  fwrite(report_data, file.path(output_dir, "Analysis_Report.csv"))
  
  # 9. 打印摘要 --------------------------------------------------------
  cat("\n[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", dataset, "分析完成\n")
  cat("总相关性计算对数:", format(nrow(results), big.mark = ","), "\n")
  cat("显著相关性对(FDR < 0.01):", format(nrow(filtered_results$FDR_0.01), big.mark = ","), "\n")
  cat("显著相关性对(FDR < 0.001):", format(nrow(filtered_results$FDR_0.001), big.mark = ","), "\n")
  cat("显著相关性对(FDR < 0.0001):", format(nrow(filtered_results$FDR_0.0001), big.mark = ","), "\n")
  cat("显著相关性对(FDR < 0.00001):", format(nrow(filtered_results$FDR_0.00001), big.mark = ","), "\n")
  cat("结果已保存至:", output_dir, "\n")
  cat("处理时间:", format(Sys.time() - start_time, digits = 2), "\n")
}

# 保存汇总结果
if(nrow(summary_results) > 0) {
  summary_file <- file.path(base_dir, "Summary_Results_All_Datasets.csv")
  fwrite(summary_results, summary_file)
  cat("\n所有数据集处理完成！汇总结果已保存至:", summary_file, "\n")
} else {
  cat("\n警告：没有成功处理任何数据集！\n")
}

cat("总运行时间:", format(Sys.time() - overall_start, digits = 2), "\n")
}  # 结束 PCD 循环
