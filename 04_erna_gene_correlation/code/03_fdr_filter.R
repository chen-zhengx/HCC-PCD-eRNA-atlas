# =========================================================
# 对应结果：按不同 FDR 阈值（1e-2 至 1e-10）筛选每个数据集的显著 eRNA–基因对并统计数量（Table 1 所用阈值为 FDR<1e-5）
# 输入：output/correlation/<PCD>/<GSE>/output/<GSE>_All_Results.csv.gz（02_spearman_correlation.R 的输出）
# 输出：output/correlation/<PCD>/<GSE>/output/filtered/<GSE>_FDR_*.csv、statistics/<GSE>_FDR_Statistics.csv
#       output/correlation/<PCD>/FDR_filter_log.txt
# 运行：在模块目录下 Rscript code/03_fdr_filter.R
# =========================================================
pcd_list <- c("Anoikis", "Apoptosis", "Ferroptosis")

# 定义要处理的数据集文件夹
dataset_folders <- c(
  "GSE148355",
  "GSE77509",
  "GSE94660",
  "GSE124535",
  "GSE214846",
  "GSE144269"
)

# 定义FDR阈值序列
fdr_thresholds <- c(0.01, 0.001, 0.0001, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
threshold_names <- c("001", "0001", "00001", "000001", "0000001", "00000001", 
                     "000000001", "0000000001", "00000000001")

# 函数：处理单个数据集
process_dataset <- function(dataset_folder) {
  # 创建路径
  dataset_dir <- file.path(base_path, dataset_folder)
  input_dir <- file.path(dataset_dir, "output")
  output_dir <- file.path(dataset_dir, "output")

  # 创建筛选文件和统计文件的专用子目录
  filter_dir <- file.path(output_dir, "filtered")
  stats_dir <- file.path(output_dir, "statistics")
  
  # 确保目录存在
  dirs_to_create <- c(input_dir, output_dir, filter_dir, stats_dir)
  for (dir_path in dirs_to_create) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      cat(paste("目录不存在，已创建:", dir_path, "\n"))
    }
  }
  
  # 获取数据集前缀（如"20GSE77509"）
  dataset_prefix <- gsub("对", "", dataset_folder)
  
  # 构建输入文件路径
  input_file <- file.path(input_dir, paste0(dataset_prefix, "_All_Results.csv.gz"))  # 相关性计算输出为 gzip 压缩
  
  # 检查输入文件是否存在
  if (!file.exists(input_file)) {
    warning(paste("输入文件不存在:", input_file, "\n跳过此数据集"))
    return(NULL)
  }
  
  # 读取数据
  data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
  
  # 检查数据列是否符合预期
  required_columns <- c("eRNA", "Gene", "Correlation", "P.Value", "FDR")
  if (!all(required_columns %in% colnames(data))) {
    missing_cols <- setdiff(required_columns, colnames(data))
    warning(paste("数据文件中缺少必要的列:", paste(missing_cols, collapse = ", "), "\n跳过此数据集"))
    return(NULL)
  }
  
  # 创建统计向量
  count_stats <- integer(length(fdr_thresholds))
  names(count_stats) <- threshold_names
  
  # 循环处理每个阈值
  for (i in seq_along(fdr_thresholds)) {
    threshold <- fdr_thresholds[i]
    label <- threshold_names[i]
    
    # 筛选数据
    filtered_data <- data[data$FDR < threshold, ]
    
    # 存储统计结果
    count_stats[i] <- nrow(filtered_data)
    
    # 保存筛选结果到筛选目录
    output_file <- file.path(filter_dir, paste0(dataset_prefix, "_FDR_", label, ".csv"))
    write.csv(filtered_data, output_file, row.names = FALSE)
  }
  
  # 生成统计文件 - 表格形式
  stats_df <- data.frame(
    FDR_Threshold = fdr_thresholds,
    Scientific_Notation = formatC(fdr_thresholds, format = "e", digits = 0),
    Number_of_Pairs = count_stats
  )
  
  # 保存统计结果到统计文件目录
  stats_csv_file <- file.path(stats_dir, paste0(dataset_prefix, "_FDR_Statistics.csv"))
  write.csv(stats_df, stats_csv_file, row.names = FALSE)
  
  # 同时保存为文本格式
  stats_txt_file <- file.path(stats_dir, paste0(dataset_prefix, "_FDR_Statistics.txt"))
  write.table(stats_df, 
              stats_txt_file, 
              sep = "\t", 
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)
  
  # 返回处理摘要
  list(
    dataset = dataset_folder,
    input_file = input_file,
    filter_dir = filter_dir,
    stats_dir = stats_dir,
    stats_df = stats_df
  )
}

for (pcd in pcd_list) {  # 按 PCD 循环
# 设置基础路径
base_path <- file.path("output", "correlation", pcd)

# 创建主日志文件
main_log_file <- file.path(base_path, "FDR_filter_log.txt")
sink(main_log_file, split = TRUE)  # 同时输出到控制台和文件

# 记录开始时间
start_time <- Sys.time()
cat(paste("开始处理时间:", start_time, "\n\n"))

# 处理所有数据集
results <- list()
for (dataset_folder in dataset_folders) {
  # 打印分隔线
  cat(paste0("\n", strrep("=", 50), "\n"))
  cat(paste("开始处理数据集:", dataset_folder, "\n"))
  
  # 处理单个数据集
  result <- process_dataset(dataset_folder)
  
  if (!is.null(result)) {
    results[[dataset_folder]] <- result
    
    # 输出成功摘要
    cat(paste("\n数据集", dataset_folder, "处理完成!\n"))
    cat(paste("输入文件:", result$input_file, "\n"))
    cat(paste("筛选文件目录:", result$filter_dir, "\n"))
    cat(paste("统计文件目录:", result$stats_dir, "\n"))
    cat("各阈值统计结果：\n")
    print(result$stats_df)
  } else {
    # 输出失败信息
    cat(paste("\n数据集", dataset_folder, "处理失败!\n"))
  }
  
  # 打印分隔线
  cat(paste0("\n", strrep("=", 50), "\n\n"))
}

# 记录结束时间
end_time <- Sys.time()
duration <- end_time - start_time
cat(paste("\n所有数据集处理完成!\n"))
cat(paste("开始时间:", start_time, "\n"))
cat(paste("结束时间:", end_time, "\n"))
cat(paste("总耗时:", round(duration, 2), attr(duration, "units"), "\n"))

# 关闭日志文件
sink()

# 显示最终摘要
cat("\n===== 最终处理摘要 =====\n")
cat(paste("处理数据集数量:", length(dataset_folders), "\n"))
cat(paste("成功处理数据集:", length(results), "\n"))
cat(paste("失败数据集:", length(dataset_folders) - length(results), "\n"))
cat(paste("主日志文件:", main_log_file, "\n"))

}  # 结束 PCD 循环
