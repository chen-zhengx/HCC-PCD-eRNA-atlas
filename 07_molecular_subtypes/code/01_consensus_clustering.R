# =========================================================
# 对应结果：分子亚型 Molecular_Cluster1/2 的一致性聚类（Figure 3D、S1D、S2D 及 Figure 5B 的上游）
# 输入：input/<pcd>_eRNA_RPM_tumor70.csv（PCD 相关 eRNA × 70 个肿瘤样本的 RPM，与 05_prognostic_model 的 LASSO-Cox 输入相同）
# 输出：output/consensus/<pcd>/cluster_assignments_k2..kmaxK.csv、consensus_results.rds、consensus*.png
# 运行：在模块目录下 Rscript code/01_consensus_clustering.R
# =========================================================

#' 运行一致性聚类分析（输入路径和文件分开）
#'
#' @param input_dir 输入文件所在目录
#' @param input_file 输入文件名（CSV格式）
#' @param output_dir 输出目录路径
#' @param maxK 最大聚类数 (默认=6)
#' @param reps 重采样次数 (默认=100)
#' @param pItem 样本抽样比例 (默认=0.8)
#' @param pFeature 特征抽样比例 (默认=0.8)
#' @param clusterAlg 聚类算法 (默认="hc")
#' @param distance 距离度量方法 (默认="pearson")
#' @param seed 随机种子 (默认=12345)
#'
#' @return 无返回值，结果保存到输出目录
#' @export
run_consensus_clustering <- function(input_dir,
                                     input_file,
                                     output_dir,
                                     maxK = 6,
                                     reps = 100,
                                     pItem = 0.8,
                                     pFeature = 0.8,
                                     clusterAlg = "hc",
                                     distance = "pearson",
                                     seed = 12345) {

  # 检查并安装必要的包
  if (!require("ConsensusClusterPlus", quietly = TRUE)) {
    install.packages("ConsensusClusterPlus")
  }
  if (!require("readr", quietly = TRUE)) {
    install.packages("readr")
  }
  library(ConsensusClusterPlus, quietly = TRUE)
  library(readr, quietly = TRUE)

  # 创建输出目录
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # 构建完整输入路径
  input_path <- file.path(input_dir, input_file)

  # 检查文件是否存在
  if (!file.exists(input_path)) {
    stop("输入文件不存在: ", input_path)
  }

  # 读取表达矩阵
  cat("读取表达矩阵...\n")
  expr_matrix <- read.csv(input_path, row.names = 1)

  # 验证数据格式
  if (nrow(expr_matrix) < 10 || ncol(expr_matrix) < 5) {
    stop("表达矩阵太小：行数=", nrow(expr_matrix), "，列数=", ncol(expr_matrix))
  }

  # 运行一致性聚类
  cat("开始一致性聚类分析...\n")
  cat("样本数: ", ncol(expr_matrix), "\n")
  cat("特征数: ", nrow(expr_matrix), "\n")

  results <- ConsensusClusterPlus(
    d = as.matrix(expr_matrix),
    maxK = maxK,
    reps = reps,
    pItem = pItem,
    pFeature = pFeature,
    clusterAlg = clusterAlg,
    distance = distance,
    title = output_dir,
    plot = "png",
    seed = seed,
    verbose = TRUE
  )

  # 保存所有聚类结果
  cat("保存聚类结果...\n")
  for(k in 2:maxK) {
    cluster_labels <- results[[k]]$consensusClass
    output_file <- file.path(output_dir, paste0("cluster_assignments_k", k, ".csv"))
    write.csv(
      data.frame(Sample = names(cluster_labels), Cluster = cluster_labels),
      output_file,
      row.names = FALSE
    )
    cat("保存: ", basename(output_file), "\n")
  }

  # 保存完整结果对象
  rds_file <- file.path(output_dir, "consensus_results.rds")
  saveRDS(results, rds_file)
  cat("保存完整结果对象: ", basename(rds_file), "\n")

  # 生成结果摘要
  result_summary <- paste(
    "\n===== 一致性聚类分析完成! =====",
    "\n输入文件: ", input_path,
    "\n输出目录: ", output_dir,
    "\n包含结果:",
    paste0("- ", maxK-1, "个聚类标签文件 (k=2-", maxK, ")"),
    paste0("- 完整结果对象: ", basename(rds_file)),
    "- 可视化图表 (PNG格式)",
    paste0("\n分析耗时: ", round(Sys.time() - start_time, 1), " 秒"),
    sep = "\n"
  )

  cat(result_summary)
}
# 定义数据集列表
datasets <- list(
  list(
    input_dir = "input",
    input_file = "anoikis_eRNA_RPM_tumor70.csv",
    output_dir = "output/consensus/anoikis"
  ),
  list(
    input_dir = "input",
    input_file = "apoptosis_eRNA_RPM_tumor70.csv",
    output_dir = "output/consensus/apoptosis"
  ),
  list(
    input_dir = "input",
    input_file = "ferroptosis_eRNA_RPM_tumor70.csv",
    output_dir = "output/consensus/ferroptosis",
    maxK = 4  # 自定义参数
  )
)

# 批量运行
start_time <- Sys.time()
lapply(datasets, function(params) {
  cat("\n===== 开始分析: ", params$input_file, " =====\n")
  tryCatch({
    do.call(run_consensus_clustering, params)
    cat("===== 分析成功完成 =====\n\n")
    return(TRUE)
  }, error = function(e) {
    cat("错误: ", e$message, "\n")
    return(FALSE)
  })
})
cat("总耗时:", round(Sys.time() - start_time, 1), "秒\n")
