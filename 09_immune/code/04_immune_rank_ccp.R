# =========================================================
# 对应结果：Figure 6C 及正文使用的免疫排名分型 ImmuneRank_Cluster1/2（n = 46/24；Pij + ConsensusClusterPlus
#           hc + average、Pearson 距离、reps = 1000、pItem = 0.8、seed = 100，取 k = 2）
# 输入：input/ssgsea_scores_28cells.csv
# 输出：output/immune_ccp/A_rank_based_CCP/cluster_assignments_k2..k10.csv、consensus_results.rds
#       stability_PAC_itemConsensus.csv、bestK_by_PAC.csv、CCP 自带图
#       output/immune_ccp/B_raw_zscore_based_CCP/*（行 Z-score 的对照分型，论文未使用）
#       其中 A_rank_based_CCP/cluster_assignments_k2.csv 即 input/immune_rank_clusters_CCP_k2.csv
#       ，B_raw_zscore_based_CCP/cluster_assignments_k2.csv 即
#       input/immune_matrix_clusters_CCP_k2.csv
# 运行：在模块目录下 Rscript code/04_immune_rank_ccp.R
# =========================================================

# ============================================================
# 免疫一致性聚类（按你run_consensus_clustering的逻辑改写）
# - 输入：12.16新作的免疫细胞表达.csv (cell x sample)
# - 两套表示：
#   A) Rank/Percentile (Pij) 结构型
#   B) Raw ssGSEA row-z      强度敏感型
# - 聚类：ConsensusClusterPlus (hc + average)
# - 距离：pearson
# - 输出：每套输出目录下保存
#   - cluster_assignments_k2..k10.csv
#   - consensus_results.rds
#   - stability_PAC_itemConsensus.csv + bestK.csv + 简单图（PAC曲线/簇大小）
# ============================================================


#' 运行免疫一致性聚类（同你原函数风格：input_dir + input_file + output_dir）
#'
#' @param input_dir 输入文件所在目录
#' @param input_file 输入文件名（CSV格式，cell x sample）
#' @param output_dir 输出目录路径
#' @param mode "rank" 或 "rawz"；rank=Pij相对结构，rawz=row-z强度敏感
#' @param maxK 最大聚类数 (默认=10)
#' @param reps 重采样次数 (默认=1000)
#' @param pItem 样本抽样比例 (默认=0.8)
#' @param pFeature 特征抽样比例 (默认=1)
#' @param clusterAlg 聚类算法 (默认="hc")
#' @param distance 距离度量方法 (默认="pearson")
#' @param seed 随机种子 (默认=100)
#' @export
run_immune_consensus_clustering <- function(input_dir,
                                            input_file,
                                            output_dir,
                                            mode = c("rank", "rawz"),
                                            maxK = 10,
                                            reps = 1000,
                                            pItem = 0.8,
                                            pFeature = 1,
                                            clusterAlg = "hc",
                                            distance = "pearson",
                                            seed = 100) {
  
  mode <- match.arg(mode)
  start_time <- Sys.time()
  
  # 1) 包
  if (!require("ConsensusClusterPlus", quietly = TRUE)) install.packages("ConsensusClusterPlus")
  library(ConsensusClusterPlus, quietly = TRUE)
  
  # 2) 创建输出目录
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 3) 构建输入路径并读取
  input_path <- file.path(input_dir, input_file)
  if (!file.exists(input_path)) stop("输入文件不存在: ", input_path)
  
  cat("读取表达矩阵...\n")
  expr <- read.csv(input_path, row.names = 1, check.names = FALSE)
  expr <- as.matrix(expr)
  mode(expr) <- "numeric"
  
  if (nrow(expr) < 5 || ncol(expr) < 5) {
    stop("表达矩阵太小：行数=", nrow(expr), "，列数=", ncol(expr))
  }
  
  cat("样本数: ", ncol(expr), "\n")
  cat("特征数(细胞类型): ", nrow(expr), "\n")
  cat("分型模式: ", mode, "\n")
  
  # 4) 构造输入到CCP的矩阵（features x samples）
  #    - rank：Pij
  #    - rawz：row-z
  if (mode == "rank") {
    n_cell <- nrow(expr)
    rank_mat <- apply(expr, 2, function(x) rank(-x, ties.method = "min"))  # 1=最高
    Pij <- (n_cell - rank_mat + 1) / n_cell * 100
    data_mat <- as.matrix(Pij)
  } else {
    row_zscore <- function(mat) {
      mu <- rowMeans(mat, na.rm = TRUE)
      sdv <- apply(mat, 1, sd, na.rm = TRUE)
      sdv[is.na(sdv) | sdv == 0] <- 1
      z <- sweep(mat, 1, mu, "-")
      z <- sweep(z, 1, sdv, "/")
      z[is.na(z)] <- 0
      z
    }
    data_mat <- row_zscore(expr)
  }
  
  # 5) 运行一致性聚类（CCP自带画图有时会报ylim问题，所以做一个兜底）
  cat("开始一致性聚类分析...\n")
  
  run_ccp_once <- function(plot_mode) {
    ConsensusClusterPlus(
      d = data_mat,
      maxK = maxK,
      reps = reps,
      pItem = pItem,
      pFeature = pFeature,
      clusterAlg = clusterAlg,
      distance = distance,
      innerLinkage = "average",
      finalLinkage = "average",
      title = output_dir,
      plot = plot_mode,     # "png" / "pdf" / "none"
      seed = seed,
      verbose = TRUE
    )
  }
  
  set.seed(seed)
  results <- tryCatch(
    run_ccp_once("png"),
    error = function(e) {
      cat("\n[警告] CCP自带画图阶段报错，改为 plot='none' 继续输出结果。\n")
      cat("报错信息：", e$message, "\n\n")
      run_ccp_once("none")
    }
  )
  
  # 6) 保存所有k的聚类标签（同你原函数习惯）
  cat("保存聚类结果...\n")
  k_min <- 2
  for (k in k_min:maxK) {
    cluster_labels <- results[[k]]$consensusClass
    output_file <- file.path(output_dir, paste0("cluster_assignments_k", k, ".csv"))
    write.csv(
      data.frame(Sample = names(cluster_labels), Cluster = as.integer(cluster_labels)),
      output_file,
      row.names = FALSE
    )
    cat("保存: ", basename(output_file), "\n")
  }
  
  # 7) 保存完整对象
  rds_file <- file.path(output_dir, "consensus_results.rds")
  saveRDS(results, rds_file)
  cat("保存完整结果对象: ", basename(rds_file), "\n")
  
  # 8) 额外：稳定性表（PAC + MeanItemConsensus）+ bestK + 简单图
  calc_PAC <- function(consensus_mat, lower = 0.1, upper = 0.9) {
    m <- consensus_mat
    diag(m) <- NA
    v <- as.numeric(m)
    v <- v[!is.na(v)]
    if (length(v) == 0) return(NA_real_)
    mean(v > lower & v < upper)
  }
  
  stab <- data.frame(
    k = k_min:maxK,
    PAC_0.1_0.9 = NA_real_,
    MeanItemConsensus = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (k in k_min:maxK) {
    cm <- results[[k]]$consensusMatrix
    stab$PAC_0.1_0.9[stab$k == k] <- calc_PAC(cm, 0.1, 0.9)
    if (!is.null(results[[k]]$itemConsensus)) {
      stab$MeanItemConsensus[stab$k == k] <- mean(results[[k]]$itemConsensus, na.rm = TRUE)
    }
  }
  
  stab_file <- file.path(output_dir, "stability_PAC_itemConsensus.csv")
  write.csv(stab, stab_file, row.names = FALSE)
  cat("保存: ", basename(stab_file), "\n")
  
  best_k <- stab$k[which.min(stab$PAC_0.1_0.9)]
  best_file <- file.path(output_dir, "bestK_by_PAC.csv")
  write.csv(data.frame(BestK = best_k), best_file, row.names = FALSE)
  cat("保存: ", basename(best_file), " (bestK=", best_k, ")\n", sep = "")
  
  # 自定义图：PAC曲线 + bestK簇大小（不依赖CCP自带图，稳）
  png(file.path(output_dir, "Plot_PAC_curve.png"), width = 1200, height = 900, res = 150)
  plot(stab$k, stab$PAC_0.1_0.9, type = "b", pch = 16,
       xlab = "k", ylab = "PAC (0.1-0.9)", main = paste0("PAC curve (", mode, ") - lower is better"))
  abline(v = best_k, lty = 2)
  dev.off()
  
  best_labels <- results[[best_k]]$consensusClass
  tab <- table(best_labels)
  png(file.path(output_dir, "Plot_bestK_cluster_sizes.png"), width = 1200, height = 900, res = 150)
  barplot(tab, xlab = "Cluster", ylab = "Sample count",
          main = paste0("bestK=", best_k, " cluster sizes (", mode, ")"))
  dev.off()
  
  # 9) 摘要输出（同你原函数的风格）
  used <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  result_summary <- paste(
    "\n===== 免疫一致性聚类完成! =====",
    paste0("\n模式: ", mode),
    paste0("\n输入文件: ", input_path),
    paste0("\n输出目录: ", output_dir),
    "\n包含结果:",
    paste0("- ", maxK-1, " 个聚类标签文件 (k=2-", maxK, ")"),
    "- consensus_results.rds",
    "- stability_PAC_itemConsensus.csv + bestK_by_PAC.csv",
    "- Plot_PAC_curve.png + Plot_bestK_cluster_sizes.png",
    paste0("\n分析耗时: ", round(used, 1), " 秒"),
    sep = "\n"
  )
  cat(result_summary)
  
  invisible(TRUE)
}


# ============================================================
# 批量运行（按你原来的datasets + lapply逻辑）
# ============================================================

base_input_dir  <- "input"
base_output_dir <- "output/immune_ccp"
input_file_name <- "ssgsea_scores_28cells.csv"

datasets <- list(
  list(
    input_dir  = base_input_dir,
    input_file = input_file_name,
    output_dir = file.path(base_output_dir, "A_rank_based_CCP"),
    mode = "rank",
    maxK = 10,
    reps = 1000,
    pItem = 0.8,
    pFeature = 1,
    clusterAlg = "hc",
    distance = "pearson",
    seed = 100
  ),
  list(
    input_dir  = base_input_dir,
    input_file = input_file_name,
    output_dir = file.path(base_output_dir, "B_raw_zscore_based_CCP"),
    mode = "rawz",
    maxK = 10,
    reps = 1000,
    pItem = 0.8,
    pFeature = 1,
    clusterAlg = "hc",
    distance = "pearson",
    seed = 100
  )
)

# 批量运行
total_start <- Sys.time()
ok <- lapply(datasets, function(params) {
  cat("\n===== 开始分析: ", params$mode, " =====\n", sep = "")
  tryCatch({
    do.call(run_immune_consensus_clustering, params)
    cat("\n===== 分析成功完成 =====\n\n")
    return(TRUE)
  }, error = function(e) {
    cat("\n错误: ", e$message, "\n\n")
    return(FALSE)
  })
})
cat("总耗时:", round(as.numeric(difftime(Sys.time(), total_start, units = "secs")), 1), "秒\n")
cat("总输出目录:", base_output_dir, "\n")

