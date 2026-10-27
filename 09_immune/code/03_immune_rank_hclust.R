# =========================================================
# 对应结果：Figure 6B 使用的免疫排名分型（ImmuneRank_Cluster；样本内排名百分位 Pij + 1−Pearson 距离 + hclust average
#           silhouette 选 k，最佳 k = 2）
# 输入：input/ssgsea_scores_28cells.csv
# 输出：output/immune_rank_hclust/A_Rank_Pearson/RankBased_Pearson_01–06_*（其中 RankBased_Pearson_04_clusters_bestK.csv
#       即 05_fig6B_rank_heatmap.R 使用的 input/RankBased_Pearson_04_clusters_bestK.csv）
#       output/immune_rank_hclust/B_Rank_Spearman/*（Spearman 距离的对照，论文未使用）
#       output/immune_rank_hclust/Summary/*
# 运行：在模块目录下 Rscript code/03_immune_rank_hclust.R
# =========================================================

## ============================================================
## 相对结构表达矩阵 Pij：Pearson vs Spearman 距离聚类（都跑）
## Pij：由 ssGSEA(cell x sample) 在“样本内”rank -> percentile 得到
## 距离：1 - cor (Pearson / Spearman)
## 聚类：hclust + average linkage
## k选择：silhouette (k=2..10)
## 输入：12.16新作的免疫细胞表达.csv (cell x sample)
## 输出：指定输出目录
## ============================================================

## ===============================
## 0. 路径与参数
## ===============================
input_dir  <- "input"
output_dir <- "output/immune_rank_hclust"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

file_in <- file.path(input_dir, "ssgsea_scores_28cells.csv")
if (!file.exists(file_in)) stop("输入文件不存在: ", file_in)

k_min <- 2
k_max <- 10

cat("输入文件:", file_in, "\n")
cat("输出目录:", output_dir, "\n\n")

## ===============================
## 1. 包
## ===============================
if (!requireNamespace("cluster", quietly = TRUE)) install.packages("cluster")
library(cluster)

## ===============================
## 2. 读取 ssGSEA (cell x sample)
## ===============================
expr <- read.csv(file_in, row.names = 1, check.names = FALSE)
expr <- as.matrix(expr)
mode(expr) <- "numeric"

n_cell   <- nrow(expr)
n_sample <- ncol(expr)

cat("读取完成：细胞类型数 =", n_cell, "; 样本数 =", n_sample, "\n\n")

## ===============================
## 3. 构建 Rij / Pij（样本内相对结构）
## ===============================
# Rij: rank(-Sij), 1=最高
rank_mat <- apply(expr, 2, function(x) rank(-x, ties.method = "min"))
rownames(rank_mat) <- rownames(expr)
colnames(rank_mat) <- colnames(expr)

# Pij: (n - Rij + 1)/n * 100
Pij <- (n_cell - rank_mat + 1) / n_cell * 100
Pij <- as.matrix(Pij)
rownames(Pij) <- rownames(expr)
colnames(Pij) <- colnames(expr)

# 样本 x 细胞
X_rank <- t(Pij)

## （可选）强度指标：原始 ssGSEA 在该样本的平均值（粗略代表“总体浸润强弱”）
intensity_mean <- colMeans(expr, na.rm = TRUE)

## ===============================
## 4. 核心函数：cor -> dist(1-cor) -> hclust -> silhouette选k -> 输出
## ===============================
cluster_rank_by_cor <- function(X_samples_by_features,
                                cor_method = c("pearson","spearman"),
                                linkage = "average",
                                k_min = 2, k_max = 10,
                                prefix = "RankBased_Pearson",
                                outdir = ".") {
  cor_method <- match.arg(cor_method)
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # X 是 samples x features
  cor_mat <- cor(t(X_samples_by_features), method = cor_method, use = "pairwise.complete.obs")
  cor_mat[is.na(cor_mat)] <- 0
  dist_mat <- as.dist(1 - cor_mat)
  
  hc <- hclust(dist_mat, method = linkage)
  
  k_max_eff <- min(k_max, nrow(X_samples_by_features) - 1)
  if (k_max_eff < k_min) stop("样本数过少，无法 silhouette 选 k。")
  
  sil_scores <- data.frame(k = k_min:k_max_eff, mean_silhouette = NA_real_)
  for (kk in k_min:k_max_eff) {
    cl <- cutree(hc, k = kk)
    sil <- silhouette(cl, dist_mat)
    sil_scores$mean_silhouette[sil_scores$k == kk] <- mean(sil[, 3])
  }
  
  best_k <- sil_scores$k[which.max(sil_scores$mean_silhouette)]
  clusters_best <- cutree(hc, k = best_k)
  
  # 输出
  write.csv(as.matrix(cor_mat),
            file.path(outdir, paste0(prefix, "_01_sample_correlation_", cor_method, ".csv")),
            row.names = TRUE)
  
  write.csv(as.matrix(dist_mat),
            file.path(outdir, paste0(prefix, "_02_sample_distance_1minus", cor_method, ".csv")),
            row.names = TRUE)
  
  write.csv(sil_scores,
            file.path(outdir, paste0(prefix, "_03_silhouette_scores.csv")),
            row.names = FALSE)
  
  write.csv(data.frame(Sample = names(clusters_best),
                       Cluster = as.integer(clusters_best),
                       stringsAsFactors = FALSE),
            file.path(outdir, paste0(prefix, "_04_clusters_bestK.csv")),
            row.names = FALSE)
  
  saveRDS(hc, file.path(outdir, paste0(prefix, "_05_hclust_object.rds")))
  
  # 输出所有k分组
  allk <- lapply(k_min:k_max_eff, function(kk) cutree(hc, k = kk))
  names(allk) <- paste0("k", k_min:k_max_eff)
  allk_df <- data.frame(Sample = rownames(X_samples_by_features),
                        do.call(cbind, allk),
                        check.names = FALSE)
  write.csv(allk_df,
            file.path(outdir, paste0(prefix, "_06_clusters_allK.csv")),
            row.names = FALSE)
  
  cat(prefix, "完成：cor =", cor_method, "; best_k =", best_k, "\n")
  return(list(best_k = best_k, clusters = clusters_best, hc = hc, dist = dist_mat, cor = cor_mat, sil = sil_scores))
}

## ===============================
## 5. 跑两套：Pearson vs Spearman（都在 Pij 上）
## ===============================
out_pearson  <- file.path(output_dir, "A_Rank_Pearson")
out_spearman <- file.path(output_dir, "B_Rank_Spearman")

res_p <- cluster_rank_by_cor(
  X_samples_by_features = X_rank,
  cor_method = "pearson",
  linkage = "average",
  k_min = k_min, k_max = k_max,
  prefix = "RankBased_Pearson",
  outdir = out_pearson
)

res_s <- cluster_rank_by_cor(
  X_samples_by_features = X_rank,
  cor_method = "spearman",
  linkage = "average",
  k_min = k_min, k_max = k_max,
  prefix = "RankBased_Spearman",
  outdir = out_spearman
)

## ===============================
## 6. 总结输出（强烈建议你保留）
## ===============================
summary_dir <- file.path(output_dir, "Summary")
if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)

# 6.1 best_k 汇总
bestk_df <- data.frame(
  Method = c("Rank(Pij) + 1-Pearson + average", "Rank(Pij) + 1-Spearman + average"),
  BestK  = c(res_p$best_k, res_s$best_k),
  stringsAsFactors = FALSE
)
write.csv(bestk_df, file.path(summary_dir, "01_bestK_summary.csv"), row.names = FALSE)

# 6.2 分群标签与强度（mean ssGSEA）合并输出
common_samples <- intersect(names(res_p$clusters), names(res_s$clusters))
df_lab <- data.frame(
  Sample = common_samples,
  RankPearson_Cluster  = as.integer(res_p$clusters[common_samples]),
  RankSpearman_Cluster = as.integer(res_s$clusters[common_samples]),
  Intensity_Mean_ssGSEA = as.numeric(intensity_mean[common_samples]),
  stringsAsFactors = FALSE
)
write.csv(df_lab, file.path(summary_dir, "02_intensity_and_clusters_rankPij.csv"), row.names = FALSE)

# 6.3 “谁更能分强弱”：用 Intensity_Mean_ssGSEA 做组间差异检验（Kruskal更稳）
test_one <- function(cluster_vec, intensity_vec) {
  # cluster_vec: named integer
  d <- data.frame(Cluster = as.factor(cluster_vec), Intensity = intensity_vec)
  # 组概况
  stat <- aggregate(Intensity ~ Cluster, d, function(x) c(mean=mean(x), median=median(x), n=length(x)))
  stat <- do.call(data.frame, stat)
  colnames(stat) <- c("Cluster","Mean","Median","N")
  # 检验
  kw <- kruskal.test(Intensity ~ Cluster, data = d)
  list(stat = stat, p = kw$p.value)
}

tp <- test_one(df_lab$RankPearson_Cluster,  df_lab$Intensity_Mean_ssGSEA)
ts <- test_one(df_lab$RankSpearman_Cluster, df_lab$Intensity_Mean_ssGSEA)

sep_sum <- rbind(
  cbind(Method="Rank(Pij)+Pearson",  tp$stat, KW_pvalue=tp$p),
  cbind(Method="Rank(Pij)+Spearman", ts$stat, KW_pvalue=ts$p)
)
write.csv(sep_sum, file.path(summary_dir, "03_strength_separation_summary_rankPij.csv"), row.names = FALSE)

cat("\n==============================\n")
cat("全部完成！输出在：", output_dir, "\n")
cat("A_Rank_Pearson：Pij 结构 + Pearson 距离\n")
cat("B_Rank_Spearman：Pij 结构 + Spearman 距离\n")
cat("Summary：best_k + 分群对照 + 强度分离检验\n")
cat("==============================\n")

