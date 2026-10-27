# =========================================================
# 对应结果：Figure 6B（基于样本内免疫细胞排名百分位 Pij 的相对免疫结构热图，按 ImmuneRank_Cluster 分组）
# 输入：input/ssgsea_scores_28cells.csv
#       input/RankBased_Pearson_04_clusters_bestK.csv
# 输出：output/fig6B/04_structure_Pij_heatmap_optimized.png，output/fig6B/pdf/*.pdf，output/fig6B/svg/*.svg
# 运行：在模块目录下 Rscript code/05_fig6B_rank_heatmap.R
# =========================================================

# ===============================
# 相对结构聚类可视化（最终版）
# 说明：
# 1）删除热图左侧 Difference 行注释及右侧对应图例
# 2）热图标题改为：Relative immune structure clustering
# 3）图例名称改为：ImmuneRank_Cluster
# 4）输出路径改为：26.4.24
# ===============================

rm(list = ls())

# ---------- 1. 加载包 ----------
pkgs <- c("data.table", "ggplot2", "pheatmap", "svglite")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

# ---------- 2. 路径 ----------
expr_file <- "input/ssgsea_scores_28cells.csv"

rank_dir <- "input"
vis_dir  <- "output/fig6B"

svg_dir <- file.path(vis_dir, "svg")
pdf_dir <- file.path(vis_dir, "pdf")

dir.create(vis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(svg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

# rank_sil_file <- file.path(rank_dir, "RankBased_Pearson_03_silhouette_scores.csv")
rank_cl_file  <- file.path(rank_dir, "RankBased_Pearson_04_clusters_bestK.csv")
# rank_hc_file  <- file.path(rank_dir, "RankBased_Pearson_05_hclust_object.rds")

# ---------- 3. 读取原始矩阵 ----------
expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
expr <- as.matrix(expr)
mode(expr) <- "numeric"

n_cell <- nrow(expr)

# ---------- 4. 构建 Rij / Pij ----------
rank_mat <- apply(expr, 2, function(x) rank(-x, ties.method = "min"))
rank_mat <- as.matrix(rank_mat)
rownames(rank_mat) <- rownames(expr)
colnames(rank_mat) <- colnames(expr)

Pij <- (n_cell - rank_mat + 1) / n_cell * 100
Pij <- as.matrix(Pij)
rownames(Pij) <- rownames(expr)
colnames(Pij) <- colnames(expr)

# write.csv(rank_mat, file.path(vis_dir, "00_Rij_rank_matrix.csv"), row.names = TRUE)
# write.csv(Pij,      file.path(vis_dir, "01_Pij_relative_structure_matrix.csv"), row.names = TRUE)

# ---------- 5. 读取聚类结果 ----------
# rank_sil <- read.csv(rank_sil_file, stringsAsFactors = FALSE)
rank_cl  <- read.csv(rank_cl_file, stringsAsFactors = FALSE)
# rank_hc  <- readRDS(rank_hc_file)

colnames(rank_cl) <- c("Sample", "ImmuneRank_Cluster")
rank_cl$ImmuneRank_Cluster <- factor(rank_cl$ImmuneRank_Cluster)

# ---------- 6. silhouette 曲线 ----------
# p_rank_sil <- ggplot(rank_sil, aes(x = k, y = mean_silhouette)) +
#   geom_line(linewidth = 0.9) +
#   geom_point(size = 3) +
#   geom_point(
#     data = rank_sil[which.max(rank_sil$mean_silhouette), , drop = FALSE],
#     aes(x = k, y = mean_silhouette),
#     size = 4
#   ) +
#   labs(
#     title = "Relative structure clustering (Pij + Pearson)",
#     x = "Number of clusters (k)",
#     y = "Mean silhouette width"
#   ) +
#   theme_bw(base_size = 14)
#
# ggsave(file.path(vis_dir, "02_structure_silhouette_curve.png"),
#        p_rank_sil, width = 8, height = 6, dpi = 300)
# ggsave(file.path(svg_dir, "02_structure_silhouette_curve.svg"),
#        p_rank_sil, width = 8, height = 6)
# ggsave(file.path(pdf_dir, "02_structure_silhouette_curve.pdf"),
#        p_rank_sil, width = 8, height = 6)

# ---------- 7. dendrogram ----------
# png(file.path(vis_dir, "03_structure_dendrogram.png"), width = 3000, height = 1600, res = 300)
# plot(rank_hc,
#      main = "Dendrogram of relative structure clustering (Pij + Pearson)",
#      xlab = "", sub = "", cex = 0.7)
# rect.hclust(rank_hc, k = 2, border = 2:3)
# dev.off()
#
# pdf(file.path(pdf_dir, "03_structure_dendrogram.pdf"), width = 12, height = 6)
# plot(rank_hc,
#      main = "Dendrogram of relative structure clustering (Pij + Pearson)",
#      xlab = "", sub = "", cex = 0.7)
# rect.hclust(rank_hc, k = 2, border = 2:3)
# dev.off()
#
# svglite::svglite(file.path(svg_dir, "03_structure_dendrogram.svg"), width = 12, height = 6)
# plot(rank_hc,
#      main = "Dendrogram of relative structure clustering (Pij + Pearson)",
#      xlab = "", sub = "", cex = 0.7)
# rect.hclust(rank_hc, k = 2, border = 2:3)
# dev.off()

# ---------- 8. 优化热图排序 ----------
samples_c1 <- rank_cl$Sample[rank_cl$ImmuneRank_Cluster == levels(rank_cl$ImmuneRank_Cluster)[1]]
samples_c2 <- rank_cl$Sample[rank_cl$ImmuneRank_Cluster == levels(rank_cl$ImmuneRank_Cluster)[2]]

get_within_order <- function(mat_sub) {
  if (ncol(mat_sub) <= 2) return(colnames(mat_sub))
  hc <- hclust(dist(t(mat_sub)), method = "average")
  return(colnames(mat_sub)[hc$order])
}

ord_c1 <- get_within_order(Pij[, samples_c1, drop = FALSE])
ord_c2 <- get_within_order(Pij[, samples_c2, drop = FALSE])

sample_order_rank <- c(ord_c1, ord_c2)

# 行顺序仍按两组均值差排序，但不再显示 Difference 注释
mean_c1 <- rowMeans(Pij[, samples_c1, drop = FALSE], na.rm = TRUE)
mean_c2 <- rowMeans(Pij[, samples_c2, drop = FALSE], na.rm = TRUE)
row_diff <- mean_c1 - mean_c2
row_order <- names(sort(row_diff, decreasing = TRUE))

Pij_ord <- Pij[row_order, sample_order_rank, drop = FALSE]

ann_col_rank <- data.frame(
  ImmuneRank_Cluster = rank_cl$ImmuneRank_Cluster[match(sample_order_rank, rank_cl$Sample)]
)
rownames(ann_col_rank) <- sample_order_rank

cluster_levels <- levels(rank_cl$ImmuneRank_Cluster)
cluster_palette <- colorRampPalette(c("#5DADE2", "#E74C3C"))(length(cluster_levels))
ann_colors <- list(
  ImmuneRank_Cluster = setNames(cluster_palette, cluster_levels)
)

heat_colors <- colorRampPalette(c("#4A90E2", "white", "#E74C3C"))(100)
heat_breaks <- seq(0, 100, length.out = length(heat_colors) + 1)

# ---------- 9. 优化版热图 ----------
draw_pij_heatmap <- function(filename = NULL, device = c("png", "pdf", "svg")) {
  device <- match.arg(device)
  
  if (device == "png") {
    png(filename, width = 3600, height = 2600, res = 300, bg = "white")
  } else if (device == "pdf") {
    pdf(filename, width = 13, height = 9, bg = "white", useDingbats = FALSE)
  } else if (device == "svg") {
    svglite::svglite(filename, width = 13, height = 9, bg = "white")
  }
  
  pheatmap(
    Pij_ord,
    color = heat_colors,
    breaks = heat_breaks,
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = ann_col_rank,
    annotation_colors = ann_colors,
    show_colnames = FALSE,
    fontsize_row = 10,
    fontsize_col = 8,
    border_color = NA,
    main = "Relative immune structure clustering",
    gaps_col = length(ord_c1),
    legend_breaks = c(0, 25, 50, 75, 100),
    legend_labels = c("0", "25", "50", "75", "100")
  )
  
  dev.off()
}

draw_pij_heatmap(file.path(vis_dir, "04_structure_Pij_heatmap_optimized.png"), "png")
draw_pij_heatmap(file.path(pdf_dir, "04_structure_Pij_heatmap_optimized.pdf"), "pdf")
draw_pij_heatmap(file.path(svg_dir, "04_structure_Pij_heatmap_optimized.svg"), "svg")

# ---------- 10. 差异免疫细胞统计 ----------
# diff_pij <- lapply(rownames(Pij), function(cell) {
#   x <- as.numeric(Pij[cell, samples_c1])
#   y <- as.numeric(Pij[cell, samples_c2])
#   wt <- wilcox.test(x, y)
#   data.frame(
#     Cell = cell,
#     Cluster1_Mean = mean(x),
#     Cluster2_Mean = mean(y),
#     Cluster1_Median = median(x),
#     Cluster2_Median = median(y),
#     MeanDiff_1_minus_2 = mean(x) - mean(y),
#     Pvalue = wt$p.value,
#     stringsAsFactors = FALSE
#   )
# })
#
# diff_pij <- rbindlist(diff_pij)
# diff_pij$FDR <- p.adjust(diff_pij$Pvalue, method = "BH")
# diff_pij <- diff_pij[order(diff_pij$FDR, -abs(diff_pij$MeanDiff_1_minus_2)), ]
#
# write.csv(diff_pij,
#           file.path(vis_dir, "05_structure_differential_immune_cells_Pij.csv"),
#           row.names = FALSE)

cat("\n==============================\n")
cat("相对结构聚类优化版可视化完成！\n")
cat("普通输出目录：", vis_dir, "\n")
cat("SVG目录：", svg_dir, "\n")
cat("PDF目录：", pdf_dir, "\n")
cat("==============================\n")
