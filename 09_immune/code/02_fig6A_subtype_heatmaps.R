rm(list = ls())

# =========================================================
# 对应结果：Figure 6A（三种分子亚型 Molecular_Cluster1/2 中 28 种免疫细胞 ssGSEA 得分热图，行 Z-score）
# 输入：input/ssgsea_scores_28cells.csv（01_ssgsea.R 的得分矩阵）
#       input/subtype_<pcd>.csv（分子亚型分组，来自 07_molecular_subtypes 模块）
# 输出：output/fig6A/<pcd>/亚型_免疫细胞热图_无聚类树_强化分组版.pdf / .svg
#       output/fig6A/<pcd>/热图用_排序后标准化表达矩阵.csv、热图用_排序后样本亚型信息.csv、
#       免疫细胞按Cluster2_vs_Cluster1差值排序结果.csv
# 运行：在模块目录下 Rscript code/02_fig6A_subtype_heatmaps.R
# =========================================================
# 失巢死亡亚型免疫细胞热图
# 要求：
# 1）不显示聚类树
# 2）样本列按 Cluster 1 / 2 手动排序
# 3）免疫细胞行按组间差值排序，让两组区分更明显
# 4）导出纯白背景 PDF + SVG
# =========================================================

# -----------------------------
# 1. 加载/安装包
# -----------------------------
pkgs <- c("data.table", "pheatmap", "svglite")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

# -----------------------------
# 2. 路径设置
# -----------------------------
input_dir  <- "input"

expr_file  <- file.path(input_dir, "ssgsea_scores_28cells.csv")

pcd_settings <- list(
  anoikis = list(
    group_file  = "subtype_anoikis.csv",
    ann_colors  = c("1" = "#5DADE2", "2" = "#E74C3C"),
    heat_colors = c("#4A90E2", "white", "#E74C3C")
  ),
  apoptosis = list(
    group_file  = "subtype_apoptosis.csv",
    ann_colors  = c("1" = "#7E57C2", "2" = "#EC407A"),
    heat_colors = c("#5B4EA3", "white", "#D81B60")
  ),
  ferroptosis = list(
    group_file  = "subtype_ferroptosis.csv",
    ann_colors  = c("1" = "#6FB39A", "2" = "#D97C6C"),
    heat_colors = c("#5FA58B", "white", "#D46A5A")
  )
)

for (pcd in names(pcd_settings)) {  # 按 PCD 循环
  cat("\n========== ", pcd, " ==========\n")

  output_dir <- file.path("output/fig6A", pcd)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  group_file <- file.path(input_dir, pcd_settings[[pcd]]$group_file)

  # -----------------------------
  # 3. 检查文件是否存在
  # -----------------------------
  if (!file.exists(expr_file)) {
    stop(paste0("表达矩阵文件不存在：", expr_file))
  }
  if (!file.exists(group_file)) {
    stop(paste0("亚型文件不存在：", group_file))
  }

  # -----------------------------
  # 4. 读取数据
  # -----------------------------
  expr <- fread(expr_file, data.table = FALSE, check.names = FALSE)
  rownames(expr) <- expr[, 1]
  expr <- expr[, -1, drop = FALSE]

  expr_mat <- as.matrix(expr)
  mode(expr_mat) <- "numeric"

  group_df <- fread(group_file, data.table = FALSE)

  if (!all(c("Sample", "Cluster") %in% colnames(group_df))) {
    stop("亚型文件中必须包含列名：Sample 和 Cluster")
  }

  # -----------------------------
  # 5. 对齐样本
  # -----------------------------
  common_samples <- intersect(colnames(expr_mat), group_df$Sample)

  if (length(common_samples) == 0) {
    stop("表达矩阵和亚型文件没有共同样本，请检查样本名是否一致。")
  }

  cat("共同样本数：", length(common_samples), "\n")

  expr_mat <- expr_mat[, common_samples, drop = FALSE]
  group_df <- group_df[group_df$Sample %in% common_samples, , drop = FALSE]
  group_df <- group_df[match(colnames(expr_mat), group_df$Sample), , drop = FALSE]

  # 固定亚型顺序：1 在前，2 在后
  group_df$Cluster <- as.character(group_df$Cluster)
  group_df$Cluster <- factor(group_df$Cluster, levels = c("1", "2"))

  # 去掉没有有效分组的样本
  valid_idx <- !is.na(group_df$Cluster)
  expr_mat  <- expr_mat[, valid_idx, drop = FALSE]
  group_df  <- group_df[valid_idx, , drop = FALSE]

  if (ncol(expr_mat) == 0) {
    stop("对齐后没有可用样本，请检查亚型信息。")
  }

  # -----------------------------
  # 6. 手动调整样本列顺序
  # -----------------------------
  ord_col <- order(group_df$Cluster, group_df$Sample)
  expr_mat <- expr_mat[, ord_col, drop = FALSE]
  group_df <- group_df[ord_col, , drop = FALSE]

  # -----------------------------
  # 7. 行标准化（Z-score）
  # -----------------------------
  expr_scaled <- t(scale(t(expr_mat)))

  # 去掉标准化后全 NA 的行
  expr_scaled <- expr_scaled[complete.cases(expr_scaled), , drop = FALSE]

  if (nrow(expr_scaled) == 0) {
    stop("标准化后没有可用免疫细胞行，可能某些行全为常数。")
  }

  # -----------------------------
  # 8. 手动调整细胞行顺序
  #    按 Cluster2 - Cluster1 的均值差排序
  # -----------------------------
  cluster1_samples <- group_df$Sample[group_df$Cluster == "1"]
  cluster2_samples <- group_df$Sample[group_df$Cluster == "2"]

  if (length(cluster1_samples) == 0 || length(cluster2_samples) == 0) {
    stop("Cluster 1 或 Cluster 2 样本数为 0，请检查 Cluster 是否正确。")
  }

  cluster1_mean <- rowMeans(expr_scaled[, cluster1_samples, drop = FALSE], na.rm = TRUE)
  cluster2_mean <- rowMeans(expr_scaled[, cluster2_samples, drop = FALSE], na.rm = TRUE)

  diff_score <- cluster2_mean - cluster1_mean

  row_info <- data.frame(
    ImmuneCell = rownames(expr_scaled),
    Cluster1Mean = cluster1_mean,
    Cluster2Mean = cluster2_mean,
    Diff = diff_score,
    stringsAsFactors = FALSE
  )

  row_info <- row_info[order(row_info$Diff, decreasing = TRUE), , drop = FALSE]
  expr_scaled <- expr_scaled[row_info$ImmuneCell, , drop = FALSE]

  # -----------------------------
  # 9. 列注释
  # -----------------------------
  annotation_col <- data.frame(
    Cluster = group_df$Cluster
  )
  rownames(annotation_col) <- group_df$Sample

  ann_colors <- list(
    Cluster = pcd_settings[[pcd]]$ann_colors
  )

  heat_colors <- colorRampPalette(pcd_settings[[pcd]]$heat_colors)(100)  # 同上

  # 在 Cluster1 和 Cluster2 之间加一个分隔位置
  gap_col_pos <- sum(group_df$Cluster == "1")

  # -----------------------------
  # 10. 导出排序信息
  # -----------------------------
  write.csv(
    data.frame(ImmuneCell = rownames(expr_scaled), expr_scaled, check.names = FALSE),
    file = file.path(output_dir, "热图用_排序后标准化表达矩阵.csv"),
    row.names = FALSE
  )

  write.csv(
    group_df,
    file = file.path(output_dir, "热图用_排序后样本亚型信息.csv"),
    row.names = FALSE
  )

  write.csv(
    row_info,
    file = file.path(output_dir, "免疫细胞按Cluster2_vs_Cluster1差值排序结果.csv"),
    row.names = FALSE
  )

  # -----------------------------
  # 11. 绘制 PDF（纯白背景）
  # -----------------------------
  pdf(
    file.path(output_dir, "亚型_免疫细胞热图_无聚类树_强化分组版.pdf"),
    width = 12,
    height = 9,
    bg = "white"
  )

  pheatmap(
    expr_scaled,
    color = heat_colors,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 10,
    fontsize_col = 8,
    border_color = NA,
    scale = "none",
    gaps_col = gap_col_pos,
    main = "Immune Cell Expression Heatmap by Cluster"
  )

  dev.off()

  # -----------------------------
  # 12. 绘制 SVG（纯白背景）
  # -----------------------------
  svglite(
    file.path(output_dir, "亚型_免疫细胞热图_无聚类树_强化分组版.svg"),
    width = 12,
    height = 9,
    bg = "white"
  )

  pheatmap(
    expr_scaled,
    color = heat_colors,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 10,
    fontsize_col = 8,
    border_color = NA,
    scale = "none",
    gaps_col = gap_col_pos,
    main = "Immune Cell Expression Heatmap by Cluster"
  )

  dev.off()

  cat("完成：热图已输出到\n", output_dir, "\n")
}
