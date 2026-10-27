# =========================================================
# 对应结果：Figure 6C（失巢凋亡分子亚型与免疫排名亚型 ImmuneRank_Cluster 的对应桑基图）
#           正文中三种 PCD 分子亚型与免疫排名亚型的对应比例（anoikis 85.29%/52.78%，apoptosis 93.10%/53.66%
#           ferroptosis 84.38%/50.00%），见各 PCD 的 *_Summary_Mapping.csv 中 Ratio_in_LeftGroup 列
# 输入：input/subtype_<pcd>.csv
#       input/immune_rank_clusters_CCP_k2.csv
#       input/risk_group_<pcd>.csv、input/immune_matrix_clusters_CCP_k2.csv（只用于已注释的 A1/A2/B2 三组图，论文未使用）
# 输出：output/fig6C/<pcd>/B1_MolecularSubtype_vs_ImmuneRankSubtype/ 下的
#       *_Sankey_Beautified.pdf/.svg、无标题/*_Sankey_NoTitle.pdf/.svg、*_Sample_Level_Mapping.csv、*_Summary_Mapping.csv
# 运行：在模块目录下 Rscript code/06_fig6C_sankey.R
# =========================================================
# Beautified Sankey Plots for Group Mapping
# 失巢死亡版本 + 统一字体版
# 输出：
# 1) PDF
# 2) SVG
# 3) 样本级对应表 CSV
# 4) 汇总对应表 CSV
# 5) 无标题 PDF + SVG
# 共生成4组：
# A1 RiskGroup vs ImmuneRankSubtype
# A2 RiskGroup vs ImmuneMatrixSubtype
# B1 MolecularSubtype vs ImmuneRankSubtype
# B2 MolecularSubtype vs ImmuneMatrixSubtype
# =========================================================

rm(list = ls())

# -----------------------------
# 1. 加载包
# -----------------------------
pkgs <- c("dplyr", "ggplot2", "ggalluvial", "svglite")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

# -----------------------------
# 2. 路径
# -----------------------------
input_dir  <- "input"
output_root <- "output/fig6C"

pcd_files <- list(
  anoikis     = c(risk = "risk_group_anoikis.csv",     subtype = "subtype_anoikis.csv"),
  apoptosis   = c(risk = "risk_group_apoptosis.csv",   subtype = "subtype_apoptosis.csv"),
  ferroptosis = c(risk = "risk_group_ferroptosis.csv", subtype = "subtype_ferroptosis.csv")
)

immune_rank_file   <- file.path(input_dir, "immune_rank_clusters_CCP_k2.csv")
immune_matrix_file <- file.path(input_dir, "immune_matrix_clusters_CCP_k2.csv")

# -----------------------------
# 3. 统一字体参数（参考前面小论文版）
# -----------------------------
title_size      <- 18
subtitle_size   <- 12
axis_title_size <- 13
axis_text_x_sz  <- 13
axis_text_y_sz  <- 12
stratum_text_sz <- 4.6

# -----------------------------
# 4. 通用读取函数
# -----------------------------
read_group_file <- function(file, group_name,
                            sample_candidates = c("Sample", "sample", "SAMPLE", "ID", "id"),
                            group_candidates = c("Group", "group", "Cluster", "cluster", "Subtype", "subtype",
                                                 "RiskGroup", "riskgroup", "Risk_Group", "risk_group",
                                                 "风险组", "亚型")) {
  dat <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(dat) <- trimws(colnames(dat))

  sample_hit <- sample_candidates[sample_candidates %in% colnames(dat)]
  if (length(sample_hit) == 0) {
    stop("文件中缺少样本列（Sample）：", file)
  }

  group_hit <- group_candidates[group_candidates %in% colnames(dat)]
  if (length(group_hit) == 0) {
    stop("文件中缺少分组列：", file)
  }

  out <- dat[, c(sample_hit[1], group_hit[1])]
  colnames(out) <- c("Sample", group_name)

  out$Sample <- trimws(as.character(out$Sample))
  out[[group_name]] <- trimws(as.character(out[[group_name]]))

  out
}

# -----------------------------
# 5. 标签标准化函数
# -----------------------------
normalize_risk <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    x %in% c("High", "HIGH", "high", "高风险", "高风险组", "Risk-High", "risk-high") ~ "High risk",
    x %in% c("Low", "LOW", "low", "低风险", "低风险组", "Risk-Low", "risk-low") ~ "Low risk",
    TRUE ~ x
  )
}

normalize_cluster <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    x %in% c("1", "Cluster1", "cluster1", "C1", "亚型1", "亚型-1", "Subtype1", "subtype1") ~ "Cluster1",
    x %in% c("2", "Cluster2", "cluster2", "C2", "亚型2", "亚型-2", "Subtype2", "subtype2") ~ "Cluster2",
    TRUE ~ x
  )
}

normalize_immune_rank <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    x %in% c("1", "Cluster1", "cluster1", "C1", "免疫1", "免疫组1", "Immune1", "immune1",
             "ImmuneRankCluster1", "ImmuneRankCluster 1") ~ "ImmuneRankCluster 1",
    x %in% c("2", "Cluster2", "cluster2", "C2", "免疫2", "免疫组2", "Immune2", "immune2",
             "ImmuneRankCluster2", "ImmuneRankCluster 2") ~ "ImmuneRankCluster 2",
    TRUE ~ x
  )
}

normalize_immune_matrix <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    x %in% c("1", "Cluster1", "cluster1", "C1", "免疫1", "免疫组1", "Immune1", "immune1",
             "ImmuneMatrixCluster1", "ImmuneMatrixCluster 1") ~ "ImmuneMatrixCluster 1",
    x %in% c("2", "Cluster2", "cluster2", "C2", "免疫2", "免疫组2", "Immune2", "immune2",
             "ImmuneMatrixCluster2", "ImmuneMatrixCluster 2") ~ "ImmuneMatrixCluster 2",
    TRUE ~ x
  )
}

# -----------------------------
# 8. 颜色函数
# -----------------------------
get_fill_colors <- function(strata_vec) {
  strata_vec <- unique(as.character(strata_vec))

  base_cols <- c(
    "Cluster1"               = "#CACFDD",
    "Cluster2"               = "#F6D3C8",
    "Low risk"               = "#9BD4CB",
    "High risk"              = "#DD7694",
    "ImmuneRankCluster 1"    = "#DD7694",
    "ImmuneRankCluster 2"    = "#9BD4CB",
    "ImmuneMatrixCluster 1"  = "#DD7694",
    "ImmuneMatrixCluster 2"  = "#9BD4CB"
  )

  missing_names <- setdiff(strata_vec, names(base_cols))
  if (length(missing_names) > 0) {
    extra_cols <- grDevices::hcl.colors(length(missing_names), palette = "Set 2")
    names(extra_cols) <- missing_names
    base_cols <- c(base_cols, extra_cols)
  }

  base_cols[strata_vec]
}

get_left_box_colors <- function(strata_vec) {
  strata_vec <- unique(as.character(strata_vec))

  base_cols <- c(
    "Cluster1"               = "#8F98BB",
    "Cluster2"               = "#DFA08A",
    "Low risk"               = "#6FAFA6",
    "High risk"              = "#BD5F7D",
    "ImmuneRankCluster 1"    = "#BD5F7D",
    "ImmuneRankCluster 2"    = "#6FAFA6",
    "ImmuneMatrixCluster 1"  = "#BD5F7D",
    "ImmuneMatrixCluster 2"  = "#6FAFA6"
  )

  missing_names <- setdiff(strata_vec, names(base_cols))
  if (length(missing_names) > 0) {
    extra_cols <- grDevices::hcl.colors(length(missing_names), palette = "Dark 3")
    names(extra_cols) <- missing_names
    base_cols <- c(base_cols, extra_cols)
  }

  base_cols[strata_vec]
}

# -----------------------------
# 9. 单组桑基图函数
# -----------------------------
plot_beautified_sankey <- function(left_df, right_df,
                                   left_col, right_col,
                                   left_axis_label, right_axis_label,
                                   prefix, out_dir,
                                   left_levels = NULL, right_levels = NULL) {

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  merged_data <- dplyr::inner_join(
    left_df[, c("Sample", left_col)],
    right_df[, c("Sample", right_col)],
    by = "Sample"
  )

  if (nrow(merged_data) == 0) {
    stop("合并后没有匹配到样本：", prefix, "，请检查 Sample 是否一致。")
  }

  cat("【", prefix, "】成功匹配样本数：", nrow(merged_data), "\n", sep = "")

  if (!is.null(left_levels)) {
    merged_data[[left_col]] <- factor(as.character(merged_data[[left_col]]), levels = left_levels)
  }
  if (!is.null(right_levels)) {
    merged_data[[right_col]] <- factor(as.character(merged_data[[right_col]]), levels = right_levels)
  }

  # 样本级对应表
  sample_map_out <- merged_data %>%
    dplyr::arrange(.data[[left_col]], .data[[right_col]], Sample) %>%
    dplyr::select(Sample, dplyr::all_of(left_col), dplyr::all_of(right_col))

  write.csv(
    sample_map_out,
    file.path(out_dir, paste0(prefix, "_Sample_Level_Mapping.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  # 汇总表
  summary_table <- merged_data %>%
    dplyr::count(.data[[left_col]], .data[[right_col]], name = "Count") %>%
    dplyr::group_by(.data[[left_col]]) %>%
    dplyr::mutate(
      LeftGroup_Total = sum(Count),
      Ratio_in_LeftGroup = round(Count / LeftGroup_Total, 4)
    ) %>%
    dplyr::ungroup()

  colnames(summary_table)[1:2] <- c(left_col, right_col)

  write.csv(
    summary_table,
    file.path(out_dir, paste0(prefix, "_Summary_Mapping.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  plot_data <- summary_table
  colnames(plot_data)[1:2] <- c("LeftGroup", "RightGroup")
  plot_data$LeftGroup <- factor(as.character(plot_data$LeftGroup), levels = left_levels)
  plot_data$RightGroup <- factor(as.character(plot_data$RightGroup), levels = right_levels)

  fill_colors <- get_fill_colors(c(as.character(plot_data$LeftGroup),
                                   as.character(plot_data$RightGroup)))

  left_box_data <- plot_data %>%
    dplyr::group_by(LeftGroup) %>%
    dplyr::summarise(Count = sum(Count), .groups = "drop") %>%
    dplyr::mutate(LeftGroup = factor(LeftGroup, levels = levels(plot_data$LeftGroup))) %>%
    dplyr::arrange(desc(as.integer(LeftGroup))) %>%
    dplyr::mutate(
      ymin = dplyr::lag(cumsum(Count), default = 0),
      ymax = cumsum(Count),
      xmin = 0.90,
      xmax = 1.10,
      box_fill = unname(get_left_box_colors(as.character(LeftGroup)))
    )

  sankey_plot <- ggplot(
    plot_data,
    aes(axis1 = LeftGroup, axis2 = RightGroup, y = Count)
  ) +
    geom_alluvium(
      aes(fill = LeftGroup),
      width = 0.20,
      alpha = 0.88,
      knot.pos = 0.42,
      color = NA,
      linewidth = 0
    ) +
    geom_stratum(
      aes(fill = after_stat(stratum)),
      width = 0.20,
      color = NA,
      linewidth = 0
    ) +
    geom_rect(
      data = left_box_data,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = left_box_data$box_fill,
      color = NA,
      linewidth = 0
    ) +
    geom_text(
      stat = "stratum",
      aes(label = paste0(after_stat(stratum), "\n(n=", after_stat(count), ")")),
      size = 4.3,
      fontface = "bold",
      lineheight = 0.92,
      color = "grey10"
    ) +
    scale_x_discrete(
      limits = c(left_axis_label, right_axis_label),
      expand = c(0.10, 0.10)
    ) +
    scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.04))
    ) +
    scale_fill_manual(values = fill_colors, guide = "none") +
    labs(
      title = paste0(left_axis_label, " vs ", right_axis_label),
      subtitle = paste0("n = ", nrow(merged_data)),
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(size = 15.5, face = "bold", hjust = 0.5, color = "grey10"),
      plot.subtitle = element_text(size = 10.5, hjust = 0.5, color = "grey35"),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 11.5, face = "bold", color = "grey10"),
      axis.text.y = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 16, 10, 16)
    )

  pdf_file <- file.path(out_dir, paste0(prefix, "_Sankey_Beautified.pdf"))
  svg_file <- file.path(out_dir, paste0(prefix, "_Sankey_Beautified.svg"))

  save_pdf_file <- function(plot, final_file, width = 9.5, height = 7) {
    tryCatch(
      ggsave(
        filename = final_file,
        plot = plot,
        device = cairo_pdf,
        width = width,
        height = height,
        units = "in",
        dpi = 300
      ),
      error = function(e) {
        message("cairo_pdf failed, falling back to grDevices::pdf: ", conditionMessage(e))
        ggsave(
          filename = final_file,
          plot = plot,
          device = grDevices::pdf,
          width = width,
          height = height,
          units = "in",
          dpi = 300,
          useDingbats = FALSE
        )
      }
    )
  }

  save_svg_file <- function(plot, final_file, width = 9.5, height = 7) {
    ggsave(
      filename = final_file,
      plot = plot,
      device = svglite::svglite,
      width = width,
      height = height,
      units = "in"
    )
  }

  save_pdf_file(sankey_plot, pdf_file)
  save_svg_file(sankey_plot, svg_file)

  no_title_dir <- file.path(out_dir, "无标题")
  if (!dir.exists(no_title_dir)) dir.create(no_title_dir, recursive = TRUE)

  sankey_plot_no_title <- sankey_plot +
    labs(title = NULL, subtitle = NULL) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )

  pdf_file_no_title <- file.path(no_title_dir, paste0(prefix, "_Sankey_NoTitle.pdf"))
  svg_file_no_title <- file.path(no_title_dir, paste0(prefix, "_Sankey_NoTitle.svg"))

  save_pdf_file(sankey_plot_no_title, pdf_file_no_title)
  save_svg_file(sankey_plot_no_title, svg_file_no_title)

  cat("✅ PDF已保存：", pdf_file, "\n")
  cat("✅ SVG已保存：", svg_file, "\n")
  cat("✅ 样本级对应表已保存：", file.path(out_dir, paste0(prefix, "_Sample_Level_Mapping.csv")), "\n")
  cat("✅ 汇总对应表已保存：", file.path(out_dir, paste0(prefix, "_Summary_Mapping.csv")), "\n")
  cat("✅ 无标题 PDF 已保存：", pdf_file_no_title, "\n")
  cat("✅ 无标题 SVG 已保存：", svg_file_no_title, "\n\n")
}

# -----------------------------
# 10. 生成四组图
# -----------------------------
for (pcd in names(pcd_files)) {  # 按 PCD 循环
  cat("\n========== ", pcd, " ==========\n")

  output_dir <- file.path(output_root, pcd)  # 按 PCD 分子目录
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  risk_file    <- file.path(input_dir, pcd_files[[pcd]][["risk"]])
  subtype_file <- file.path(input_dir, pcd_files[[pcd]][["subtype"]])

  cat("输入文件1：", risk_file, "\n")
  cat("输入文件2：", subtype_file, "\n")
  cat("输入文件3：", immune_rank_file, "\n")
  cat("输入文件4：", immune_matrix_file, "\n")
  cat("输出目录：", output_dir, "\n")

  for (f in c(risk_file, subtype_file, immune_rank_file, immune_matrix_file)) {
    if (!file.exists(f)) stop("找不到文件：", f)
  }

  # -----------------------------
  # 6. 读取四类数据
  # -----------------------------
  risk_data <- read_group_file(
    file = risk_file,
    group_name = "RiskGroup",
    group_candidates = c("RiskGroup", "riskgroup", "Risk_Group", "risk_group", "Group", "group", "Risk", "风险组")
  )

  subtype_data <- read_group_file(
    file = subtype_file,
    group_name = "Subtype",
    group_candidates = c("Cluster", "cluster", "Subtype", "subtype", "Type", "type", "亚型")
  )

  immune_rank_data <- read_group_file(
    file = immune_rank_file,
    group_name = "ImmuneRankCluster",
    group_candidates = c("Cluster", "cluster", "Group", "group", "Subtype", "subtype", "Type", "type", "亚型")
  )

  immune_matrix_data <- read_group_file(
    file = immune_matrix_file,
    group_name = "ImmuneMatrixCluster",
    group_candidates = c("Cluster", "cluster", "Group", "group", "Subtype", "subtype", "Type", "type", "亚型")
  )

  # -----------------------------
  # 7. 标签标准化
  # -----------------------------
  risk_data$RiskGroup <- normalize_risk(risk_data$RiskGroup)
  subtype_data$Subtype <- normalize_cluster(subtype_data$Subtype)
  immune_rank_data$ImmuneRankCluster <- normalize_immune_rank(immune_rank_data$ImmuneRankCluster)
  immune_matrix_data$ImmuneMatrixCluster <- normalize_immune_matrix(immune_matrix_data$ImmuneMatrixCluster)

  # A1 Risk group vs Immune rank subtype
  # plot_beautified_sankey(
  #   left_df = risk_data,
  #   right_df = immune_rank_data,
  #   left_col = "RiskGroup",
  #   right_col = "ImmuneRankCluster",
  #   left_axis_label = "Risk group",
  #   right_axis_label = "Immune rank subtype",
  #   prefix = "A1_RiskGroup_vs_ImmuneRankSubtype",
  #   out_dir = file.path(output_dir, "A1_RiskGroup_vs_ImmuneRankSubtype"),
  #   left_levels = c("Low risk", "High risk"),
  #   right_levels = c("ImmuneRankCluster 1", "ImmuneRankCluster 2")
  # )

  # A2 Risk group vs Immune matrix subtype
  # plot_beautified_sankey(
  #   left_df = risk_data,
  #   right_df = immune_matrix_data,
  #   left_col = "RiskGroup",
  #   right_col = "ImmuneMatrixCluster",
  #   left_axis_label = "Risk group",
  #   right_axis_label = "Immune matrix subtype",
  #   prefix = "A2_RiskGroup_vs_ImmuneMatrixSubtype",
  #   out_dir = file.path(output_dir, "A2_RiskGroup_vs_ImmuneMatrixSubtype"),
  #   left_levels = c("Low risk", "High risk"),
  #   right_levels = c("ImmuneMatrixCluster 1", "ImmuneMatrixCluster 2")
  # )

  # B1 Molecular subtype vs Immune rank subtype
  plot_beautified_sankey(
    left_df = subtype_data,
    right_df = immune_rank_data,
    left_col = "Subtype",
    right_col = "ImmuneRankCluster",
    left_axis_label = "Molecular subtype",
    right_axis_label = "Immune rank subtype",
    prefix = "B1_MolecularSubtype_vs_ImmuneRankSubtype",
    out_dir = file.path(output_dir, "B1_MolecularSubtype_vs_ImmuneRankSubtype"),
    left_levels = c("Cluster1", "Cluster2"),
    right_levels = c("ImmuneRankCluster 1", "ImmuneRankCluster 2")
  )

  # B2 Molecular subtype vs Immune matrix subtype
  # plot_beautified_sankey(
  #   left_df = subtype_data,
  #   right_df = immune_matrix_data,
  #   left_col = "Subtype",
  #   right_col = "ImmuneMatrixCluster",
  #   left_axis_label = "Molecular subtype",
  #   right_axis_label = "Immune matrix subtype",
  #   prefix = "B2_MolecularSubtype_vs_ImmuneMatrixSubtype",
  #   out_dir = file.path(output_dir, "B2_MolecularSubtype_vs_ImmuneMatrixSubtype"),
  #   left_levels = c("Cluster1", "Cluster2"),
  #   right_levels = c("ImmuneMatrixCluster 1", "ImmuneMatrixCluster 2")
  # )
}

# -----------------------------
# 11. 全部完成
# -----------------------------
cat("========================================\n")
cat("B1 sankey plot and mapping tables exported.\n")
if (FALSE)
cat("✅ 全部4组桑基图及对应表已输出完成！\n")
cat("输出总目录：", output_root, "\n")
cat("========================================\n")
