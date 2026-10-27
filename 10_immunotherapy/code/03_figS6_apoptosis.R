# =========================================================
# 对应结果：Figure S6（凋亡（apoptosis）相关分子亚型 Molecular_Cluster1 vs Cluster2 的 tigeR 免疫治疗反应标志物比较：
#           双侧 Wilcoxon 秩和检验，展示 P < 0.05 的标志物，图中为原始 P 值）
# 输入：input/tigeR_23signatures_scores.csv
#       input/subtype_apoptosis.csv
# 输出：output/figS6_apoptosis/tigeR_23signatures_diff_clusters_statistics_byP.csv（23 种标志物的检验结果）
#       output/figS6_apoptosis/SignificantSignatures_ViolinBox_Facet_Plt0.05.png/.pdf（P < 0.05 的标志物分面图，即论文图
#       论文图中删去了总标题，横轴标签 Cluster1/2 改为 Molecular Cluster1/2，均在作图软件中完成）
#       output/figS6_apoptosis/<标志物>_violin_boxplot_withP.png/.pdf（单个标志物图，论文未使用）
# 运行：在模块目录下 Rscript code/03_figS6_apoptosis.R
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
})

# =========================
# 0) 路径
# =========================
in_dir  <- "input"
out_dir <- "output/figS6_apoptosis"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

f_scores <- file.path(in_dir, "tigeR_23signatures_scores.csv")
f_group  <- file.path(in_dir, "subtype_apoptosis.csv")

# =========================
# 1) 读取数据
# =========================
scores <- fread(f_scores)
group_info <- fread(f_group)

# =========================
# 2) 检查并统一列名
# =========================
if (!("Sample" %in% colnames(scores))) {
  setnames(scores, 1, "Sample")
}
if (!("Sample" %in% colnames(group_info))) {
  setnames(group_info, 1, "Sample")
}
if (!("Cluster" %in% colnames(group_info))) {
  setnames(group_info, 2, "Cluster")
}

scores$Sample <- trimws(as.character(scores$Sample))
group_info$Sample <- trimws(as.character(group_info$Sample))
group_info$Cluster <- trimws(as.character(group_info$Cluster))

# 统一亚型标签
group_info$Cluster <- case_when(
  group_info$Cluster %in% c("1", "Cluster1", "cluster1", "C1", "亚型1", "亚型-1") ~ "Cluster1",
  group_info$Cluster %in% c("2", "Cluster2", "cluster2", "C2", "亚型2", "亚型-2") ~ "Cluster2",
  TRUE ~ group_info$Cluster
)
group_info$Cluster <- factor(group_info$Cluster, levels = c("Cluster1", "Cluster2"))

# =========================
# 3) 合并数据
# =========================
merged <- merge(scores, group_info, by = "Sample")
merged$Cluster <- factor(as.character(merged$Cluster), levels = c("Cluster1", "Cluster2"))

sig_cols <- setdiff(colnames(merged), c("Sample", "Cluster"))

for (cc in sig_cols) {
  merged[[cc]] <- as.numeric(merged[[cc]])
}

# =========================
# 4) 定义统计函数
# =========================
get_stats <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]
  g <- g[ok]
  
  if (length(unique(g)) < 2) {
    return(data.frame(
      p_value = NA_real_,
      Cluster1_median = NA_real_,
      Cluster2_median = NA_real_,
      Difference = NA_real_,
      HigherGroup = NA_character_
    ))
  }
  
  c1_vals <- x[g == "Cluster1"]
  c2_vals <- x[g == "Cluster2"]
  
  p <- tryCatch({
    wilcox.test(x ~ g, exact = FALSE)$p.value
  }, error = function(e) NA_real_)
  
  c1_med <- if (length(c1_vals) > 0) median(c1_vals, na.rm = TRUE) else NA_real_
  c2_med <- if (length(c2_vals) > 0) median(c2_vals, na.rm = TRUE) else NA_real_
  diff_med <- c2_med - c1_med
  
  higher <- ifelse(
    is.na(diff_med), NA_character_,
    ifelse(diff_med > 0, "Cluster2",
           ifelse(diff_med < 0, "Cluster1", "Equal"))
  )
  
  data.frame(
    p_value = p,
    Cluster1_median = c1_med,
    Cluster2_median = c2_med,
    Difference = diff_med,
    HigherGroup = higher
  )
}

# =========================
# 5) 批量统计
# =========================
stats_list <- lapply(sig_cols, function(s) {
  tmp <- get_stats(merged[[s]], merged$Cluster)
  tmp$Signature <- s
  tmp
})

stats_df <- bind_rows(stats_list) %>%
  dplyr::select(Signature, everything()) %>%
  arrange(p_value)

stats_df <- stats_df %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  mutate(Significance = case_when(
    is.na(p_value) ~ NA_character_,
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

write.csv(
  stats_df,
  file.path(out_dir, "tigeR_23signatures_diff_clusters_statistics_byP.csv"),
  row.names = FALSE
)

journal_fill_cols <- c("Cluster1" = "#E64B35", "Cluster2" = "#4DBBD5")
journal_point_cols <- c("Cluster1" = "#8F2D24", "Cluster2" = "#267A96")

theme_journal <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "#222222"),
      plot.title = element_text(hjust = 0.5, color = "#1F1F1F"),
      axis.title = element_text(color = "#222222"),
      axis.text = element_text(color = "#222222"),
      axis.line = element_line(color = "#222222", linewidth = 0.45),
      axis.ticks = element_line(color = "#222222", linewidth = 0.45),
      strip.background = element_rect(fill = "#F3F4F4", color = NA),
      strip.text = element_text(color = "#222222"),
      panel.grid = element_blank(),
      panel.spacing = grid::unit(0.9, "lines"),
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

save_plot <- function(filename, plot, width, height, dpi = 300) {
  ext <- tolower(tools::file_ext(filename))
  
  if (ext == "pdf") {
    tmp_pdf <- tempfile(fileext = ".pdf")
    on.exit(if (file.exists(tmp_pdf)) unlink(tmp_pdf), add = TRUE)
    ggsave(
      tmp_pdf, plot,
      width = width, height = height,
      device = grDevices::pdf,
      useDingbats = FALSE,
      bg = "white"
    )
    file.copy(tmp_pdf, filename, overwrite = TRUE)
  } else {
    ggsave(
      filename, plot,
      width = width, height = height,
      dpi = dpi,
      bg = "white"
    )
  }
}

# =========================
# 6) 按 p 值筛选显著 signature
# =========================
sig_keep <- stats_df %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  pull(Signature)

if (length(sig_keep) == 0) {
  message("没有任何 signature 达到 p < 0.05。统计表已输出。")
} else {
  
  long_df <- merged %>%
    as.data.frame() %>%
    dplyr::select(Sample, Cluster, all_of(sig_keep)) %>%
    pivot_longer(
      cols = all_of(sig_keep),
      names_to = "Signature",
      values_to = "Score"
    ) %>%
    filter(is.finite(Score), !is.na(Score)) %>%
    mutate(Signature = factor(Signature, levels = sig_keep))
  
  y_pos <- long_df %>%
    group_by(Signature) %>%
    summarise(
      ymax = max(Score, na.rm = TRUE),
      ymin = min(Score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(stats_df %>% dplyr::select(Signature, p_value), by = "Signature") %>%
    mutate(
      label = ifelse(p_value < 0.001, "P < 0.001", paste0("P = ", signif(p_value, 3))),
      y = ymax + pmax(ymax - ymin, .Machine$double.eps) * 0.18
    )
  
  # =========================
  # 7) 分面总图：小提琴 + 细箱线图 + 散点
  # =========================
  p_all <- ggplot(long_df, aes(x = Cluster, y = Score)) +
    geom_violin(
      aes(fill = Cluster),
      width = 0.78,
      alpha = 0.34,
      color = NA,
      trim = FALSE
    ) +
    geom_boxplot(
      aes(fill = Cluster),
      width = 0.16,
      outlier.shape = NA,
      alpha = 0.76,
      color = "#2B2B2B",
      linewidth = 0.45
    ) +
    geom_jitter(
      aes(color = Cluster),
      width = 0.07,
      size = 1.05,
      alpha = 0.58
    ) +
    facet_wrap(~Signature, scales = "free_y") +
    geom_text(
      data = y_pos,
      aes(x = 1.5, y = y, label = label),
      inherit.aes = FALSE,
      size = 3.0,
      color = "#222222"
    ) +
    scale_fill_manual(values = journal_fill_cols) +
    scale_color_manual(values = journal_point_cols) +
    labs(
      title = "tigeR signatures in apoptosis clusters",
      x = "Cluster",
      y = "Signature Score"
    ) +
    theme_journal(base_size = 12) +
    theme(
      strip.text = element_text(size = 9),
      plot.title = element_text(size = 15),
      axis.text.x = element_text()
    )
  
  save_plot(
    file.path(out_dir, "SignificantSignatures_ViolinBox_Facet_Plt0.05.png"),
    p_all, width = 14, height = 10, dpi = 300
  )
  save_plot(
    file.path(out_dir, "SignificantSignatures_ViolinBox_Facet_Plt0.05.pdf"),
    p_all, width = 14, height = 10
  )
  
  # =========================
  # 8) 单个 signature 图
  # =========================
  for (s in sig_keep) {
    df1 <- long_df %>% filter(Signature == s)
    
    p1 <- ggplot(df1, aes(x = Cluster, y = Score)) +
      geom_violin(
        aes(fill = Cluster),
        width = 0.74,
        alpha = 0.34,
        color = NA,
        trim = FALSE
      ) +
      geom_boxplot(
        aes(fill = Cluster),
        width = 0.14,
        outlier.shape = NA,
        alpha = 0.78,
        color = "#2B2B2B",
        linewidth = 0.5
      ) +
      geom_jitter(
        aes(color = Cluster),
        width = 0.07,
        size = 1.45,
        alpha = 0.62
      ) +
      stat_compare_means(
        method = "wilcox.test",
        label = "p.format",
        comparisons = list(c("Cluster1", "Cluster2")),
        tip.length = 0.01,
        bracket.size = 0.4,
        size = 4.2,
        color = "#222222"
      ) +
      scale_fill_manual(values = journal_fill_cols) +
      scale_color_manual(values = journal_point_cols) +
      labs(
        title = s,
        x = "Cluster",
        y = "Signature Score"
      ) +
      theme_journal(base_size = 14) +
      theme(
        plot.title = element_text(size = 16),
        axis.text.x = element_text()
      )
    
    safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", s)
    
    save_plot(
      file.path(out_dir, paste0(safe_name, "_violin_boxplot_withP.png")),
      p1, width = 5.8, height = 5.6, dpi = 300
    )
    save_plot(
      file.path(out_dir, paste0(safe_name, "_violin_boxplot_withP.pdf")),
      p1, width = 5.8, height = 5.6
    )
  }
}

message("分析完成！输出文件位于：", out_dir)

