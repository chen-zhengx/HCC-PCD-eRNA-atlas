# =========================================================
# 对应结果：Figure 5B（三种分子亚型分组的两两一致率热图，90.00%–95.71%）
#           同时输出三种亚型分组的桑基图、组合模式柱状图及统计表（未用于论文图）
# 输入：input/subtype_anoikis.csv、input/subtype_apoptosis.csv、input/subtype_ferroptosis.csv
#       （论文使用的亚型分组，列 Sample, Cluster；保留的聚类代码无法重新生成这些分组，故作为数据提供）
# 输出：output/figure5B_subtype_concordance/ 下：
#       三种死亡亚型_两两一致率热图_浅蓝版.pdf/.png（Figure 5B）、三种死亡亚型_两两一致率矩阵.csv
#       三种死亡亚型_样本对应表.csv、三种死亡亚型_组合模式统计.csv、三种死亡亚型_桑基图.*
#       三种死亡亚型_组合模式柱状图.*、三种死亡亚型_结果摘要.txt
# 运行：在模块目录下 Rscript code/04_subtype_concordance.R
# =========================================================
# =========================================================
# 三种死亡亚型横向对比可视化
# 输出：
# 1) 桑基图 / Alluvial plot
# 2) 组合模式柱状图
# 3) 两两一致率热图（浅蓝色版）
# 4) 合并后的样本亚型对应表
# 5) 组合统计表、一致率矩阵
# 6) 结果摘要 txt
# =========================================================

rm(list = ls())

# -----------------------------
# 1. 加载包
# -----------------------------
pkgs <- c("dplyr", "ggplot2", "ggalluvial", "reshape2")
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
output_dir <- "output/figure5B_subtype_concordance"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Journal-style colors and shared legend settings
journal_colors <- c("Cluster1" = "#00A087", "Cluster2" = "#3C5488")
heatmap_colors <- c("#c1e1ed", "#b0d4e3", "#4DBBD5", "#3C5488")
legend_limits <- c(0.9, 1)
legend_breaks <- seq(0.9, 1, by = 0.025)
percent_labels <- function(x) {
  pct <- round(x * 100, 1)
  ifelse(pct == round(pct), sprintf("%.0f%%", pct), sprintf("%.1f%%", pct))
}

theme_journal <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        hjust = 0.5, face = "bold", color = "#222222",
        size = base_size + 2
      ),
      axis.title = element_text(face = "bold", color = "#222222"),
      axis.text = element_text(color = "#222222"),
      axis.line = element_line(color = "#222222", linewidth = 0.5),
      axis.ticks = element_line(color = "#222222", linewidth = 0.4),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(color = "#222222"),
      panel.grid = element_blank(),
      plot.margin = margin(10, 14, 10, 10)
    )
}

file_anoikis     <- file.path(input_dir, "subtype_anoikis.csv")
file_apoptosis   <- file.path(input_dir, "subtype_apoptosis.csv")
file_ferroptosis <- file.path(input_dir, "subtype_ferroptosis.csv")

# -----------------------------
# 3. 读取数据函数
# -----------------------------
read_cluster_file <- function(file, group_name) {
  df <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  
  # 默认读取前两列：Sample + Cluster
  df <- df[, 1:2]
  colnames(df) <- c("Sample", group_name)
  
  # 去空格
  df$Sample <- trimws(df$Sample)
  df[[group_name]] <- trimws(as.character(df[[group_name]]))
  
  # 统一写法：1/2 -> Cluster1/Cluster2
  df[[group_name]] <- ifelse(df[[group_name]] %in% c("1", "Cluster1", "cluster1", "C1", "c1"),
                             "Cluster1",
                             ifelse(df[[group_name]] %in% c("2", "Cluster2", "cluster2", "C2", "c2"),
                                    "Cluster2",
                                    df[[group_name]]))
  
  return(df)
}

# -----------------------------
# 4. 读取三组亚型数据
# -----------------------------
anoikis_df     <- read_cluster_file(file_anoikis, "Anoikis")
apoptosis_df   <- read_cluster_file(file_apoptosis, "Apoptosis")
ferroptosis_df <- read_cluster_file(file_ferroptosis, "Ferroptosis")

# -----------------------------
# 5. 合并数据
# -----------------------------
merged_df <- anoikis_df %>%
  inner_join(apoptosis_df, by = "Sample") %>%
  inner_join(ferroptosis_df, by = "Sample") %>%
  distinct()

# 因子顺序：Cluster1 在上，Cluster2 在下
merged_df$Anoikis     <- factor(merged_df$Anoikis, levels = c("Cluster1", "Cluster2"))
merged_df$Apoptosis   <- factor(merged_df$Apoptosis, levels = c("Cluster1", "Cluster2"))
merged_df$Ferroptosis <- factor(merged_df$Ferroptosis, levels = c("Cluster1", "Cluster2"))

# 保存样本对应表
write.csv(
  merged_df,
  file.path(output_dir, "三种死亡亚型_样本对应表.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -----------------------------
# 6. 组合模式统计
# -----------------------------
combo_df <- merged_df %>%
  mutate(Combination = paste(Anoikis, Apoptosis, Ferroptosis, sep = "-")) %>%
  count(Combination, sort = TRUE)

write.csv(
  combo_df,
  file.path(output_dir, "三种死亡亚型_组合模式统计.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -----------------------------
# 7. 两两一致率矩阵
# -----------------------------
group_names <- c("Anoikis", "Apoptosis", "Ferroptosis")

agree_mat <- matrix(NA, nrow = 3, ncol = 3)
rownames(agree_mat) <- group_names
colnames(agree_mat) <- group_names

for (i in 1:3) {
  for (j in 1:3) {
    agree_mat[i, j] <- mean(merged_df[[group_names[i]]] == merged_df[[group_names[j]]])
  }
}

write.csv(
  round(agree_mat, 4),
  file.path(output_dir, "三种死亡亚型_两两一致率矩阵.csv"),
  row.names = TRUE,
  fileEncoding = "UTF-8"
)

agree_long <- reshape2::melt(agree_mat)
colnames(agree_long) <- c("Group1", "Group2", "Agreement")

# -----------------------------
# 8. 桑基图 / Alluvial plot
# -----------------------------
p_alluvial <- ggplot(
  merged_df,
  aes(axis1 = Anoikis, axis2 = Apoptosis, axis3 = Ferroptosis, y = 1)
) +
  geom_alluvium(
    aes(fill = Anoikis),
    width = 0.18, alpha = 0.82, color = "white", linewidth = 0.15
  ) +
  geom_stratum(width = 0.18, fill = "#F7F7F7", color = "#3A3A3A", linewidth = 0.45) +
  geom_text(
    stat = "stratum", aes(label = after_stat(stratum)),
    size = 4.2, fontface = "bold", color = "#222222"
  ) +
  scale_x_discrete(
    limits = c("Anoikis", "Apoptosis", "Ferroptosis"),
    expand = c(0.08, 0.08)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  scale_fill_manual(values = journal_colors) +
  labs(
    title = "Alluvial Plot of Molecular Subtypes Across Three Programmed Cell Death Types",
    x = NULL,
    y = "Sample Count",
    fill = NULL
  ) +
  guides(fill = guide_legend(title = NULL, override.aes = list(alpha = 0.9))) +
  theme_journal(base_size = 13) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(margin = margin(r = 8))
  )

ggsave(
  file.path(output_dir, "三种死亡亚型_桑基图.pdf"),
  p_alluvial, width = 10, height = 6
)
ggsave(
  file.path(output_dir, "三种死亡亚型_桑基图.png"),
  p_alluvial, width = 10, height = 6, dpi = 300
)

# -----------------------------
# 9. 组合模式柱状图
# -----------------------------
combo_df$Combination <- factor(
  combo_df$Combination,
  levels = combo_df$Combination[order(combo_df$n, decreasing = TRUE)]
)

p_bar <- ggplot(combo_df, aes(x = Combination, y = n)) +
  geom_col(width = 0.68, fill = "#3C5488", color = "#222222", linewidth = 0.25) +
  geom_text(
    aes(label = n),
    vjust = -0.35, size = 4, fontface = "bold", color = "#222222"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Combination Pattern Counts of Three Molecular Subtype Classifications",
    x = "Subtype Combination",
    y = "Sample Count"
  ) +
  theme_journal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, color = "#222222"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8))
  )

ggsave(
  file.path(output_dir, "三种死亡亚型_组合模式柱状图.pdf"),
  p_bar, width = 10, height = 6
)
ggsave(
  file.path(output_dir, "三种死亡亚型_组合模式柱状图.png"),
  p_bar, width = 10, height = 6, dpi = 300
)

# -----------------------------
# 10. 两两一致率热图（浅蓝色版）
# -----------------------------
p_heat <- ggplot(agree_long, aes(x = Group1, y = Group2, fill = Agreement)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(
    aes(label = sprintf("%.1f%%", Agreement * 100), color = Agreement >= 0.65),
    size = 4.7, fontface = "bold", show.legend = FALSE
  ) +
  scale_fill_gradientn(
    colors = heatmap_colors,
    limits = legend_limits,
    breaks = legend_breaks,
    labels = percent_labels,
    oob = scales::squish
  ) +
  scale_color_manual(
    values = c("FALSE" = "#222222", "TRUE" = "white")
  ) +
  coord_fixed() +
  labs(
    title = "Pairwise Agreement of Molecular Subtypes",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  guides(
    fill = guide_colorbar(
      title = NULL,
      barheight = grid::unit(45, "mm"),
      barwidth = grid::unit(4, "mm")
    )
  ) +
  theme_journal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "#222222", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    axis.text = element_text(face = "bold", color = "#222222"),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  )

ggsave(
  file.path(output_dir, "三种死亡亚型_两两一致率热图_浅蓝版.pdf"),
  p_heat, width = 7, height = 6
)
ggsave(
  file.path(output_dir, "三种死亡亚型_两两一致率热图_浅蓝版.png"),
  p_heat, width = 7, height = 6, dpi = 300
)

# -----------------------------
# 11. 输出文字摘要
# -----------------------------
combo_text <- paste0(combo_df$Combination, ": ", combo_df$n, collapse = "\n")

sink(file.path(output_dir, "三种死亡亚型_结果摘要.txt"))

cat("三种死亡亚型横向对比结果摘要\n")
cat("====================================\n\n")
cat("总样本数：", nrow(merged_df), "\n\n")

cat("组合模式统计：\n")
cat(combo_text, "\n\n")

cat("两两一致率矩阵：\n")
print(round(agree_mat, 4))
cat("\n")

cat("说明：\n")
cat("1. 桑基图用于展示同一样本在三种死亡亚型之间的流向关系。\n")
cat("2. 组合模式柱状图用于展示不同 Cluster 组合的样本数。\n")
cat("3. 热图用于展示三种亚型分组之间的两两一致率。\n")

sink()

# -----------------------------
# 12. 控制台提示
# -----------------------------
cat("分析完成，结果已输出到：\n", output_dir, "\n")

