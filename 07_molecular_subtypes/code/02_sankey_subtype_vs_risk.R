# =========================================================
# 对应结果：Figure 3D（失巢凋亡 anoikis）、Figure S1D（凋亡 apoptosis）、Figure S2D（铁死亡 ferroptosis）
#           分子亚型（Cluster1/2）与高/低风险组的对应桑基图；正文比例（如 anoikis 91.18%/88.89%）
#           来自输出的 Summary_Cluster_RiskGroup_Mapping.csv
# 输入：input/subtype_<pcd>.csv（论文使用的亚型分组，列 Sample, Cluster）
#       input/risk_group_<pcd>.csv（中位数风险评分分组，列 Sample, RiskGroup；来自 05_prognostic_model）
# 输出：output/sankey/<pcd>/Cluster_RiskGroup_Sankey_Beautified.pdf/.svg
#       Sample_Level_Cluster_RiskGroup_Mapping.csv、Summary_Cluster_RiskGroup_Mapping.csv
# 运行：在模块目录下 Rscript code/02_sankey_subtype_vs_risk.R
# =========================================================

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

for (pcd in c("anoikis", "apoptosis", "ferroptosis")) {
cat("\n==========", pcd, "==========\n")

# -----------------------------
# 2. 路径
# -----------------------------
input_dir  <- "input"
output_dir <- file.path("output", "sankey", pcd)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cluster_file <- file.path(input_dir, paste0("subtype_", pcd, ".csv"))
risk_file    <- file.path(input_dir, paste0("risk_group_", pcd, ".csv"))

cat("输入文件1：", cluster_file, "\n")
cat("输入文件2：", risk_file, "\n")
cat("输出目录：", output_dir, "\n")

if (!file.exists(cluster_file)) stop("找不到文件：", cluster_file)
if (!file.exists(risk_file)) stop("找不到文件：", risk_file)

# -----------------------------
# 3. 读取数据
# -----------------------------
cluster_data <- read.csv(cluster_file, stringsAsFactors = FALSE, check.names = FALSE)
risk_data    <- read.csv(risk_file, stringsAsFactors = FALSE, check.names = FALSE)

# -----------------------------
# 4. 检查列名
# -----------------------------
if (!"Sample" %in% colnames(cluster_data)) {
  cand <- c("sample", "SAMPLE", "ID", "id")
  hit <- cand[cand %in% colnames(cluster_data)]
  if (length(hit) == 0) stop("cluster_data 中缺少 Sample 列")
  cluster_data$Sample <- cluster_data[[hit[1]]]
}
if (!"Sample" %in% colnames(risk_data)) {
  cand <- c("sample", "SAMPLE", "ID", "id")
  hit <- cand[cand %in% colnames(risk_data)]
  if (length(hit) == 0) stop("risk_data 中缺少 Sample 列")
  risk_data$Sample <- risk_data[[hit[1]]]
}

if (!"Cluster" %in% colnames(cluster_data)) {
  cand <- c("cluster", "Subtype", "subtype", "Type", "type", "亚型")
  hit <- cand[cand %in% colnames(cluster_data)]
  if (length(hit) == 0) stop("cluster_data 中缺少 Cluster 列")
  cluster_data$Cluster <- cluster_data[[hit[1]]]
}

if (!"RiskGroup" %in% colnames(risk_data)) {
  cand <- c("riskgroup", "Risk_Group", "risk_group", "Group", "group", "Risk", "风险组")
  hit <- cand[cand %in% colnames(risk_data)]
  if (length(hit) == 0) stop("risk_data 中缺少 RiskGroup 列")
  risk_data$RiskGroup <- risk_data[[hit[1]]]
}

cluster_data$Sample   <- trimws(as.character(cluster_data$Sample))
cluster_data$Cluster  <- trimws(as.character(cluster_data$Cluster))
risk_data$Sample      <- trimws(as.character(risk_data$Sample))
risk_data$RiskGroup   <- trimws(as.character(risk_data$RiskGroup))

# -----------------------------
# 5. 合并
# -----------------------------
merged_data <- dplyr::inner_join(
  cluster_data[, c("Sample", "Cluster")],
  risk_data[, c("Sample", "RiskGroup")],
  by = "Sample"
)

if (nrow(merged_data) == 0) {
  stop("合并后没有匹配到样本，请检查 Sample 是否一致。")
}

cat("成功匹配样本数：", nrow(merged_data), "\n")

# -----------------------------
# 6. 标签标准化
# Cluster1 / Cluster2
# Low risk / High risk
# -----------------------------
merged_data <- merged_data %>%
  mutate(
    Cluster = case_when(
      Cluster %in% c("1", "Cluster1", "cluster1", "C1", "亚型1", "亚型-1") ~ "Cluster1",
      Cluster %in% c("2", "Cluster2", "cluster2", "C2", "亚型2", "亚型-2") ~ "Cluster2",
      TRUE ~ as.character(Cluster)
    ),
    RiskGroup = case_when(
      RiskGroup %in% c("High", "HIGH", "high", "高风险", "高风险组", "Risk-High") ~ "High risk",
      RiskGroup %in% c("Low", "LOW", "low", "低风险", "低风险组", "Risk-Low") ~ "Low risk",
      TRUE ~ as.character(RiskGroup)
    )
  )

merged_data$Cluster <- factor(merged_data$Cluster, levels = c("Cluster1", "Cluster2"))
# Keep the right-side Sankey strata ordered from top to bottom as:
# High risk, then Low risk.
merged_data$RiskGroup <- factor(merged_data$RiskGroup, levels = c("High risk", "Low risk"))

# -----------------------------
# 7. 导出样本级对应表
# -----------------------------
sample_map <- merged_data %>%
  arrange(Cluster, RiskGroup, Sample)

write.csv(
  sample_map,
  file.path(output_dir, "Sample_Level_Cluster_RiskGroup_Mapping.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -----------------------------
# 8. 导出汇总表
# -----------------------------
summary_table <- merged_data %>%
  count(Cluster, RiskGroup, name = "Count") %>%
  group_by(Cluster) %>%
  mutate(
    Cluster_Total = sum(Count),
    Ratio_in_Cluster = round(Count / Cluster_Total, 4)
  ) %>%
  ungroup()

write.csv(
  summary_table,
  file.path(output_dir, "Summary_Cluster_RiskGroup_Mapping.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -----------------------------
# 9. 配色
# -----------------------------
fill_colors <- c(
  "Cluster1"  = "#CACFDD",
  "Cluster2"  = "#F6D3C8",
  "High risk" = "#dd7694",
  "Low risk"  = "#9BD4CB"
)

cluster_box_colors <- c(
  "Cluster1" = "#8F98BB",
  "Cluster2" = "#DFA08A"
)

# -----------------------------
# 10. 作图
# -----------------------------
plot_data <- summary_table

cluster_box_data <- summary_table %>%
  group_by(Cluster) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  mutate(
    Cluster = factor(Cluster, levels = levels(merged_data$Cluster))
  ) %>%
  arrange(desc(as.integer(Cluster))) %>%
  mutate(
    ymin = dplyr::lag(cumsum(Count), default = 0),
    ymax = cumsum(Count),
    xmin = 0.90,
    xmax = 1.10
  )

sankey_plot <- ggplot(
  plot_data,
  aes(axis1 = Cluster, axis2 = RiskGroup, y = Count)
) +
  geom_alluvium(
    aes(fill = Cluster),
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
    data = subset(cluster_box_data, Cluster == "Cluster1"),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = cluster_box_colors["Cluster1"],
    color = NA,
    linewidth = 0
  ) +
  geom_rect(
    data = subset(cluster_box_data, Cluster == "Cluster2"),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = cluster_box_colors["Cluster2"],
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
    limits = c("Cluster subtype", "Risk group"),
    expand = c(0.10, 0.10)
  ) +
  scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.02, 0.04))
  ) +
  scale_fill_manual(values = fill_colors, guide = "none") +
  labs(
    title = "Association Between Cluster Subtypes and Risk Groups",
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

# -----------------------------
# 11. 保存 PDF 和 SVG
# -----------------------------
pdf_file <- file.path(output_dir, "Cluster_RiskGroup_Sankey_Beautified.pdf")
svg_file <- file.path(output_dir, "Cluster_RiskGroup_Sankey_Beautified.svg")

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

# -----------------------------
# 12. 输出提示
# -----------------------------
cat("✅ PDF已保存：", pdf_file, "\n")
cat("✅ SVG已保存：", svg_file, "\n")
cat("✅ 样本级对应表已保存：", file.path(output_dir, "Sample_Level_Cluster_RiskGroup_Mapping.csv"), "\n")
cat("✅ 汇总对应表已保存：", file.path(output_dir, "Summary_Cluster_RiskGroup_Mapping.csv"), "\n")
}  # 循环结束
