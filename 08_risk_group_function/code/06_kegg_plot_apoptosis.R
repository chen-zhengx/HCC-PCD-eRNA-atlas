# =========================================================
# 对应结果：Figure S4C（高风险组上调基因 KEGG 富集条形图）、Figure S4D（低风险组上调基因 KEGG 富集条形图）——凋亡 apoptosis
# 输入：output/figureS4/KEGG_ORA/High_up_KEGG_Plt0.05.csv、Low_up_KEGG_Plt0.05.csv（05_kegg_enrich_apoptosis.R 的输出）
# 输出：output/figureS4/KEGG/KEGG_High_up_Top_barplot_apoptosis.pdf/.svg（S4C）、KEGG_Low_up_Top_barplot_apoptosis.pdf/.svg（S4D）
#       KEGG_High_up_Top_for_plot.csv、KEGG_Low_up_Top_for_plot.csv
# 运行：在模块目录下 Rscript code/06_kegg_plot_apoptosis.R（需先运行 05）
# =========================================================
# ============================================================
# KEGG ORA 可视化（统一风格版，凋亡，PDF + SVG）
# 固定使用：
#   High_up_KEGG_Plt0.05.csv
#   Low_up_KEGG_Plt0.05.csv
# 全部按 pvalue 作图
# 使用通路名称 Description
# 图标题：显示
# 高低组统一显示前15个
# barplot 横轴为 Count
# barplot 填色为连续 -log10(Pvalue) 渐变，与 Figure 4 KEGG 一致
# dotplot 和合并发散图保留代码，但默认不输出
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

# ---------- 0) 路径 ----------
in_dir  <- "output/figureS4/KEGG_ORA"
out_dir <- "output/figureS4/KEGG"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

high_file <- file.path(in_dir, "High_up_KEGG_Plt0.05.csv")
low_file  <- file.path(in_dir, "Low_up_KEGG_Plt0.05.csv")

if (!file.exists(high_file)) stop("找不到文件：", high_file)
if (!file.exists(low_file))  stop("找不到文件：", low_file)

# ---------- 1) 参数 ----------
topN            <- 15
base_size       <- 12
axis_text_y_sz  <- 11
axis_text_x_sz  <- 13
legend_text_sz  <- 12
wrap_width      <- 28
plot_width      <- 12
plot_height     <- 9.5
div_width       <- 12
div_height      <- 9.5
left_margin_pt  <- 90

# 保持凋亡配色
# Plot colors follow the Figure 4 KEGG gradient.
low_color  <- "#9ECBF4"
high_color <- "#D84B3A"
col_high <- high_color
col_low  <- low_color

# ---------- 2) p值图例标签格式 ----------
format_p_legend <- function(pmax) {
  txt <- sprintf("%.2e", pmax)
  txt <- gsub("-", "−", txt, fixed = TRUE)
  paste0("p_value<", txt)
}

# ---------- 3) 统一列名 ----------
normalize_cols <- function(dt){
  if (!("Description" %in% names(dt))) {
    if ("Term" %in% names(dt)) {
      setnames(dt, "Term", "Description")
    } else if ("Pathway" %in% names(dt)) {
      setnames(dt, "Pathway", "Description")
    } else if ("Name" %in% names(dt)) {
      setnames(dt, "Name", "Description")
    }
  }
  
  if (!("pvalue" %in% names(dt))) {
    if ("p.value" %in% names(dt)) setnames(dt, "p.value", "pvalue")
    if ("PValue" %in% names(dt))  setnames(dt, "PValue", "pvalue")
  }
  
  if (!("p.adjust" %in% names(dt))) {
    if ("padj" %in% names(dt)) setnames(dt, "padj", "p.adjust")
    if ("FDR" %in% names(dt))  setnames(dt, "FDR", "p.adjust")
  }
  
  if (!("Count" %in% names(dt)) && ("count" %in% names(dt))) {
    setnames(dt, "count", "Count")
  }
  
  if (!("GeneRatio" %in% names(dt)) && ("GeneRatio_num" %in% names(dt))) {
    dt[, GeneRatio := as.character(GeneRatio_num)]
  }
  
  dt
}

# ---------- 4) 保存：PDF + SVG ----------
save_pdf_svg <- function(p, file_base, w, h){
  # Output: PDF first, using grDevices::pdf for stable Windows export.
  ggsave(
    filename = file.path(out_dir, paste0(file_base, ".pdf")),
    plot = p, width = w, height = h,
    device = grDevices::pdf,
    useDingbats = FALSE
  )
  
  # Output: SVG is optional and only written when svglite is installed.
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(
      filename = file.path(out_dir, paste0(file_base, ".svg")),
      plot = p, width = w, height = h, device = svglite::svglite
    )
  } else {
    message("未安装 svglite：已跳过 SVG 输出")
  }
}

# ---------- 5) 统一主题 ----------
theme_kegg_clean <- function(){
  theme_bw() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = axis_text_y_sz, colour = "black", lineheight = 0.88),
      axis.text.x = element_text(size = axis_text_x_sz, colour = "black"),
      axis.title.x = element_text(size = 16, colour = "black"),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = legend_text_sz),
      plot.margin = margin(t = 15, r = 20, b = 15, l = left_margin_pt)
    )
}

# ---------- 6) 标签处理 ----------
wrap_desc <- function(x, width = wrap_width){
  stringr::str_wrap(as.character(x), width = width)
}

parse_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(z) {
    if (length(z) == 2) {
      as.numeric(z[1]) / as.numeric(z[2])
    } else {
      NA_real_
    }
  })
}

prep_top_dt <- function(dt, topN = 15){
  dt <- copy(dt)
  dt <- dt[!is.na(pvalue)]
  dt <- dt[order(pvalue)]
  dt <- dt[1:min(topN, .N)]
  
  # Plot color uses the same continuous -log10(Pvalue) scale as Figure 4.
  dt[, neg_log10_p := -log10(pmax(pvalue, .Machine$double.xmin))]
  dt[, sig := neg_log10_p]
  dt[, Description_raw := Description]
  dt[, Description := wrap_desc(Description_raw)]
  dt[, Description2 := factor(Description, levels = rev(Description))]
  
  if ("GeneRatio" %in% names(dt)) {
    dt[, GeneRatio_num := parse_ratio(GeneRatio)]
  } else {
    dt[, GeneRatio_num := NA_real_]
  }
  
  if (!("Count" %in% names(dt))) {
    dt[, Count := 1]
  }
  
  dt[, legend_label := format_p_legend(max(pvalue, na.rm = TRUE))]
  
  dt
}

# ---------- 7) 读入数据 ----------
hi <- normalize_cols(fread(high_file))
lo <- normalize_cols(fread(low_file))

req <- c("Description", "pvalue")
if (!all(req %in% names(hi))) stop("High_up 文件缺少列：", paste(setdiff(req, names(hi)), collapse = ", "))
if (!all(req %in% names(lo))) stop("Low_up 文件缺少列：", paste(setdiff(req, names(lo)), collapse = ", "))

for (col in intersect(c("pvalue", "p.adjust", "Count"), names(hi))) {
  hi[, (col) := as.numeric(get(col))]
}
for (col in intersect(c("pvalue", "p.adjust", "Count"), names(lo))) {
  lo[, (col) := as.numeric(get(col))]
}

hi_top <- prep_top_dt(hi, topN = topN)
lo_top <- prep_top_dt(lo, topN = topN)

# ---------- 8) barplot：横轴 Count ----------
make_barplot_manual <- function(dt, fill_color, title_text){
  x_break_max <- max(as.numeric(dt$Count), na.rm = TRUE)
  if (!is.finite(x_break_max)) x_break_max <- 0
  x_breaks <- seq(0, ceiling(x_break_max / 2) * 2, by = 2)
  
  ggplot(dt, aes(x = Count, y = Description2, fill = neg_log10_p)) +
    geom_col(width = 0.7) +
    scale_x_continuous(breaks = x_breaks) +
    scale_fill_gradient(
      low = low_color,
      high = high_color,
      name = expression(-log[10](Pvalue))
    ) +
    labs(
      title = title_text,
      x = "Count",
      y = NULL
    ) +
    theme_kegg_clean()
}

# ---------- 9) dotplot ----------
make_dotplot_manual <- function(dt, point_color, title_text){
  ggplot(dt, aes(x = GeneRatio_num, y = Description2)) +
    geom_point(aes(size = Count, colour = neg_log10_p)) +
    scale_colour_gradient(
      low = low_color,
      high = high_color,
      name = expression(-log[10](Pvalue))
    ) +
    labs(
      title = title_text,
      x = "GeneRatio",
      y = NULL,
      size = "Count"
    ) +
    guides(
      colour = guide_colorbar(order = 1),
      size = guide_legend(order = 2)
    ) +
    theme_kegg_clean()
}

# ---------- 10) High_up 条形图 ----------
p_hi_bar <- make_barplot_manual(
  hi_top,
  fill_color = col_high,
  title_text = "KEGG enrichment of upregulated genes"
)
# Active output: barplot, matching the Figure 4 KEGG script.
save_pdf_svg(p_hi_bar, "KEGG_High_up_Top_barplot_apoptosis", w = plot_width, h = plot_height)

# ---------- 11) Low_up 条形图 ----------
p_lo_bar <- make_barplot_manual(
  lo_top,
  fill_color = col_low,
  title_text = "KEGG enrichment of downregulated genes"
)
# Active output: barplot, matching the Figure 4 KEGG script.
save_pdf_svg(p_lo_bar, "KEGG_Low_up_Top_barplot_apoptosis", w = plot_width, h = plot_height)

# ---------- 12) High_up dotplot ----------
p_hi_dot <- make_dotplot_manual(
  hi_top,
  point_color = col_high,
  title_text = "KEGG enrichment of upregulated genes"
)

# ---------- 13) Low_up dotplot ----------
p_lo_dot <- make_dotplot_manual(
  lo_top,
  point_color = col_low,
  title_text = "KEGG enrichment of downregulated genes"
)

# Disabled output: dotplot code is retained for optional reuse, but not written.
if (FALSE) save_pdf_svg(p_hi_dot, "KEGG_High_up_Top_dotplot_apoptosis", w = plot_width, h = plot_height)
if (FALSE) save_pdf_svg(p_lo_dot, "KEGG_Low_up_Top_dotplot_apoptosis",  w = plot_width, h = plot_height)

# ---------- 14) 合并发散条形图 ----------
hi_merge <- copy(hi_top)[, .(Group = "High_up", Description_raw, score = sig)]
lo_merge <- copy(lo_top)[, .(Group = "Low_up",  Description_raw, score = -sig)]
merge_df <- rbindlist(list(hi_merge, lo_merge), use.names = TRUE, fill = TRUE)

ord <- merge_df[, .(maxabs = max(abs(score), na.rm = TRUE)), by = Description_raw][order(maxabs)]
merge_df[, Description := wrap_desc(Description_raw)]
merge_df[, Description2 := factor(Description, levels = wrap_desc(ord$Description_raw))]

p_div <- ggplot(merge_df, aes(x = score, y = Description2, fill = Group)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.6, color = "grey40") +
  scale_fill_manual(values = c("High_up" = col_high, "Low_up" = col_low)) +
  labs(
    title = "KEGG enrichment of upregulated and downregulated genes",
    x = "-log10(pvalue) (Low shown as negative)",
    y = NULL
  ) +
  theme_kegg_clean()

# Disabled output: the combined diverging plot is not part of the Figure 4-style output.
if (FALSE) save_pdf_svg(p_div, "KEGG_Diverging_barplot_High_vs_Low_apoptosis", w = div_width, h = div_height)

# ---------- 15) 导出 Top 表 ----------
fwrite(hi_top, file.path(out_dir, "KEGG_High_up_Top_for_plot.csv"))
fwrite(lo_top, file.path(out_dir, "KEGG_Low_up_Top_for_plot.csv"))

cat("Done!\n输出目录：", out_dir, "\n", sep = "")

