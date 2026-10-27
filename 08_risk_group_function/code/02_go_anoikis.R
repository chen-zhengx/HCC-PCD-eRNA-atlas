# =========================================================
# 对应结果：Figure 4A（高风险组上调基因 GO BP 富集）、Figure 4B（低风险组上调基因 GO BP 富集）——失巢凋亡 anoikis
# 输入：output/limma/anoikis/High_vs_Low_limma_all.csv（背景基因）、High_up_sig_FDR0.01_logFC1.csv
#       Low_up_sig_FDR0.01_logFC1.csv
# 输出：output/figure4/GO/High_up_GO_BP_all.csv、High_up_GO_BP_Plt0.05.csv、High_up_GO_BP_dotplot_p_lt_0.05.pdf/.svg（4A）
#       Low_up_GO_BP_all.csv、Low_up_GO_BP_Plt0.05.csv、Low_up_GO_BP_dotplot_p_lt_0.05.pdf/.svg（4B）
# 运行：在模块目录下 Rscript code/02_go_anoikis.R（需先运行 01）
# =========================================================
# =========================================================
# GO ORA (enrichGO) for High_up / Low_up
# 仅保留 BP
# 输出：全量结果 + pvalue筛选结果 + dotplot 的 PDF 和 SVG
# 高低组尺寸/字体统一
# 高低组都只显示前15个条目
# 图例改为单一色块，不再显示连续渐变
# 图例标签固定为：p_value<0.05
# Count 图例强制使用整数刻度
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(stringr)
})

# ---------- 0) 路径 ----------
limma_all_file <- "output/limma/anoikis/High_vs_Low_limma_all.csv"
high_up_file   <- "output/limma/anoikis/High_up_sig_FDR0.01_logFC1.csv"
low_up_file    <- "output/limma/anoikis/Low_up_sig_FDR0.01_logFC1.csv"

out_dir <- "output/figure4/GO"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 1) 阈值 ----------
P_cut <- 0.05

# ---------- 2) 固定颜色 ----------
low_color  <- "#9ECBF4"
high_color <- "#D84B3A"
fixed_color <- high_color

# ---------- 3) p值图例标签自动格式化 ----------
format_p_legend <- function(pmax) {
  txt <- sprintf("%.2e", pmax)      # 例如 1.23e-02 / 7.88e-16
  txt <- gsub("-", "−", txt, fixed = TRUE)
  paste0("p_value<", txt)
}

# ---------- 4) Count 图例使用整数 ----------
count_breaks_arithmetic <- function(x, n = 4) {
  x <- as.numeric(x[!is.na(x)])

  if (length(x) == 0) {
    return(numeric(0))
  }

  if (length(unique(x)) == 1) {
    return(as.integer(round(unique(x))))
  }

  rng <- range(x)
  brks <- pretty(rng, n = n)
  brks <- brks[brks >= rng[1] & brks <= rng[2]]
  brks <- brks[abs(brks - round(brks)) < 1e-8]
  brks <- unique(as.integer(round(brks)))

  if (length(brks) < 2) {
    brks <- seq(from = floor(rng[1]), to = ceiling(rng[2]), by = 1)
    brks <- brks[brks >= rng[1] & brks <= rng[2]]
  }

  if (length(brks) < 2) {
    brks <- unique(as.integer(c(floor(rng[1]), ceiling(rng[2]))))
  }

  brks
}

# ---------- 5) 读入基因 ----------
res_all <- fread(limma_all_file)
high_up <- fread(high_up_file)
low_up  <- fread(low_up_file)

stopifnot("Gene" %in% names(res_all), "Gene" %in% names(high_up), "Gene" %in% names(low_up))

genes_high <- unique(na.omit(high_up$Gene))
genes_low  <- unique(na.omit(low_up$Gene))
bg_genes   <- unique(na.omit(res_all$Gene))

# ---------- 6) SYMBOL -> ENTREZID ----------
sym2ent <- function(sym_vec) {
  sym_vec <- unique(sym_vec)
  m <- bitr(sym_vec, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  m <- unique(m[, c("SYMBOL", "ENTREZID")])
  list(
    map    = m,
    entrez = unique(m$ENTREZID),
    n_in   = length(sym_vec),
    n_map  = length(unique(m$ENTREZID))
  )
}

m_high <- sym2ent(genes_high)
m_low  <- sym2ent(genes_low)
m_bg   <- sym2ent(bg_genes)

cat("Mapping summary:\n")
cat("High: in=", m_high$n_in, " mapped_entrez=", m_high$n_map, "\n", sep = "")
cat("Low : in=", m_low$n_in,  " mapped_entrez=", m_low$n_map,  "\n", sep = "")
cat("BG  : in=", m_bg$n_in,   " mapped_entrez=", m_bg$n_map,   "\n", sep = "")

# ---------- 7) 更短术语缩写 ----------
abbr_go_terms <- function(x) {
  x <- gsub("^mitotic sister chromatid segregation$", "mitotic sister chr. seg.", x)
  x <- gsub("^sister chromatid segregation$", "sister chr. seg.", x)
  x <- gsub("^nuclear chromosome segregation$", "nuclear chr. seg.", x)
  x <- gsub("^chromosome segregation$", "chromosome seg.", x)
  
  x <- gsub("^microtubule cytoskeleton organization involved in mitosis$",
            "microtubule org. in mitosis", x)
  
  x <- gsub("^regulation of chromosome segregation$", "reg. chr. seg.", x)
  x <- gsub("^regulation of nuclear division$", "reg. nuclear division", x)
  x <- gsub("^regulation of mitotic nuclear division$", "reg. mitotic nuclear division", x)
  
  x <- gsub("^regulation of mitotic metaphase/anaphase transition$",
            "reg. mitotic M/A trans.", x)
  x <- gsub("^metaphase/anaphase transition of mitotic cell cycle$",
            "M/A trans. of mitotic cycle", x)
  x <- gsub("^metaphase/anaphase transition of cell cycle$",
            "M/A trans. of cell cycle", x)
  
  x <- gsub("^regulation of mitotic sister chromatid separation$",
            "reg. mitotic sister chr. sep.", x)
  x <- gsub("^mitotic sister chromatid separation$", "mitotic sister chr. sep.", x)
  
  x <- gsub("^mitotic spindle organization$", "mitotic spindle org.", x)
  x <- gsub("^spindle organization$", "spindle organization", x)
  x <- gsub("^mitotic spindle assembly$", "mitotic spindle assembly", x)
  
  return(x)
}

# ---------- 8) 自动换行 ----------
wrap_terms <- function(x, width = 24) {
  stringr::str_wrap(x, width = width)
}

# ---------- 9) GeneRatio 转数值 ----------
parse_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(z) {
    if (length(z) == 2) {
      as.numeric(z[1]) / as.numeric(z[2])
    } else {
      NA_real_
    }
  })
}

# ---------- 10) 统一作图参数 ----------
plot_width   <- 12
plot_height  <- 9
y_text_size  <- 11
y_lineheight <- 0.88
left_margin  <- 90
wrap_width   <- 24

# ---------- 11) 手动画 dotplot ----------
make_dotplot_manual <- function(df_plot, title_text) {
  ggplot(df_plot, aes(x = GeneRatio_num, y = Description)) +
    geom_point(aes(size = Count, colour = neg_log10_p)) +
    scale_size_continuous(
      breaks = count_breaks_arithmetic(df_plot$Count),
      labels = function(x) as.integer(round(x)),
      range = c(4, 12)
    ) +
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
      colour = guide_colorbar(
        order = 1
      ),
      size = guide_legend(order = 2)
    ) +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_blank(),
      axis.text.x  = element_text(size = 13, color = "black"),
      axis.text.y  = element_text(size = y_text_size, color = "black", lineheight = y_lineheight),
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12),
      plot.margin  = margin(t = 15, r = 20, b = 15, l = left_margin)
    )
}

# ---------- 12) 手动画 barplot ----------
make_barplot_manual <- function(df_plot, title_text) {
  this_label <- unique(df_plot$legend_label)
  
  ggplot(df_plot, aes(x = Count, y = Description)) +
    geom_col(aes(fill = legend_label), width = 0.7) +
    scale_fill_manual(
      values = setNames(fixed_color, this_label),
      breaks = this_label,
      labels = this_label,
      name = NULL
    ) +
    labs(
      title = title_text,
      x = "Count",
      y = NULL
    ) +
    guides(
      fill = guide_legend(order = 1)
    ) +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_blank(),
      axis.text.x  = element_text(size = 13, color = "black"),
      axis.text.y  = element_text(size = y_text_size, color = "black", lineheight = y_lineheight),
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12),
      plot.margin  = margin(t = 15, r = 20, b = 15, l = left_margin)
    )
}

# ---------- 13) 只跑 BP ----------
run_enrichGO_BP <- function(entrez_vec, universe_entrez, tag, title_text) {
  
  eg <- enrichGO(
    gene          = entrez_vec,
    universe      = universe_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
  
  df <- as.data.frame(eg)
  
  # 全量输出
  fwrite(df, file.path(out_dir, sprintf("%s_GO_BP_all.csv", tag)))
  
  if (nrow(df) == 0) {
    cat(sprintf("[%s-BP] No terms.\n", tag))
    return(invisible(NULL))
  }
  
  # 按 pvalue 筛选并排序
  df_p <- df[!is.na(df$pvalue) & df$pvalue < P_cut, ]
  df_p <- df_p[order(df_p$pvalue), ]
  
  # 输出筛选结果
  fwrite(df_p, file.path(out_dir, sprintf("%s_GO_BP_Plt%s.csv", tag, P_cut)))
  
  if (nrow(df_p) > 0) {
    showN <- min(15, nrow(df_p))
    df_plot <- df_p[1:showN, ]
    
    # 缩写 + 换行
    df_plot$Description <- abbr_go_terms(df_plot$Description)
    df_plot$Description <- wrap_terms(df_plot$Description, width = wrap_width)
    
    # GeneRatio 转数值
    df_plot$GeneRatio_num <- parse_ratio(df_plot$GeneRatio)
    
    # Map p values to a continuous color legend.
    df_plot$neg_log10_p <- -log10(pmax(df_plot$pvalue, .Machine$double.xmin))
    
    # 图例标签固定显示筛选阈值
    df_plot$legend_label <- format_p_legend(P_cut)
    
    # 因子顺序：显著性最高的放最上面
    df_plot$Description <- factor(df_plot$Description, levels = rev(df_plot$Description))
    
    # 手动画图
    p_dot <- make_dotplot_manual(df_plot, title_text)
    # p_bar <- make_barplot_manual(df_plot, title_text)
    
    # dotplot PDF
    ggsave(
      filename = file.path(out_dir, sprintf("%s_GO_BP_dotplot_p_lt_%s.pdf", tag, P_cut)),
      plot = p_dot,
      width = plot_width,
      height = plot_height
    )
    
    # dotplot SVG
    if (requireNamespace("svglite", quietly = TRUE)) {
      ggsave(
        filename = file.path(out_dir, sprintf("%s_GO_BP_dotplot_p_lt_%s.svg", tag, P_cut)),
        plot = p_dot,
        width = plot_width,
        height = plot_height,
        device = svglite::svglite
      )
    } else {
      message("未安装 svglite：已跳过 SVG 输出（可运行 install.packages('svglite')）")
    }
    
    # # barplot PDF
    # ggsave(
    #   filename = file.path(out_dir, sprintf("%s_GO_BP_barplot_p_lt_%s.pdf", tag, P_cut)),
    #   plot = p_bar,
    #   width = plot_width,
    #   height = plot_height
    # )
    #
    # # barplot SVG
    # if (requireNamespace("svglite", quietly = TRUE)) {
    #   ggsave(
    #     filename = file.path(out_dir, sprintf("%s_GO_BP_barplot_p_lt_%s.svg", tag, P_cut)),
    #     plot = p_bar,
    #     width = plot_width,
    #     height = plot_height,
    #     device = svglite::svglite
    #   )
    # } else {
    #   message("未安装 svglite：已跳过 SVG 输出（可运行 install.packages('svglite')）")
    # }
  }
  
  cat(sprintf("[%s-BP] all=%d | p<%.3g: %d\n", tag, nrow(df), P_cut, nrow(df_p)))
  invisible(eg)
}

# ---------- 14) 分别跑 High / Low ----------
ego_high_bp <- run_enrichGO_BP(
  entrez_vec = m_high$entrez,
  universe_entrez = m_bg$entrez,
  tag = "High_up",
  title_text = "GO Biological Process enrichment of upregulated genes"
)

ego_low_bp <- run_enrichGO_BP(
  entrez_vec = m_low$entrez,
  universe_entrez = m_bg$entrez,
  tag = "Low_up",
  title_text = "GO Biological Process enrichment of downregulated genes"
)

cat("\n✅ GO BP ORA 完成。输出目录：", out_dir, "\n", sep = "")
cat("输出文件说明：\n")
cat("  *_all.csv       = 全量条目\n")
cat("  *_Plt0.05.csv   = 按 pvalue < 0.05 筛选后的条目\n")
cat("  dotplot         = 单一颜色图例（1个色块），同时输出 PDF 和 SVG\n")
