# =========================================================
# 对应结果：Figure 4C（高风险组上调基因 KEGG 富集）、Figure 4D（低风险组上调基因 KEGG 富集）——失巢凋亡 anoikis
# 输入：output/limma/anoikis/High_vs_Low_limma_all.csv（背景基因）、High_up_sig_FDR0.01_logFC1.csv
#       Low_up_sig_FDR0.01_logFC1.csv
# 输出：output/figure4/KEGG/High_up_KEGG_all.csv、High_up_KEGG_Plt0.05.csv、High_up_KEGG_barplot_p_lt_0.05.pdf/.svg（4C）
#       Low_up_KEGG_all.csv、Low_up_KEGG_Plt0.05.csv、Low_up_KEGG_barplot_p_lt_0.05.pdf/.svg（4D），SYMBOL_to_ENTREZ_*.csv
# 运行：在模块目录下 Rscript code/03_kegg_anoikis.R（需先运行 01）
# =========================================================
# =========================================================
# KEGG ORA (clusterProfiler::enrichKEGG)
# 仅保留 pvalue < 0.05 的绘图结果
# 输出：全量结果 + pvalue筛选结果 + barplot 的 PDF 和 SVG
# 高低组尺寸/字体统一
# 高低组都只显示前15个条目
# 失巢死亡版本
# 保留 -log10(Pvalue) 图例，但颜色统一成一个颜色（非渐变）
# 使用 ggplot2 手动画图，避免 enrichplot::dotplot 渐变色残留
# =========================================================

# # 手动下载 KEGG.db 包
# download.file("https://bioconductor.org/packages/3.6/data/annotation/src/contrib/KEGG.db_3.2.3.tar.gz",
#               destfile = "KEGG.db_3.2.3.tar.gz")
# # 本地安装
# install.packages("KEGG.db_3.2.3.tar.gz", repos = NULL, type = "source")

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(stringr)
  library(KEGG.db)
})

# ---------- 0) 路径 ----------
limma_all_file <- "output/limma/anoikis/High_vs_Low_limma_all.csv"
high_up_file   <- "output/limma/anoikis/High_up_sig_FDR0.01_logFC1.csv"
low_up_file    <- "output/limma/anoikis/Low_up_sig_FDR0.01_logFC1.csv"

out_dir <- "output/figure4/KEGG"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 1) 阈值 ----------
P_cut <- 0.05
low_color  <- "#9ECBF4"
high_color <- "#D84B3A"
fixed_color <- high_color

# ---------- 2) 工具函数 ----------
normalize_cols <- function(dt, tag) {
  if (!("pvalue" %in% names(dt))) {
    if ("p.value" %in% names(dt)) setnames(dt, "p.value", "pvalue")
    if ("PValue" %in% names(dt))  setnames(dt, "PValue", "pvalue")
  }
  
  if (!("p.adjust" %in% names(dt))) {
    if ("padj" %in% names(dt))       setnames(dt, "padj", "p.adjust")
    if ("p.adjusted" %in% names(dt)) setnames(dt, "p.adjusted", "p.adjust")
    if ("FDR" %in% names(dt))        setnames(dt, "FDR", "p.adjust")
  }
  
  if (!("pvalue" %in% names(dt))) {
    stop(sprintf("[%s] missing pvalue. cols=%s", tag, paste(names(dt), collapse = ", ")))
  }
  if (!("p.adjust" %in% names(dt))) {
    stop(sprintf("[%s] missing p.adjust. cols=%s", tag, paste(names(dt), collapse = ", ")))
  }
  
  dt
}

# KEGG ID -> 通路名
kegg_pathid2name_map <- function() {
  envs <- list()
  
  if ("package:KEGG.db" %in% search()) {
    envs <- c(envs, list(as.environment("package:KEGG.db")))
  }
  
  if (requireNamespace("KEGG.db", quietly = TRUE)) {
    envs <- c(envs, list(asNamespace("KEGG.db")))
  }
  
  for (env in envs) {
    if (exists("KEGGPATHID2NAME", envir = env, inherits = FALSE)) {
      return(get("KEGGPATHID2NAME", envir = env, inherits = FALSE))
    }
  }
  
  data_env <- new.env(parent = emptyenv())
  data_loaded <- tryCatch({
    suppressWarnings(data("KEGGPATHID2NAME", package = "KEGG.db", envir = data_env))
  }, error = function(e) character())
  
  if (length(data_loaded) > 0 &&
      exists("KEGGPATHID2NAME", envir = data_env, inherits = FALSE)) {
    return(get("KEGGPATHID2NAME", envir = data_env, inherits = FALSE))
  }
  
  NULL
}

kegg_id_to_name <- function(id_vec) {
  id_vec <- as.character(id_vec)
  map_obj <- kegg_pathid2name_map()
  id2name <- tryCatch(as.list(map_obj), error = function(e) NULL)
  
  if (is.null(id2name)) {
    warning("KEGGPATHID2NAME not found in KEGG.db; keep KEGG IDs as Description.")
    return(id_vec)
  }
  
  out <- rep(NA_character_, length(id_vec))
  
  for (i in seq_along(id_vec)) {
    x <- id_vec[i]
    if (is.na(x) || !nzchar(x)) {
      out[i] <- x
      next
    }
    
    keys <- unique(c(gsub("^hsa", "", x), x))
    hit <- NULL
    
    for (key in keys) {
      if (!is.null(id2name[[key]])) {
        hit <- id2name[[key]]
        break
      }
    }
    
    if (!is.null(hit) && length(hit) > 0 && !is.na(hit[1])) {
      out[i] <- as.character(hit[1])
    } else {
      out[i] <- x
    }
  }
  
  out
}

sym2ent <- function(sym_vec) {
  sym_vec <- unique(na.omit(sym_vec))
  m <- bitr(sym_vec, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  m <- unique(m[, c("SYMBOL", "ENTREZID")])
  
  list(
    map    = as.data.table(m),
    entrez = unique(m$ENTREZID),
    n_in   = length(sym_vec),
    n_map  = length(unique(m$ENTREZID))
  )
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

# ---------- 3) 名称缩写 ----------
abbr_kegg_terms <- function(x) {
  x <- gsub("^ECM-receptor interaction$", "ECM-receptor interaction", x)
  x <- gsub("^PI3K-Akt signaling pathway$", "PI3K-Akt signaling pathway", x)
  x <- gsub("^Cytokine-cytokine receptor interaction$", "Cytokine-receptor interaction", x)
  x <- gsub("^Cell adhesion molecules$", "Cell adhesion molecules", x)
  x <- gsub("^Protein digestion and absorption$", "Protein digestion and absorption", x)
  x <- gsub("^Neuroactive ligand-receptor interaction$", "Neuroactive ligand-receptor interaction", x)
  x <- gsub("^Complement and coagulation cascades$", "Complement/coagulation cascades", x)
  x <- gsub("^Metabolism of xenobiotics by cytochrome P450$", "Xenobiotic metabolism by CYP", x)
  x <- gsub("^Drug metabolism - cytochrome P450$", "Drug metabolism - CYP", x)
  x <- gsub("^Chemical carcinogenesis - DNA adducts$", "Chemical carcinogenesis - DNA adducts", x)
  x <- gsub("^Chemical carcinogenesis - receptor activation$", "Chemical carcinogenesis - receptor act.", x)
  x <- gsub("^Chemical carcinogenesis - reactive oxygen species$", "Chemical carcinogenesis - ROS", x)
  x <- gsub("^Arachidonic acid metabolism$", "Arachidonic acid metabolism", x)
  x <- gsub("^Retinol metabolism$", "Retinol metabolism", x)
  x <- gsub("^Drug metabolism - other enzymes$", "Drug metabolism - other enzymes", x)
  x <- gsub("^AGE-RAGE signaling pathway in diabetic complications$", "AGE-RAGE signaling in diabetic complications", x)
  x <- gsub("^cAMP signaling pathway$", "cAMP signaling pathway", x)
  x <- gsub("^cGMP-PKG signaling pathway$", "cGMP-PKG signaling pathway", x)
  x <- gsub("^Calcium signaling pathway$", "Calcium signaling pathway", x)
  x <- gsub("^Focal adhesion$", "Focal adhesion", x)
  x <- gsub("^MAPK signaling pathway$", "MAPK signaling pathway", x)
  x <- gsub("^Rap1 signaling pathway$", "Rap1 signaling pathway", x)
  x <- gsub("^Regulation of actin cytoskeleton$", "Regulation of actin cytoskeleton", x)
  x <- gsub("^Cell cycle$", "Cell cycle", x)
  x <- gsub("^Oocyte meiosis$", "Oocyte meiosis", x)
  x <- gsub("^Ubiquitin mediated proteolysis$", "Ubiquitin-mediated proteolysis", x)
  x <- gsub("^Inositol phosphate metabolism$", "Inositol phosphate metabolism", x)
  return(x)
}

# ---------- 4) 自动换行 ----------
wrap_terms <- function(x, width = 28) {
  stringr::str_wrap(x, width = width)
}

# ---------- 5) 统一作图参数 ----------
plot_width   <- 12
plot_height  <- 9.5
y_text_size  <- 11
y_lineheight <- 0.88
left_margin  <- 90
wrap_width   <- 28

# ---------- 6) 读入基因 ----------
res_all <- fread(limma_all_file)
high_up <- fread(high_up_file)
low_up  <- fread(low_up_file)

stopifnot(
  "Gene" %in% names(res_all),
  "Gene" %in% names(high_up),
  "Gene" %in% names(low_up)
)

genes_high <- unique(na.omit(high_up$Gene))
genes_low  <- unique(na.omit(low_up$Gene))
bg_genes   <- unique(na.omit(res_all$Gene))

# ---------- 7) SYMBOL -> ENTREZ ----------
m_high <- sym2ent(genes_high)
m_low  <- sym2ent(genes_low)
m_bg   <- sym2ent(bg_genes)

fwrite(m_high$map, file.path(out_dir, "SYMBOL_to_ENTREZ_High_up.csv"))
fwrite(m_low$map,  file.path(out_dir, "SYMBOL_to_ENTREZ_Low_up.csv"))
fwrite(m_bg$map,   file.path(out_dir, "SYMBOL_to_ENTREZ_Background.csv"))

cat("Mapping summary:\n")
cat("High: in=", m_high$n_in, " mapped_entrez=", m_high$n_map, "\n", sep = "")
cat("Low : in=", m_low$n_in,  " mapped_entrez=", m_low$n_map,  "\n", sep = "")
cat("BG  : in=", m_bg$n_in,   " mapped_entrez=", m_bg$n_map,   "\n", sep = "")

# ---------- 8) 手动画 dotplot ----------
make_dotplot_manual <- function(df_plot, title_text) {
  ggplot(df_plot, aes(x = GeneRatio_num, y = Description)) +
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

# ---------- 9) 手动画 barplot ----------
make_barplot_manual <- function(df_plot, title_text) {
  x_break_max <- max(as.numeric(df_plot$Count), na.rm = TRUE)
  if (!is.finite(x_break_max)) x_break_max <- 0
  x_breaks <- seq(0, ceiling(x_break_max / 2) * 2, by = 2)
  
  ggplot(df_plot, aes(x = Count, y = Description, fill = neg_log10_p)) +
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

# ---------- 10) KEGG ORA + 作图 ----------
run_kegg_with_plot <- function(entrez_vec, universe_entrez, tag, title_text) {
  
  ek <- enrichKEGG(
    gene              = entrez_vec,
    universe          = universe_entrez,
    organism          = "hsa",
    keyType           = "kegg",
    pAdjustMethod     = "BH",
    pvalueCutoff      = 1,
    qvalueCutoff      = 1,
    use_internal_data = TRUE
  )
  
  df <- as.data.table(as.data.frame(ek))
  
  if (nrow(df) == 0) {
    fwrite(df, file.path(out_dir, sprintf("%s_KEGG_all.csv", tag)))
    cat(sprintf("[%s-KEGG] No pathways.\n", tag))
    return(invisible(NULL))
  }
  
  df <- normalize_cols(df, paste0(tag, "-KEGG"))
  
  # 修正 Description
  if ("ID" %in% names(df)) {
    need_fix <- FALSE
    
    if (!("Description" %in% names(df))) {
      need_fix <- TRUE
    } else {
      desc_chr <- as.character(df$Description)
      if (all(is.na(desc_chr)) || all(grepl("^hsa\\d+$", desc_chr))) {
        need_fix <- TRUE
      }
    }
    
    if (need_fix) {
      df[, Description := kegg_id_to_name(ID)]
    }
  }
  
  # 保存全量结果
  fwrite(df, file.path(out_dir, sprintf("%s_KEGG_all.csv", tag)))
  
  # pvalue筛选
  df_p <- df[!is.na(pvalue) & pvalue < P_cut]
  df_p <- df_p[order(pvalue)]
  
  fwrite(df_p, file.path(out_dir, sprintf("%s_KEGG_Plt%s.csv", tag, P_cut)))
  
  if (nrow(df_p) > 0) {
    showN <- min(15, nrow(df_p))
    df_plot <- copy(df_p[1:showN])
    
    # 缩写 + 换行
    df_plot[, Description := abbr_kegg_terms(Description)]
    df_plot[, Description := wrap_terms(Description, width = wrap_width)]
    
    # GeneRatio 转数值
    df_plot[, GeneRatio_num := parse_ratio(GeneRatio)]
    
    # p values for continuous legend.
    df_plot[, neg_log10_p := -log10(pmax(pvalue, .Machine$double.xmin))]
    
    # 因子顺序：最显著的放最上面
    df_plot[, Description := factor(Description, levels = rev(Description))]
    
    # p_dot <- make_dotplot_manual(df_plot, title_text)
    p_bar <- make_barplot_manual(df_plot, title_text)
    
    # # dotplot PDF
    # ggsave(
    #   filename = file.path(out_dir, sprintf("%s_KEGG_dotplot_p_lt_%s.pdf", tag, P_cut)),
    #   plot = p_dot,
    #   width = plot_width,
    #   height = plot_height
    # )
    
    # # dotplot SVG
    # if (requireNamespace("svglite", quietly = TRUE)) {
    #   ggsave(
    #     filename = file.path(out_dir, sprintf("%s_KEGG_dotplot_p_lt_%s.svg", tag, P_cut)),
    #     plot = p_dot,
    #     width = plot_width,
    #     height = plot_height,
    #     device = svglite::svglite
    #   )
    # } else {
    #   message("未安装 svglite：已跳过 SVG 输出（可运行 install.packages('svglite')）")
    # }
    
    # barplot PDF
    ggsave(
      filename = file.path(out_dir, sprintf("%s_KEGG_barplot_p_lt_%s.pdf", tag, P_cut)),
      plot = p_bar,
      width = plot_width,
      height = plot_height
    )
    
    # barplot SVG
    if (requireNamespace("svglite", quietly = TRUE)) {
      ggsave(
        filename = file.path(out_dir, sprintf("%s_KEGG_barplot_p_lt_%s.svg", tag, P_cut)),
        plot = p_bar,
        width = plot_width,
        height = plot_height,
        device = svglite::svglite
      )
    } else {
      message("未安装 svglite：已跳过 SVG 输出（可运行 install.packages('svglite')）")
    }
  }
  
  cat(sprintf("[%s-KEGG] all=%d | p<%.3g: %d\n", tag, nrow(df), P_cut, nrow(df_p)))
  invisible(ek)
}

# ---------- 11) 分别跑 High / Low ----------
ek_high <- run_kegg_with_plot(
  entrez_vec = m_high$entrez,
  universe_entrez = m_bg$entrez,
  tag = "High_up",
  title_text = "KEGG enrichment of upregulated genes"
)

ek_low <- run_kegg_with_plot(
  entrez_vec = m_low$entrez,
  universe_entrez = m_bg$entrez,
  tag = "Low_up",
  title_text = "KEGG enrichment of downregulated genes"
)

cat("\n✅ KEGG ORA 完成。输出目录：", out_dir, "\n", sep = "")
cat("输出文件说明：\n")
cat("  *_KEGG_all.csv               = 全量通路（已尽量修正 Description 为真正通路名）\n")
cat("  *_KEGG_Plt0.05.csv           = pvalue < 0.05\n")
# cat("  *_KEGG_dotplot_p_lt_0.05.*   = dotplot（PDF + SVG）\n")
cat("  *_KEGG_barplot_p_lt_0.05.*   = barplot（PDF + SVG）\n")
cat("  SYMBOL_to_ENTREZ_*.csv       = SYMBOL->ENTREZ 映射表\n")

