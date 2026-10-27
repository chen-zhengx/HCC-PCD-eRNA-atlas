# =========================================================
# 对应结果：Figure S5C/S5D 的 KEGG 富集计算（铁死亡 ferroptosis；高/低风险组上调基因的 enrichKEGG）
# 输入：output/limma/ferroptosis/High_vs_Low_limma_all.csv（背景基因）、High_up_sig_FDR0.01_logFC1.csv、Low_up_sig_FDR0.01_logFC1.csv（01 的输出）
# 输出：output/figureS5/KEGG_ORA/High_up_KEGG_*.csv、Low_up_KEGG_*.csv（all / Plt0.05 / FDRlt0.05 / Plt0.05_FDRlt0.05）
#       SYMBOL_to_ENTREZ_*.csv、KEGG_run_summary.csv；svg/*_KEGG_dotplot_top*.svg
# 运行：在模块目录下 Rscript code/08_kegg_enrich_ferroptosis.R（需先运行 01）
# =========================================================
# =========================================================
# KEGG ORA (clusterProfiler::enrichKEGG) - OFFLINE NO-PLOT
# High_up vs Low_up
# 失巢死亡风格版，输入输出为铁死亡
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(KEGG.db)
  library(ggplot2)
})

# ---------- 0) 路径 ----------
limma_all_file <- "output/limma/ferroptosis/High_vs_Low_limma_all.csv"
high_up_file   <- "output/limma/ferroptosis/High_up_sig_FDR0.01_logFC1.csv"
low_up_file    <- "output/limma/ferroptosis/Low_up_sig_FDR0.01_logFC1.csv"

out_dir <- "output/figureS5/KEGG_ORA"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 1) 阈值 ----------
P_cut   <- 0.05
FDR_cut <- 0.05

# ---------- 2) 工具函数 ----------
normalize_cols <- function(dt, tag){
  # pvalue
  if (!("pvalue" %in% names(dt))) {
    if ("p.value" %in% names(dt)) setnames(dt, "p.value", "pvalue")
    if ("PValue" %in% names(dt))  setnames(dt, "PValue", "pvalue")
  }
  # p.adjust
  if (!("p.adjust" %in% names(dt))) {
    if ("padj" %in% names(dt))       setnames(dt, "padj", "p.adjust")
    if ("p.adjusted" %in% names(dt)) setnames(dt, "p.adjusted", "p.adjust")
    if ("FDR" %in% names(dt))        setnames(dt, "FDR", "p.adjust")
  }
  dt
}

# KEGG ID -> 通路名称
kegg_id_to_name <- function(id_vec){
  id_vec <- as.character(id_vec)
  id2name <- as.list(KEGGPATHID2NAME)
  sapply(id_vec, function(x){
    x2 <- gsub("^hsa", "", x)
    if (!is.null(id2name[[x2]])) id2name[[x2]] else x
  })
}

sym2ent <- function(sym_vec){
  sym_vec <- unique(na.omit(sym_vec))
  m <- bitr(sym_vec, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  m <- unique(m[, c("SYMBOL","ENTREZID")])
  list(
    map    = as.data.table(m),
    entrez = unique(m$ENTREZID),
    n_in   = length(sym_vec),
    n_map  = length(unique(m$ENTREZID))
  )
}

# ---------- 3) 读取基因 ----------
res_all <- fread(limma_all_file)
high_up <- fread(high_up_file)
low_up  <- fread(low_up_file)

stopifnot("Gene" %in% names(res_all),
          "Gene" %in% names(high_up),
          "Gene" %in% names(low_up))

genes_high <- unique(na.omit(high_up$Gene))
genes_low  <- unique(na.omit(low_up$Gene))
bg_genes   <- unique(na.omit(res_all$Gene))

# ---------- 4) SYMBOL -> ENTREZID ----------
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

# ---------- 5) KEGG ORA（离线，无画图） ----------
run_kegg_offline_noplot <- function(entrez_vec, universe_entrez, tag){
  
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
    return(list(n_all=0, n_p=0, n_f=0, n_both=0))
  }
  
  df <- normalize_cols(df, paste0(tag, "-KEGG"))
  
  # 自动映射 KEGG ID -> 通路名称
  if ("ID" %in% names(df)) {
    df[, Description := kegg_id_to_name(ID)]
  }
  
  # 保存全量结果
  fwrite(df, file.path(out_dir, sprintf("%s_KEGG_all.csv", tag)))
  
  # 各种筛选
  df_p    <- df[!is.na(pvalue)   & pvalue   < P_cut]
  df_f    <- df[!is.na(p.adjust) & p.adjust < FDR_cut]
  df_both <- df[!is.na(pvalue) & !is.na(p.adjust) & pvalue < P_cut & p.adjust < FDR_cut]
  
  fwrite(df_p,    file.path(out_dir, sprintf("%s_KEGG_Plt%s.csv", tag, P_cut)))
  fwrite(df_f,    file.path(out_dir, sprintf("%s_KEGG_FDRlt%s.csv", tag, FDR_cut)))
  fwrite(df_both, file.path(out_dir, sprintf("%s_KEGG_Plt%s_FDRlt%s.csv", tag, P_cut, FDR_cut)))
  
  cat(sprintf("[%s-KEGG] all=%d | p<%.3g: %d | FDR<%.3g: %d | both: %d\n",
              tag, nrow(df), P_cut, nrow(df_p), FDR_cut, nrow(df_f), nrow(df_both)))
  
  list(n_all=nrow(df), n_p=nrow(df_p), n_f=nrow(df_f), n_both=nrow(df_both))
}

# ---------- 6) 分别跑 High / Low ----------
ret_high <- run_kegg_offline_noplot(m_high$entrez, m_bg$entrez, tag="High_up")
ret_low  <- run_kegg_offline_noplot(m_low$entrez,  m_bg$entrez, tag="Low_up")

# ---------- 7) 汇总统计 ----------
summary_dt <- data.table(
  Group  = c("High_up","Low_up"),
  n_all  = c(ret_high$n_all,  ret_low$n_all),
  n_p    = c(ret_high$n_p,    ret_low$n_p),
  n_FDR  = c(ret_high$n_f,    ret_low$n_f),
  n_both = c(ret_high$n_both, ret_low$n_both)
)
fwrite(summary_dt, file.path(out_dir, "KEGG_run_summary.csv"))

cat("\n✅ KEGG ORA（离线、无画图，自动映射 KEGG 名称）完成。输出目录：", out_dir, "\n", sep="")
cat("输出文件说明：\n")
cat("  *_KEGG_all.csv                   = 全量通路（已映射通路名称）\n")
cat("  *_KEGG_Plt0.05.csv               = pvalue < 0.05\n")
cat("  *_KEGG_FDRlt0.05.csv             = BH-FDR(p.adjust) < 0.05\n")
cat("  *_KEGG_Plt0.05_FDRlt0.05.csv     = 同时满足 pvalue<0.05 且 FDR<0.05\n")
cat("  SYMBOL_to_ENTREZ_*.csv           = SYMBOL->ENTREZ 映射表\n")
cat("  KEGG_run_summary.csv             = 两组富集条目数汇总\n")

# ---------- 8) 追加 SVG 输出 ----------
if (!requireNamespace("svglite", quietly = TRUE)) {
  message("未安装 svglite：SVG 输出将被跳过")
} else {
  svg_dir <- file.path(out_dir, "svg")
  dir.create(svg_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (tag in c("High_up","Low_up")) {
    for (ont in c("KEGG")) {
      f_all <- file.path(out_dir, sprintf("%s_KEGG_all.csv", tag))
      if (!file.exists(f_all)) next
      df <- fread(f_all)
      if (nrow(df) == 0) next
      
      df_p <- df[!is.na(pvalue) & pvalue < P_cut]
      if (nrow(df_p) == 0) next
      
      topN <- min(20, nrow(df_p))
      df_plot <- df_p[1:topN]
      df_plot[, Description2 := factor(Description, levels = rev(Description))]
      
      p_dot <- ggplot(df_plot, aes(x = Count, y = Description2)) +
        geom_point(aes(size = Count, alpha = -log10(pvalue)), color = ifelse(tag=="High_up","#B2182B","#1B9E77")) +
        theme_bw(base_size = 12) +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank(),
              legend.title = element_blank(),
              plot.title = element_blank()) +
        guides(alpha = "none") +
        labs(x="Count", y=NULL)
      
      ggsave(file.path(svg_dir, sprintf("%s_KEGG_dotplot_top%s.svg", tag, topN)),
             plot=p_dot, width=10, height=7, device=svglite::svglite)
    }
  }
  
  cat("\n✅ SVG 图已额外输出到：", svg_dir, "\n", sep="")
}

