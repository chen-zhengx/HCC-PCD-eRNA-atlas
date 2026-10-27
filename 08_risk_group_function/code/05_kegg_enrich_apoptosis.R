# =========================================================
# 对应结果：Figure S4C/S4D 的 KEGG 富集计算（凋亡 apoptosis；高/低风险组上调基因的 enrichKEGG，无作图）
# 输入：output/limma/apoptosis/High_vs_Low_limma_all.csv（背景基因）、High_up_sig_FDR0.01_logFC1.csv、Low_up_sig_FDR0.01_logFC1.csv（01 的输出）
# 输出：output/figureS4/KEGG_ORA/High_up_KEGG_*.csv、Low_up_KEGG_*.csv（all / Plt0.05 / FDRlt0.05 / Plt0.05_FDRlt0.05）
#       SYMBOL_to_ENTREZ_*.csv、KEGG_run_summary.csv
# 运行：在模块目录下 Rscript code/05_kegg_enrich_apoptosis.R（需先运行 01）
# =========================================================
# =========================================================
# KEGG ORA (clusterProfiler::enrichKEGG) - OFFLINE NO-PLOT
# 凋亡：High_up vs Low_up
# 修正版：自动把 KEGG ID（如 hsa04110）映射为通路名称
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(KEGG.db)
})

# ---------- 0) 路径 ----------
limma_all_file <- "output/limma/apoptosis/High_vs_Low_limma_all.csv"
high_up_file   <- "output/limma/apoptosis/High_up_sig_FDR0.01_logFC1.csv"
low_up_file    <- "output/limma/apoptosis/Low_up_sig_FDR0.01_logFC1.csv"

out_dir <- "output/figureS4/KEGG_ORA"
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
  
  if (!("pvalue" %in% names(dt))) {
    stop(sprintf("[%s] missing pvalue. cols=%s", tag, paste(names(dt), collapse = ", ")))
  }
  if (!("p.adjust" %in% names(dt))) {
    stop(sprintf("[%s] missing p.adjust. cols=%s", tag, paste(names(dt), collapse = ", ")))
  }
  
  dt
}

# KEGG ID -> 通路名
# 例如 hsa04110 -> Cell cycle
kegg_id_to_name <- function(id_vec){
  id_vec <- as.character(id_vec)
  id2name <- as.list(KEGGPATHID2NAME)
  
  out <- rep(NA_character_, length(id_vec))
  
  for (i in seq_along(id_vec)) {
    x <- id_vec[i]
    x2 <- gsub("^hsa", "", x)  # hsa04110 -> 04110
    
    if (!is.null(id2name[[x2]])) {
      out[i] <- id2name[[x2]]
    } else {
      out[i] <- x
    }
  }
  
  out
}

sym2ent <- function(sym_vec){
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

# ---------- 3) 读入基因 ----------
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

# ---------- 4) SYMBOL -> ENTREZ ----------
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
    return(list(n_all = 0, n_p = 0, n_f = 0, n_both = 0))
  }
  
  df <- normalize_cols(df, paste0(tag, "-KEGG"))
  
  # ---------- 修正 Description ----------
  # 如果 Description 缺失，或其实还是 hsa04110 这种编号，则自动映射通路名称
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
  
  # 保存修正后的全量结果
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
  
  list(
    n_all  = nrow(df),
    n_p    = nrow(df_p),
    n_f    = nrow(df_f),
    n_both = nrow(df_both)
  )
}

# ---------- 6) 分别跑 High / Low ----------
ret_high <- run_kegg_offline_noplot(m_high$entrez, m_bg$entrez, tag = "High_up")
ret_low  <- run_kegg_offline_noplot(m_low$entrez,  m_bg$entrez, tag = "Low_up")

# ---------- 7) 汇总统计 ----------
summary_dt <- data.table(
  Group  = c("High_up", "Low_up"),
  n_all  = c(ret_high$n_all,  ret_low$n_all),
  n_p    = c(ret_high$n_p,    ret_low$n_p),
  n_FDR  = c(ret_high$n_f,    ret_low$n_f),
  n_both = c(ret_high$n_both, ret_low$n_both)
)
fwrite(summary_dt, file.path(out_dir, "KEGG_run_summary.csv"))

cat("\n✅ KEGG ORA（离线、无画图）完成。输出目录：", out_dir, "\n", sep = "")
cat("输出文件说明：\n")
cat("  *_KEGG_all.csv                   = 全量通路（已尽量修正 Description 为真正通路名）\n")
cat("  *_KEGG_Plt0.05.csv               = pvalue < 0.05\n")
cat("  *_KEGG_FDRlt0.05.csv             = BH-FDR(p.adjust) < 0.05\n")
cat("  *_KEGG_Plt0.05_FDRlt0.05.csv     = 同时满足 pvalue<0.05 且 FDR<0.05\n")
cat("  SYMBOL_to_ENTREZ_*.csv           = SYMBOL->ENTREZ 映射表\n")
cat("  KEGG_run_summary.csv             = 两组富集条目数汇总\n")

