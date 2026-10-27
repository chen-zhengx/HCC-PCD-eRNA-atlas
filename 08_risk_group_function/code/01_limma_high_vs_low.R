# =========================================================
# 对应结果：高/低风险组差异基因（limma，High vs Low，70 个 GSE144269 肿瘤样本）
#           代码阈值 FDR(BH) < 0.01 且 |log2FC| > 1，得到 高风险上调/低风险上调 基因数：
#           失巢凋亡 anoikis 92/111，凋亡 apoptosis 139/129，铁死亡 ferroptosis 105/142
#           是 Figure 4、Figure S4、Figure S5 GO/KEGG 富集（脚本 02–09）的上游
# 输入：../data/expression/gene_TPM/GSE144269_gene_TPM.csv.gz（19,416 个编码基因 TPM，140 列）
#       ../data/samples/GSE144269_samples.csv（用于取 70 个肿瘤样本列）
#       input/risk_group_<pcd>.csv（Sample, RiskGroup = High/Low；中位数风险评分分组，来自 05_prognostic_model）
# 输出：output/limma/<pcd>/High_vs_Low_limma_all.csv、High_vs_Low_limma_sig_FDR0.01_logFC1.csv
#       High_up_sig_FDR0.01_logFC1.csv、Low_up_sig_FDR0.01_logFC1.csv、High_vs_Low_limma_threshold_stats.csv
#       High_vs_Low_limma_results.xlsx
# 运行：在模块目录下 Rscript code/01_limma_high_vs_low.R
# =========================================================

# =========================================================
# High vs Low 全流程分析（TPM表达矩阵）
# 包含：
# 1) limma差异分析（过滤零方差基因）
# 2) 输出完整表 / 显著表 / 分方向表 / 阈值统计 / Excel
# 3) 输出 GSEA 排序文件 .rnk
# 4) fgsea 运行 Hallmark GSEA
# 5) 输出 High / Low top20 通路
# 6) 差异基因 ∩ leading-edge 交集分析
# =========================================================

suppressPackageStartupMessages({
  library(limma)
  library(openxlsx)
})

for (pcd in c("anoikis", "apoptosis", "ferroptosis")) {
cat("\n==========", pcd, "==========\n")

# =========================================================
# 0. 路径与参数
# =========================================================
limma_dir <- file.path("output", "limma", pcd)

dir.create(limma_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- "../data/expression/gene_TPM/GSE144269_gene_TPM.csv.gz"   # gene x sample, TPM
grp_file  <- file.path("input", paste0("risk_group_", pcd, ".csv"))     # Sample, RiskGroup(High/Low)

# 差异阈值
FDR_cut   <- 0.01
logFC_cut <- 1.0

# =========================================================
# 1. limma差异分析
# =========================================================
cat("========================================\n")
cat("Step 1/3: limma差异分析开始\n")
cat("========================================\n")

# ---------- 1.1 读取表达矩阵 ----------
expr <- read.csv(expr_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# 取 70 个肿瘤样本列（第 1 列为基因名）
sample_info <- read.csv("../data/samples/GSE144269_samples.csv", stringsAsFactors = FALSE)
expr <- expr[, c(colnames(expr)[1], sample_info$name[sample_info$Treatment == "tumor"])]
gene_col <- colnames(expr)[1]
rownames(expr) <- expr[[gene_col]]
expr[[gene_col]] <- NULL
expr <- as.matrix(expr)

# 转数值矩阵
mode(expr) <- "numeric"

# log2(TPM + 1)
E <- log2(expr + 1)

# ---------- 1.2 读取分组并对齐样本 ----------
grp <- read.csv(grp_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

need <- c("Sample", "RiskGroup")
miss <- setdiff(need, colnames(grp))
if (length(miss) > 0) stop("分组文件缺少列：", paste(miss, collapse = ", "))

grp$Sample <- trimws(grp$Sample)
grp$RiskGroup <- trimws(grp$RiskGroup)

# 仅保留表达矩阵中存在的样本
keep <- grp$Sample[grp$Sample %in% colnames(E)]
missing_samples <- setdiff(grp$Sample, colnames(E))
if (length(missing_samples) > 0) {
  message("注意：以下分组样本不在表达矩阵列名中，将跳过：\n",
          paste(missing_samples, collapse = ", "))
}
if (length(keep) < 2) stop("匹配到的样本太少，无法分析。请检查样本名一致性。")

grp <- grp[grp$Sample %in% keep, , drop = FALSE]
E   <- E[, grp$Sample, drop = FALSE]  # 严格按分组表顺序排列

grp$RiskGroup <- factor(grp$RiskGroup, levels = c("Low", "High"))
if (any(is.na(grp$RiskGroup))) stop("RiskGroup 存在非 High/Low 值，请检查分组文件。")

cat("样本数：", nrow(grp), "（Low=", sum(grp$RiskGroup=="Low"),
    ", High=", sum(grp$RiskGroup=="High"), "）\n", sep = "")

# ---------- 1.3 过滤零方差基因 ----------
gene_variance <- apply(E, 1, var, na.rm = TRUE)
E_filtered <- E[gene_variance > 0, , drop = FALSE]

cat("过滤零方差基因后剩余：", nrow(E_filtered), "个基因\n", sep = "")

if (nrow(E_filtered) < 10) stop("过滤后基因数过少，请检查表达矩阵。")

# ---------- 1.4 limma建模 ----------
design <- model.matrix(~ 0 + RiskGroup, data = grp)
colnames(design) <- levels(grp$RiskGroup)  # Low / High

fit <- lmFit(E_filtered, design)
fit <- contrasts.fit(fit, makeContrasts(High_vs_Low = High - Low, levels = design))
fit <- eBayes(fit)

res <- topTable(
  fit,
  coef = "High_vs_Low",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

res$Gene <- rownames(res)
res <- res[, c("Gene", setdiff(colnames(res), "Gene"))]

# ---------- 1.5 输出limma结果 ----------
f_all <- file.path(limma_dir, "High_vs_Low_limma_all.csv")
f_sig <- file.path(limma_dir, sprintf("High_vs_Low_limma_sig_FDR%s_logFC%s.csv", FDR_cut, logFC_cut))
f_hup <- file.path(limma_dir, sprintf("High_up_sig_FDR%s_logFC%s.csv", FDR_cut, logFC_cut))
f_lup <- file.path(limma_dir, sprintf("Low_up_sig_FDR%s_logFC%s.csv", FDR_cut, logFC_cut))
f_stats <- file.path(limma_dir, "High_vs_Low_limma_threshold_stats.csv")
f_xlsx <- file.path(limma_dir, "High_vs_Low_limma_results.xlsx")

write.csv(res, f_all, row.names = FALSE)

sig <- res[res$adj.P.Val < FDR_cut & abs(res$logFC) > logFC_cut, , drop = FALSE]
high_up <- sig[sig$logFC > 0, , drop = FALSE]  # High > Low
low_up  <- sig[sig$logFC < 0, , drop = FALSE]  # Low > High

write.csv(sig,     f_sig, row.names = FALSE)
write.csv(high_up, f_hup, row.names = FALSE)
write.csv(low_up,  f_lup, row.names = FALSE)

# 阈值统计
fc_thresholds <- c(0.5, 0.6, 0.8, 1.0)
stats <- do.call(rbind, lapply(fc_thresholds, function(fc){
  s <- res[res$adj.P.Val < FDR_cut & abs(res$logFC) > fc, , drop = FALSE]
  data.frame(
    Threshold = paste0("FDR<", FDR_cut, " & |logFC|>", fc),
    Total = nrow(s),
    High_up = sum(s$logFC > 0),
    Low_up  = sum(s$logFC < 0),
    stringsAsFactors = FALSE
  )
}))
write.csv(stats, f_stats, row.names = FALSE)

# Excel 汇总
xlsx_list <- list(
  All = res,
  Stats = stats
)
xlsx_list[[sprintf("Sig_FDR%s_logFC%s", FDR_cut, logFC_cut)]] <- sig
xlsx_list[[sprintf("High_up_FDR%s_logFC%s", FDR_cut, logFC_cut)]] <- high_up
xlsx_list[[sprintf("Low_up_FDR%s_logFC%s", FDR_cut, logFC_cut)]] <- low_up

write.xlsx(xlsx_list, file = f_xlsx, rowNames = FALSE)

cat("\n✅ limma差异分析完成\n")
cat("完整结果：", f_all, "\n", sep = "")
cat("显著结果：", f_sig, "\n", sep = "")
cat("高风险上调：", f_hup, "\n", sep = "")
cat("低风险上调：", f_lup, "\n", sep = "")
cat("阈值统计：", f_stats, "\n", sep = "")
cat("Excel汇总：", f_xlsx, "\n", sep = "")
}  # 循环结束
