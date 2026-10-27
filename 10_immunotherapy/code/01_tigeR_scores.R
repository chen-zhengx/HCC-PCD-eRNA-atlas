# =========================================================
# 对应结果：70 个肿瘤样本的 23 种免疫治疗反应标志物得分（tigeR::score_biomk；Figure 7、S6、S7 的上游
#           Methods "Immunotherapy-related molecular biomarkers"）
# 输入：../09_immune/input/GSE144269_tumor70_allgenes_TPM.csv
#       ../data/clinical/GSE144269_survival.csv
# 输出：output/tigeR/tigeR_23signatures_scores.csv
#       output/tigeR/tigeR_23signatures_NA_rate.csv、tigeR_23signatures_univariateCox.csv（单因素 Cox，论文未使用）
# 运行：在模块目录下 Rscript code/01_tigeR_scores.R
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(tigeR)
  library(survival)
})

# =========================
# 0) 路径
# =========================
input_dir  <- "input"
output_dir <- "output/tigeR"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

f_expr <- "../09_immune/input/GSE144269_tumor70_allgenes_TPM.csv"
f_surv <- "../data/clinical/GSE144269_survival.csv"

# =========================
# 1) 读表达矩阵（基因×样本）
# =========================
dt <- fread(f_expr, header = TRUE, data.table = TRUE, check.names = FALSE)
gene_col <- colnames(dt)[1]
genes <- dt[[gene_col]]

mat <- as.matrix(dt[, -1, with = FALSE])
rownames(mat) <- genes

# 强制数值化
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- genes

# 去掉重复基因（保留第一次出现）
dup <- duplicated(rownames(mat))
if (any(dup)) mat <- mat[!dup, , drop = FALSE]

# 如果你已是log2过的数据，把下一行注释掉
mat <- log2(mat + 1)

# =========================
# 2) 读生存信息
# =========================
surv <- fread(f_surv, header = TRUE, data.table = TRUE)
surv[, Sample := as.character(Sample)]
surv[, months := as.numeric(months)]
surv[, status := as.integer(status)]  # 1=Dead, 0=Alive

# =========================
# 3) 对齐样本：取交集并按表达矩阵列顺序排列
# =========================
expr_samples <- colnames(mat)
common <- intersect(expr_samples, surv$Sample)

if (length(common) < 5) {
  stop("表达矩阵与生存表匹配到的样本太少（<5）。请检查Sample命名是否一致。")
}

mat <- mat[, expr_samples[expr_samples %in% common], drop = FALSE]
surv <- surv[match(colnames(mat), Sample)]
stopifnot(all(colnames(mat) == surv$Sample))

# =========================
# 4) 构建 SummarizedExperiment（按 tigeR 期望字段名）
#    关键：Treatment 必须是 PRE（治疗前/未治疗用于生存分析）
# =========================
meta <- DataFrame(
  sample_id = surv$Sample,
  Treatment = rep("PRE", ncol(mat)),                     # ✅关键
  overall_survival_days = surv$months * 30.4375,         # months -> days（近似）
  vital_status = ifelse(surv$status == 1, "Dead", "Alive")
)
rownames(meta) <- meta$sample_id

se <- SummarizedExperiment(
  assays = list(expr = mat),
  colData = meta
)

# =========================
# 5) 计算 23 种免疫评分
# =========================
sig_res <- score_biomk(se)

# =========================
# 6) 输出：评分矩阵 + NA比例
# =========================
out_scores <- file.path(output_dir, "tigeR_23signatures_scores.csv")
fwrite(as.data.table(sig_res, keep.rownames = "Sample"), out_scores)

na_rate <- colMeans(is.na(sig_res))
out_na <- file.path(output_dir, "tigeR_23signatures_NA_rate.csv")
fwrite(data.table(Signature = names(na_rate), NA_rate = as.numeric(na_rate)), out_na)

# =========================
# 7) （可选）每个signature做单因素 Cox（基于你提供的生存信息）
# =========================
cox_list <- lapply(colnames(sig_res), function(sig) {
  x <- sig_res[, sig]
  time <- meta$overall_survival_days
  status <- ifelse(meta$vital_status == "Dead", 1L, 0L)
  
  ok <- is.finite(x) & !is.na(x) & is.finite(time) & !is.na(time) & !is.na(status)
  if (sum(ok) < 10) return(NULL)
  
  fit <- coxph(Surv(time[ok], status[ok]) ~ x[ok])
  s <- summary(fit)
  
  data.table(
    Signature = sig,
    HR = unname(s$coef[1, "exp(coef)"]),
    HR_lower95 = unname(s$conf.int[1, "lower .95"]),
    HR_upper95 = unname(s$conf.int[1, "upper .95"]),
    P = unname(s$coef[1, "Pr(>|z|)"])
  )
})

cox_dt <- rbindlist(cox_list, fill = TRUE)
if (nrow(cox_dt) > 0) {
  cox_dt[, FDR := p.adjust(P, method = "fdr")]
  out_cox <- file.path(output_dir, "tigeR_23signatures_univariateCox.csv")
  fwrite(cox_dt[order(P)], out_cox)
}

cat("✅ 完成！\n",
    "评分矩阵: ", out_scores, "\n",
    "NA比例: ", out_na, "\n",
    if (exists("out_cox")) paste0("Cox结果: ", out_cox, "\n") else "",
    sep = "")
# 看哪些signature NA比例最高
na <- fread(file.path(output_dir, "tigeR_23signatures_NA_rate.csv"))
na[order(-NA_rate)][1:10]

# 看评分矩阵维度是否合理
sc <- fread(file.path(output_dir, "tigeR_23signatures_scores.csv"))
dim(sc)

