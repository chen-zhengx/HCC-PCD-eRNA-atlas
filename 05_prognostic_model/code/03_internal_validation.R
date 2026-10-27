# =========================================================
# 对应结果：Table S7（三个模型的风险 eRNA 与 LASSO 系数：apparent_model_coefficients.csv）
#           Table S8（apparent 与重复 90:10 划分内部验证的 C-index 和 1/2/3 年 AUC：summary_performance.csv）
#           正文中的队列描述（70 例、23 例死亡、中位随访 39.3 个月、λ 值：summary_cohort.csv）
#           以及 06_clinical 使用的每例风险分数和高/低风险组（patient_scores_apparent.csv）
# 输入：input/<module>_eRNA_RPM_tumor70.csv、../data/clinical/GSE144269_survival.csv
# 输出：output/<module>/ 下的 summary_cohort.csv、apparent_model_coefficients.csv、patient_scores_apparent.csv
#       summary_performance.csv、summary_risk_groups.csv、split90_per_repeat.csv、sessionInfo.txt
# 运行：在模块目录下依次运行（每次一个模型）
#       Rscript code/03_internal_validation.R anoikis
#       Rscript code/03_internal_validation.R apoptosis
#       Rscript code/03_internal_validation.R ferroptosis
# =========================================================

# =========================================================
# PCD-eRNA 预后模型的内部验证（随机 90:10 划分，重复 10 次；特征选择包含在每次划分的训练集内）
# 用法：Rscript code/03_internal_validation.R <module>
#   module = anoikis / apoptosis / ferroptosis
# 建模流程与 02_lasso_km_roc.R 完全一致：
#   log2(x+1) -> 1%/99% winsorize -> 单因素 Cox P<0.2 ->
#   LASSO-Cox（固定 lambda 网格、10 折、standardize=TRUE、lambda.1se，无特征时退回 lambda.min）->
#   风险分数 = 线性预测值；以中位数分高低风险组
# LASSO 在 lambda.1se 和 lambda.min 都未选中任何 eRNA 时：取单因素 Cox P 值最小的 5 个 eRNA，以单因素 log(HR) 为权重
# 内部验证：按生存状态分层，随机取 90% 为训练集、10% 为测试集，重复 10 次；
#   winsorize 分位数、单因素筛选、LASSO、系数和中位数阈值全部只用训练集估计，再用于测试集
# 输出（<module>/）：
#   summary_cohort.csv            队列描述、全样本模型的筛选结果、90:10 训练集的 LASSO 兜底比例
#   apparent_model_coefficients.csv  全样本模型（全部 70 例）的 eRNA 和系数
#   patient_scores_apparent.csv   每个病人的全样本模型风险分数和高低风险组
#   summary_performance.csv       C 指数和 1/2/3 年 AUC：apparent（95% CI）与 90:10 测试集（均值、最小值、最大值）
#   summary_risk_groups.csv       全样本模型高低风险组的 HR 和 log-rank P
#   split90_per_repeat.csv        10 次划分的逐次结果
# =========================================================

suppressPackageStartupMessages({
  library(survival); library(glmnet); library(timeROC)
})

args   <- commandArgs(trailingOnly = TRUE)
module <- if (length(args) >= 1) args[1] else "anoikis"
n_rep  <- 10
p_train <- 0.9
times  <- c(12, 24, 36)   # 月

base  <- "."
files <- list(
  anoikis     = c("input/anoikis_eRNA_RPM_tumor70.csv",
                  "../data/clinical/GSE144269_survival.csv"),
  apoptosis   = c("input/apoptosis_eRNA_RPM_tumor70.csv",
                  "../data/clinical/GSE144269_survival.csv"),
  ferroptosis = c("input/ferroptosis_eRNA_RPM_tumor70.csv",
                  "../data/clinical/GSE144269_survival.csv")
)
out_dir <- file.path("output", module)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush.console() }

# ---------- 数据 ----------
expr <- read.csv(file.path(base, files[[module]][1]), row.names = 1,
                 check.names = FALSE, fileEncoding = "UTF-8-BOM")
surv <- read.csv(file.path(base, files[[module]][2]), fileEncoding = "UTF-8-BOM")
X    <- log2(t(as.matrix(expr)) + 1)
surv <- surv[match(rownames(X), surv$Sample), ]
stopifnot(identical(as.character(surv$Sample), rownames(X)))
time <- surv$months; status <- surv$status; n <- nrow(X)

# ---------- 建模流程 ----------
lambda_seq <- exp(seq(log(0.001), log(1), length.out = 100))
winsor <- function(M, lo, hi) sweep(sweep(M, 2, lo, pmax), 2, hi, pmin)

fit_pipeline <- function(Xtr, ttr, str, seed = 123, p_uni = 0.2) {
  lo <- apply(Xtr, 2, quantile, 0.01, na.rm = TRUE)
  hi <- apply(Xtr, 2, quantile, 0.99, na.rm = TRUE)
  Xw <- winsor(Xtr, lo, hi)
  y  <- Surv(ttr, str)
  uni <- t(apply(Xw, 2, function(x) {
    f <- tryCatch(suppressWarnings(coxph(y ~ x)), error = function(e) NULL)
    if (is.null(f)) c(1, 0) else summary(f)$coefficients[1, c(5, 1)]
  }))
  pv <- uni[, 1]; uni_beta <- uni[, 2]
  pv[is.na(pv)] <- 1; uni_beta[is.na(uni_beta)] <- 0
  keep <- names(which(pv < p_uni))
  if (length(keep) == 0) keep <- colnames(Xw)
  if (length(keep) < 2) keep <- names(sort(pv))[1:2]   # glmnet 至少需要 2 列
  set.seed(seed)
  cv <- cv.glmnet(Xw[, keep], y, family = "cox", alpha = 1, nfolds = 10,
                  standardize = TRUE, lambda = lambda_seq)
  b <- coef(cv, s = "lambda.1se"); s_use <- "lambda.1se"
  if (sum(b != 0) == 0) { b <- coef(cv, s = "lambda.min"); s_use <- "lambda.min" }
  nz <- which(as.vector(b) != 0)
  if (length(nz) == 0) {                       # LASSO 未选中任何 eRNA 时的兜底
    sel <- names(sort(pv))[1:min(5, length(pv))]
    cf  <- unname(uni_beta[sel])
    s_use <- "top5_signed"
  } else {
    sel <- rownames(b)[nz]; cf <- as.vector(b)[nz]
  }
  lp_tr <- as.vector(Xw[, sel, drop = FALSE] %*% cf)
  list(lo = lo, hi = hi, sel = sel, coef = cf, cutoff = median(lp_tr), s_use = s_use,
       n_uni = length(keep), lambda = if (s_use == "lambda.1se") cv$lambda.1se else cv$lambda.min)
}
predict_lp <- function(m, Xnew) {
  Xw <- winsor(Xnew[, names(m$lo), drop = FALSE], m$lo, m$hi)
  as.vector(Xw[, m$sel, drop = FALSE] %*% m$coef)
}
# 带重试的拟合（个别训练集中 glmnet 可能报错）
safe_fit <- function(..., seed) {
  for (tr in 0:5) {
    m <- tryCatch(fit_pipeline(..., seed = seed + tr * 7919), error = function(e) NULL)
    if (!is.null(m)) return(m)
  }
  NULL
}

# ---------- 性能指标 ----------
c_index <- function(lp, t, s) {
  cc <- concordance(Surv(t, s) ~ lp, reverse = TRUE)
  c(C = unname(cc$concordance), se = unname(sqrt(cc$var)))
}
td_auc <- function(lp, t, s) {
  r <- tryCatch(timeROC(T = t, delta = s, marker = lp, cause = 1, times = times, iid = FALSE),
                error = function(e) NULL)
  if (is.null(r)) return(setNames(rep(NA, length(times)), paste0("AUC", times)))
  setNames(unname(r$AUC)[seq_along(times)], paste0("AUC", times))
}
group_stats <- function(g, t, s) {
  f  <- coxph(Surv(t, s) ~ factor(g, levels = c(FALSE, TRUE)))
  ci <- exp(confint(f)); lr <- survdiff(Surv(t, s) ~ g)
  c(HR = unname(exp(coef(f))), HR_low = ci[1, 1], HR_high = ci[1, 2],
    logrank_P = 1 - pchisq(lr$chisq, df = 1))
}

# ---------- 0. 队列描述 ----------
fu <- survfit(Surv(time, 1 - status) ~ 1)                # reverse KM 计算中位随访
fu_q <- quantile(fu, probs = c(0.25, 0.5, 0.75))$quantile
cohort <- data.frame(module = module, N = n, deaths = sum(status), censored = sum(status == 0),
                     median_followup_months = unname(fu_q[2]),
                     followup_IQR_low = unname(fu_q[1]), followup_IQR_high = unname(fu_q[3]),
                     max_time_months = max(time), n_candidate_eRNAs = ncol(X),
                     EPV_candidates = sum(status) / ncol(X))
log_msg(module, ": N =", n, "deaths =", sum(status), "candidates =", ncol(X))

# ---------- 1. Apparent（开发队列的性能）----------
m0  <- fit_pipeline(X, time, status, seed = 123)
lp0 <- predict_lp(m0, X)
app_c <- c_index(lp0, time, status)
roc0 <- timeROC(T = time, delta = status, marker = lp0, cause = 1, times = times, iid = TRUE)
roc0_ci <- confint(roc0)$CI_AUC / 100
app_grp <- group_stats(lp0 > m0$cutoff, time, status)
cohort$n_univariate_pass <- m0$n_uni
cohort$n_final_eRNAs <- length(m0$sel)
cohort$EPV_final_model <- sum(status) / length(m0$sel)
cohort$lambda_rule <- m0$s_use; cohort$lambda <- m0$lambda
log_msg("Apparent model:", length(m0$sel), "eRNAs; C =", round(app_c["C"], 3))
write.csv(data.frame(eRNA = m0$sel, coefficient = m0$coef),
          file.path(out_dir, "apparent_model_coefficients.csv"), row.names = FALSE)
write.csv(data.frame(Sample = rownames(X), months = time, status = status,
                     apparent_lp = lp0, apparent_high = lp0 > m0$cutoff),
          file.path(out_dir, "patient_scores_apparent.csv"), row.names = FALSE)

# ---------- 2. 随机 90:10 划分，重复 10 次 ----------
set.seed(2026)
res <- list()
for (r in seq_len(n_rep)) {
  train <- logical(n)
  for (g in 0:1) { idx <- which(status == g); train[sample(idx, round(length(idx) * p_train))] <- TRUE }
  m  <- safe_fit(X[train, ], time[train], status[train], seed = r * 100)
  lp <- predict_lp(m, X[!train, , drop = FALSE]); tt <- time[!train]; ss <- status[!train]
  cc <- c_index(lp, tt, ss)
  grp <- lp > m$cutoff
  gs <- tryCatch(suppressWarnings(group_stats(grp, tt, ss)),
                 error = function(e) c(HR = NA, HR_low = NA, HR_high = NA, logrank_P = NA))
  res[[r]] <- data.frame(
    repeat_id = r, n_train = sum(train), deaths_train = sum(status[train]),
    n_test = sum(!train), deaths_test = sum(ss),
    test_deaths_by_12m = sum(ss == 1 & tt <= 12), test_deaths_by_24m = sum(ss == 1 & tt <= 24),
    test_deaths_by_36m = sum(ss == 1 & tt <= 36),
    n_eRNA = length(m$sel), lambda_rule = m$s_use, eRNAs = paste(m$sel, collapse = ";"),
    C_test = cc[["C"]], C_low = cc[["C"]] - 1.96 * cc[["se"]], C_high = cc[["C"]] + 1.96 * cc[["se"]],
    t(td_auc(lp, tt, ss)),
    n_high_test = sum(grp), deaths_high_test = sum(ss[grp]),
    HR_test = gs[["HR"]], logrank_P = gs[["logrank_P"]], row.names = NULL)
}
res <- do.call(rbind, res)
# 测试集低风险组没有死亡时，Cox 的 HR 会发散，只作记录，不报告
res$HR_estimable <- res$deaths_high_test > 0 & res$deaths_high_test < res$deaths_test
write.csv(res, file.path(out_dir, "split90_per_repeat.csv"), row.names = FALSE)
cohort$split_n_repeats <- n_rep
cohort$split_lasso_fallback <- sum(grepl("^top5", res$lambda_rule))
cohort$split_median_n_eRNA <- median(res$n_eRNA)

# ---------- 3. 汇总 ----------
sm <- function(x) c(mean = mean(x, na.rm = TRUE), min = min(x, na.rm = TRUE),
                    max = max(x, na.rm = TRUE), n_evaluable = sum(!is.na(x)))
split_cols <- c("C_test", paste0("AUC", times))
perf <- data.frame(
  module = module, metric = c("C-index", paste0(times / 12, "-year AUC")),
  apparent = c(app_c[["C"]], roc0$AUC),
  apparent_low = c(app_c[["C"]] - 1.96 * app_c[["se"]], roc0_ci[, 1]),
  apparent_high = c(app_c[["C"]] + 1.96 * app_c[["se"]], roc0_ci[, 2]),
  t(sapply(split_cols, function(k) sm(res[[k]]))), n_repeats = n_rep, row.names = NULL)
names(perf)[6:9] <- c("split90_mean", "split90_min", "split90_max", "split90_n_evaluable")
grp_tab <- data.frame(module = module, grouping = "Apparent (development cohort)",
                      n_high = sum(lp0 > m0$cutoff), deaths_high = sum(status[lp0 > m0$cutoff]), t(app_grp))
write.csv(perf,    file.path(out_dir, "summary_performance.csv"), row.names = FALSE)
write.csv(grp_tab, file.path(out_dir, "summary_risk_groups.csv"), row.names = FALSE)
write.csv(cohort,  file.path(out_dir, "summary_cohort.csv"), row.names = FALSE)

print(cohort); print(perf, digits = 3); print(grp_tab, digits = 3)
print(res[, c("repeat_id", "deaths_test", "n_eRNA", "lambda_rule", "C_test", "AUC12", "AUC24", "AUC36", "logrank_P")],
      digits = 3, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
log_msg("done:", module)
