# =========================================================
# 对应结果：Table S4（临床特征及单因素 Cox：table_clinical_characteristics.csv）
#           Table S9（风险分组、TNM 分期、肝硬化的单因素/多因素 Cox：cox_results.csv，及表注中的敏感性分析、PH 检验）
#           混杂因素筛选表（confounder_screening.csv）
# 输入：output/clinical_analysis_data.csv（由 code/01_build_clinical_data.py 生成）
# 输出：output/table_clinical_characteristics.csv、confounder_screening.csv、cox_results.csv
#       cox_PH_test.csv、sessionInfo.txt
# 运行：在模块目录下 Rscript code/02_cox_analysis.R
# =========================================================

# =========================================================
# 临床特征，以及风险分组在校正临床混杂因素后与总生存的关系（Cox 回归）
# 数据：clinical_analysis_data.csv（由 01_build_clinical_data.py 生成；70 例，23 例死亡）
# 风险分组：论文 Figure 3B/S1B/S2B 的高/低风险组（开发队列风险评分中位数切分）
# 混杂因素：TNM 分期（III–IV vs I–II）与肝硬化 —— 本队列中与总生存显著相关、HCC 公认预后因素、且在两组间不平衡
#   主分析：风险分组 + TNM 分期 + 肝硬化（3 个参数，EPV = 23/3 ≈ 7.7），缺失值多重插补（m=20，Rubin 规则合并）
#   敏感性分析：(1) 完整病例；(2) Firth 惩罚 Cox，额外加入年龄、性别、AFP（6 个参数）
#   比例风险假设：Schoenfeld 残差检验（cox.zph）
# =========================================================
suppressPackageStartupMessages({ library(survival); library(mice); library(coxphf) })
root <- "output"
d <- read.csv(file.path(root, "clinical_analysis_data.csv"))
d$age10 <- d$age / 10; d$bmi5 <- d$bmi / 5
models <- c(anoikis = "Anoikis", apoptosis = "Apoptosis", ferroptosis = "Ferroptosis")
for (m in names(models)) d[[paste0("high_", m)]] <- as.integer(as.logical(d[[paste0("app_high_", m)]]))
vars <- c(age10 = "Age (per 10 years)", male = "Male sex", stage_III_IV = "TNM stage III-IV (vs. I-II)",
          tumor_large = "Tumor size >5 cm", multinodular = "Multinodular tumor", cirrhosis = "Cirrhosis",
          afp_abn = "Abnormal AFP", alb_abn = "Abnormal albumin", bil_abn = "Abnormal total bilirubin",
          alt_abn = "Abnormal ALT", hbv = "HBV infection", hcv = "HCV infection", hdv = "HDV infection",
          bmi5 = "BMI (per 5 kg/m2)", alcohol = "Alcohol consumption", smoker = "Ever smoker",
          fhx_lc = "Family history of liver cancer")
fmt_p  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001",
                      ifelse(p < 0.01 | (p >= 0.04 & p < 0.06), sprintf("%.3f", p), sprintf("%.2f", p))))
fmt_hr <- function(h, l, u) sprintf("%.2f (%.2f-%.2f)", h, l, u)
cox_row <- function(fit, term) {
  s <- summary(fit)$conf.int; p <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  list(HR = s[term, "exp(coef)"], lo = s[term, "lower .95"], hi = s[term, "upper .95"], P = p)
}

# ---------- 1. 临床特征 + 单因素 Cox（Table 2）----------
char <- list(); uni <- list()
add <- function(...) char[[length(char) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
pct <- function(k) sprintf("%d (%.1f)", k, 100 * k / nrow(d))
add(Variable = "Patients", Category = "", `n (%)` = "70 (100)", Deaths = sum(d$status), `HR (95% CI)` = "", P = "", check.names = FALSE)
for (v in names(vars)) {
  x <- d[[v]]; cc <- !is.na(x)
  f <- coxph(Surv(months, status) ~ x, data = data.frame(months = d$months[cc], status = d$status[cc], x = x[cc]))
  r <- cox_row(f, "x"); uni[[v]] <- r
  hr <- fmt_hr(r$HR, r$lo, r$hi); p <- fmt_p(r$P)
  if (v %in% c("age10", "bmi5")) {
    raw <- if (v == "age10") d$age else d$bmi
    q <- quantile(raw, c(.5, .25, .75), na.rm = TRUE)
    add(Variable = vars[v], Category = "Median (IQR)", `n (%)` = sprintf("%.1f (%.1f-%.1f)", q[1], q[2], q[3]),
        Deaths = "", `HR (95% CI)` = hr, P = p, check.names = FALSE)
  } else if (v == "stage_III_IV") {
    for (s in 1:4) add(Variable = "TNM stage", Category = as.character(as.roman(s)), `n (%)` = pct(sum(d$tnm == s, na.rm = TRUE)),
                       Deaths = sum(d$status[d$tnm %in% s]), `HR (95% CI)` = "", P = "", check.names = FALSE)
    add(Variable = vars[v], Category = "III-IV vs. I-II", `n (%)` = "", Deaths = "", `HR (95% CI)` = hr, P = p, check.names = FALSE)
  } else {
    add(Variable = vars[v], Category = "Yes", `n (%)` = pct(sum(x == 1, na.rm = TRUE)), Deaths = sum(d$status[x %in% 1]),
        `HR (95% CI)` = hr, P = p, check.names = FALSE)
    add(Variable = vars[v], Category = "No", `n (%)` = pct(sum(x == 0, na.rm = TRUE)), Deaths = sum(d$status[x %in% 0]),
        `HR (95% CI)` = "Reference", P = "", check.names = FALSE)
  }
  if (any(!cc)) add(Variable = vars[v], Category = "Missing", `n (%)` = pct(sum(!cc)), Deaths = sum(d$status[!cc]),
                    `HR (95% CI)` = "", P = "", check.names = FALSE)
}
write.csv(do.call(rbind, char), file.path(root, "table_clinical_characteristics.csv"), row.names = FALSE)

# ---------- 2. 混杂因素筛选：与生存的关系、缺失、两组间不平衡（SMD）----------
smd <- function(x, g) { ok <- !is.na(x); a <- x[ok & g == 1]; b <- x[ok & g == 0]
                        (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2) }
decision <- c(stage_III_IV = "Included", cirrhosis = "Included",
              age10 = "Not included: not associated with survival and balanced between risk groups",
              male = "Not included: not associated with survival and balanced between risk groups",
              afp_abn = "Not included in the main model: not associated with survival; evaluated in a sensitivity analysis",
              tumor_large = "Not included: component of TNM stage (T category); not associated with survival",
              multinodular = "Not included: component of TNM stage (T category)",
              alb_abn = "Not included: liver function overlaps with cirrhosis; not associated with survival",
              bil_abn = "Not included: liver function overlaps with cirrhosis; not associated with survival",
              alt_abn = "Not included: not associated with survival",
              hbv = "Not included: etiological factor not associated with survival",
              hcv = "Not included: etiological factor not associated with survival",
              hdv = "Not included: etiological factor not associated with survival",
              bmi5 = "Not included: risk factor for HCC development; not associated with survival",
              alcohol = "Not included: risk factor for HCC development; not associated with survival",
              smoker = "Not included: risk factor for HCC development; not associated with survival",
              fhx_lc = "Not included: risk factor for HCC development; not associated with survival")
scr <- do.call(rbind, lapply(names(vars), function(v) data.frame(
  variable = vars[v], missing_n = sum(is.na(d[[v]])), univariable_HR = fmt_hr(uni[[v]]$HR, uni[[v]]$lo, uni[[v]]$hi),
  univariable_P = fmt_p(uni[[v]]$P),
  SMD_anoikis = round(smd(d[[v]], d$high_anoikis), 2), SMD_apoptosis = round(smd(d[[v]], d$high_apoptosis), 2),
  SMD_ferroptosis = round(smd(d[[v]], d$high_ferroptosis), 2), decision = decision[v])))
write.csv(scr, file.path(root, "confounder_screening.csv"), row.names = FALSE)

# ---------- 3. Cox 回归 ----------
conf <- c("stage_III_IV", "cirrhosis"); extra <- c("age10", "male", "afp_abn")
aux  <- c("multinodular", "tumor_large", "alb_abn", "bil_abn", "hbv", "hcv", "hdv")
M <- 20
pool <- function(B, V) {                                   # Rubin 规则：B、V 为 M × p 的系数与方差
  q <- colMeans(B); u <- colMeans(V); b <- apply(B, 2, var); t <- u + (1 + 1 / M) * b
  dfr <- (M - 1) * (1 + u / ((1 + 1 / M) * b))^2
  data.frame(term = colnames(B), HR = exp(q), lo = exp(q - qt(.975, dfr) * sqrt(t)), hi = exp(q + qt(.975, dfr) * sqrt(t)),
             P = 2 * pt(-abs(q / sqrt(t)), dfr), row.names = NULL)
}
res <- list(); ph <- list()
for (m in names(models)) {
  g <- paste0("high_", m)
  fml <- function(cov) as.formula(paste("Surv(months, status) ~", paste(c(g, cov), collapse = " + ")))
  # 单因素
  for (v in c(g, conf)) {
    cc <- d[!is.na(d[[v]]), ]
    r <- cox_row(coxph(as.formula(paste("Surv(months, status) ~", v)), data = cc), v)
    res[[length(res) + 1]] <- data.frame(model = models[m], analysis = "Univariable", term = v, n = nrow(cc),
                                         events = sum(cc$status), HR = r$HR, lo = r$lo, hi = r$hi, P = r$P)
  }
  # 多重插补
  dd <- d[, c("months", "status", g, conf, extra, aux)]; dd$na_haz <- nelsonaalen(dd, months, status)
  for (b in c("male", "stage_III_IV", "cirrhosis", "afp_abn", "multinodular", "tumor_large", "alb_abn", "bil_abn")) dd[[b]] <- factor(dd[[b]])
  meth <- make.method(dd); meth[c("months", "status", g)] <- ""
  pred <- make.predictorMatrix(dd); pred[, "months"] <- 0
  imp <- mice(dd, m = M, method = meth, predictorMatrix = pred, printFlag = FALSE, seed = 2026)
  comp <- lapply(1:M, function(i) { di <- complete(imp, i)
    for (b in c(conf, extra)) di[[b]] <- as.numeric(as.character(di[[b]])); di })
  # 主分析：分组 + TNM + 肝硬化（标准 Cox，多重插补）
  fits <- lapply(comp, function(di) { f <- coxph(fml(conf), data = di, x = TRUE, model = TRUE)
    list(coef = coef(f), var = diag(vcov(f)), zph = cox.zph(f)$table) })      # 比例风险检验须在数据所在的环境中计算
  pm <- pool(t(sapply(fits, `[[`, "coef")), t(sapply(fits, `[[`, "var")))
  res[[length(res) + 1]] <- data.frame(model = models[m], analysis = "Multivariable (multiple imputation)", pm[, "term", drop = FALSE],
                                       n = 70, events = 23, pm[, c("HR", "lo", "hi", "P")])
  ph[[m]] <- data.frame(model = models[m], data = "multiple imputation (median of 20)",
                        global_P = median(sapply(fits, function(f) f$zph["GLOBAL", "p"])),
                        group_P = median(sapply(fits, function(f) f$zph[g, "p"])))
  # 敏感性 1：完整病例
  cc <- d[complete.cases(d[, conf]), ]; fc <- coxph(fml(conf), data = cc, x = TRUE, model = TRUE)
  s <- summary(fc)
  res[[length(res) + 1]] <- data.frame(model = models[m], analysis = "Multivariable (complete cases)", term = rownames(s$conf.int),
                                       n = nrow(cc), events = sum(cc$status), HR = s$conf.int[, 1], lo = s$conf.int[, 3],
                                       hi = s$conf.int[, 4], P = s$coefficients[, 5], row.names = NULL)
  z <- cox.zph(fc)$table
  ph[[paste0(m, "_cc")]] <- data.frame(model = models[m], data = "complete cases", global_P = z["GLOBAL", "p"], group_P = z[g, "p"])
  # 敏感性 2：Firth 惩罚 Cox，额外加入年龄、性别、AFP（多重插补）
  ff <- lapply(comp, function(di) coxphf(fml(c(conf, extra)), data = di, pl = FALSE))
  pf <- pool(t(sapply(ff, function(f) f$coefficients)), t(sapply(ff, function(f) diag(f$var))))
  res[[length(res) + 1]] <- data.frame(model = models[m], analysis = "Firth, additionally adjusted for age, sex and AFP (multiple imputation)",
                                       pf[, "term", drop = FALSE], n = 70, events = 23, pf[, c("HR", "lo", "hi", "P")])
}
res <- do.call(rbind, res); rownames(res) <- NULL
res$term <- sub("^high_.*", "High-risk group", res$term)
res$HR_CI <- fmt_hr(res$HR, res$lo, res$hi); res$P_fmt <- fmt_p(res$P)
write.csv(res, file.path(root, "cox_results.csv"), row.names = FALSE)
ph <- do.call(rbind, ph); rownames(ph) <- NULL
write.csv(ph, file.path(root, "cox_PH_test.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(root, "sessionInfo.txt"))
options(width = 200)
print(scr[, c("variable", "missing_n", "univariable_HR", "univariable_P", "SMD_anoikis", "SMD_apoptosis", "SMD_ferroptosis")], row.names = FALSE)
print(res[, c("model", "analysis", "term", "n", "events", "HR_CI", "P_fmt")], row.names = FALSE)
print(ph, digits = 3, row.names = FALSE)
