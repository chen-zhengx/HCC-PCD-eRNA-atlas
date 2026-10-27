# =========================================================
# 对应结果：Table S11 及其扩展版：198 种药物的预测 IC50 中位数（IQR）、配对差值的 Hodges–Lehmann
#           估计（95% CI）、配对秩二列相关 r、肿瘤低于癌旁的配对数、单侧 P 和 BH FDR，并标出 19 种候选药物
# 输入：input/drug_predictions_GDSC2_140samples.csv
#       output/screening/AllDrugs_StatisticalScreening.csv
# 输出：output/tableS11/Table_S11_all_198_drugs.csv
# 运行：在模块目录下先运行 02，再 Rscript code/03_tableS11_all_drugs.R
# =========================================================

# =========================================================
# 药物敏感性预测：预测反应值、配对效应量和校正后 P 值
# 输入：oncoPredict（GDSC2）对 70 对肿瘤/癌旁样本的预测 IC50（DrugPredictions9.4整理.csv）
# 检验：单侧配对 Wilcoxon 符号秩检验（肿瘤 < 癌旁），198 种药物 BH 校正；
#   候选药物 = FDR<0.01 且平均差值（癌旁 − 肿瘤）> 0
# 新增：中位数（IQR）、配对差值的 Hodges–Lehmann 估计（95% CI）、配对秩二列相关 r、肿瘤低于癌旁的配对比例
# =========================================================
root <- "output/tableS11"
dir.create(root, recursive = TRUE, showWarnings = FALSE)
d <- read.csv("input/drug_predictions_GDSC2_140samples.csv", check.names = FALSE)
tum <- d[d$Treatment == "tumor", ]; nor <- d[d$Treatment == "nontumor", ]
tum <- tum[order(tum$PatientID), ]; nor <- nor[match(tum$PatientID, nor$PatientID), ]
stopifnot(identical(tum$PatientID, nor$PatientID), nrow(tum) == 70)
drugs <- colnames(d)[4:ncol(d)]

rank_biserial <- function(diff) {                     # 配对秩二列相关：(正差秩和 − 负差秩和) / 总秩和
  diff <- diff[diff != 0]; r <- rank(abs(diff))
  (sum(r[diff > 0]) - sum(r[diff < 0])) / sum(r)
}
res <- do.call(rbind, lapply(drugs, function(g) {
  t <- tum[[g]]; n <- nor[[g]]; dif <- t - n
  w1 <- wilcox.test(t, n, paired = TRUE, alternative = "less")                 # 筛选所用的检验
  w2 <- suppressWarnings(wilcox.test(t, n, paired = TRUE, conf.int = TRUE))   # Hodges–Lehmann 估计与 95% CI
  data.frame(drug_raw = g,
             tumor_median = median(t), tumor_q1 = quantile(t, .25), tumor_q3 = quantile(t, .75),
             adjacent_median = median(n), adjacent_q1 = quantile(n, .25), adjacent_q3 = quantile(n, .75),
             mean_difference_adjacent_minus_tumor = mean(n) - mean(t),
             HL_difference_tumor_minus_adjacent = unname(w2$estimate),
             HL_CI_low = w2$conf.int[1], HL_CI_high = w2$conf.int[2],
             rank_biserial_r = rank_biserial(dif),
             pairs_lower_in_tumor = sum(dif < 0), P_one_sided = w1$p.value)
}))
res$FDR <- p.adjust(res$P_one_sided, "BH")
res$candidate <- res$FDR < 0.01 & res$mean_difference_adjacent_minus_tumor > 0
res$GDSC2_ID <- sub(".*_", "", res$drug_raw)
nm <- sub("_[0-9]+$", "", res$drug_raw)
nm[nm == "Sepantronium.bromide"] <- "Sepantronium bromide"; nm[nm == "GDC0810"] <- "GDC-0810"
res$drug <- gsub("\\.", "-", nm)
res <- res[order(res$FDR, -res$mean_difference_adjacent_minus_tumor), ]
rownames(res) <- NULL

# ---------- 与 02 的筛选结果核对 ----------
orig <- read.csv("output/screening/AllDrugs_StatisticalScreening.csv")
m <- match(orig$Drug, make.names(res$drug_raw)); m[is.na(m)] <- match(orig$Drug[is.na(m)], res$drug_raw)
cat("matched drugs:", sum(!is.na(m)), "/ 198\n")
cat("max |P diff|:", max(abs(orig$Wilcoxon_Pvalue - res$P_one_sided[m])), "| max |FDR diff|:",
    max(abs(orig$Wilcoxon_FDR - res$FDR[m])), "\n")
cat("candidates:", sum(res$candidate), "| identical to 02 screening:",
    setequal(make.names(res$drug_raw[res$candidate]), orig$Drug[orig$Wilcoxon_FDR < 0.01 & orig$Mean_Difference > 0]), "\n")

out <- res[, c("drug", "GDSC2_ID", "tumor_median", "tumor_q1", "tumor_q3", "adjacent_median", "adjacent_q1",
               "adjacent_q3", "HL_difference_tumor_minus_adjacent", "HL_CI_low", "HL_CI_high", "rank_biserial_r",
               "pairs_lower_in_tumor", "P_one_sided", "FDR", "candidate")]
write.csv(out, file.path(root, "Table_S11_all_198_drugs.csv"), row.names = FALSE)
options(width = 220)
print(head(out[out$candidate, ], 19), digits = 3, row.names = FALSE)
cat("\nnon-candidate drugs with FDR<0.01 but mean difference <=0:", sum(res$FDR < 0.01 & !res$candidate), "\n")
