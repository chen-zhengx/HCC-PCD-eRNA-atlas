# =========================================================
# 对应结果：Table S10（28 种免疫细胞 ssGSEA 得分在 Molecular_Cluster2 vs Cluster1 间的比较：秩二列相关 r、P、BH FDR）
# 输入：input/ssgsea_scores_28cells.csv
#       input/subtype_<pcd>.csv
# 输出：output/tableS10/TableS10_immune_tests_FDR.csv
# 运行：在模块目录下 Rscript code/07_tableS10_immune_tests.R
# =========================================================

# =========================================================
# 28 种免疫细胞 ssGSEA 富集得分在分子亚型间比较的多重检验校正（Figure 6A）
# 每种分子亚型分类（失巢凋亡、凋亡、铁死亡）内：Molecular_Cluster2 vs Cluster1，双侧 Wilcoxon 秩和检验
# （与 tigeR 标志物的比较相同：wilcox.test, exact = FALSE），再对 28 个检验做 BH 校正
# 效应量：秩二列相关 r = 2U / (n1 n2) − 1，正值表示 Cluster2 更高
# tigeR 标志物（Figure 7、S6、S7）不做校正，沿用原始 P 值
# =========================================================
out  <- "output/tableS10"
dir.create(out, recursive = TRUE, showWarnings = FALSE)  # 新建输出目录
mods <- c(Anoikis = "失巢死亡", Apoptosis = "凋亡", Ferroptosis = "铁死亡")

ss <- read.csv("input/ssgsea_scores_28cells.csv",
               row.names = 1, check.names = FALSE)                      # 28 细胞 × 70 肿瘤
ss <- as.data.frame(t(ss))

test_set <- function(X, grp, model) {
  X <- X[names(grp), , drop = FALSE]
  r <- do.call(rbind, lapply(colnames(X), function(v) {
    x <- as.numeric(X[[v]]); ok <- is.finite(x)
    a <- x[ok & grp == 2]; b <- x[ok & grp == 1]
    w <- wilcox.test(a, b, exact = FALSE)                               # W 即 Cluster2 的 U 统计量
    data.frame(model = model, feature = v, n1 = length(b), n2 = length(a),
               median_C1 = median(b), median_C2 = median(a),
               r_rb = unname(2 * w$statistic / (length(a) * length(b)) - 1), P = w$p.value)
  }))
  r$FDR <- p.adjust(r$P, "BH")
  r
}

res <- list()
for (m in names(mods)) {
  g <- read.csv(file.path("input", paste0("subtype_", tolower(m), ".csv")))
  res[[m]] <- test_set(ss, setNames(as.integer(g[[2]]), g[[1]]), m)
}
res <- do.call(rbind, res); rownames(res) <- NULL
res$higher <- ifelse(res$r_rb > 0, "Cluster2", "Cluster1")
write.csv(res, file.path(out, "TableS10_immune_tests_FDR.csv"), row.names = FALSE)

# 汇总：未校正 P<0.05 与校正后 FDR<0.05 的数目
options(width = 200)
print(aggregate(cbind(P_lt_0.05 = P < 0.05, FDR_lt_0.05 = FDR < 0.05) ~ model, data = res, FUN = sum), row.names = FALSE)
show <- res[res$P < 0.05, c("model", "feature", "median_C1", "median_C2", "r_rb", "P", "FDR", "higher")]
show <- show[order(show$model, show$P), ]
print(transform(show, median_C1 = signif(median_C1, 3), median_C2 = signif(median_C2, 3), r_rb = round(r_rb, 2),
                P = signif(P, 3), FDR = signif(FDR, 3)), row.names = FALSE)
