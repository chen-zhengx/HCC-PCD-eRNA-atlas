# =========================================================
# 对应结果：20 个 PCD 基因集的超几何富集检验（筛出失巢凋亡、凋亡、铁死亡）
# 输入：input/pcd_genesets_TableS2.csv
#       input/common_DE_genes_649.csv（6 个数据集共有差异基因，= 02_differential_expression 的 common_DE_genes.csv）
#       input/background_coding_genes_19416.csv（背景：19,416 个编码基因）
# 输出：output/TableS3_hypergeometric_20categories.csv
# 运行：在模块目录下 Rscript code/01_hypergeometric_test.R
# =========================================================
dir.create("output", showWarnings = FALSE)
genesets <- read.csv("input/pcd_genesets_TableS2.csv", stringsAsFactors = FALSE)
de_genes <- read.csv("input/common_DE_genes_649.csv")$gene
background <- read.csv("input/background_coding_genes_19416.csv")$gene

res <- data.frame()
for (type in unique(genesets$Type)) {
  genes <- unique(genesets$Symbol[genesets$Type == type])
  a <- length(background)                 # 总编码基因 a = 19416
  b <- length(de_genes)                   # 差异表达基因 b = 649
  c <- length(genes)                      # 程序性死亡基因 c
  d <- length(intersect(genes, de_genes)) # 差异程序性死亡基因 d

  # 计算右侧p值（富集检验）
  p_value <- 1 - phyper(d-1, c, a-c, b)
  cat(type, "p值 =", format.pval(p_value, digits=3), "\n")
  res <- rbind(res, data.frame(Type = type, a = a, b = b, c = c, d = d, P = p_value))
}
# Methods：超几何 P 值在各 PCD 类别间做 Benjamini–Hochberg 校正，FDR<0.05 为显著富集
res$FDR <- p.adjust(res$P, method = "BH")
res <- res[order(tolower(res$Type)), ]
write.csv(res, "output/TableS3_hypergeometric_20categories.csv", row.names = FALSE)
print(res)
