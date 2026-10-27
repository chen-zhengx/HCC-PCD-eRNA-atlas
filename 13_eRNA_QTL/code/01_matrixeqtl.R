# =========================================================
# 对应结果：eRNA-QTL 分析——4 个 cis-eRNA-QTL（FDR < 0.05，肿瘤基因型主分析）、157 个 cis / 109,943 个 trans 检验
#           以及癌旁组织基因型的敏感性分析（Results "Exploratory eRNA-QTL analysis"、Figure 9、Supplementary Table 13）
# 输入：input/genotype_tumor/（或 genotype_normal/）中的 erna_expression.csv（12 个风险 eRNA 在 70 个肿瘤样本中的 RPM
#       前 3 列为染色体、起点、终点）和 genotype_matrix.csv（SNP × 70 个肿瘤样本，ALT 等位基因剂量 0/1/2，NA 为缺失）
# 输出：output/tumor/、output/normal/：QTL_results_cis.csv、QTL_results_trans.csv
#       cis_QTLs.csv、trans_QTLs.csv、FDR_0.05/cis_QTLs_FDR_0.05.csv、MatrixQTL_results.RData 等
# 运行：在模块目录下
#   Rscript code/01_matrixeqtl.R input/genotype_tumor/erna_expression.csv input/genotype_tumor/genotype_matrix.csv output/tumor
#   Rscript code/01_matrixeqtl.R input/genotype_normal/erna_expression.csv input/genotype_normal/genotype_matrix.csv output/normal
# =========================================================
# Original MatrixEQTL script (eRNA-QTL分析/qtl代码.txt). Only the input/output paths
# are now command-line arguments, MatrixEQTL is loaded from tools/Rlib and the
# unused ggplot2 import is commented out.
# Usage: Rscript code/01_matrixeqtl.R ERNA_CSV GENOTYPE_CSV OUTPUT_DIR
# 加载必要的包
library(MatrixEQTL)
library(data.table)

# 设定新的文件路径
args <- commandArgs(trailingOnly = TRUE)
eRNA_file <- args[1]
genotype_file <- args[2]
output_dir <- args[3]

# 创建主输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 创建FDR阈值子目录
fdr_0.05_dir <- file.path(output_dir, "FDR_0.05")
fdr_0.01_dir <- file.path(output_dir, "FDR_0.01")
dir.create(fdr_0.05_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fdr_0.01_dir, recursive = TRUE, showWarnings = FALSE)

cat("创建FDR阈值文件夹:\n")
cat("FDR < 0.05:", fdr_0.05_dir, "\n")
cat("FDR < 0.01:", fdr_0.01_dir, "\n")

# 1. 读取和检查原始数据
cat("=== 数据读取和检查 ===\n")

# 读取eRNA表达矩阵
cat("读取eRNA表达矩阵...\n")
eRNA_data <- read.csv(eRNA_file, header = TRUE, row.names = 1)
cat("eRNA数据维度:", dim(eRNA_data), "\n")
cat("前几列名称:", colnames(eRNA_data)[1:10], "\n")

# 提取eRNA位置信息和表达矩阵
# 根据你的描述，前三列可能是染色体、起始位置、终止位置
# 第四列开始是表达矩阵
eRNA_names <- rownames(eRNA_data)

# 检查eRNA数据的结构
cat("\neRNA数据结构:\n")
cat("总列数:", ncol(eRNA_data), "\n")
cat("前三列名称:", colnames(eRNA_data)[1:3], "\n")

# 假设前三列是染色体、起始、终止位置
eRNA_chr <- eRNA_data[, 1]
eRNA_start <- eRNA_data[, 2]
eRNA_end <- eRNA_data[, 3]

# 提取表达矩阵（从第4列开始）
expression_matrix <- as.matrix(eRNA_data[, 4:ncol(eRNA_data)])

cat("表达矩阵维度:", dim(expression_matrix), "\n")
cat("表达矩阵样本数:", ncol(expression_matrix), "\n")

# 创建基因位置信息
gene_pos <- data.frame(
  geneid = eRNA_names,
  chr = eRNA_chr,
  left = eRNA_start,
  right = eRNA_end,
  stringsAsFactors = FALSE
)

# 读取基因型矩阵
cat("\n读取基因型矩阵...\n")
genotype_data <- fread(genotype_file, header = TRUE)
cat("基因型数据维度:", dim(genotype_data), "\n")
cat("基因型数据列名:", colnames(genotype_data), "\n")

# 检查基因型数据结构
cat("基因型数据前几行:\n")
print(head(genotype_data[, 1:5]))

# 2. 检查SNP ID格式
cat("\n=== 检查SNP ID格式 ===\n")
# 新文件已经有CHROM, POS, ID列，可以直接使用
cat("SNP ID示例:\n")
print(head(genotype_data$ID))

# 3. 样本对齐
cat("\n=== 样本对齐 ===\n")
eRNA_samples <- colnames(expression_matrix)
genotype_samples <- colnames(genotype_data)[4:ncol(genotype_data)]  # 从第4列开始是样本

cat("eRNA样本数:", length(eRNA_samples), "\n")
cat("基因型样本数:", length(genotype_samples), "\n")

# 获取公共样本
common_samples <- intersect(eRNA_samples, genotype_samples)
cat("公共样本数量:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  cat("警告：没有公共样本！请检查样本名称是否匹配。\n")
  cat("eRNA样本示例:", head(eRNA_samples, 5), "\n")
  cat("基因型样本示例:", head(genotype_samples, 5), "\n")
} else {
  # 筛选公共样本
  expression_matrix <- expression_matrix[, common_samples]

  # 保留基因型数据的前三列（CHROM, POS, ID）和公共样本
  genotype_cols <- c("CHROM", "POS", "ID", common_samples)
  genotype_data <- genotype_data[, ..genotype_cols]

  cat("筛选后eRNA数据维度:", dim(expression_matrix), "\n")
  cat("筛选后基因型数据维度:", dim(genotype_data), "\n")
}

# 4. 处理缺失值 - 不删除任何SNP
cat("\n=== 处理缺失值（保留所有SNP） ===\n")

# 移除有缺失值的eRNA行
original_eRNA_count <- nrow(expression_matrix)
expression_matrix <- expression_matrix[complete.cases(expression_matrix), ]
gene_pos <- gene_pos[gene_pos$geneid %in% rownames(expression_matrix), ]

cat("原始eRNA数量:", original_eRNA_count, "\n")
cat("去除缺失值后eRNA数量:", nrow(expression_matrix), "\n")

# 5. 转换基因型数据 - 保留缺失值
cat("\n=== 转换基因型数据（保留缺失值） ===\n")

# 提取样本列（从第4列开始）
genotype_matrix <- as.matrix(genotype_data[, 4:ncol(genotype_data), with = FALSE])

# 检查基因型编码格式
cat("基因型编码示例:\n")
print(table(as.vector(genotype_matrix[1:5, 1:5])))

# 转换基因型为0,1,2编码，保留缺失值为NA
# 假设基因型编码格式可能是0/0, 0/1, 1/1等
genotype_numeric <- matrix(NA, nrow = nrow(genotype_matrix), ncol = ncol(genotype_matrix))

for (i in 1:nrow(genotype_matrix)) {
  for (j in 1:ncol(genotype_matrix)) {
    gt <- as.character(genotype_matrix[i, j])
    if (is.na(gt) || gt == "" || gt == "./." || gt == ".") {
      genotype_numeric[i, j] <- NA
    } else if (gt == "0/0" || gt == "0|0") {
      genotype_numeric[i, j] <- 0
    } else if (gt == "1/1" || gt == "1|1") {
      genotype_numeric[i, j] <- 2
    } else if (gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0") {
      genotype_numeric[i, j] <- 1
    } else {
      # 尝试直接转换为数值
      genotype_numeric[i, j] <- as.numeric(gt)
    }
  }
}

# 检查转换后的基因型
cat("转换后基因型值分布:\n")
print(table(genotype_numeric, useNA = "always"))

# 添加SNP信息作为行名
rownames(genotype_numeric) <- genotype_data$ID

cat("最终保留的SNP数量:", nrow(genotype_numeric), "\n")

# 6. 转换表达数据
cat("\n=== 转换表达数据 ===\n")
# 检查表达数据是否已经是log转换过的
expression_summary <- summary(as.vector(expression_matrix))
cat("表达数据统计摘要:\n")
print(expression_summary)

# 如果表达数据值很大，可能需要log转换
if (max(expression_matrix, na.rm = TRUE) > 100) {
  cat("执行log2转换...\n")
  expression_matrix <- log2(expression_matrix + 1)
} else {
  cat("表达数据可能已经过标准化，跳过log转换。\n")
}

# 7. 准备SNP位置信息
cat("\n=== 准备SNP位置信息 ===\n")
snp_pos <- data.frame(
  snpid = genotype_data$ID,
  chr = paste0("chr", genotype_data$CHROM),  # 添加chr前缀以匹配基因位置格式
  pos = genotype_data$POS,
  stringsAsFactors = FALSE
)

cat("SNP位置信息示例:\n")
print(head(snp_pos))

# 8. 创建SlicedData对象
cat("\n=== 创建SlicedData对象 ===\n")
snps = SlicedData$new(genotype_numeric)
gene = SlicedData$new(expression_matrix)

# 9. 运行QTL分析
cat("\n=== 运行QTL分析 ===\n")

output_file_name = file.path(output_dir, "QTL_results.csv")
output_file_name.cis = file.path(output_dir, "QTL_results_cis.csv")
output_file_name.trans = file.path(output_dir, "QTL_results_trans.csv")

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = SlicedData$new(),
  output_file_name = output_file_name.trans,
  pvOutputThreshold = 0.05,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = output_file_name.cis,
  pvOutputThreshold.cis = 0.05,
  snpspos = snp_pos,
  genepos = gene_pos,
  cisDist = 1e6,
  pvalue.hist = FALSE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)

# 10. 保存完整结果
cat("\n=== 保存完整结果 ===\n")
save(me, file = file.path(output_dir, "MatrixQTL_results.RData"))

# 保存所有显著结果
if (!is.null(me$all$eqtls) && nrow(me$all$eqtls) > 0) {
  write.csv(me$all$eqtls, file.path(output_dir, "all_significant_QTLs.csv"), row.names = FALSE)
  cat("已保存所有显著QTL到 all_significant_QTLs.csv\n")
} else {
  cat("没有发现达到阈值的显著QTL\n")
}

# 分别保存cis和trans结果
if (!is.null(me$cis$eqtls) && nrow(me$cis$eqtls) > 0) {
  write.csv(me$cis$eqtls, file.path(output_dir, "cis_QTLs.csv"), row.names = FALSE)
  cat("已保存", nrow(me$cis$eqtls), "个cis-QTL到 cis_QTLs.csv\n")
} else {
  cat("没有发现cis-QTL\n")
}

if (!is.null(me$trans$eqtls) && nrow(me$trans$eqtls) > 0) {
  write.csv(me$trans$eqtls, file.path(output_dir, "trans_QTLs.csv"), row.names = FALSE)
  cat("已保存", nrow(me$trans$eqtls), "个trans-QTL到 trans_QTLs.csv\n")
} else {
  cat("没有发现trans-QTL\n")
}

# 11. 按FDR阈值筛选和保存结果
cat("\n=== 按FDR阈值筛选结果 ===\n")

filter_by_fdr <- function(eqtl_data, fdr_threshold) {
  if (is.null(eqtl_data) || nrow(eqtl_data) == 0) {
    return(NULL)
  }
  return(eqtl_data[eqtl_data$FDR < fdr_threshold, ])
}

# FDR < 0.05 的结果
cat("筛选FDR < 0.05的结果...\n")
cis_fdr_0.05 <- filter_by_fdr(me$cis$eqtls, 0.05)
trans_fdr_0.05 <- filter_by_fdr(me$trans$eqtls, 0.05)
all_fdr_0.05 <- filter_by_fdr(me$all$eqtls, 0.05)

# 保存FDR < 0.05 的结果
if (!is.null(cis_fdr_0.05) && nrow(cis_fdr_0.05) > 0) {
  write.csv(cis_fdr_0.05, file.path(fdr_0.05_dir, "cis_QTLs_FDR_0.05.csv"), row.names = FALSE)
  cat("已保存", nrow(cis_fdr_0.05), "个cis-QTL (FDR<0.05)\n")
}

if (!is.null(trans_fdr_0.05) && nrow(trans_fdr_0.05) > 0) {
  write.csv(trans_fdr_0.05, file.path(fdr_0.05_dir, "trans_QTLs_FDR_0.05.csv"), row.names = FALSE)
  cat("已保存", nrow(trans_fdr_0.05), "个trans-QTL (FDR<0.05)\n")
}

if (!is.null(all_fdr_0.05) && nrow(all_fdr_0.05) > 0) {
  write.csv(all_fdr_0.05, file.path(fdr_0.05_dir, "all_QTLs_FDR_0.05.csv"), row.names = FALSE)
  cat("已保存", nrow(all_fdr_0.05), "个所有QTL (FDR<0.05)\n")
}

# FDR < 0.01 的结果
cat("筛选FDR < 0.01的结果...\n")
cis_fdr_0.01 <- filter_by_fdr(me$cis$eqtls, 0.01)
trans_fdr_0.01 <- filter_by_fdr(me$trans$eqtls, 0.01)
all_fdr_0.01 <- filter_by_fdr(me$all$eqtls, 0.01)

# 保存FDR < 0.01 的结果
if (!is.null(cis_fdr_0.01) && nrow(cis_fdr_0.01) > 0) {
  write.csv(cis_fdr_0.01, file.path(fdr_0.01_dir, "cis_QTLs_FDR_0.01.csv"), row.names = FALSE)
  cat("已保存", nrow(cis_fdr_0.01), "个cis-QTL (FDR<0.01)\n")
}

if (!is.null(trans_fdr_0.01) && nrow(trans_fdr_0.01) > 0) {
  write.csv(trans_fdr_0.01, file.path(fdr_0.01_dir, "trans_QTLs_FDR_0.01.csv"), row.names = FALSE)
  cat("已保存", nrow(trans_fdr_0.01), "个trans-QTL (FDR<0.01)\n")
}

if (!is.null(all_fdr_0.01) && nrow(all_fdr_0.01) > 0) {
  write.csv(all_fdr_0.01, file.path(fdr_0.01_dir, "all_QTLs_FDR_0.01.csv"), row.names = FALSE)
  cat("已保存", nrow(all_fdr_0.01), "个所有QTL (FDR<0.01)\n")
}

# 12. 生成分析报告
cat("\n=== QTL分析报告 ===\n")
cat("🎉 分析成功完成！\n")
cat("📊 总共测试了", me$all$ntests, "个SNP-基因对\n")
cat("🧬 使用的SNP数量:", nrow(genotype_numeric), "\n")
cat("🔬 使用的eRNA数量:", nrow(expression_matrix), "\n")
cat("👥 样本数量:", ncol(expression_matrix), "\n")

if (!is.null(me$all$eqtls) && nrow(me$all$eqtls) > 0) {
  cat("🔍 发现", nrow(me$all$eqtls), "个显著的QTL (p < 0.05)\n")
} else {
  cat("🔍 发现 0 个显著的QTL (p < 0.05)\n")
}

if (!is.null(me$cis$eqtls) && nrow(me$cis$eqtls) > 0) {
  cat("📍 其中cis-QTL:", nrow(me$cis$eqtls), "个\n")
} else {
  cat("📍 其中cis-QTL: 0 个\n")
}

if (!is.null(me$trans$eqtls) && nrow(me$trans$eqtls) > 0) {
  cat("🌐 其中trans-QTL:", nrow(me$trans$eqtls), "个\n")
} else {
  cat("🌐 其中trans-QTL: 0 个\n")
}

# FDR结果汇总
cat("\n=== FDR校正结果 ===\n")
cat("FDR < 0.05:\n")
cat("  - cis-QTL:", ifelse(!is.null(cis_fdr_0.05), nrow(cis_fdr_0.05), 0), "\n")
cat("  - trans-QTL:", ifelse(!is.null(trans_fdr_0.05), nrow(trans_fdr_0.05), 0), "\n")

cat("FDR < 0.01:\n")
cat("  - cis-QTL:", ifelse(!is.null(cis_fdr_0.01), nrow(cis_fdr_0.01), 0), "\n")
cat("  - trans-QTL:", ifelse(!is.null(trans_fdr_0.01), nrow(trans_fdr_0.01), 0), "\n")

# 13. 显示最显著的关联
cat("\n=== 最显著的关联 ===\n")
if (!is.null(me$all$eqtls) && nrow(me$all$eqtls) > 0) {
  cat("最显著的5个关联:\n")
  print(head(me$all$eqtls[order(me$all$eqtls$pvalue), c("snps", "gene", "pvalue", "FDR", "beta")], 5))
} else if (!is.null(me$trans$eqtls) && nrow(me$trans$eqtls) > 0) {
  cat("最显著的5个trans关联:\n")
  print(head(me$trans$eqtls[order(me$trans$eqtls$pvalue), c("snps", "gene", "pvalue", "FDR", "beta")], 5))
} else if (!is.null(me$cis$eqtls) && nrow(me$cis$eqtls) > 0) {
  cat("最显著的5个cis关联:\n")
  print(head(me$cis$eqtls[order(me$cis$eqtls$pvalue), c("snps", "gene", "pvalue", "FDR", "beta")], 5))
}

cat("\n💾 结果文件保存在:", output_dir, "\n")
cat("✨ 所有分析步骤完成！\n")

# 14. 额外检查：数据质量报告
cat("\n=== 数据质量报告 ===\n")

# 检查基因型缺失率
missing_rate <- sum(is.na(genotype_numeric)) / (nrow(genotype_numeric) * ncol(genotype_numeric))
cat("基因型缺失率:", round(missing_rate * 100, 2), "%\n")

# 检查表达数据分布
cat("表达数据统计:\n")
print(summary(as.vector(expression_matrix)))

# 检查染色体分布
cat("\nSNP染色体分布:\n")
print(table(genotype_data$CHROM))

cat("\neRNA染色体分布:\n")
print(table(gene_pos$chr))

# 保存运行参数
run_info <- data.frame(
  参数 = c("eRNA文件", "基因型文件", "输出目录", "cis距离", "p值阈值", "分析时间"),
  值 = c(eRNA_file, genotype_file, output_dir, "1e6", "0.05", Sys.time())
)

write.csv(run_info, file.path(output_dir, "运行参数.csv"), row.names = FALSE)
cat("\n运行参数已保存到 运行参数.csv\n")





####第二部分
# 14. 额外检查：数据质量报告
cat("\n=== 数据质量报告 ===\n")

# 检查基因型缺失率
missing_rate <- sum(is.na(genotype_numeric)) / (nrow(genotype_numeric) * ncol(genotype_numeric))
cat("基因型缺失率:", round(missing_rate * 100, 2), "%\n")

# 检查表达数据分布
cat("表达数据统计:\n")
print(summary(as.vector(expression_matrix)))

# 检查染色体分布
cat("\nSNP染色体分布:\n")
print(table(genotype_data$CHROM))

cat("\neRNA染色体分布:\n")
print(table(gene_pos$chr))

# 保存运行参数
run_info <- data.frame(
  参数 = c("eRNA文件", "基因型文件", "输出目录", "cis距离", "p值阈值", "分析时间"),
  值 = c(eRNA_file, genotype_file, output_dir, "1e6", "0.05", Sys.time())
)

write.csv(run_info, file.path(output_dir, "运行参数.csv"), row.names = FALSE)
cat("\n运行参数已保存到 运行参数.csv\n")






####第二部分
# 15. 按p值阈值筛选和保存结果（新添加部分）
cat("\n=== 按p值阈值筛选和保存结果 ===\n")

# 创建p值阈值子目录
p_0.05_dir <- file.path(output_dir, "p_value_0.05")
p_0.01_dir <- file.path(output_dir, "p_value_0.01")
dir.create(p_0.05_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(p_0.01_dir, recursive = TRUE, showWarnings = FALSE)

cat("创建p值阈值文件夹:\n")
cat("p < 0.05:", p_0.05_dir, "\n")
cat("p < 0.01:", p_0.01_dir, "\n")

# 定义按p值筛选的函数
filter_by_pvalue <- function(eqtl_data, p_threshold) {
  if (is.null(eqtl_data) || nrow(eqtl_data) == 0) {
    return(NULL)
  }
  return(eqtl_data[eqtl_data$pvalue < p_threshold, ])
}

# p < 0.05 的结果
cat("筛选p < 0.05的结果...\n")
cis_p_0.05 <- filter_by_pvalue(me$cis$eqtls, 0.05)
trans_p_0.05 <- filter_by_pvalue(me$trans$eqtls, 0.05)
all_p_0.05 <- filter_by_pvalue(me$all$eqtls, 0.05)

# 保存p < 0.05 的结果
if (!is.null(cis_p_0.05) && nrow(cis_p_0.05) > 0) {
  write.csv(cis_p_0.05, file.path(p_0.05_dir, "cis_QTLs_p_0.05.csv"), row.names = FALSE)
  cat("已保存", nrow(cis_p_0.05), "个cis-QTL (p<0.05)\n")
} else {
  cat("没有发现cis-QTL (p<0.05)\n")
}

if (!is.null(trans_p_0.05) && nrow(trans_p_0.05) > 0) {
  write.csv(trans_p_05, file.path(p_0.05_dir, "trans_QTLs_p_0.05.csv"), row.names = FALSE)
  cat("已保存", nrow(trans_p_0.05), "个trans-QTL (p<0.05)\n")
} else {
  cat("没有发现trans-QTL (p<0.05)\n")
}

if (!is.null(all_p_0.05) && nrow(all_p_0.05) > 0) {
  write.csv(all_p_0.05, file.path(p_0.05_dir, "all_QTLs_p_0.05.csv"), row.names = FALSE)
  cat("已保存", nrow(all_p_0.05), "个所有QTL (p<0.05)\n")
} else {
  cat("没有发现任何QTL (p<0.05)\n")
}

# p < 0.01 的结果
cat("筛选p < 0.01的结果...\n")
cis_p_0.01 <- filter_by_pvalue(me$cis$eqtls, 0.01)
trans_p_0.01 <- filter_by_pvalue(me$trans$eqtls, 0.01)
all_p_0.01 <- filter_by_pvalue(me$all$eqtls, 0.01)

# 保存p < 0.01 的结果
if (!is.null(cis_p_0.01) && nrow(cis_p_0.01) > 0) {
  write.csv(cis_p_0.01, file.path(p_0.01_dir, "cis_QTLs_p_0.01.csv"), row.names = FALSE)
  cat("已保存", nrow(cis_p_0.01), "个cis-QTL (p<0.01)\n")
} else {
  cat("没有发现cis-QTL (p<0.01)\n")
}

if (!is.null(trans_p_0.01) && nrow(trans_p_0.01) > 0) {
  write.csv(trans_p_0.01, file.path(p_0.01_dir, "trans_QTLs_p_0.01.csv"), row.names = FALSE)
  cat("已保存", nrow(trans_p_0.01), "个trans-QTL (p<0.01)\n")
} else {
  cat("没有发现trans-QTL (p<0.01)\n")
}

if (!is.null(all_p_0.01) && nrow(all_p_0.01) > 0) {
  write.csv(all_p_0.01, file.path(p_0.01_dir, "all_QTLs_p_0.01.csv"), row.names = FALSE)
  cat("已保存", nrow(all_p_0.01), "个所有QTL (p<0.01)\n")
} else {
  cat("没有发现任何QTL (p<0.01)\n")
}

# 16. 更新分析报告，包含p值结果汇总
cat("\n=== 更新分析报告（包含p值结果） ===\n")

cat("📊 总共测试了", me$all$ntests, "个SNP-基因对\n")
cat("🧬 使用的SNP数量:", nrow(genotype_numeric), "\n")
cat("🔬 使用的eRNA数量:", nrow(expression_matrix), "\n")
cat("👥 样本数量:", ncol(expression_matrix), "\n")

# p值结果汇总
cat("\n=== p值阈值结果汇总 ===\n")
cat("p < 0.05:\n")
cat("  - cis-QTL:", ifelse(!is.null(cis_p_0.05), nrow(cis_p_0.05), 0), "个\n")
cat("  - trans-QTL:", ifelse(!is.null(trans_p_0.05), nrow(trans_p_0.05), 0), "个\n")
cat("  - 总计:", ifelse(!is.null(all_p_0.05), nrow(all_p_0.05), 0), "个\n")

cat("p < 0.01:\n")
cat("  - cis-QTL:", ifelse(!is.null(cis_p_0.01), nrow(cis_p_0.01), 0), "个\n")
cat("  - trans-QTL:", ifelse(!is.null(trans_p_0.01), nrow(trans_p_0.01), 0), "个\n")
cat("  - 总计:", ifelse(!is.null(all_p_0.01), nrow(all_p_0.01), 0), "个\n")

# FDR结果汇总（已存在，再次显示以便对比）
cat("\n=== FDR阈值结果汇总 ===\n")
cat("FDR < 0.05:\n")
cat("  - cis-QTL:", ifelse(!is.null(cis_fdr_0.05), nrow(cis_fdr_0.05), 0), "个\n")
cat("  - trans-QTL:", ifelse(!is.null(trans_fdr_0.05), nrow(trans_fdr_0.05), 0), "个\n")
cat("  - 总计:", ifelse(!is.null(all_fdr_0.05), nrow(all_fdr_0.05), 0), "个\n")

cat("FDR < 0.01:\n")
cat("  - cis-QTL:", ifelse(!is.null(cis_fdr_0.01), nrow(cis_fdr_0.01), 0), "个\n")
cat("  - trans-QTL:", ifelse(!is.null(trans_fdr_0.01), nrow(trans_fdr_0.01), 0), "个\n")
cat("  - 总计:", ifelse(!is.null(all_fdr_0.01), nrow(all_fdr_0.01), 0), "个\n")

# 17. 创建汇总统计报告
cat("\n=== 创建汇总统计报告 ===\n")

summary_stats <- data.frame(
  类别 = c("总测试对数", "SNP数量", "eRNA数量", "样本数量",
         "p<0.05总数", "p<0.01总数",
         "cis-QTL总数", "trans-QTL总数",
         "FDR<0.05总数", "FDR<0.01总数",
         "基因型缺失率"),
  数值 = c(
    me$all$ntests,
    nrow(genotype_numeric),
    nrow(expression_matrix),
    ncol(expression_matrix),
    ifelse(!is.null(all_p_0.05), nrow(all_p_0.05), 0),
    ifelse(!is.null(all_p_0.01), nrow(all_p_0.01), 0),
    ifelse(!is.null(me$cis$eqtls), nrow(me$cis$eqtls), 0),
    ifelse(!is.null(me$trans$eqtls), nrow(me$trans$eqtls), 0),
    ifelse(!is.null(all_fdr_0.05), nrow(all_fdr_0.05), 0),
    ifelse(!is.null(all_fdr_0.01), nrow(all_fdr_0.01), 0),
    paste0(round(missing_rate * 100, 2), "%")
  )
)

write.csv(summary_stats, file.path(output_dir, "分析汇总统计.csv"), row.names = FALSE)
cat("汇总统计已保存到 分析汇总统计.csv\n")

# 18. 生成目录结构报告
cat("\n=== 输出文件目录结构 ===\n")
output_files <- list.files(output_dir, recursive = TRUE, full.names = FALSE)
cat("生成的文件:\n")
for (file in output_files) {
  cat("  ", file, "\n")
}

cat("\n🎉 所有分析步骤完成！结果文件已保存在:", output_dir, "\n")
cat("✨ 分析已完成，请查看上述目录中的结果文件。\n")


###第三部分
##############
# 14. 额外检查：数据质量报告
cat("\n=== 数据质量报告 ===\n")

# 检查基因型缺失率
missing_rate <- sum(is.na(genotype_numeric)) / (nrow(genotype_numeric) * ncol(genotype_numeric))
cat("基因型缺失率:", round(missing_rate * 100, 2), "%\n")

# 检查表达数据分布
cat("表达数据统计:\n")
print(summary(as.vector(expression_matrix)))

# 检查染色体分布
cat("\nSNP染色体分布:\n")
print(table(genotype_data$CHROM))

cat("\neRNA染色体分布:\n")
print(table(gene_pos$chr))

# 保存运行参数
run_info <- data.frame(
  参数 = c("eRNA文件", "基因型文件", "输出目录", "cis距离", "p值阈值", "分析时间"),
  值 = c(eRNA_file, genotype_file, output_dir, "1e6", "0.05", Sys.time())
)

write.csv(run_info, file.path(output_dir, "运行参数.csv"), row.names = FALSE)
cat("\n运行参数已保存到 运行参数.csv\n")

# 15. 按p值阈值筛选和保存结果
cat("\n=== 按p值阈值筛选和保存结果 ===\n")

# 创建p值阈值子目录
p_0.05_dir <- file.path(output_dir, "p_value_0.05")
p_0.01_dir <- file.path(output_dir, "p_value_0.01")
dir.create(p_0.05_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(p_0.01_dir, recursive = TRUE, showWarnings = FALSE)

cat("创建p值阈值文件夹:\n")
cat("p < 0.05:", p_0.05_dir, "\n")
cat("p < 0.01:", p_0.01_dir, "\n")

# 定义按p值筛选的函数
filter_by_pvalue <- function(eqtl_data, p_threshold) {
  if (is.null(eqtl_data) || nrow(eqtl_data) == 0) {
    return(NULL)
  }
  return(eqtl_data[eqtl_data$pvalue < p_threshold, ])
}

# p < 0.05 的结果
cat("筛选p < 0.05的结果...\n")
cis_p_0.05 <- filter_by_pvalue(me$cis$eqtls, 0.05)
trans_p_0.05 <- filter_by_pvalue(me$trans$eqtls, 0.05)
all_p_0.05 <- filter_by_pvalue(me$all$eqtls, 0.05)

# 保存p < 0.05 的结果
if (!is.null(cis_p_0.05) && nrow(cis_p_0.05) > 0) {
  write.csv(cis_p_0.05, file.path(p_0.05_dir, "cis_QTLs_p_0.05.csv"), row.names = FALSE)
  cat("已保存", nrow(cis_p_0.05), "个cis-QTL (p<0.05)\n")
} else {
  cat("没有发现cis-QTL (p<0.05)\n")
}

if (!is.null(trans_p_0.05) && nrow(trans_p_0.05) > 0) {
  write.csv(trans_p_0.05, file.path(p_0.05_dir, "trans_QTLs_p_0.05.csv"), row.names = FALSE)
  cat("已保存", nrow(trans_p_0.05), "个trans-QTL (p<0.05)\n")
} else {
  cat("没有发现trans-QTL (p<0.05)\n")
}

if (!is.null(all_p_0.05) && nrow(all_p_0.05) > 0) {
  write.csv(all_p_0.05, file.path(p_0.05_dir, "all_QTLs_p_0.05.csv"), row.names = FALSE)
  cat("已保存", nrow(all_p_0.05), "个所有QTL (p<0.05)\n")
} else {
  cat("没有发现任何QTL (p<0.05)\n")
}

# p < 0.01 的结果
cat("筛选p < 0.01的结果...\n")
cis_p_0.01 <- filter_by_pvalue(me$cis$eqtls, 0.01)
trans_p_0.01 <- filter_by_pvalue(me$trans$eqtls, 0.01)
all_p_0.01 <- filter_by_pvalue(me$all$eqtls, 0.01)

# 保存p < 0.01 的结果
if (!is.null(cis_p_0.01) && nrow(cis_p_0.01) > 0) {
  write.csv(cis_p_0.01, file.path(p_0.01_dir, "cis_QTLs_p_0.01.csv"), row.names = FALSE)
  cat("已保存", nrow(cis_p_0.01), "个cis-QTL (p<0.01)\n")
} else {
  cat("没有发现cis-QTL (p<0.01)\n")
}

if (!is.null(trans_p_0.01) && nrow(trans_p_0.01) > 0) {
  write.csv(trans_p_0.01, file.path(p_0.01_dir, "trans_QTLs_p_0.01.csv"), row.names = FALSE)
  cat("已保存", nrow(trans_p_0.01), "个trans-QTL (p<0.01)\n")
} else {
  cat("没有发现trans-QTL (p<0.01)\n")
}

if (!is.null(all_p_0.01) && nrow(all_p_0.01) > 0) {
  write.csv(all_p_0.01, file.path(p_0.01_dir, "all_QTLs_p_0.01.csv"), row.names = FALSE)
  cat("已保存", nrow(all_p_0.01), "个所有QTL (p<0.01)\n")
} else {
  cat("没有发现任何QTL (p<0.01)\n")
}

# 16. 创建详细的统计总结文件
cat("\n=== 创建统计总结文件 ===\n")

# 获取当前时间
current_time <- Sys.time()

# 创建统计总结的文本文件
summary_file <- file.path(output_dir, "QTL分析统计总结.txt")
file_conn <- file(summary_file, "w")

# 写入统计总结
writeLines("==========================================", file_conn)
writeLines("        QTL分析统计总结报告", file_conn)
writeLines("==========================================", file_conn)
writeLines(paste("分析时间:", current_time), file_conn)
writeLines(paste("分析人员:", Sys.info()["user"]), file_conn)
writeLines("", file_conn)

# 第一部分：数据基本信息
writeLines("一、数据基本信息", file_conn)
writeLines(paste("1. eRNA文件:", eRNA_file), file_conn)
writeLines(paste("   数据维度:", nrow(eRNA_data), "行 ×", ncol(eRNA_data), "列"), file_conn)
writeLines(paste("   前三列为:", paste(colnames(eRNA_data)[1:3], collapse = ", ")), file_conn)
writeLines(paste("2. 基因型文件:", genotype_file), file_conn)
writeLines(paste("   数据维度:", nrow(genotype_data), "行 ×", ncol(genotype_data), "列"), file_conn)
writeLines(paste("   前三列为:", paste(colnames(genotype_data)[1:3], collapse = ", ")), file_conn)
writeLines("", file_conn)

# 第二部分：数据预处理信息
writeLines("二、数据预处理信息", file_conn)
writeLines(paste("1. 公共样本数量:", length(common_samples)), file_conn)
writeLines(paste("2. 处理后eRNA数量:", nrow(expression_matrix)), file_conn)
writeLines(paste("3. 处理后SNP数量:", nrow(genotype_numeric)), file_conn)
writeLines(paste("4. 基因型缺失率:", round(missing_rate * 100, 2), "%"), file_conn)
writeLines(paste("5. 表达数据是否log2转换:", ifelse(max(expression_matrix, na.rm = TRUE) > 100, "是", "否")), file_conn)
writeLines("", file_conn)

# 第三部分：分析参数
writeLines("三、分析参数", file_conn)
writeLines(paste("1. cis距离阈值:", "1 Mb (1,000,000 bp)"), file_conn)
writeLines(paste("2. p值阈值(cis):", "0.05"), file_conn)
writeLines(paste("3. p值阈值(trans):", "0.05"), file_conn)
writeLines(paste("4. 使用的模型:", "modelLINEAR"), file_conn)
writeLines(paste("5. 输出目录:", output_dir), file_conn)
writeLines("", file_conn)

# 第四部分：QTL分析总体统计
writeLines("四、QTL分析总体统计", file_conn)
writeLines(paste("1. 总共测试的SNP-基因对:", format(me$all$ntests, big.mark = ",")), file_conn)
writeLines(paste("2. 测试的SNP数量:", nrow(genotype_numeric)), file_conn)
writeLines(paste("3. 测试的eRNA数量:", nrow(expression_matrix)), file_conn)
writeLines(paste("4. 样本数量:", ncol(expression_matrix)), file_conn)
writeLines("", file_conn)

# 第五部分：p值阈值结果统计
writeLines("五、p值阈值结果统计", file_conn)
writeLines("1. p < 0.05:", file_conn)
writeLines(paste("   - cis-QTL数量:", ifelse(!is.null(cis_p_0.05), nrow(cis_p_0.05), 0)), file_conn)
writeLines(paste("   - trans-QTL数量:", ifelse(!is.null(trans_p_0.05), nrow(trans_p_0.05), 0)), file_conn)
writeLines(paste("   - 总计:", ifelse(!is.null(all_p_0.05), nrow(all_p_0.05), 0)), file_conn)
writeLines("", file_conn)

writeLines("2. p < 0.01:", file_conn)
writeLines(paste("   - cis-QTL数量:", ifelse(!is.null(cis_p_0.01), nrow(cis_p_0.01), 0)), file_conn)
writeLines(paste("   - trans-QTL数量:", ifelse(!is.null(trans_p_0.01), nrow(trans_p_0.01), 0)), file_conn)
writeLines(paste("   - 总计:", ifelse(!is.null(all_p_0.01), nrow(all_p_0.01), 0)), file_conn)
writeLines("", file_conn)

# 第六部分：FDR阈值结果统计
writeLines("六、FDR阈值结果统计", file_conn)
writeLines("1. FDR < 0.05:", file_conn)
writeLines(paste("   - cis-QTL数量:", ifelse(!is.null(cis_fdr_0.05), nrow(cis_fdr_0.05), 0)), file_conn)
writeLines(paste("   - trans-QTL数量:", ifelse(!is.null(trans_fdr_0.05), nrow(trans_fdr_0.05), 0)), file_conn)
writeLines(paste("   - 总计:", ifelse(!is.null(all_fdr_0.05), nrow(all_fdr_0.05), 0)), file_conn)
writeLines("", file_conn)

writeLines("2. FDR < 0.01:", file_conn)
writeLines(paste("   - cis-QTL数量:", ifelse(!is.null(cis_fdr_0.01), nrow(cis_fdr_0.01), 0)), file_conn)
writeLines(paste("   - trans-QTL数量:", ifelse(!is.null(trans_fdr_0.01), nrow(trans_fdr_0.01), 0)), file_conn)
writeLines(paste("   - 总计:", ifelse(!is.null(all_fdr_0.01), nrow(all_fdr_0.01), 0)), file_conn)
writeLines("", file_conn)

# 第七部分：染色体分布统计
writeLines("七、染色体分布统计", file_conn)
writeLines("1. SNP染色体分布:", file_conn)
snp_chr_table <- table(genotype_data$CHROM)
for (chr in names(snp_chr_table)) {
  writeLines(paste("   - 染色体", chr, ":", snp_chr_table[chr], "个SNP"), file_conn)
}
writeLines("", file_conn)

writeLines("2. eRNA染色体分布:", file_conn)
gene_chr_table <- table(gene_pos$chr)
for (chr in names(gene_chr_table)) {
  writeLines(paste("   - 染色体", chr, ":", gene_chr_table[chr], "个eRNA"), file_conn)
}
writeLines("", file_conn)

# 第八部分：最显著关联
writeLines("八、最显著的QTL关联", file_conn)
if (!is.null(me$all$eqtls) && nrow(me$all$eqtls) > 0) {
  top5 <- head(me$all$eqtls[order(me$all$eqtls$pvalue), ], 5)
  writeLines("最显著的5个关联:", file_conn)
  for (i in 1:nrow(top5)) {
    writeLines(paste("   ", i, ". SNP:", top5$snps[i], " 基因:", top5$gene[i],
                     " p值:", format(top5$pvalue[i], scientific = TRUE, digits = 3),
                     " FDR:", format(top5$FDR[i], scientific = TRUE, digits = 3),
                     " beta:", round(top5$beta[i], 4)), file_conn)
  }
} else if (!is.null(me$cis$eqtls) && nrow(me$cis$eqtls) > 0) {
  top5 <- head(me$cis$eqtls[order(me$cis$eqtls$pvalue), ], 5)
  writeLines("最显著的5个cis关联:", file_conn)
  for (i in 1:nrow(top5)) {
    writeLines(paste("   ", i, ". SNP:", top5$snps[i], " 基因:", top5$gene[i],
                     " p值:", format(top5$pvalue[i], scientific = TRUE, digits = 3),
                     " FDR:", format(top5$FDR[i], scientific = TRUE, digits = 3),
                     " beta:", round(top5$beta[i], 4)), file_conn)
  }
} else if (!is.null(me$trans$eqtls) && nrow(me$trans$eqtls) > 0) {
  top5 <- head(me$trans$eqtls[order(me$trans$eqtls$pvalue), ], 5)
  writeLines("最显著的5个trans关联:", file_conn)
  for (i in 1:nrow(top5)) {
    writeLines(paste("   ", i, ". SNP:", top5$snps[i], " 基因:", top5$gene[i],
                     " p值:", format(top5$pvalue[i], scientific = TRUE, digits = 3),
                     " FDR:", format(top5$FDR[i], scientific = TRUE, digits = 3),
                     " beta:", round(top5$beta[i], 4)), file_conn)
  }
} else {
  writeLines("没有发现显著的QTL关联", file_conn)
}
writeLines("", file_conn)

# 第九部分：输出文件清单
writeLines("九、生成的输出文件", file_conn)
output_files <- list.files(output_dir, recursive = TRUE, full.names = FALSE)
if (length(output_files) > 0) {
  for (file in output_files) {
    file_size <- file.size(file.path(output_dir, file))
    writeLines(paste("   - ", file, " (", round(file_size/1024, 2), " KB)", sep = ""), file_conn)
  }
} else {
  writeLines("   没有生成任何输出文件", file_conn)
}
writeLines("", file_conn)

# 第十部分：分析总结
writeLines("十、分析总结", file_conn)
total_significant <- ifelse(!is.null(me$all$eqtls), nrow(me$all$eqtls), 0)
total_cis <- ifelse(!is.null(me$cis$eqtls), nrow(me$cis$eqtls), 0)
total_trans <- ifelse(!is.null(me$trans$eqtls), nrow(me$trans$eqtls), 0)

writeLines(paste("本次QTL分析共发现了", total_significant, "个显著的SNP-eRNA关联。"), file_conn)
writeLines(paste("其中，cis关联有", total_cis, "个，trans关联有", total_trans, "个。"), file_conn)

if (total_significant > 0) {
  if (total_cis > total_trans) {
    writeLines("分析结果表明，cis关联占主导地位，这与大多数eQTL研究的结果一致。", file_conn)
  } else if (total_trans > total_cis) {
    writeLines("分析结果表明，trans关联占主导地位，这可能表明存在广泛的调控网络。", file_conn)
  } else {
    writeLines("分析结果表明，cis和trans关联数量相当。", file_conn)
  }
} else {
  writeLines("本次分析未发现显著的SNP-eRNA关联，可能需要调整分析参数或增加样本量。", file_conn)
}
writeLines("", file_conn)

writeLines("==========================================", file_conn)
writeLines("          QTL分析统计总结报告结束", file_conn)
writeLines("==========================================", file_conn)

# 关闭文件连接
close(file_conn)

cat("统计总结文件已保存到:", summary_file, "\n")

# 17. 在控制台显示关键统计信息
cat("\n=== QTL分析关键统计信息 ===\n")
cat("📊 总共测试了", format(me$all$ntests, big.mark = ","), "个SNP-基因对\n")
cat("🧬 使用的SNP数量:", nrow(genotype_numeric), "\n")
cat("🔬 使用的eRNA数量:", nrow(expression_matrix), "\n")
cat("👥 样本数量:", ncol(expression_matrix), "\n")
cat("🔍 发现", ifelse(!is.null(me$all$eqtls), nrow(me$all$eqtls), 0), "个显著的QTL (p < 0.05)\n")
cat("📍 其中cis-QTL:", ifelse(!is.null(me$cis$eqtls), nrow(me$cis$eqtls), 0), "个\n")
cat("🌐 其中trans-QTL:", ifelse(!is.null(me$trans$eqtls), nrow(me$trans$eqtls), 0), "个\n")
cat("📈 基因型缺失率:", round(missing_rate * 100, 2), "%\n")

# p值结果汇总
cat("\n=== p值阈值结果 ===\n")
cat("p < 0.05:\n")
cat("  - cis-QTL:", ifelse(!is.null(cis_p_0.05), nrow(cis_p_0.05), 0), "个\n")
cat("  - trans-QTL:", ifelse(!is.null(trans_p_0.05), nrow(trans_p_0.05), 0), "个\n")
cat("  - 总计:", ifelse(!is.null(all_p_0.05), nrow(all_p_0.05), 0), "个\n")

cat("p < 0.01:\n")
cat("  - cis-QTL:", ifelse(!is.null(cis_p_0.01), nrow(cis_p_0.01), 0), "个\n")
cat("  - trans-QTL:", ifelse(!is.null(trans_p_0.01), nrow(trans_p_0.01), 0), "个\n")
cat("  - 总计:", ifelse(!is.null(all_p_0.01), nrow(all_p_0.01), 0), "个\n")

# FDR结果汇总
cat("\n=== FDR阈值结果 ===\n")
cat("FDR < 0.05:\n")
cat("  - cis-QTL:", ifelse(!is.null(cis_fdr_0.05), nrow(cis_fdr_0.05), 0), "个\n")
cat("  - trans-QTL:", ifelse(!is.null(trans_fdr_0.05), nrow(trans_fdr_0.05), 0), "个\n")
cat("  - 总计:", ifelse(!is.null(all_fdr_0.05), nrow(all_fdr_0.05), 0), "个\n")

cat("FDR < 0.01:\n")
cat("  - cis-QTL:", ifelse(!is.null(cis_fdr_0.01), nrow(cis_fdr_0.01), 0), "个\n")
cat("  - trans-QTL:", ifelse(!is.null(trans_fdr_0.01), nrow(trans_fdr_0.01), 0), "个\n")
cat("  - 总计:", ifelse(!is.null(all_fdr_0.01), nrow(all_fdr_0.01), 0), "个\n")

cat("\n💾 统计总结文件已保存到:", summary_file, "\n")
cat("💾 结果文件保存在:", output_dir, "\n")
cat("✨ 所有分析步骤完成！请查看输出目录中的结果文件。\n")
yn
