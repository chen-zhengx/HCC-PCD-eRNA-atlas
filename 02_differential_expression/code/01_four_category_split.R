# =========================================================
# 对应结果：差异分析第 1 步（Figure 2A/2B 的上游）——按配对样本的表达状态把 eRNA/基因分为 a、b、c、d 四类
#           再按三个条件拆分表达矩阵：条件 1 d > n/2（剔除）；条件 2 a > n/2（limma）；条件 3 a+b+c > n/2 且 a ≤ n/2（McNemar）
# 输入：../data/expression/eRNA_RPM/<GSE>_eRNA_RPM.csv.gz、../data/expression/gene_TPM/<GSE>_gene_TPM.csv.gz
#       ../data/samples/<GSE>_samples.csv（name, PatientID, Treatment）、../data/annotation/protein_coding_genes_19416.csv
# 输出：output/<eRNA|gene>/<GSE>/pair_type_counts.csv、Condition1/2/3_*.csv 及对应表达矩阵
# 运行：在模块目录下 Rscript code/01_four_category_split.R
# =========================================================
# 加载包
library(limma)
library(openxlsx)
library(ggplot2)  # 用于验证绘图

datasets <- c("GSE148355", "GSE77509", "GSE94660", "GSE124535", "GSE214846", "GSE144269")
features <- list(
  eRNA = "../data/expression/eRNA_RPM/%s_eRNA_RPM.csv.gz",   # eRNA：RPM，RPM>=1 为表达
  gene = "../data/expression/gene_TPM/%s_gene_TPM.csv.gz"    # 基因：TPM，TPM>=1 为表达
)
gene_ids <- read.csv("../data/annotation/protein_coding_genes_19416.csv")

for (feature in names(features)) {
for (dataset in datasets) {
cat("\n==========", feature, dataset, "==========\n")
output_dir <- file.path("output", feature, dataset)
# 创建输出目录
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------- 数据加载与预处理 ----------
# 读取表达矩阵
expr_file <- sprintf(features[[feature]], dataset)
expr_data <- read.csv(expr_file,
                      row.names = 1, check.names = FALSE)
# GSE214846 的基因 TPM 矩阵行名为 NCBI GeneID，换成基因符号（与其余数据集一致）
if (feature == "gene" && dataset == "GSE214846") {
  rownames(expr_data) <- gene_ids$Symbol[match(rownames(expr_data), gene_ids$GeneID)]
}
expr_matrix <- as.matrix(expr_data)
#读取样本信息（包含name, PatientID, Treatment）
sample_info <- read.csv(sprintf("../data/samples/%s_samples.csv", dataset),
                        stringsAsFactors = FALSE)
# ---------- 关键修改1：重构样本配对逻辑 ----------
# 定义分组标签
sample_info$Group <- factor(
  sample_info$Treatment,
  levels = c("nontumor", "tumor"),
  labels = c("Normal", "Tumor")
)
n_pairs <- length(unique(sample_info$PatientID))
if (nrow(sample_info) != 2 * n_pairs) {
  stop("样本信息必须为成对样本！")
}
# 按PatientID排序确保成对顺序
sample_info <- sample_info[order(sample_info$PatientID), ]
expr_matrix <- expr_matrix[, sample_info$name]  # 对齐表达矩阵列顺序
# ---------- 关键修改2：动态提取正常/肿瘤样本列 ----------
# 获取Normal和Tumor的列索引
normal_cols <- which(sample_info$Group == "Normal")
tumor_cols <- which(sample_info$Group == "Tumor")
# 分割表达矩阵（按配对的Normal和Tumor）
norm_part <- expr_matrix[, normal_cols]  # 正常组
tumor_part <- expr_matrix[, tumor_cols]  # 肿瘤组
# 计算四种情况的配对数
a <- rowSums(norm_part >= 1 & tumor_part >= 1)  # 两者均>=1
b <- rowSums(norm_part >= 1 & tumor_part < 1)   # 正常>=1且肿瘤<1
c <- rowSums(norm_part < 1 & tumor_part >= 1)   # 正常<1且肿瘤>=1
d <- rowSums(norm_part < 1 & tumor_part < 1)    # 两者均<1
# 合并统计结果
pair_counts <- data.frame(
  Gene = rownames(expr_matrix),
  a = a,
  b = b,
  c = c,
  d = d,
  row.names = NULL
)
# 保存结果
write.csv(pair_counts,
          file.path(output_dir, "pair_type_counts.csv"),
          row.names = FALSE)
# ---------- 定义三个筛选条件 ----------
half <- n_pairs / 2
# 条件1：d > n/2
condition1 <- pair_counts$d > half
# 条件2：a > n/2
condition2 <- pair_counts$a > half
# 条件3：a + b + c > n/2 且 a <= n/2
condition3 <- (pair_counts$a + pair_counts$b + pair_counts$c) > half &
  pair_counts$a <= half
# ---------- 分别筛选并保存 ----------
# 条件1结果
genes_cond1 <- pair_counts[condition1, ]
write.csv(genes_cond1,
          file.path(output_dir, "Condition1_d_gt_half.csv"),
          row.names = FALSE)
# 条件2结果
genes_cond2 <- pair_counts[condition2, ]
write.csv(genes_cond2,
          file.path(output_dir, "Condition2_a_gt_half.csv"),
          row.names = FALSE)
# 条件3结果
genes_cond3 <- pair_counts[condition3, ]
write.csv(genes_cond3,
          file.path(output_dir, "Condition3_abc_gt_half.csv"),
          row.names = FALSE)
# 统计各条件的基因数
count_stats <- data.frame(
  Condition = c("Condition1 (d > n/2)",
                "Condition2 (a > n/2)",
                "Condition3 (a <= n/2 & a+b+c > n/2)"),
  Gene_Count = c(nrow(genes_cond1), nrow(genes_cond2), nrow(genes_cond3))
)
write.csv(count_stats, file.path(output_dir, "各条件基因数量统计.csv"), row.names = FALSE)
print(count_stats)
# ---------- 检查是否覆盖所有基因 ----------
total_genes <- nrow(pair_counts)
selected_genes <- unique(rbind(genes_cond1, genes_cond2, genes_cond3))
if (nrow(selected_genes) == total_genes) {
  cat("✅ 所有基因均被覆盖，无遗漏。\n")
} else {
  cat("❌ 存在未覆盖基因，需检查条件逻辑！\n")
  uncovered <- pair_counts[!(condition1 | condition2 | condition3), ]
  write.csv(uncovered, file.path(output_dir, "未覆盖基因列表.csv"), row.names = FALSE)
}
# ---------- 从原始数据中提取基因表达值并保存 ----------
# 确保 expr_data 的行名与 pair_counts$Gene 一致
# 提取并保存条件1的基因表达数据
expr_cond1 <- expr_data[genes_cond1$Gene, ]
write.csv(expr_cond1,
          file.path(output_dir, "Condition1_d_gt_half_表达矩阵.csv"),
          row.names = TRUE)
# 提取并保存条件2的基因表达数据
expr_cond2 <- expr_data[genes_cond2$Gene, ]
write.csv(expr_cond2,
          file.path(output_dir, "Condition2_a_gt_half_表达矩阵.csv"),
          row.names = TRUE)
# 提取并保存条件3的基因表达数据
expr_cond3 <- expr_data[genes_cond3$Gene, ]
write.csv(expr_cond3,
          file.path(output_dir, "Condition3_abc_gt_half_表达矩阵.csv"),
          row.names = TRUE)
}
}
