# =========================================================
# 对应结果：Figure 3A–C、Figure S1A–C、Figure S2A–C、Table S7、Table S8 的建模输入
#           （三种 PCD 相关 eRNA 在 GSE144269 70 个肿瘤样本中的 RPM 表达矩阵）
# 输入：../data/expression/eRNA_RPM/GSE144269_eRNA_RPM.csv.gz
#       ../data/samples/GSE144269_samples.csv
#       ../data/clinical/GSE144269_survival.csv
#       input/<pcd>_eRNA_RPM_tumor70.csv
# 输出：output/extracted/<pcd>_eRNA_RPM_tumor70_from_RPM.csv
# 运行：在模块目录下 Rscript code/01_extract_eRNA_matrix.R
# =========================================================


# ====== 设置路径 ======
rpm_file   <- "../data/expression/eRNA_RPM/GSE144269_eRNA_RPM.csv.gz"
samp_file  <- "../data/samples/GSE144269_samples.csv"
surv_file  <- "../data/clinical/GSE144269_survival.csv"
input_dir  <- "input"
output_dir <- "output/extracted"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

pcd_list <- c("anoikis", "apoptosis", "ferroptosis")   # 失巢凋亡、凋亡、铁死亡

# ====== 读取数据 ======
cat("📖 读取 GSE144269 eRNA RPM 矩阵（约 21 万行 × 140 列，需要 1–2 分钟）...\n")
rpm_data  <- read.csv(rpm_file, row.names = 1, check.names = FALSE)
samples   <- read.csv(samp_file, stringsAsFactors = FALSE)
surv_data <- read.csv(surv_file, stringsAsFactors = FALSE)
cat("✅ RPM 矩阵:", nrow(rpm_data), "个 eRNA ×", ncol(rpm_data), "个样本\n")

# 70 个肿瘤样本（Treatment == "tumor"）；
# 列顺序与生存数据（建模输入）一致，即按病人编号 pat 1 – pat 70 排列。
# 注意：cv.glmnet 按行位置随机分折，样本顺序会影响 LASSO 结果，因此必须保持这一顺序
tumor_samples <- samples$name[samples$Treatment == "tumor"]
stopifnot(length(tumor_samples) == 70, setequal(tumor_samples, surv_data$Sample))
tumor_samples <- as.character(surv_data$Sample)

# ====== 按 PCD 提取 ======
for (pcd in pcd_list) {
  erna_file <- file.path(input_dir, paste0(pcd, "_eRNA_RPM_tumor70.csv"))
  erna_data <- read.csv(erna_file, row.names = 1, check.names = FALSE)

  # PCD 相关 eRNA：上游 eRNA–PCD 基因相关性分析中六个数据集共有的 eRNA（196 / 240 / 276 个），
  # 顺序与 RPM 矩阵相同
  erna_ids <- rownames(erna_data)
  stopifnot(all(erna_ids %in% rownames(rpm_data)))

  exp_sub <- rpm_data[erna_ids, tumor_samples, drop = FALSE]
  write.csv(exp_sub, file.path(output_dir, paste0(pcd, "_eRNA_RPM_tumor70_from_RPM.csv")))
  cat(sprintf("%-12s %d 个 eRNA × %d 个样本\n", pcd, nrow(exp_sub), ncol(exp_sub)))
}

cat("\n🎉 完成，输出位于:", output_dir, "\n")
