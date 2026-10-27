# =========================================================
# 对应结果：Figure 8、Figure S8、Figure S9 与 Supplementary Table 12 的多层共表达关联网络——网络各层的输入矩阵
#           （GSE144269 的 70 个肿瘤样本：风险 eRNA 表达、PCD 相关基因、±1 Mb 候选靶基因、免疫治疗相关基因）
# 输入：../data/samples/GSE144269_samples.csv、../data/expression/eRNA_RPM/GSE144269_eRNA_RPM.csv.gz
#       ../data/expression/gene_TPM/GSE144269_gene_TPM.csv.gz、input/risk_eRNAs_by_PCD.csv、input/PCD_genes_<PCD>.csv
#       input/candidate_target_genes_1Mb.csv、input/immunotherapy_gene_set.csv
# 输出：output/inputs/eRNA_expression/<eRNA>表达.csv、output/inputs/GSE144269_tumor70_gene_TPM.csv
#       output/inputs/PCD_genes/<PCD>N个高度相关基因表达矩阵.csv、output/inputs/target_genes/<eRNA>_1MB靶基因矩阵.csv
#       output/inputs/免疫治疗基因表达矩阵.csv（及 提取统计.txt、免疫治疗基因_未匹配到Symbol的清单.csv）
#       免疫细胞 ssGSEA 得分与 19 种候选药物的预测 IC50 不在这里生成，直接使用 input/ 中的矩阵
# 运行：在模块目录下 Rscript code/02_prepare_network_inputs.R
# =========================================================

# ---------------------------------------------------------
# 0. 读入共享表达矩阵，取 GSE144269 的 70 个肿瘤样本（列顺序与样本表相同，所有矩阵一致）
# ---------------------------------------------------------
out_root <- file.path("output", "inputs")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
samples <- read.csv("../data/samples/GSE144269_samples.csv", stringsAsFactors = FALSE)
tumor_samples <- samples$name[samples$Treatment == "tumor"]
rpm <- read.csv(gzfile("../data/expression/eRNA_RPM/GSE144269_eRNA_RPM.csv.gz"), check.names = FALSE, stringsAsFactors = FALSE)
tpm <- read.csv(gzfile("../data/expression/gene_TPM/GSE144269_gene_TPM.csv.gz"), check.names = FALSE, stringsAsFactors = FALSE)
risk <- read.csv("input/risk_eRNAs_by_PCD.csv", stringsAsFactors = FALSE)
cat("肿瘤样本数:", length(tumor_samples), "\n")

# 70 个肿瘤样本的编码基因 TPM（第一列为 Symbol）
tpm_tumor <- data.frame(Symbol = tpm[[1]], tpm[, tumor_samples], check.names = FALSE)
tpm_file <- file.path(out_root, "GSE144269_tumor70_gene_TPM.csv")
write.csv(tpm_tumor, tpm_file, row.names = FALSE)

# ---------------------------------------------------------
# ---------------------------------------------------------
dir.create(file.path(out_root, "eRNA_expression"), showWarnings = FALSE)
for (erna in unique(risk$eRNA)) {
  v <- log2(as.numeric(rpm[rpm[[1]] == erna, tumor_samples]) + 1)
  df <- data.frame(eRNA = erna, t(v), check.names = FALSE)
  colnames(df) <- c("eRNA", tumor_samples)
  write.csv(df, file.path(out_root, "eRNA_expression", paste0(gsub(":", "", erna), "表达.csv")), row.names = FALSE)
}

# ---------------------------------------------------------
# 2. PCD 相关基因（6 个数据集中至少 5 个显著的 eRNA–基因对所含基因：失巢凋亡 20 个、凋亡 40 个、铁死亡 50 个）
# ---------------------------------------------------------
pcd_list <- unique(risk[, c("PCD", "PCD_cn")])
n_genes <- sapply(pcd_list$PCD, function(p) nrow(read.csv(file.path("input", paste0("PCD_genes_", p, ".csv")))))
for (k in seq_len(nrow(pcd_list))) {
  pcd <- pcd_list$PCD[k]; pcd_cn <- pcd_list$PCD_cn[k]
  # 设置工作路径（可选，非必须）

  # 定义输入和输出路径（注意：Windows 路径中的反斜杠需转为正斜杠 / 或双反斜杠 \\）
  input_dir  <- "input"
  output_dir <- file.path(out_root, "PCD_genes")

  # 创建输出目录（如果不存在）
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # 文件路径
  expr_file      <- tpm_file
  gene_list_file <- file.path(input_dir, paste0("PCD_genes_", pcd, ".csv"))
  output_file    <- file.path(output_dir, paste0(pcd_cn, n_genes[pcd], "个高度相关基因表达矩阵.csv"))

  # 读取基因列表（只有一列，列名为 'Gene'）
  gene_list_df <- read.csv(gene_list_file, header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8")
  target_genes <- trimws(gene_list_df$Gene)  # 去除首尾空格
  target_genes <- target_genes[target_genes != ""]  # 去除空字符串

  # 读取表达矩阵（第一列为 Symbol）
  expr_df <- read.csv(expr_file, header = TRUE, row.names = NULL, stringsAsFactors = FALSE, encoding = "UTF-8")

  # 检查是否存在 'Symbol' 列
  if (!"Symbol" %in% colnames(expr_df)) {
    stop("错误：表达矩阵中未找到 'Symbol' 列！")
  }

  # 提取匹配的行（保留原始顺序）
  matched_rows <- expr_df$Symbol %in% target_genes
  filtered_df <- expr_df[matched_rows, , drop = FALSE]

  # 可选：报告缺失基因
  found_genes <- expr_df$Symbol[matched_rows]
  missing_genes <- setdiff(target_genes, found_genes)
  if (length(missing_genes) > 0) {
    cat("⚠️ 注意：以下", length(missing_genes), "个基因在表达矩阵中未找到：\n")
    print(sort(missing_genes))
  }

  # 写入输出文件（UTF-8 with BOM，确保 Excel 正确显示中文）
  write.csv(filtered_df, output_file, row.names = FALSE, fileEncoding = "UTF-8")

  cat("\n✅ 成功提取", nrow(filtered_df), "个基因的表达数据！\n")
  cat("📁 已保存至：", output_file, "\n")
}

# ---------------------------------------------------------
# ---------------------------------------------------------
tg <- read.csv("input/candidate_target_genes_1Mb.csv", stringsAsFactors = FALSE)
dir.create(file.path(out_root, "target_genes"), showWarnings = FALSE)
for (erna in unique(risk$eRNA)) {
  target_genes <- tg$Symbol[tg$eRNA == erna]
  filtered_df <- tpm_tumor[tpm_tumor$Symbol %in% target_genes, , drop = FALSE]
  missing_genes <- setdiff(target_genes, filtered_df$Symbol)
  if (length(missing_genes) > 0) cat("⚠️ 注意：", erna, "的以下靶基因在表达矩阵中未找到：", missing_genes, "\n")
  write.csv(filtered_df, file.path(out_root, "target_genes", paste0(gsub(":", "", erna), "_1MB靶基因矩阵.csv")), row.names = FALSE)
  cat(erna, ": ±1 Mb 候选靶基因", nrow(filtered_df), "个\n")
}

# ---------------------------------------------------------
# ---------------------------------------------------------
gene_set <- read.csv("input/immunotherapy_gene_set.csv", stringsAsFactors = FALSE)
set_genes <- unique(trimws(gene_set$Gene))
expr_unique <- tpm_tumor[!duplicated(tpm_tumor$Symbol), ]
matched <- set_genes[set_genes %in% expr_unique$Symbol]
not_found <- setdiff(set_genes, expr_unique$Symbol)
immuno_df <- expr_unique[match(matched, expr_unique$Symbol), ]
colnames(immuno_df)[1] <- "Gene"
write.csv(immuno_df, file.path(out_root, "免疫治疗基因表达矩阵.csv"), row.names = FALSE)
write.csv(data.frame(Gene = not_found), file.path(out_root, "免疫治疗基因_未匹配到Symbol的清单.csv"), row.names = FALSE)
writeLines(c(paste("Expression rows (unique Symbol kept):", nrow(expr_unique)),
             paste("Gene set size (unique Gene kept):", length(set_genes)),
             paste("Matched genes:", length(matched)),
             paste("Not found genes:", length(not_found))),
           file.path(out_root, "提取统计.txt"))
cat("免疫治疗基因：", length(matched), "/", length(set_genes), "个匹配\n")
