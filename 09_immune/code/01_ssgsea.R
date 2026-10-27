# =========================================================
# 对应结果：70 个肿瘤样本 28 种免疫细胞的 ssGSEA 富集得分（Figure 6A–C、Table S10 的上游；Methods "Immune signature analysis"）
# 输入：input/GSE144269_tumor70_allgenes_TPM.csv（GSE144269 的 70 个肿瘤样本，39,373 个基因 × 70 列，TPM）
#       input/immune_28celltypes.gmt
# 输出：output/ssgsea/ssgsea_immune_enrichment_scores.csv
#       output/ssgsea/ssgsea_results.RData
#       output/ssgsea/ssgsea_scores_28cells.csv（ 行名整理后的得分矩阵，与 input/ssgsea_scores_28cells.csv 相同）
# 运行：在模块目录下 Rscript code/01_ssgsea.R
# =========================================================

# ===============================
# 0. 设置路径（修正版）
# ===============================
input_path  <- "input"
output_path <- "output/ssgsea"
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# ===============================
# 1. 加载必要的包
# ===============================
cat("\n===== 步骤1: 加载必要的R包 =====\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GSVA", quietly = TRUE)) BiocManager::install("GSVA")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(GSVA)

# ===============================
# 2. 准备表达数据（Gene列做行名）
# ===============================
cat("\n===== 步骤2: 准备表达数据 =====\n")
expression_file <- file.path(input_path, "GSE144269_tumor70_allgenes_TPM.csv")
cat("正在读取表达矩阵文件:", expression_file, "\n")

expression_df <- read.csv(expression_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

if (!("Gene" %in% colnames(expression_df))) {
  stop("表达矩阵缺少 Gene 列。当前列名前5个：", paste(colnames(expression_df)[1:5], collapse = ", "))
}

gene_names <- toupper(trimws(as.character(expression_df$Gene)))
expression_df$Gene <- NULL

expression_data <- as.matrix(expression_df)
mode(expression_data) <- "numeric"
rownames(expression_data) <- gene_names

keep <- !is.na(rownames(expression_data)) & rownames(expression_data) != ""
expression_data <- expression_data[keep, , drop = FALSE]

cat("表达矩阵维度:", dim(expression_data), "\n")
cat("样本数量:", ncol(expression_data), "\n")
cat("基因数量:", nrow(expression_data), "\n")
cat("基因名示例:", paste(head(rownames(expression_data), 5), collapse = ", "), "\n")

dup_n <- sum(duplicated(rownames(expression_data)))
cat("重复基因名数量:", dup_n, "\n")
if (dup_n > 0) {
  cat("检测到重复基因名，正在合并（取max）...\n")
  df_tmp <- as.data.frame(expression_data, check.names = FALSE)
  df_tmp$Gene <- rownames(expression_data)
  df_tmp <- df_tmp |>
    dplyr::group_by(Gene) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~max(.x, na.rm = TRUE)), .groups = "drop")
  expression_data <- as.matrix(df_tmp[, -1, drop = FALSE])
  rownames(expression_data) <- df_tmp$Gene
  mode(expression_data) <- "numeric"
  cat("合并后表达矩阵维度:", dim(expression_data), "\n")
}

# ===============================
# 3. 准备基因集（读GMT）
# ===============================
cat("\n===== 步骤3: 准备基因集 =====\n")
gmt_file <- file.path(input_path, "immune_28celltypes.gmt")
cat("正在读取GMT文件:", gmt_file, "\n")

read_gmt <- function(file) {
  gmt <- list()
  lines <- readLines(file, warn = FALSE)
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    name  <- parts[1]
    genes <- parts[3:length(parts)]
    genes <- genes[genes != ""]
    genes <- toupper(genes)
    gmt[[name]] <- genes
  }
  gmt
}

gene_sets <- read_gmt(gmt_file)
cat("加载的基因集数量:", length(gene_sets), "\n")
cat("前5个基因集:", paste(names(gene_sets)[1:min(5, length(gene_sets))], collapse = ", "), "\n")

# ===============================
# 4. 基因匹配诊断
# ===============================
cat("\n===== 基因匹配诊断 =====\n")

check_gene_overlap <- function(expr_data, gene_sets) {
  expr_genes <- rownames(expr_data)
  all_geneset_genes <- unique(unlist(gene_sets))

  overlap <- intersect(expr_genes, all_geneset_genes)
  overlap_percent <- round(length(overlap) / length(all_geneset_genes) * 100, 1)

  cat("表达矩阵基因数:", length(expr_genes), "\n")
  cat("基因集总基因数:", length(all_geneset_genes), "\n")
  cat("匹配基因数:", length(overlap), "(", overlap_percent, "%)\n")

  set_overlap <- sapply(gene_sets, function(geneset) mean(geneset %in% expr_genes))
  cat("\n基因集匹配度分布:\n")
  print(summary(set_overlap))

  low_overlap_sets <- names(which(set_overlap < 0.3))
  if (length(low_overlap_sets) > 0) {
    cat("\n警告: 以下基因集匹配度低于30%:\n")
    for (set_name in low_overlap_sets) {
      total_genes <- length(gene_sets[[set_name]])
      matched_genes <- sum(gene_sets[[set_name]] %in% expr_genes)
      cat(set_name, ": ", matched_genes, "/", total_genes,
          "(", round(matched_genes / total_genes * 100, 1), "%)\n", sep = "")
    }
  }

  list(matched_genes = overlap, set_overlap = set_overlap)
}

match_results <- check_gene_overlap(expression_data, gene_sets)

if (mean(match_results$set_overlap) < 0.5) {
  cat("\n警告: 平均基因匹配度低于50%，分析结果可能不可靠！\n")
} else if (min(match_results$set_overlap) < 0.1) {
  cat("\n警告: 某些基因集匹配度极低，考虑移除或修正这些基因集\n")
}

# ===============================
# 5. 运行 ssGSEA（兼容新旧GSVA）
# ===============================
cat("\n===== 步骤5: 使用GSVA包进行ssGSEA分析 =====\n")
start_time <- Sys.time()

if (exists("ssgseaParam")) {
  param <- ssgseaParam(
    exprData = expression_data,
    geneSets = gene_sets,
    alpha = 0.25,
    normalize = TRUE
  )
  ssgsea_results <- gsva(param, verbose = TRUE)
} else {
  ssgsea_results <- gsva(
    expr = expression_data,
    gset.idx.list = gene_sets,
    method = "ssgsea",
    ssgsea.norm = TRUE,
    verbose = TRUE
  )
}

end_time <- Sys.time()
duration <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 1)
cat("\nssGSEA分析成功完成! 耗时:", duration, "分钟\n")

# ===============================
# 6. 保存结果
# ===============================
cat("\n===== 步骤6: 保存结果 =====\n")
output_file <- file.path(output_path, "ssgsea_immune_enrichment_scores.csv")
write.csv(ssgsea_results, output_file, row.names = TRUE)
cat("富集分数矩阵已保存到:", output_file, "\n")

rdata_file <- file.path(output_path, "ssgsea_results.RData")
save(ssgsea_results, file = rdata_file)
cat("完整结果已保存到:", rdata_file, "\n")

# ===============================
# 7. 整理行名（作图与后续分析用的得分矩阵）
# ===============================
cell_names <- sub("__(Adaptive|Innate)$", "", rownames(ssgsea_results))
keep_suffix <- rownames(ssgsea_results) %in% c("Eosinophil__Innate", "Neutrophil__Innate")
cell_names[keep_suffix] <- rownames(ssgsea_results)[keep_suffix]
ssgsea_scores_28 <- ssgsea_results
rownames(ssgsea_scores_28) <- cell_names
write.csv(ssgsea_scores_28, file.path(output_path, "ssgsea_scores_28cells.csv"), row.names = TRUE)
cat("行名整理后的得分矩阵已保存到:", file.path(output_path, "ssgsea_scores_28cells.csv"), "\n")

cat("\n===== 分析成功完成! =====\n")
