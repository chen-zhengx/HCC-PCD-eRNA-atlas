# =========================================================
# 对应结果：GSE144269 70 对肿瘤/癌旁样本的药物敏感性预测（oncoPredict::calcPhenotype，GDSC2 训练数据，198 种药物的
#           预测 IC50；Table S11 及 02、03 的上游；Methods "Drug sensitivity prediction"）
# 输入：input/GDSC2/GDSC2_Expr (RMA Normalized and Log Transformed).rds、input/GDSC2/GDSC2_Res.rds
#       （GDSC2 训练数据，包中不提供，需从 oncoPredict 官方数据仓库 https://osf.io/c6tfx/ 的 DataFiles/Training Data
#        下载后放入 input/GDSC2/）
#       ../data/expression/gene_TPM/GSE144269_gene_TPM.csv.gz（GSE144269 140 个样本的编码基因 TPM）
#       ../data/samples/GSE144269_samples.csv（用于加入配对编号和组织类型）
# 输出：output/oncoPredict/DrugPredictions.csv/.rds（140 个样本 × 198 种药物的预测 IC50）
#       output/oncoPredict/DrugPredictions_annotated.csv（加入 PatientID、Treatment 两列，即 input/drug_predictions_GDSC2_140samples.csv）
#       ./calcPhenotype_Output/（calcPhenotype 默认写在工作目录下的结果）
# 运行：在模块目录下 Rscript code/01_oncoPredict.R（需要先下载 GDSC2 训练数据；耗时较长）
# =========================================================

# 清理环境
rm(list = ls())
options(stringsAsFactors = F)

# 加载必要包
library(oncoPredict)
library(biomaRt)
library(dplyr)

# 设置训练数据目录
dir <- "input/GDSC2/"

# 1. 加载GDSC2表达数据
GDSC2_Expr <- readRDS(file.path(dir, 'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
cat("GDSC2基因数:", nrow(GDSC2_Expr), "\n")

# 2. 加载GDSC2药物反应数据
GDSC2_Res <- readRDS(file.path(dir, "GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res)  # 指数转换
cat("GDSC2药物数:", ncol(GDSC2_Res), "\n")

# 3. 加载并预处理您的肿瘤样本数据
your_data_path <- "../data/expression/gene_TPM/GSE144269_gene_TPM.csv.gz"
your_data <- read.csv(your_data_path, row.names = 1, check.names = FALSE)
cat("\n=== 您的数据信息 ===\n")
cat("原始数据维度:", dim(your_data), "\n")

# 4. 基因标识符转换（如果必要）
# 检查GDSC使用的基因标识符类型
cat("GDSC示例基因:", head(rownames(GDSC2_Expr)), "\n")
cat("您的示例基因:", head(rownames(your_data)), "\n")

# 5. 数据标准化
# 如果数据范围较大（如TPM值），进行log2转换
if (max(your_data, na.rm = TRUE) > 100) {
  cat("进行log2转换...\n")
  your_data <- log2(your_data + 1)
}

# 6. 基因匹配
common_genes <- intersect(rownames(GDSC2_Expr), rownames(your_data))
cat("\n共同基因数:", length(common_genes), "\n")

# 确保数据是矩阵格式
GDSC2_Expr_common <- as.matrix(GDSC2_Expr[common_genes, ])  # 转换为矩阵
testExpr <- as.matrix(your_data[common_genes, ])            # 转换为矩阵

cat("处理后训练数据维度:", dim(GDSC2_Expr_common), "\n")
cat("处理后测试数据维度:", dim(testExpr), "\n")
cat("测试数据类型:", class(testExpr), "\n")  # 添加类型检查

# 7. 创建输出目录
output_dir <- "output/oncoPredict"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 8. 运行药物敏感性预测
cat("\n开始药物敏感性预测...\n")
start_time <- Sys.time()

# 确保testExpr是矩阵
if (!is.matrix(testExpr)) {
  testExpr <- as.matrix(testExpr)
  cat("已将testExpr转换为矩阵\n")
}

predictedPtype <- calcPhenotype(trainingExprData = GDSC2_Expr,
             trainingPtype = GDSC2_Res,
             testExprData = testExpr,
             batchCorrect = 'eb',  #   "eb" for ComBat  
             powerTransformPhenotype = TRUE,
             removeLowVaryingGenes = 0.2,
             minNumSamples = 10, 
             printOutput = TRUE, 
             removeLowVaringGenesFrom = 'rawData' )


end_time <- Sys.time()
cat("预测完成! 耗时:", round(end_time - start_time, 1), "分钟\n")

# 9. 保存结果
saveRDS(predictedPtype, file.path(output_dir, "DrugPredictions.rds"))
write.csv(predictedPtype, file.path(output_dir, "DrugPredictions.csv"))
cat("预测结果已保存至:", output_dir, "\n")

# 10. 输出简要预览
cat("\n预测结果预览:\n")
cat("样本数:", nrow(predictedPtype), "\n")
cat("药物数:", ncol(predictedPtype), "\n")
cat("前5个样本:", head(rownames(predictedPtype)), "\n")
cat("前5种药物:", head(colnames(predictedPtype)), "\n")
cat("敏感性范围:", range(predictedPtype, na.rm = TRUE), "\n")

cat("\n=== 核心预测步骤完成! ===\n")



# 确保正确获取预测结果
if (exists("predictedPtype")) {
  # 保存结果
  saveRDS(predictedPtype, file.path(output_dir, "DrugPredictions.rds"))
  write.csv(predictedPtype, file.path(output_dir, "DrugPredictions.csv"))
  cat("预测结果已保存至:", output_dir, "\n")
} else {
  # 检查默认输出位置
  default_output <- "./calcPhenotype_Output"
  if (dir.exists(default_output)) {
    cat("在默认目录找到预测结果:", default_output, "\n")
    
    # 从默认位置读取结果
    predictedPtype <- readRDS(file.path(default_output, "DrugPredictions.rds"))
    
    # 保存到您指定的目录
    saveRDS(predictedPtype, file.path(output_dir, "DrugPredictions.rds"))
    write.csv(predictedPtype, file.path(output_dir, "DrugPredictions.csv"))
    cat("结果已复制到:", output_dir, "\n")
  } else {
    cat("警告: 未能找到预测结果\n")
    cat("请检查函数文档: ?calcPhenotype\n")
  }
}

# 如果仍然找不到，尝试直接读取默认文件
if (!exists("predictedPtype")) {
  if (file.exists("./calcPhenotype_Output/DrugPredictions.rds")) {
    predictedPtype <- readRDS("./calcPhenotype_Output/DrugPredictions.rds")
    cat("从默认位置加载预测结果\n")
  }
}

# 确保结果存在
if (exists("predictedPtype")) {
  # 输出简要预览
  cat("\n预测结果预览:\n")
  cat("样本数:", nrow(predictedPtype), "\n")
  cat("药物数:", ncol(predictedPtype), "\n")
  cat("前5个样本:", head(rownames(predictedPtype)), "\n")
  cat("前5种药物:", head(colnames(predictedPtype)), "\n")
  cat("敏感性范围:", range(predictedPtype, na.rm = TRUE), "\n")
} else {
  cat("\n!!! 未能获取预测结果 !!!\n")
  cat("请尝试以下解决方案:\n")
  cat("1. 检查默认目录: ", getwd(), "/calcPhenotype_Output\n")
  cat("2. 重新运行预测函数\n")
}




# =========================================================
# 为预测结果加入 PatientID（70 对肿瘤/癌旁依次编号 1–70）和 Treatment（tumor / nontumor）两列
# =========================================================
sample_info <- read.csv("../data/samples/GSE144269_samples.csv", stringsAsFactors = FALSE)
pred <- read.csv(file.path(output_dir, "DrugPredictions.csv"), row.names = 1, check.names = FALSE)
sample_info <- sample_info[match(trimws(rownames(pred)), sample_info$name), ]
stopifnot(!anyNA(sample_info$name))
pred_annot <- data.frame(
  PatientID = as.integer(factor(sample_info$PatientID, levels = unique(sample_info$PatientID))),
  Treatment = sample_info$Treatment,
  pred,
  check.names = FALSE
)
write.csv(pred_annot, file.path(output_dir, "DrugPredictions_annotated.csv"))
cat("整理后的预测结果已保存至:", file.path(output_dir, "DrugPredictions_annotated.csv"), "\n")
