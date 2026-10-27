# =========================================================
# 对应结果：eRNA count 矩阵与 RPM 矩阵（../data/expression/eRNA_counts/、../data/expression/eRNA_RPM/）
# 输入：<计数目录>/<SRR>_bd.txt（multicov 输出：chr, start, end, eRNA, type, count）
#       <flagstat 目录>/<SRR>_aligned_sorted.txt（../batch_scripts/03_flagstat_multicov_batch.sh 的输出）
# 输出：<输出前缀>_eRNA_counts.csv、<输出前缀>_eRNA_RPM.csv（行 = eRNA，列 = SRR）
# 运行：Rscript code/01_erna_rpm.R <计数目录> <flagstat 目录> <输出前缀>
# =========================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("用法：Rscript code/01_erna_rpm.R <计数目录> <flagstat 目录> <输出前缀>")
count_dir <- args[1]; flagstat_dir <- args[2]; out_prefix <- args[3]

count_files <- sort(list.files(count_dir, pattern = "_bd\\.txt$", full.names = TRUE))
samples <- sub("_bd\\.txt$", "", basename(count_files))
cat("找到", length(samples), "个样本\n")

# 读取每个样本的 multicov 计数（第 4 列为 eRNA 名称，最后一列为计数）
counts <- NULL
for (i in seq_along(count_files)) {
  bd <- read.delim(count_files[i], header = FALSE, stringsAsFactors = FALSE)
  if (is.null(counts)) {
    counts <- matrix(0, nrow(bd), length(samples), dimnames = list(bd[[4]], samples))
  }
  stopifnot(identical(rownames(counts), bd[[4]]))
  counts[, i] <- bd[[ncol(bd)]]
}

# 比对 reads 数：flagstat 中第一个 “mapped (” 行
mapped <- sapply(samples, function(s) {
  fs <- readLines(file.path(flagstat_dir, paste0(s, "_aligned_sorted.txt")))
  as.numeric(strsplit(fs[grep(" mapped \\(", fs)[1]], " ")[[1]][1])
})
print(data.frame(Sample = samples, Mapped = mapped))

rpm <- sweep(counts, 2, mapped, "/") * 1e6
write.csv(data.frame(name = rownames(counts), counts, check.names = FALSE),
          paste0(out_prefix, "_eRNA_counts.csv"), row.names = FALSE)
write.csv(data.frame(name = rownames(rpm), rpm, check.names = FALSE),
          paste0(out_prefix, "_eRNA_RPM.csv"), row.names = FALSE)
cat("已输出:", paste0(out_prefix, c("_eRNA_counts.csv", "_eRNA_RPM.csv")), "\n")
