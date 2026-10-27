# =========================================================
# 对应结果：Figure 9（A–D，4 个 cis-eRNA-QTL 的肿瘤 eRNA 表达按基因型分组的箱线图）及图注中的各基因型人数
#           source data（各 SNP 的 0/1/2 人数、携带者与非携带者的表达中位数、β 及 95% CI、P、FDR）
#           癌旁组织基因型敏感性分析的对应数值（rs4239538 β = 1.61、P = 8.3e-4、FDR = 0.23；rs10858200、rs149423305
#           在癌旁组织缺失率 14%、20%，未检验）
# 输入：input/genotype_tumor/genotype_matrix.csv、input/genotype_normal/genotype_matrix.csv
#       input/genotype_tumor/erna_expression.csv、input/selected_snps_masked_gt.tsv
#       output/tumor/QTL_results_cis.csv、output/normal/QTL_results_cis.csv（01_matrixeqtl.R 的输出）
#       input/GSE144269_runs.tsv（run、gsm、patient、tissue）
# 输出：output/source_data_genotype_expression.csv、output/source_data_effects.csv
#       output/FigS_eRNA_QTL_genotypes.svg/.pdf/.tiff、output/FigS_eRNA_QTL_genotypes_preview.png
# 运行：在模块目录下先运行 01_matrixeqtl.R（肿瘤、癌旁两套输入），再 Rscript code/02_fig9_genotype_boxplots.R
# =========================================================
# eRNA expression by genotype for the four cis-eRNA-QTLs.
#
#   A-D  tumor eRNA RPM by tumor-derived genotype for the cis-eRNA-QTLs with FDR < 0.05
#
# Statistics are refitted with lm(), which is the MatrixEQTL linear model when no
# covariates are used (missing genotypes set to the SNP mean, as MatrixEQTL does), and
# checked against the MatrixEQTL output. FDR values are taken from the MatrixEQTL
# output (BH within cis tests). The same statistics with adjacent-tissue genotypes are
# written to the source data for the sensitivity analysis described in the text.
#
# Run from the project root:
#   Rscript code/02_fig9_genotype_boxplots.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggbeeswarm)
})

in_dir <- "input"
res_dir <- "output"  # MatrixEQTL 结果由 01_matrixeqtl.R 写到 output/tumor、output/normal
out_dir <- "output"
samples_tsv <- "input/GSE144269_runs.tsv"

# ---- Associations shown -----------------------------------------------------------

assoc <- data.frame(
  panel = c("A", "B", "C", "D"),
  snp = c("rs10858200", "rs4842071", "rs149423305", "rs4239538"),
  erna = c("ENSR00000338188", "ENSR00000338188", "ENSR00000324014", "ENSR00000110509"),
  stringsAsFactors = FALSE
)

# ---- Data -----------------------------------------------------------------------------

read_genotypes <- function(set) {
  g <- read.csv(file.path(in_dir, paste0("genotype_", set), "genotype_matrix.csv"),
                check.names = FALSE, stringsAsFactors = FALSE)
  g[g$ID %in% assoc$snp, ]
}

expr <- read.csv(file.path(in_dir, "genotype_tumor/erna_expression.csv"), check.names = FALSE)
rownames(expr) <- expr$erna_id
runs <- colnames(expr)[-(1:4)]
geno <- list(tumor = read_genotypes("tumor"), adjacent = read_genotypes("normal"))
stopifnot(all(runs %in% colnames(geno$tumor)), all(runs %in% colnames(geno$adjacent)))

samples <- read.delim(samples_tsv, stringsAsFactors = FALSE)
patient_of <- setNames(samples$patient, samples$run)

# MatrixEQTL cis output, used for FDR and to check the refitted statistics
read_meqtl <- function(set) {
  r <- read.delim(file.path(res_dir, set, "QTL_results_cis.csv"), check.names = FALSE)
  setNames(r[, c("SNP", "gene", "beta", "p-value", "FDR")], c("snp", "erna", "beta", "p", "fdr"))
}
meqtl <- rbind(cbind(source = "tumor", read_meqtl("tumor")),
               cbind(source = "adjacent", read_meqtl("normal")))

# Adjacent-tissue missing rate for SNPs absent from the adjacent-tissue set
# (computed from the masked joint genotypes before site filtering)
raw_gt <- read.delim(file.path(in_dir, "selected_snps_masked_gt.tsv"), check.names = FALSE)
normal_runs <- samples$run[samples$tissue == "normal"]
adjacent_missing <- setNames(
  apply(raw_gt[, normal_runs], 1, function(x) mean(!x %in% c("0/0", "0/1", "1/1"))),
  raw_gt$ID
)

# ---- Long table and statistics ----------------------------------------------------

long <- list()
effects <- list()
for (i in seq_len(nrow(assoc))) {
  a <- assoc[i, ]
  y <- unlist(expr[a$erna, runs])
  for (src in names(geno)) {
    row <- geno[[src]][geno[[src]]$ID == a$snp, ]
    if (nrow(row) == 0) {
      effects[[length(effects) + 1]] <- data.frame(
        panel = a$panel, snp = a$snp, erna = a$erna, source = src, tested = FALSE,
        n_0 = NA, n_1 = NA, n_2 = NA, median_noncarrier = NA, median_carrier = NA,
        beta = NA, beta_lo = NA, beta_hi = NA, p = NA, fdr = NA,
        adjacent_missing = adjacent_missing[[a$snp]])
      next
    }
    d <- as.numeric(unlist(row[, runs]))
    ok <- !is.na(d)
    # MatrixEQTL replaces missing genotypes with the SNP mean before fitting
    d_imputed <- ifelse(ok, d, mean(d, na.rm = TRUE))
    fit <- lm(y ~ d_imputed)
    est <- coef(summary(fit))[2, ]
    ci <- confint(fit)[2, ]
    m <- meqtl[meqtl$source == src & meqtl$snp == a$snp & meqtl$erna == a$erna, ]
    stopifnot(nrow(m) == 1,
              abs(est[["Estimate"]] - m$beta) < 1e-6 * abs(m$beta),
              abs(est[["Pr(>|t|)"]] - m$p) < 1e-6 * m$p)
    counts <- tabulate(d[ok] + 1, 3)
    effects[[length(effects) + 1]] <- data.frame(
      panel = a$panel, snp = a$snp, erna = a$erna, source = src, tested = TRUE,
      n_0 = counts[1], n_1 = counts[2], n_2 = counts[3],
      median_noncarrier = median(y[ok & d == 0]), median_carrier = median(y[ok & d > 0]),
      beta = est[["Estimate"]], beta_lo = ci[[1]], beta_hi = ci[[2]], p = m$p, fdr = m$fdr,
      adjacent_missing = adjacent_missing[[a$snp]])
    alleles <- c(row$REF, row$ALT)
    gt_label <- c(strrep(alleles[1], 2), paste0(alleles[1], alleles[2]), strrep(alleles[2], 2))
    long[[length(long) + 1]] <- data.frame(
      panel = a$panel, snp = a$snp, erna = a$erna, genotype_source = src,
      patient = patient_of[runs], tumor_run = runs, dosage = d,
      genotype = ifelse(ok, gt_label[d + 1], NA), rpm = y, row.names = NULL)
  }
}
long <- do.call(rbind, long)
effects <- do.call(rbind, effects)

write.csv(long, file.path(out_dir, "source_data_genotype_expression.csv"), row.names = FALSE)
write.csv(effects, file.path(out_dir, "source_data_effects.csv"), row.names = FALSE)

# ---- Style ----------------------------------------------------------------------------

col_dark <- "#272727"
dosage_fill <- c("0" = "#C9DCEE", "1" = "#79A9D6", "2" = "#2F6FAE")

theme_set(
  theme_classic(base_size = 6.5, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(size = 6, colour = "black"),
      axis.title = element_text(size = 6.5),
      plot.title = element_text(size = 6.5, face = "bold", margin = margin(b = 3)),
      plot.title.position = "plot",
      plot.margin = margin(4, 6, 2, 4),
      panel.grid = element_blank()
    )
)

fmt_p <- function(p) {
  if (p >= 0.001) return(sprintf("'%s'", formatC(p, digits = 2, format = "fg", flag = "#")))
  e <- floor(log10(p))
  sprintf("%.1f %%*%% 10^%d", p / 10^e, e)
}
fmt_fdr <- function(q) sprintf("'%s'", formatC(q, digits = 2, format = "fg", flag = "#"))

# ---- Panels A-D: expression by genotype -----------------------------------------------

genotype_panel <- function(a, show_y_title) {
  d <- long[long$panel == a$panel & long$genotype_source == "tumor" & !is.na(long$dosage), ]
  e <- effects[effects$panel == a$panel & effects$source == "tumor", ]
  row <- geno$tumor[geno$tumor$ID == a$snp, ]
  levels_gt <- c(strrep(row$REF, 2), paste0(row$REF, row$ALT), strrep(row$ALT, 2))
  counts <- table(factor(d$dosage, levels = 0:2))
  d$x <- factor(d$dosage, levels = 0:2)

  y_top <- max(d$rpm)
  stat_line <- sprintf("italic(P) == %s*','~~FDR == %s", fmt_p(e$p), fmt_fdr(e$fdr))

  ggplot(d, aes(x, rpm)) +
    geom_boxplot(data = d[d$dosage %in% names(counts)[counts >= 3], ],
                 aes(fill = x), width = 0.55, outlier.shape = NA, linewidth = 0.3,
                 colour = col_dark, alpha = 0.55) +
    geom_quasirandom(aes(fill = x), shape = 21, size = 0.9, stroke = 0.2,
                     colour = col_dark, width = 0.18) +
    annotate("text", x = 0.45, y = y_top * 1.17, label = stat_line, parse = TRUE,
             hjust = 0, vjust = 1, size = 5.6 / .pt, family = "Arial") +
    scale_x_discrete(drop = FALSE, labels = levels_gt) +
    scale_y_continuous(limits = c(0, y_top * 1.18), expand = expansion(mult = c(0.02, 0))) +
    scale_fill_manual(values = dosage_fill, guide = "none") +
    labs(title = sprintf("%s – %s", a$snp, a$erna),
         x = NULL, y = if (show_y_title) "Tumor eRNA expression (RPM)" else NULL)
}

p_geno <- lapply(seq_len(nrow(assoc)), function(i) genotype_panel(assoc[i, ], i == 1))

# ---- Assemble and export ----------------------------------------------------------------

fig <- wrap_plots(p_geno, nrow = 1) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = "Arial"))

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 60, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi,
                 compression = "lzw")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(filename, "_preview.png"), width = w, height = h, units = "in", res = 300)
  print(plot)
  dev.off()
}

save_pub_r(fig, file.path(out_dir, "FigS_eRNA_QTL_genotypes"))
print(effects[, c("panel", "snp", "source", "tested", "n_0", "n_1", "n_2",
                  "median_noncarrier", "median_carrier", "beta", "p", "fdr")], digits = 3)
