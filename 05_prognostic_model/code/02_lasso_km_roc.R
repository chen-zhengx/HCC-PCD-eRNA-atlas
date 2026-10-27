# =========================================================
# 对应结果：Figure 3A–C（失巢凋亡 anoikis）、Figure S1A–C（凋亡 apoptosis）、Figure S2A–C（铁死亡 ferroptosis）
#   A：LASSO-Cox 交叉验证曲线（1_Lambda_Selection.pdf）
#   B：高/低风险组 Kaplan–Meier 曲线（7_KM_Curve.pdf）
#   C：1/2/3 年时间依赖 ROC 曲线及 AUC（8_Time_ROC.pdf）
#   另在屏幕打印 λmin / λ1se 及入选 eRNA 数（λmin = 0.0811 / 0.115 / 0.107，见 Table S7 表注）
# 输入：input/<pcd>_eRNA_RPM_tumor70.csv（PCD 相关 eRNA × 70 个肿瘤样本的 RPM）
#       ../data/clinical/GSE144269_survival.csv（Sample, months, status）
# 输出：output/<pcd>/1_Lambda_Selection.pdf、7_KM_Curve.pdf、8_Time_ROC.pdf（pcd = anoikis / apoptosis / ferroptosis）
# 运行：在模块目录下 Rscript code/02_lasso_km_roc.R
# =========================================================

# ====== 用户配置区（请先修改这里）======
use_english <- TRUE  # TRUE = 英文输出；FALSE = 中文输出

# 加载所需包（自动安装缺失包）
required_packages <- c("survival", "glmnet", "ggplot2", "survminer", "timeROC", "Cairo", "svglite")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
  library(pkg, character.only = TRUE)
}

# 设置路径
pcd_list <- c(anoikis = "Figure 3", apoptosis = "Figure S1", ferroptosis = "Figure S2")

cat("🚀 开始 eRNA 预后模型分析...\n")
if (use_english) {
  cat("📁 Language: English\n")
} else {
  cat("📁 语言：中文\n")
}

# ====== 语言标签定义 ======
if (use_english) {
  risk_dist_title <- "Risk Score Distribution"
  x_patients <- "Patients"
  y_risk_score <- "Risk Score"

  surv_status_title <- "Survival Status"
  y_surv_time <- "Survival Time (months)"

  km_title <- "Kaplan-Meier Survival Curve"
  roc_title <- "Time-dependent ROC Curves"
  forest_title <- "Forest Plot of Multivariate Cox Regression"
  forest_subtitle <- function(n) paste("Based on", n, "selected eRNAs")
  forest_xlab <- "Hazard Ratio (HR, log scale)"
  forest_caption <- "*P<0.05, **P<0.01, ***P<0.001"
  protective_label <- "Protective"
  risky_label <- "Risky"

  status_labels <- c(Deceased = "Deceased", Alive = "Alive")
  risk_groups <- c(High = "High", Low = "Low")

  output_msg <- "✅ Figure saved to:"
} else {
  risk_dist_title <- "风险评分分布"
  x_patients <- "患者"
  y_risk_score <- "风险评分"

  surv_status_title <- "生存状态"
  y_surv_time <- "生存时间（月）"

  km_title <- "Kaplan-Meier 生存曲线"
  roc_title <- "时间依赖性 ROC 曲线"
  forest_title <- "多因素 Cox 回归森林图"
  forest_subtitle <- function(n) paste(n, "个关键 eRNA 的森林图")
  forest_xlab <- "风险比 HR（对数尺度）"
  forest_caption <- "*P<0.05, **P<0.01, ***P<0.001"
  protective_label <- "保护性"
  risky_label <- "危险性"

  status_labels <- c(Deceased = "死亡", Alive = "存活")
  risk_groups <- c(High = "高风险", Low = "低风险")

  output_msg <- "✅ 图形已成功保存至:"
}

# 复制文件到临时目录

# 设置编码
Sys.setlocale("LC_ALL", "en_US.UTF-8")
options(encoding = "UTF-8")

for (pcd in names(pcd_list)) {
cat("\n====== ", pcd, " (", pcd_list[[pcd]], ") ======\n", sep = "")
input_dir <- "input"
output_dir <- file.path("output", pcd)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ====== 主分析流程 ======
tryCatch({

  # ====== 数据读取与预处理 ======
  cat("\n📖 读取数据文件...\n")

  exp_file <- file.path(input_dir, paste0(pcd, "_eRNA_RPM_tumor70.csv"))
  surv_file <- "../data/clinical/GSE144269_survival.csv"

  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!file.exists(exp_file)) {
    exp_patterns <- c("eRNA", "表达", "matrix", "Matrix")
    for (p in exp_patterns) {
      m <- grep(p, basename(csv_files), ignore.case = TRUE)
      if (length(m) > 0) {
        exp_file <- csv_files[m[1]]
        break
      }
    }
  }
  if (!file.exists(surv_file)) {
    surv_patterns <- c("生存", "survival", "status", "clinical")
    for (p in surv_patterns) {
      m <- grep(p, basename(csv_files), ignore.case = TRUE)
      if (length(m) > 0) {
        surv_file <- csv_files[m[1]]
        break
      }
    }
  }

  if (!file.exists(exp_file)) stop("❌ 表达文件未找到")
  if (!file.exists(surv_file)) stop("❌ 生存文件未找到")

  exp_data <- read.csv(exp_file, row.names = 1, check.names = FALSE, fileEncoding = "UTF-8")
  surv_data <- read.csv(surv_file, fileEncoding = "UTF-8")

  # 列名标准化
  required_cols <- c("Sample", "months", "status")
  col_map <- list(
    Sample = c("sample", "id", "ID", "patient", "Patient"),
    months = c("time", "Time", "OS.time", "survival.time", "month", "Month"),
    status = c("event", "Event", "OS", "survival.status", "death", "Death")
  )

  for (col in required_cols) {
    if (!(col %in% colnames(surv_data))) {
      for (cand in col_map[[col]]) {
        if (cand %in% colnames(surv_data)) {
          surv_data[[col]] <- surv_data[[cand]]
          break
        }
      }
    }
  }
  if (!all(required_cols %in% colnames(surv_data))) stop("❌ 生存数据列不全")

  # 数据处理
  exp_matrix <- t(exp_data)
  exp_matrix <- log2(exp_matrix + 1)

  winsorize_vec <- function(x, lower_p = 0.01, upper_p = 0.99) {
    lower <- quantile(x, lower_p, na.rm = TRUE)
    upper <- quantile(x, upper_p, na.rm = TRUE)
    x[x < lower] <- lower
    x[x > upper] <- upper
    x
  }
  exp_matrix <- apply(exp_matrix, 2, winsorize_vec)

  # 样本匹配
  common_samples <- intersect(rownames(exp_matrix), as.character(surv_data$Sample))
  if (length(common_samples) == 0) {
    exp_clean <- gsub("[[:space:]]", "", rownames(exp_matrix))
    surv_clean <- gsub("[[:space:]]", "", as.character(surv_data$Sample))
    common_samples <- intersect(exp_clean, surv_clean)
    if (length(common_samples) == 0) stop("❌ 无共同样本")
    rownames(exp_matrix) <- exp_clean
    surv_data$Sample <- surv_clean
  }

  exp_matrix <- exp_matrix[common_samples, , drop = FALSE]
  surv_data <- surv_data[match(common_samples, surv_data$Sample), ]

  # 移除缺失
  valid_rows <- complete.cases(as.data.frame(exp_matrix)) &
    complete.cases(surv_data$months, surv_data$status)

  x <- exp_matrix[valid_rows, , drop = FALSE]
  surv_data <- surv_data[valid_rows, ]
  y <- Surv(surv_data$months, surv_data$status)

  cat("✅ 最终样本量:", nrow(x), "个样本；eRNA数量:", ncol(x), "\n")

  # ====== 单变量 Cox 筛选 ======
  p_values <- numeric(ncol(x))
  names(p_values) <- colnames(x)

  for (i in seq_len(ncol(x))) {
    fit <- tryCatch(coxph(y ~ x[, i]), error = function(e) NULL)
    p_values[i] <- if (!is.null(fit)) summary(fit)$coefficients[1, 5] else 1
  }

  keep_genes <- names(which(p_values < 0.2))
  if (length(keep_genes) == 0) keep_genes <- colnames(x)
  x_filtered <- x[, keep_genes, drop = FALSE]

  # ====== LASSO-Cox ======
  set.seed(123)
  lambda_seq <- exp(seq(log(0.001), log(1), length.out = 100))
  cv_fit <- cv.glmnet(
    x_filtered, y,
    family = "cox",
    alpha = 1,
    nfolds = min(10, nrow(x_filtered)),
    standardize = TRUE,
    lambda = lambda_seq
  )

  # 统一使用 lambda.1se，若无特征则用 lambda.min
  coef_1se <- coef(cv_fit, s = "lambda.1se")
  selected_idx <- which(as.vector(coef_1se) != 0)

  fallback_mode <- FALSE

  if (length(selected_idx) == 0) {
    cat("⚠️ lambda.1se 未选中特征，改用 lambda.min\n")
    coef_use <- coef(cv_fit, s = "lambda.min")
    s_use <- "lambda.min"
  } else {
    coef_use <- coef_1se
    s_use <- "lambda.1se"
  }

  selected_features <- rownames(coef_use)[which(as.vector(coef_use) != 0)]

  if (length(selected_features) == 0) {
    top_genes <- names(sort(p_values))[1:min(5, length(p_values))]
    selected_features <- top_genes
    coef_use <- matrix(1, nrow = length(top_genes), ncol = 1)
    rownames(coef_use) <- top_genes
    s_use <- "top5_by_pvalue"
    fallback_mode <- TRUE
  }

  # 保存 Lambda 信息
  lambda_info <- data.frame(
    lambda_min = cv_fit$lambda.min,
    lambda_1se = cv_fit$lambda.1se,
    n_features_lambda_min = sum(coef(cv_fit, s = "lambda.min") != 0),
    n_features_lambda_1se = sum(coef(cv_fit, s = "lambda.1se") != 0),
    lambda_used = s_use,
    n_features_used = length(selected_features)
  )

  # write.csv(lambda_info, file.path(output_dir, "0_Lambda_Selection_Details.csv"), row.names = FALSE)
  # file.copy(
  #   file.path(output_dir, "0_Lambda_Selection_Details.csv"),
  #   file.path(output_dir, "0_Lambda_Selection_Details.csv"),
  #   overwrite = TRUE
  # )

  cat("\n🔍 Lambda 选择详情:\n")
  print(lambda_info)

  # 绘制 Lambda 图
  cairo_pdf(file.path(output_dir, "1_Lambda_Selection.pdf"), width = 8, height = 6,
            symbolfamily = cairoSymbolFont("Arial", usePUA = FALSE))               #   使 λ 与其他文字同为 Arial 字体

  plot(
    cv_fit,
    main = "",
    xlab = "−log(λ)",  # plot.cv.glmnet 默认 sign.lambda = -1，横轴实际为 -log(λ)
    ylab = ifelse(use_english, "Partial Likelihood Deviance", "部分似然偏差")
  )

  lambda_sel <- if (s_use == "lambda.min") cv_fit$lambda.min else cv_fit$lambda.1se
  abline(v = -log(lambda_sel), col = "red", lty = 2, lwd = 1.5)
  usr <- par("usr")                                # 图例放在两条竖线之间，不遮挡竖线
  legend(x = min(-log(cv_fit$lambda.min), -log(cv_fit$lambda.1se)) + 0.03 * diff(usr[1:2]),
         y = usr[4] - 0.02 * diff(usr[3:4]),
         legend = expression(lambda[min], lambda["1se"]),
         col = c(ifelse(s_use == "lambda.min", "red", "black"), ifelse(s_use == "lambda.1se", "red", "black")),
         lty = c(ifelse(s_use == "lambda.min", 2, 3), ifelse(s_use == "lambda.1se", 2, 3)), bty = "n")

  dev.off()

  # 计算风险评分
  if (!fallback_mode) {
    risk_score_raw <- predict(cv_fit, newx = x_filtered, s = s_use, type = "link")
  } else {
    risk_score_raw <- rowSums(x_filtered[, selected_features, drop = FALSE])
  }

  risk_score <- as.numeric(scale(risk_score_raw))
  risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")

  # 输出数据框
  surv_output <- data.frame(
    Sample = rownames(x),
    RawRiskScore = as.numeric(risk_score_raw),
    ScaledRiskScore = risk_score,
    RiskGroup = risk_group,
    Months = surv_data$months,
    Status = surv_data$status,
    stringsAsFactors = FALSE
  )

  # 保存结果
  coef_vector <- as.vector(coef_use)
  names(coef_vector) <- rownames(coef_use)
  coef_selected <- coef_vector[selected_features]
  coef_df <- data.frame(eRNA = selected_features, Coefficient = coef_selected, row.names = NULL)

  # write.csv(coef_df, file.path(output_dir, "2_Selected_Features.csv"), row.names = FALSE)
  # write.csv(surv_output, file.path(output_dir, "3_Patient_Risk_Score.csv"), row.names = FALSE)
  # write.csv(x_filtered[, selected_features, drop = FALSE], file.path(output_dir, "4_Key_eRNA_Expression.csv"), row.names = TRUE)

  # ====== 可视化 ======
  surv_output <- surv_output[order(surv_output$ScaledRiskScore), ]
  surv_output$x_id <- seq_len(nrow(surv_output))
  median_score <- median(surv_output$ScaledRiskScore)
  low_count <- sum(surv_output$RiskGroup == "Low")

  # 1. 风险评分分布
  p1 <- ggplot(surv_output, aes(x = x_id, y = ScaledRiskScore, fill = RiskGroup)) +
    geom_bar(stat = "identity", width = 0.8, alpha = 0.8) +
    geom_hline(yintercept = median_score, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("rect", xmin = -Inf, xmax = low_count + 0.5, ymin = -Inf, ymax = Inf, alpha = 0.05, fill = "#4DBBD5") +
    annotate("rect", xmin = low_count + 0.5, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.05, fill = "#E64B35") +
    scale_fill_manual(values = c("High" = "#E64B35", "Low" = "#4DBBD5")) +
    labs(x = x_patients, y = y_risk_score, title = risk_dist_title) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  # ggsave(
  #   file.path(output_dir, "5_Risk_Score_Distribution.pdf"),
  #   p1,
  #   device = cairo_pdf,
  #   width = 10,
  #   height = 6,
  #   dpi = 300
  # )

  # 2. 生存状态图
  surv_output$Status_Label <- ifelse(surv_output$Status == 1, "Deceased", "Alive")
  p2 <- ggplot(surv_output, aes(x = x_id, y = Months, color = Status_Label, shape = RiskGroup)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_vline(xintercept = low_count + 0.5, linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_color_manual(values = c("Deceased" = "#E64B35", "Alive" = "#4DBBD5")) +
    scale_shape_manual(values = c("High" = 17, "Low" = 16)) +
    labs(x = x_patients, y = y_surv_time, title = surv_status_title) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  # ggsave(
  #   file.path(output_dir, "6_Survival_Status.pdf"),
  #   p2,
  #   device = cairo_pdf,
  #   width = 10,
  #   height = 6,
  #   dpi = 300
  # )

  # 3. KM 曲线
  surv_output$RiskGroup <- factor(surv_output$RiskGroup, levels = c("High", "Low"))
  fit_km <- survfit(Surv(Months, Status) ~ RiskGroup, data = surv_output)

  if (!is.null(fit_km$strata)) {
    names(fit_km$strata) <- c("High Group", "Low Group")
  }

  p_km <- ggsurvplot(
    fit_km,
    data = surv_output,
    pval = TRUE,
    risk.table = TRUE,
    palette = c("#E64B35", "#4DBBD5"),
    legend.title = "Risk group",
    legend.labs = c("High Group", "Low Group"),
    xlab = "Time (months)"
  )

  # 保存 PDF
  pdf(file.path(output_dir, "7_KM_Curve.pdf"), width = 7, height = 6)
  print(p_km)
  dev.off()

  # 保存 SVG
  # svglite::svglite(file.path(output_dir, "7_KM_Curve.svg"), width = 7, height = 6)
  # print(p_km)
  # dev.off()

  # 4. 时间依赖 ROC
  if (sum(surv_output$Status == 1) >= 5) {
    max_followup <- max(surv_output$Months, na.rm = TRUE)
    times_try <- c(12, 24, 36)
    times_use <- times_try[times_try < max_followup]

    if (length(times_use) > 0) {
      roc <- timeROC(
        T = surv_output$Months,
        delta = surv_output$Status,
        marker = surv_output$RawRiskScore,
        cause = 1,
        times = times_use,
        iid = TRUE
      )

      auc_df <- data.frame(
        Time_months = times_use,
        Time_years = times_use / 12,
        AUC = roc$AUC
      )

      roc_cols <- c("#B24D00", "#005A8C", "#007A5E")
      roc_cols_use <- roc_cols[seq_along(times_use)]

      draw_time_roc <- function() {
        plot(
          roc,
          time = times_use[1],
          col = roc_cols_use[1],
          lwd = 1,
          title = FALSE,
          xlab = "1-Specificity",
          ylab = "Sensitivity"
        )

        if (length(times_use) > 1) {
          for (i in 2:length(times_use)) {
            plot(
              roc,
              time = times_use[i],
              col = roc_cols_use[i],
              lwd = 1,
              add = TRUE
            )
          }
        }

        abline(0, 1, lty = 2, col = "grey70", lwd = 1.5)
        title(main = "Time-dependent ROC")

        legend_text <- paste0(
          auc_df$Time_years, "-year AUC = ",
          formatC(auc_df$AUC, digits = 3, format = "f")
        )

        legend(
          "bottomright",
          legend = legend_text,
          col = roc_cols_use,
          lty = 1,
          lwd = 1,
          bty = "n",
          cex = 0.95,
          text.col = roc_cols_use,
          seg.len = 2.2,
          y.intersp = 1.1
        )
      }

      Cairo::CairoPDF(
        file.path(output_dir, "8_Time_ROC.pdf"),
        width = 7,
        height = 6
      )
      draw_time_roc()
      dev.off()

      # write.csv(auc_df, file.path(output_dir, "9_ROC_AUC.csv"), row.names = FALSE)
    }
  }

  # ====== 森林图 ======
  # forest_dir <- file.path(output_dir, "森林图")
  # dir.create(forest_dir, recursive = TRUE, showWarnings = FALSE)
  #
  # if (exists("selected_features") && length(selected_features) > 0) {
  #   tryCatch({
  #     cox_data <- cbind(
  #       surv_output[, c("Months", "Status")],
  #       x_filtered[, selected_features, drop = FALSE]
  #     )
  #
  #     formula_str <- paste("Surv(Months, Status) ~", paste(selected_features, collapse = " + "))
  #     multi_cox <- coxph(as.formula(formula_str), data = cox_data)
  #
  #     cox_summary <- summary(multi_cox)
  #     cox_results <- data.frame(
  #       eRNA = rownames(cox_summary$coefficients),
  #       HR = exp(cox_summary$coefficients[, "coef"]),
  #       lower_95 = exp(cox_summary$coefficients[, "coef"] - 1.96 * cox_summary$coefficients[, "se(coef)"]),
  #       upper_95 = exp(cox_summary$coefficients[, "coef"] + 1.96 * cox_summary$coefficients[, "se(coef)"]),
  #       p_value = cox_summary$coefficients[, "Pr(>|z|)"],
  #       stringsAsFactors = FALSE
  #     )
  #
  #     cox_results$significance <- ifelse(
  #       cox_results$p_value < 0.001, "***",
  #       ifelse(cox_results$p_value < 0.01, "**",
  #              ifelse(cox_results$p_value < 0.05, "*", "NS"))
  #     )
  #
  #     # write.csv(cox_results, file.path(forest_dir, "10_Multivariate_Cox_Results.csv"), row.names = FALSE)
  #
  #     forest_data <- cox_results[order(cox_results$HR), ]
  #     forest_data$eRNA <- factor(forest_data$eRNA, levels = forest_data$eRNA)
  #
  #     p_forest <- ggplot(forest_data, aes(x = HR, y = eRNA)) +
  #       geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  #       geom_errorbar(aes(xmin = lower_95, xmax = upper_95), height = 0.2, linewidth = 0.8, color = "gray50") +
  #       geom_point(aes(color = significance), size = 3) +
  #       scale_color_manual(values = c("***" = "#E64B35", "**" = "#F39B7F", "*" = "#8491B4", "NS" = "gray70")) +
  #       scale_x_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10)) +
  #       labs(
  #         x = forest_xlab,
  #         y = "eRNA",
  #         title = forest_title,
  #         subtitle = forest_subtitle(length(selected_features)),
  #         color = "Significance",
  #         caption = forest_caption
  #       ) +
  #       theme_minimal() +
  #       theme(
  #         panel.grid.major = element_line(color = "gray90"),
  #         axis.text = element_text(size = 10),
  #         axis.title = element_text(size = 12, face = "bold"),
  #         plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  #         plot.subtitle = element_text(size = 12, hjust = 0.5),
  #         legend.position = "top"
  #       ) +
  #       annotate("text", x = 0.5, y = nrow(forest_data) + 0.5, label = protective_label, hjust = 0.5, color = "blue", size = 4, fontface = "bold") +
  #       annotate("text", x = 2, y = nrow(forest_data) + 0.5, label = risky_label, hjust = 0.5, color = "red", size = 4, fontface = "bold")
  #
  #     # ggsave(
  #     #   file.path(forest_dir, "12_Custom_Forest_Plot.pdf"),
  #     #   p_forest,
  #     #   device = cairo_pdf,
  #     #   width = 10,
  #     #   height = 8,
  #     #   dpi = 300
  #     # )
  #
  #   }, error = function(e) {
  #     cat("⚠️ 多因素 Cox 失败，使用 LASSO 系数绘制森林图...\n")
  #
  #     lasso_forest <- coef_df
  #     lasso_forest$HR <- exp(lasso_forest$Coefficient)
  #     lasso_forest$lower_95 <- exp(lasso_forest$Coefficient - 1.96 * abs(lasso_forest$Coefficient) / 2)
  #     lasso_forest$upper_95 <- exp(lasso_forest$Coefficient + 1.96 * abs(lasso_forest$Coefficient) / 2)
  #     lasso_forest$p_value <- ifelse(abs(lasso_forest$Coefficient) > 0.1, 0.05, 0.5)
  #
  #     # write.csv(lasso_forest, file.path(forest_dir, "10_LASSO_Forest_Results.csv"), row.names = FALSE)
  #
  #     lasso_forest <- lasso_forest[order(lasso_forest$HR), ]
  #     lasso_forest$eRNA <- factor(lasso_forest$eRNA, levels = lasso_forest$eRNA)
  #
  #     p_lasso <- ggplot(lasso_forest, aes(x = HR, y = eRNA)) +
  #       geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  #       geom_errorbar(aes(xmin = lower_95, xmax = upper_95), height = 0.2, color = "gray50") +
  #       geom_point(aes(color = Coefficient > 0), size = 3) +
  #       scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5")) +
  #       scale_x_log10() +
  #       labs(
  #         x = "Approximate HR (from LASSO)",
  #         y = "eRNA",
  #         title = ifelse(use_english, "Forest Plot Based on LASSO Coefficients", "基于 LASSO 系数的森林图")
  #       ) +
  #       theme_minimal() +
  #       theme(legend.position = "none")
  #
  #     # ggsave(
  #     #   file.path(forest_dir, "14_LASSO_Forest_Plot.pdf"),
  #     #   p_lasso,
  #     #   device = cairo_pdf,
  #     #   width = 10,
  #     #   height = 8,
  #     #   dpi = 300
  #     # )
  #   })
  # }


  cat("\n🎉 分析完成！\n")
  cat(output_msg, "\n")
  cat(file.path(output_dir, "1_Lambda_Selection.pdf"), "\n")
  cat(file.path(output_dir, "7_KM_Curve.pdf"), "\n")
  cat(file.path(output_dir, "8_Time_ROC.pdf"), "\n")

}, error = function(e) {
  cat("\n❌ 错误:", e$message, "\n")
})

}  # PCD 循环结束

# 清理提示
