# =========================================================
# 对应结果：Figure 8、Figure S8、Figure S9 与 Supplementary Table 12 的多层共表达关联网络——各风险 eRNA 的网络文件与
#           三个整合网络（整合并去除重复边后，失巢凋亡、凋亡、铁死亡网络分别为 219、169、220 个节点，1,837、1,998、3,130 条边）
# 输入：output/correlations/<PCD>/<eRNA>/ 中 9 类相关的 FDR_0.01 结果（03、05 脚本）、input/risk_eRNAs_by_PCD.csv
# 输出：output/networks/<PCD>/<eRNA>_网络文件.csv、<eRNA>_节点属性.csv；output/networks/<PCD>_network_edges.csv
#       <PCD>_network_nodes.csv（整合网络）
# 运行：在模块目录下 Rscript code/06_build_networks.R
# =========================================================

library(data.table)
library(dplyr)
library(stringr)

risk <- read.csv("input/risk_eRNAs_by_PCD.csv", stringsAsFactors = FALSE)
out_dir <- file.path("output", "networks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. 9 类相关结果：所在子文件夹、文件名（03、05 脚本的输出），第 1、2 列节点的类型
# -----------------------------
layer_files <- function(erna_id, pcd_cn) {
  data.frame(
    dir = c("eRNA_vs_target_genes", "eRNA_vs_PCD_genes", "eRNA_vs_immune_cells", "eRNA_vs_drugs",
            "eRNA_vs_immunotherapy_genes", "target_genes_vs_PCD_genes", "immune_cells_vs_PCD_genes",
            "drugs_vs_PCD_genes", "immunotherapy_genes_vs_PCD_genes"),
    file = c(paste0(erna_id, "_vs_1MB_靶基因_Spearman_FDR_0.01.csv"),
             paste0(erna_id, "_vs_", pcd_cn, "基因_Spearman_FDR_0.01.csv"),
             paste0(erna_id, "表达_vs_免疫细胞_Spearman_FDR_0.01.csv"),
             paste0(erna_id, "表达_vs_敏感药物_Spearman_FDR_0.01.csv"),
             paste0(erna_id, "_vs_免疫基因_Spearman_FDR_0.01.csv"),
             paste0("所有靶基因_vs_", pcd_cn, "基因_Spearman_FDR_0.01.csv"),
             paste0("所有免疫细胞_vs_", pcd_cn, "基因_Spearman_FDR_0.01.csv"),
             paste0("所有药物_vs_", pcd_cn, "基因_Spearman_FDR_0.01.csv"),
             paste0("免疫治疗基因_vs_", pcd_cn, "基因_Spearman_FDR_0.01.csv")),
    type1 = c("eRNA", "eRNA", "eRNA", "eRNA", "eRNA", "Target_Gene", "Immune_Cell", "Drug", "Immunotherapy_Gene"),
    type2 = c("Target_Gene", "Death_Gene", "Immune_Cell", "Drug", "Immunotherapy_Gene",
              "Death_Gene", "Death_Gene", "Death_Gene", "Death_Gene"),
    stringsAsFactors = FALSE
  )
}
# 节点类型的顺序（同一名称出现在多层时，后面的类型优先，与 画图2/代码.txt 中“靶基因 -> 死亡基因 -> eRNA”的赋值顺序一致）
type_levels <- c("eRNA", "Target_Gene", "Death_Gene", "Immune_Cell", "Drug", "Immunotherapy_Gene")
type_priority <- c("Target_Gene", "Immunotherapy_Gene", "Immune_Cell", "Drug", "Death_Gene", "eRNA")

# 简单清洗函数（同 画图2/代码.txt）
clean_vec <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA
  return(x)
}
# 无方向边的键：用于去除重复边（同 07_network_backbone.R 的 read_edges）
edge_key <- function(a, b) ifelse(a < b, paste(a, b, sep = "|"), paste(b, a, sep = "|"))

# -----------------------------
# 2. 每个风险 eRNA 的网络文件与节点属性
# -----------------------------
pcd_edges <- list()
pcd_nodes <- list()
for (r in seq_len(nrow(risk))) {
  pcd <- risk$PCD[r]; pcd_cn <- risk$PCD_cn[r]; erna <- risk$eRNA[r]
  erna_id <- gsub(":", "", erna)
  run_dir <- file.path("output", "correlations", pcd, erna_id)
  lf <- layer_files(erna_id, pcd_cn)

  net_list <- list()
  node_list <- list()
  for (i in seq_len(nrow(lf))) {
    f <- file.path(run_dir, lf$dir[i], lf$file[i])
    if (!file.exists(f)) next                       # 该层没有 FDR < 0.01 的边时不写出文件
    df <- fread(f, check.names = FALSE, data.table = FALSE)
    if (nrow(df) == 0) next
    df <- df[, 1:3]                                  # 只取前三列：节点 1、节点 2、相关系数
    colnames(df) <- c("node1", "node2", "Correlation")
    df$node1 <- clean_vec(df$node1)
    df$node2 <- clean_vec(df$node2)
    df$Correlation <- as.numeric(df$Correlation)
    df <- df[!is.na(df$node1) & !is.na(df$node2) & !is.na(df$Correlation), ]
    # 对 node1 和 node2 去掉末尾 _数字（只影响类似 BMS-345541_1249 的药物名）
    df$node1 <- str_replace(df$node1, "_\\d+$", "")
    df$node2 <- str_replace(df$node2, "_\\d+$", "")
    net_list[[length(net_list) + 1]] <- df
    node_list[[length(node_list) + 1]] <- rbind(
      data.frame(Nodename = df$node1, type = lf$type1[i], stringsAsFactors = FALSE),
      data.frame(Nodename = df$node2, type = lf$type2[i], stringsAsFactors = FALSE))
  }
  network_file <- bind_rows(net_list)

  # 节点属性：按类型优先级赋值，排序：eRNA -> Target_Gene -> Death_Gene -> 其他层
  nodes <- bind_rows(node_list)
  node_attr <- data.frame(Nodename = unique(nodes$Nodename), type = NA_character_, stringsAsFactors = FALSE)
  for (tp in type_priority) {
    node_attr$type[node_attr$Nodename %in% nodes$Nodename[nodes$type == tp]] <- tp
  }
  node_attr$type <- factor(node_attr$type, levels = type_levels)
  node_attr <- node_attr %>% arrange(type, Nodename)
  node_attr$type <- as.character(node_attr$type)

  dir.create(file.path(out_dir, pcd), showWarnings = FALSE)
  fwrite(network_file, file.path(out_dir, pcd, paste0(erna_id, "_网络文件.csv")), bom = TRUE)
  fwrite(node_attr, file.path(out_dir, pcd, paste0(erna_id, "_节点属性.csv")), bom = TRUE)
  cat(pcd, erna, ": 网络文件", nrow(network_file), "行，节点", nrow(node_attr), "个\n")

  pcd_edges[[pcd]] <- rbind(pcd_edges[[pcd]], network_file)
  pcd_nodes[[pcd]] <- rbind(pcd_nodes[[pcd]], node_attr)
}

# -----------------------------
# 3. 整合网络：合并同一 PCD 的各 eRNA 网络，去除重复边
# -----------------------------
for (pcd in names(pcd_edges)) {
  e <- pcd_edges[[pcd]]
  e <- e[!duplicated(edge_key(e$node1, e$node2)), ]
  nd <- pcd_nodes[[pcd]]
  nd <- nd[!duplicated(nd$Nodename), ]
  colnames(nd) <- c("node", "type")
  fwrite(e, file.path(out_dir, paste0(pcd, "_network_edges.csv")), bom = TRUE)
  fwrite(nd, file.path(out_dir, paste0(pcd, "_network_nodes.csv")), bom = TRUE)
}
