# =========================================================
# 对应结果：Figure 8、Figure S8、Figure S9（以 eRNA 为中心的骨架网络：每个风险 eRNA 在每层保留 |ρ| 最大的 2 条边，另加枢纽 PCD
#           基因与 eRNA 之间最强的 2 条边；图中 88、56、83 条边）；网络规模与连接度（失巢凋亡/凋亡/铁死亡 219/169/220 个节点
#           1,837/1,998/3,130 条边；各网络的枢纽 eRNA 与枢纽 PCD 基因，三个网络共同的枢纽 ENSR00000110509、PARP1
#           Results “Multilayer network analysis”）；Supplementary Table 12 的数据（图中所示边与完整网络的边表）
# 输入：input/<PCD>_network_edges.csv、input/<PCD>_network_nodes.csv
# 输出：output/network_summary.csv（网络规模、枢纽节点）、output/full_network_edge_list.csv（去重后的完整网络）
#       output/network_backbone_edges_shown.csv（图中所示边）、output/cytoscape/<PCD>_backbone_*.csv 与 <PCD>_full_*.csv
#       （08 脚本在 Cytoscape 中作图的输入）、output/Figure 8|Supplementary Figure 8|9_network_backbone.svg/.pdf/.tif/.png
#       （ggplot2 版本的骨架网络图；Figure 8、S8、S9 为 08 脚本绘制的 Cytoscape 版本）
# 运行：在模块目录下 Rscript code/07_network_backbone.R
# =========================================================
# =========================================================
# 多层共表达关联网络（Figure 8、S8、S9）的骨架网络图
# 图中为“以 eRNA 为中心的骨架网络”（完整网络见边表）：
#   (1) 每个风险 eRNA 在 5 个层（邻近基因 ±1 Mb、PCD 基因、免疫细胞特征、药物、免疫治疗相关基因）中
#       各保留 |ρ| 最大的 2 条边；
#   (2) 枢纽 PCD 基因（全网络连接度前 3，同连接度按 Σ|ρ| 排序）额外保留其与 eRNA 之间 |ρ| 最大的 2 条边；
#   (3) 节点全部为圆形，颜色 = 节点类型，大小 = 全网络中的连接度；图中只标注节点名称（eRNA 名称在节点上方，
#       其余节点名称在外侧），不放层标题，不单独高亮枢纽节点；图例只列出图中出现的类别；
#   (4) 边颜色 = 正/负相关，边宽 = |ρ|；完整网络以边列表形式另附
# 布局：二分（三栏）布局。风险 eRNA 居中一列；左栏为邻近基因与免疫治疗相关基因，右栏为 PCD 相关基因、
#       免疫细胞特征与候选药物；边为水平 S 形曲线。三张图的层顺序相同；层内节点与 eRNA 用重心法反复排序以减少交叉
# 规格：双栏宽 183 mm，Arial，节点名称 5.5 pt、图例 5.8–6.2 pt；导出可编辑 SVG/PDF 与 600 dpi TIFF，并附图中所示边的源数据
# =========================================================
suppressPackageStartupMessages({ library(ggplot2); library(ggforce) })
src <- "input"
out <- "output"
dir.create(out, showWarnings = FALSE)  # 输出目录可能还不存在
mods <- list(
  list(key = "Anoikis", fig = "Figure 8",
       e = "anoikis_network_edges.csv", n = "anoikis_network_nodes.csv"),
  list(key = "Apoptosis", fig = "Supplementary Figure 8",
       e = "apoptosis_network_edges.csv", n = "apoptosis_network_nodes.csv"),
  list(key = "Ferroptosis", fig = "Supplementary Figure 9",
       e = "ferroptosis_network_edges.csv", n = "ferroptosis_network_nodes.csv"))
layers <- c(Target_Gene = "Nearby genes (\u00b11 Mb)", Death_Gene = "PCD-related genes", Immune_Cell = "Immune cell-type signatures",
            Drug = "Candidate drugs", Immunotherapy_Gene = "Immunotherapy-related genes")
# 三张图使用相同的层顺序（自上而下），便于对照；层内节点与 eRNA 的顺序仍按交叉最少排列
side_layers <- list(left = c("Immunotherapy_Gene", "Target_Gene"), right = c("Drug", "Death_Gene", "Immune_Cell"))
type_lab <- c(eRNA = "Risk eRNA", layers)
# 低饱和、色盲友好的配色
type_col <- c(eRNA = "#3F3F7A", Target_Gene = "#E9A03B", Death_Gene = "#2A9D8F", Immune_Cell = "#CC79A7",
              Drug = "#56B4E9", Immunotherapy_Gene = "#8C6BB1")
edge_col <- c(Positive = "#C8553D", Negative = "#2F6DB5")
FONT <- "Arial"; pt <- ggplot2::.pt                                 # 字号换算：文字大小(mm) = 字号(pt) / .pt
W_MM <- 183; LEGEND_MM <- 50; MM_PER_UNIT <- 3.3                    # 1 个坐标单位 = 1 行 ≈ 3.3 mm
CHAR <- 1.0 / MM_PER_UNIT                                           # 5.5 pt Arial 每个字符约 1.0 mm
MAX_SIZE <- 3.4; BLOCK_GAP <- 1.6                                   # 层间留空
K <- 2; N_HUB <- 3
full_list <- list(); shown_list <- list(); summ <- list()

save_pub_r <- function(plot, filename, width_mm, height_mm, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = FONT); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tif"), width = w, height = h, units = "in", res = dpi, compression = "lzw"); print(plot); dev.off()
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = 200); print(plot); dev.off()
}
disp_name <- function(x) { y <- gsub("_", " ", sub("__Innate$", "", sub("memeory", "memory", x))); y[y == "MDSC"] <- "Myeloid-derived suppressor cell"; y }
# 读入边表并去重：失巢凋亡网络的边表中有 1,316 条完全重复的边（PCD 基因与药物、免疫治疗基因、免疫细胞之间），
# 按无序节点对去重后为 1,837 条；凋亡与铁死亡网络无重复
read_edges <- function(f) {
  e <- read.csv(f, check.names = FALSE)[, 1:3]; names(e) <- c("a", "b", "rho")
  if (basename(f) == "apoptosis_network_edges.csv") {
    add <- read.csv(file.path(dirname(f), "apoptosis_missing_edges_ENSR00000110509.csv"), check.names = FALSE)[, 1:3]
    names(add) <- c("a", "b", "rho"); e <- rbind(e, add)
  }
  key <- ifelse(e$a < e$b, paste(e$a, e$b, sep = "|"), paste(e$b, e$a, sep = "|"))
  e[!duplicated(key), ]
}
perms <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))

# 跨网络共同枢纽：先读入三个网络，按连接度排名（并列取最小名次）
rk <- lapply(mods, function(m) {
  e <- read_edges(file.path(src, m$e)); nd <- read.csv(file.path(src, m$n), check.names = FALSE)
  d <- table(c(e[[1]], e[[2]])); ty <- setNames(nd[[2]], nd[[1]])[names(d)]
  r <- ave(-as.numeric(d), ty, FUN = function(v) rank(v, ties.method = "min")); setNames(r, names(d))
})
rec_erna <- Reduce(intersect, lapply(rk, function(r) names(r)[r <= 3 & grepl("^ENSR|^[0-9]+:", names(r))]))
rec_pcd  <- Reduce(intersect, lapply(seq_along(mods), function(i) {
  nd <- read.csv(file.path(src, mods[[i]]$n), check.names = FALSE); pcd <- nd[[1]][nd[[2]] == "Death_Gene"]
  names(rk[[i]])[rk[[i]] <= 10 & names(rk[[i]]) %in% pcd] }))
rec_hub <- c(rec_erna, rec_pcd); cat("Shared hubs:", rec_hub, "\n")

for (m in mods) {
  e <- read_edges(file.path(src, m$e))
  nd <- read.csv(file.path(src, m$n), check.names = FALSE); names(nd)[1:2] <- c("node", "type")
  ty <- setNames(nd$type, nd$node)
  deg <- table(c(e$a, e$b)); deg <- setNames(as.integer(deg), names(deg))
  # 统一方向：一端为 eRNA 时放在 a 列
  sw <- ty[e$b] == "eRNA" & ty[e$a] != "eRNA"; tmp <- e$a[sw]; e$a[sw] <- e$b[sw]; e$b[sw] <- tmp
  full_list[[m$key]] <- data.frame(Network = m$key, Node1 = e$a, Type1 = type_lab[ty[e$a]], Node2 = e$b,
                                   Type2 = type_lab[ty[e$b]], Spearman_rho = round(e$rho, 3))
  er <- e[ty[e$a] == "eRNA", ]; er$layer <- ty[er$b]
  # (1) 每个 eRNA × 每层取 |ρ| 前 K
  er <- er[order(-abs(er$rho)), ]
  keep <- do.call(rbind, lapply(split(er, list(er$a, er$layer), drop = TRUE), head, K))
  # (2) 枢纽 PCD 基因及其与 eRNA 的最强边
  strength <- tapply(abs(c(e$rho, e$rho)), c(e$a, e$b), sum)
  top_by <- function(v) v[order(-deg[v], -strength[v])][1:N_HUB]
  hub_pcd  <- top_by(nd$node[nd$type == "Death_Gene"])
  de <- deg[names(deg) %in% nd$node[nd$type == "eRNA"]]; hub_erna <- top_by(names(de))
  shown_pcd <- union(hub_pcd, intersect(rec_pcd, nd$node))
  keep <- rbind(keep, do.call(rbind, lapply(shown_pcd, function(h) head(er[er$b == h, ], K))))
  keep <- unique(keep)
  keep$side <- ifelse(keep$layer %in% side_layers$left, "left", "right")

  # ---------- 布局：二分布局 + 交叉最少化 ----------
  erna <- names(sort(de, decreasing = TRUE)); ne <- length(erna)
  partners <- unique(keep[, c("b", "layer")])
  # 给定各栏的层顺序与层内节点顺序，计算每个外栏节点的 y（自上而下，层间留空）
  col_y <- function(order_list) {
    y <- c(); for (s in names(order_list)) {
      cur <- 0
      for (L in order_list[[s]]$layers) { nodes <- order_list[[s]]$nodes[[L]]
        y <- c(y, setNames(-(cur + seq_along(nodes) - 1), nodes)); cur <- cur + length(nodes) + BLOCK_GAP }
      h <- cur - BLOCK_GAP - 1; y[unlist(order_list[[s]]$nodes)] <- y[unlist(order_list[[s]]$nodes)] + h / 2  # 各栏上下居中
    }
    y
  }
  erna_y <- function(ord) { sp <- 2.0; setNames((rev(seq_along(ord)) - (length(ord) + 1) / 2) * sp, ord) }
  crossings <- function(ye, yp) sum(sapply(c("left", "right"), function(s) {
    k <- keep[keep$side == s, ]; a <- ye[k$a]; b <- yp[k$b]
    sum(outer(a, a, "-") * outer(b, b, "-") < 0) / 2 }))
  best <- NULL
  for (lo in list(side_layers$left)) for (ro in list(side_layers$right)) {
    lo2 <- lo[lo %in% partners$layer]; ro2 <- ro[ro %in% partners$layer]
    ord <- list(left = list(layers = lo2, nodes = lapply(setNames(lo2, lo2), function(L) partners$b[partners$layer == L])),
                right = list(layers = ro2, nodes = lapply(setNames(ro2, ro2), function(L) partners$b[partners$layer == L])))
    eo <- erna
    for (it in 1:15) {                                             # 重心法：交替按对侧平均位置排序
      yp <- col_y(ord); ye <- erna_y(eo)
      for (s in c("left", "right")) for (L in ord[[s]]$layers) {
        v <- ord[[s]]$nodes[[L]]; ord[[s]]$nodes[[L]] <- v[order(-sapply(v, function(x) mean(ye[keep$a[keep$b == x]])))] }
      yp <- col_y(ord)
      eo <- eo[order(-sapply(eo, function(x) mean(yp[keep$b[keep$a == x]])))]
    }
    yp <- col_y(ord); ye <- erna_y(eo); cr <- crossings(ye, yp)
    if (is.null(best) || cr < best$cr) best <- list(cr = cr, ord = ord, ye = ye, yp = yp)
  }
  ye <- best$ye; yp <- best$yp; ord <- best$ord

  # 节点坐标：eRNA 在 x = 0，外栏在 x = ±X；X 按外侧标签长度计算，使整图宽度 = 183 mm
  lab_l <- CHAR * max(nchar(disp_name(unlist(ord$left$nodes)))); lab_r <- CHAR * max(nchar(disp_name(unlist(ord$right$nodes))))
  X <- max(10, ((W_MM - LEGEND_MM) / MM_PER_UNIT - (0.75 + lab_l) - (0.75 + lab_r) - 0.6) / 2)
  pos <- data.frame(node = names(ye), x = 0, y = unname(ye), type = "eRNA", side = "center")
  for (s in c("left", "right")) { v <- unlist(ord[[s]]$nodes)
    pos <- rbind(pos, data.frame(node = v, x = ifelse(s == "left", -X, X), y = unname(yp[v]), type = ty[v], side = s)) }
  pos$deg <- deg[pos$node]
  pos$type <- factor(pos$type, levels = names(type_col))
  maxdeg <- max(pos$deg)

  # 边：水平 S 形三次贝塞尔曲线；弱边先画，强边在上
  keep <- keep[order(abs(keep$rho)), ]
  keep$sign <- factor(ifelse(keep$rho > 0, "Positive", "Negative"), levels = c("Positive", "Negative"))
  bez <- do.call(rbind, lapply(seq_len(nrow(keep)), function(i) {
    s <- ifelse(keep$side[i] == "left", -1, 1); y0 <- ye[keep$a[i]]; y1 <- yp[keep$b[i]]
    x0 <- 0; x1 <- s * X; xm <- (x0 + x1) / 2                      # 从节点中心到节点中心（节点画在边之上）
    data.frame(id = i, x = c(x0, xm, xm, x1), y = c(y0, y0, y1, y1), sign = keep$sign[i], arho = abs(keep$rho[i]))
  }))

  # 层的竖向色带（标示各层节点的范围）
  blk <- do.call(rbind, lapply(c("left", "right"), function(s) do.call(rbind, lapply(ord[[s]]$layers, function(L) {
    v <- ord[[s]]$nodes[[L]]; data.frame(layer = L, side = s, top = max(yp[v]), bottom = min(yp[v])) }))))
  blk$x <- ifelse(blk$side == "left", -X, X)
  ylim <- c(min(pos$y) - 1.2, max(pos$y) + 1.2)
  xlim <- c(-(X + 0.75 + lab_l + 0.3), X + 0.75 + lab_r + 0.3)
  # 节点名称：外栏在节点外侧；eRNA 在节点正上方（白色半透明底衬，避免被连线遮挡）
  lab <- pos[pos$type != "eRNA", ]; lab$disp <- disp_name(lab$node)
  lab$lx <- ifelse(lab$side == "left", -(X + 0.75), X + 0.75); lab$hj <- ifelse(lab$side == "left", 1, 0)
  el <- pos[pos$type == "eRNA", ]

  # 图例只列出图中实际出现的类别与数值范围
  types_shown <- names(type_col)[names(type_col) %in% as.character(pos$type)]
  rng <- range(bez$arho)
  rho_br <- c(0.4, 0.6, 0.8); rho_br <- rho_br[rho_br >= rng[1] & rho_br <= rng[2]]
  if (length(rho_br) < 2) rho_br <- signif(seq(rng[1], rng[2], length.out = 3), 2)
  deg_br <- pretty(c(0, maxdeg), 4); deg_br <- deg_br[deg_br >= min(pos$deg) & deg_br <= maxdeg & deg_br > 0]

  node_type_text <- theme(legend.text = element_text(size = 5.8, face = "bold"))   # Node type 图例文字加粗
  g <- ggplot() +
    geom_segment(data = blk, aes(x = x, xend = x, y = top + 0.45, yend = bottom - 0.45, colour = I(type_col[layer])),
                 linewidth = 2.4, alpha = 0.25, lineend = "round") +
    # 边
    geom_bezier(data = bez, aes(x, y, group = id, colour = sign, linewidth = arho), alpha = 0.55, lineend = "round") +
    scale_colour_manual(values = edge_col, name = "Spearman correlation", drop = TRUE,
                        guide = guide_legend(order = 2, override.aes = list(linewidth = 0.9, alpha = 1))) +
    scale_linewidth(range = c(0.15, 1.0), limits = c(0.3, 1), breaks = rho_br, name = "|ρ|",
                    guide = guide_legend(order = 3, override.aes = list(colour = "grey35", alpha = 1))) +
    # 节点：全部为圆形，颜色表示类型，大小表示全网络连接度
    geom_point(data = pos, aes(x, y, fill = type, size = deg), shape = 21, colour = "grey30", stroke = 0.2) +
    geom_text(data = lab, aes(lx, y, label = disp, hjust = hj), size = 5.5 / pt, family = FONT, colour = "grey10") +
    geom_label(data = el, aes(x, y + 0.85, label = node), size = 5.5 / pt, family = FONT, colour = "#2B2B5A",
               fill = scales::alpha("white", 0.85), linewidth = 0, label.padding = unit(0.06, "lines")) +
    scale_fill_manual(values = type_col, labels = type_lab, breaks = types_shown, name = "Node type",
                      guide = guide_legend(order = 1, theme = node_type_text,
                                           override.aes = list(size = 2.6, colour = "grey30", stroke = 0.3))) +
    scale_size_area(max_size = MAX_SIZE, limits = c(0, maxdeg), name = "Degree in the\nfull network", breaks = deg_br,
                    guide = guide_legend(order = 4, override.aes = list(shape = 21, fill = "grey75", colour = "grey30", stroke = 0.3))) +
    coord_equal(xlim = xlim, ylim = ylim, clip = "off", expand = FALSE) +
    theme_void(base_size = 6, base_family = FONT) +
    theme(legend.position = "right", legend.justification = "center",
          legend.title = element_text(size = 6.2, face = "bold"), legend.text = element_text(size = 5.8),
          legend.key.size = unit(3.2, "mm"), legend.spacing.y = unit(1.2, "mm"),
          legend.box.spacing = unit(2, "mm"), plot.margin = margin(2, 2, 2, 2, "mm"),
          plot.background = element_rect(fill = "white", colour = NA))
  # 画布尺寸：每个坐标单位约 3.3 mm，右侧留出图例宽度
  w_mm <- round(diff(xlim) * MM_PER_UNIT + LEGEND_MM); h_mm <- round(diff(ylim) * MM_PER_UNIT + 5)
  save_pub_r(g, file.path(out, paste0(m$fig, "_network_backbone")), w_mm, h_mm)
  cat(sprintf("%s: canvas %d x %d mm; edge crossings %d; node types shown: %s; signs: %s\n", m$fig, w_mm, h_mm, best$cr,
              paste(types_shown, collapse = ", "), paste(unique(as.character(keep$sign)), collapse = ", ")))

  shown_list[[m$key]] <- data.frame(Network = m$key, eRNA = keep$a, Partner = keep$b, Layer = layers[keep$layer],
                                    Spearman_rho = round(keep$rho, 3))
  # 供 Cytoscape 使用的输入：骨架网络的节点坐标与边表、完整网络的节点与边表
  cy_dir <- file.path(out, "cytoscape"); dir.create(cy_dir, showWarnings = FALSE)
  write.csv(data.frame(node = pos$node, label = ifelse(pos$type == "eRNA", pos$node, disp_name(pos$node)),
                       type = as.character(pos$type), side = pos$side, degree = pos$deg, x = pos$x, y = pos$y),
            file.path(cy_dir, paste0(m$key, "_backbone_nodes.csv")), row.names = FALSE)
  write.csv(data.frame(source = keep$a, target = keep$b, rho = keep$rho, abs_rho = abs(keep$rho),
                       sign = as.character(keep$sign), layer = keep$layer),
            file.path(cy_dir, paste0(m$key, "_backbone_edges.csv")), row.names = FALSE)
  write.csv(data.frame(node = nd$node, label = ifelse(nd$type == "eRNA", nd$node, disp_name(nd$node)), type = nd$type,
                       degree = deg[nd$node]), file.path(cy_dir, paste0(m$key, "_full_nodes.csv")), row.names = FALSE)
  write.csv(data.frame(source = e$a, target = e$b, rho = e$rho, abs_rho = abs(e$rho),
                       sign = ifelse(e$rho > 0, "Positive", "Negative")),
            file.path(cy_dir, paste0(m$key, "_full_edges.csv")), row.names = FALSE)
  summ[[m$key]] <- data.frame(network = m$key, full_nodes = nrow(nd), full_edges = nrow(e), shown_nodes = nrow(pos),
                              shown_edges = nrow(keep), hub_eRNAs = paste(hub_erna, collapse = ", "),
                              hub_PCD_genes = paste(hub_pcd, collapse = ", "), shared_hubs = paste(rec_hub, collapse = ", "))
}
summ <- do.call(rbind, summ); print(summ, row.names = FALSE)
write.csv(summ, file.path(out, "network_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, full_list), file.path(out, "full_network_edge_list.csv"), row.names = FALSE)
write.csv(do.call(rbind, shown_list), file.path(out, "network_backbone_edges_shown.csv"), row.names = FALSE)   # 图中所示边的源数据
