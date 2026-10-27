# =========================================================
# 对应结果：Figure 8、Figure S8、Figure S9（在 Cytoscape 中绘制的骨架网络图）
# 输入：output/cytoscape/<PCD>_backbone_nodes.csv、<PCD>_backbone_edges.csv（07_network_backbone.R 的输出）
# 输出：output/cytoscape/Figure 8|Supplementary Figure 8|Supplementary Figure 9_cytoscape_polished.pdf/.svg/.tif/.png
#       <图名>_polished.sif、Cytoscape 会话 Figure 8_S8_S9_cytoscape.cys
# 运行：先打开 Cytoscape（CyREST 端口 1234），在模块目录下 python code/08_cytoscape_networks.py
# =========================================================
# 在 Cytoscape 中绘制三个 PCD 相关风险 eRNA 骨架网络（Figure 8、Supplementary Figure 8、Supplementary Figure 9）
# 数据：07_network_backbone.R 导出的 cytoscape/<网络>_backbone_*.csv（每个风险 eRNA 在每层保留 |ρ| 最大的 2 条边，
#       另加枢纽 PCD 基因与 eRNA 之间最强的 2 条边）；三张图使用同一版式：
#   - 布局：风险 eRNA 排成椭圆环（共享伴侣节点的 eRNA 尽量相邻；彼此共享 ≥3 个伴侣节点的 eRNA 构成核心模块，放在左侧）；
#           只与一个 eRNA 相连的伴侣节点排在外圈，紧挨各自的 eRNA 呈扇形展开，名称沿径向朝外；
#           与多个 eRNA 相连的伴侣节点放在圈内：与核心模块相连的在中央排成一列，其余放在所连 eRNA 的重心附近
#   - eRNA 画成深色圆角矩形，名称写在节点内；其余节点为圆形，颜色 = 节点类型，大小 = 全网络连接度
#   - 边颜色 = 正/负相关，边宽与 |ρ| 成正比；低饱和、色盲友好配色；字体 Arial（PDF 中映射为同字宽的 Helvetica，文字可编辑）
#   - 图例用注释画在网络左上角与右下角（三张图位置一致），只列出图中出现的类别
# 输出：按内容裁边，三张图比例尺相同（节点名称 6 pt），宽度不超过 183 mm；可编辑 PDF/SVG，
#       600 dpi TIFF 与 300 dpi PNG 由 PDF 栅格化得到；Cytoscape 会话另存为 Figure 8_S8_S9_cytoscape.cys
# 需要 Cytoscape 已打开（CyREST 端口 1234）
import itertools, math, os, random, re, sys, time
import fitz
import pandas as pd
import requests
from PIL import Image
from pypdf import PdfReader, PdfWriter, Transformation

sys.stdout.reconfigure(encoding="utf-8")
BASE = "http://127.0.0.1:1234/v1"
OUT = os.path.abspath(os.path.join("output", "cytoscape"))
NETS = [("Anoikis", "Figure 8"), ("Apoptosis", "Supplementary Figure 8"), ("Ferroptosis", "Supplementary Figure 9")]
TYPE_LAB = {"eRNA": "Risk eRNAs", "Target_Gene": "Candidate target genes (\u00b11 Mb)", "Death_Gene": "PCD-related genes",
            "Immune_Cell": "Immune cell-type signatures", "Drug": "Candidate drugs", "Immunotherapy_Gene": "Immunotherapy-related genes"}
TYPE_COL = {"eRNA": "#34416B", "Target_Gene": "#E69F00", "Death_Gene": "#009E73", "Immune_Cell": "#CC79A7",
            "Drug": "#56B4E9", "Immunotherapy_Gene": "#8C6BB1"}
TYPE_ORDER = ["Target_Gene", "Death_Gene", "Immune_Cell", "Drug", "Immunotherapy_Gene"]
EDGE_COL = {"positive": "#C8553D", "negative": "#2F6DB5"}
FONT, FS = "Arial", 30                                               # 节点名称字号（px）
W_LO, W_HI = 1.5, 7.0                                                # 边宽（px），|ρ| 0.4 → 1.0 线性映射
A1, B1 = 640, 470                                                    # 中圈（eRNA）椭圆半轴
A2, B2 = 1020, 790                                                   # 外圈（单连接伴侣节点）椭圆半轴
DELTA, G0 = 5.6, 8.0                                                 # 外圈相邻节点的角度间隔、相邻 eRNA 扇形之间的间隔（度）
COL_X, ROW = -60, 62                                                 # 中央核心枢纽列的横坐标与行距（px）
SHRINK = 0.62                                                        # 外围共享节点向圆心收缩的比例
LABEL_PT, MAX_W_MM = 6.0, 183                                        # 成图中节点名称字号（pt）与最大宽度（mm）
LF, TF, RH = 30, 32, 46                                              # 图例文字、标题字号与行高（px）
STYLE = "PCD_network_polished"


def req(method, path, **kw):
    r = requests.request(method, BASE + path, timeout=180, **kw)
    if r.status_code >= 400:
        raise RuntimeError(f"{method} {path}: {r.status_code} {r.text[:300]}")
    return r.json() if r.text and r.headers.get("content-type", "").startswith("application/json") else r.text


def cmd(ns_cmd, **args):
    ns, c = ns_cmd.split(" ", 1)
    out = req("POST", f"/commands/{ns}/{c.replace(' ', '%20')}", json={k: str(v) for k, v in args.items()})
    if out.get("errors"):
        raise RuntimeError(f"{ns_cmd}: {out['errors']}")
    return out["data"]


def text_w(s, fs=FS, bold=False):                                    # Arial 文字宽度估计（px）
    w = sum(0.667 if c.isupper() else 0.556 if c.isdigit() else 0.278 if c in " .:-()" else 0.5 for c in s)
    return w * fs * (1.06 if bold else 1.0)


def edge_width(r):
    return W_LO + (W_HI - W_LO) * (r - 0.4) / 0.6


def node_size(deg):
    return 16 + 3.2 * math.sqrt(deg)


def on_ellipse(a, b, deg_):                                          # 椭圆上极角 deg_（度，屏幕坐标顺时针）处的点
    t = math.radians(deg_); r = a * b / math.hypot(b * math.cos(t), a * math.sin(t))
    return r * math.cos(t), r * math.sin(t)


def seg_dist(p, a, b):                                               # 点 p 到线段 ab 的距离
    (px, py), (ax, ay), (bx, by) = p, a, b
    L2 = (bx - ax) ** 2 + (by - ay) ** 2
    s = 0 if L2 == 0 else max(0, min(1, ((px - ax) * (bx - ax) + (py - ay) * (by - ay)) / L2))
    return math.hypot(px - ax - s * (bx - ax), py - ay - s * (by - ay))


def overlap(a, b):                                                   # 矩形（左、右、上、下）是否重叠
    return min(a[1], b[1]) > max(a[0], b[0]) and min(a[3], b[3]) > max(a[2], b[2])


def layout(bn, be):
    """返回节点坐标、外观与分组。"""
    ty = dict(zip(bn.node, bn.type)); lab = dict(zip(bn.node, bn.label)); deg = dict(zip(bn.node, bn.degree))
    pairs = [(a, b) if ty[a] == "eRNA" else (b, a) for a, b in zip(be.source, be.target)]
    assert all(ty[a] == "eRNA" and ty[b] != "eRNA" for a, b in pairs)   # 骨架网络中每条边都是 eRNA–伴侣节点
    erna = list(bn.node[bn.type == "eRNA"])
    links = {}
    for a, b in pairs:
        links.setdefault(b, []).append(a)
    shared = {p: v for p, v in links.items() if len(v) > 1}
    leaf_of = {p: v[0] for p, v in links.items() if len(v) == 1}

    # eRNA 在环上的顺序：共享同一伴侣节点的 eRNA 尽量相邻（随机起点 + 两两交换的局部搜索）
    def ring_cost(o):
        k, idx = len(o), {x: i for i, x in enumerate(o)}
        return sum(min(abs(idx[a] - idx[b]), k - abs(idx[a] - idx[b])) ** 1.5
                   for v in shared.values() for a, b in itertools.combinations(v, 2))
    random.seed(1); best = None
    for _ in range(200):
        o = erna[:]; random.shuffle(o); c = ring_cost(o); better = True
        while better:
            better = False
            for i, j in itertools.combinations(range(len(o)), 2):
                o2 = o[:]; o2[i], o2[j] = o2[j], o2[i]; c2 = ring_cost(o2)
                if c2 < c - 1e-9:
                    o, c, better = o2, c2, True
        if best is None or c < best[0]:
            best = (c, o)
    order = best[1]

    # 核心模块：两两共享 ≥3 个伴侣节点的 eRNA 连成的最大连通块
    w2 = {(a, b): sum(a in v and b in v for v in shared.values()) for a, b in itertools.combinations(order, 2)}
    comp = {x: {x} for x in order}
    for (a, b), w in w2.items():
        if w >= 3 and comp[a] is not comp[b]:
            u = comp[a] | comp[b]
            for x in u:
                comp[x] = u
    blocks = {frozenset(s) for s in comp.values()}
    core = max(blocks, key=lambda s: (len(s), sum(w for (a, b), w in w2.items() if a in s and b in s)))

    leaves = {x: sorted([p for p, v in leaf_of.items() if v == x], key=lambda p: (TYPE_ORDER.index(ty[p]), lab[p])) for x in order}
    n_sh = {x: sum(x in v for v in shared.values()) for x in order}
    # 中圈：每个 eRNA 分到的角度 = 其外圈节点数 × DELTA + 固定间隔 G0 + 按共享节点数分配的剩余角度
    base = {x: len(leaves[x]) * DELTA + G0 for x in order}
    extra = (360 - sum(base.values())) / sum(n_sh.values())
    ang, a0 = {}, 0.0
    for x in order:
        w = base[x] + extra * n_sh[x]; ang[x] = a0 + w / 2; a0 += w
    c = math.degrees(math.atan2(sum(math.sin(math.radians(ang[x])) for x in core), sum(math.cos(math.radians(ang[x])) for x in core)))
    ang = {x: a + 180 - c for x, a in ang.items()}                     # 核心模块转到左侧
    pos = {}
    def pill_rect(x):                                                # eRNA 节点（含留白）所占矩形
        w = text_w(lab[x], 30, True) + 34
        return pos[x][0] - w / 2 - 20, pos[x][0] + w / 2 + 20, pos[x][1] - 46, pos[x][1] + 46
    for it in range(400):                                            # 相邻 eRNA 节点重叠时把二者沿环推开
        for x in order:
            pos[x] = on_ellipse(A1, B1, ang[x])
        hit = [(order[i], order[(i + 1) % len(order)]) for i in range(len(order))
               if overlap(pill_rect(order[i]), pill_rect(order[(i + 1) % len(order)]))]
        if not hit:
            break
        for x, y in hit:
            ang[x] -= 0.4; ang[y] += 0.4
    for x in order:                                                  # 外圈：以该 eRNA 的方向为中心呈扇形展开
        k = len(leaves[x])
        for i, p in enumerate(leaves[x]):
            pos[p] = on_ellipse(A2, B2, ang[x] + (i - (k - 1) / 2) * DELTA)
    fan = [(ang[x] - (len(leaves[x]) - 1) / 2 * DELTA, ang[x] + (len(leaves[x]) - 1) / 2 * DELTA) for x in order]
    for i in range(len(order)):                                      # 检查相邻扇形是否挤在一起
        g = (fan[(i + 1) % len(order)][0] - fan[i][1]) % 360
        if g < DELTA:
            print("  fan gap too small:", lab[order[i]], round(g, 1))

    # 圈内：与核心模块相连的共享节点在中央排成一列（按所连 eRNA 的平均纵坐标排序、等间距，名称在右）；
    # 其余以所连 eRNA 的重心向圆心收缩后的位置为目标（名称在左，朝向圆心），从连接多的节点开始，
    # 依次放到离目标最近、且不与已放置节点及名称重叠的位置
    col = [p for p, v in shared.items() if set(v) & core]
    col2 = [p for p in shared if p not in col]
    bary = {p: (sum(pos[x][0] for x in v) / len(v), sum(pos[x][1] for x in v) / len(v)) for p, v in shared.items()}
    ys = sorted(col, key=lambda p: bary[p][1]); y0 = sum(bary[p][1] for p in col) / len(col)
    for i, p in enumerate(ys):
        pos[p] = (COL_X, y0 + (i - (len(ys) - 1) / 2) * ROW)
    def rect(v):
        if ty[v] == "eRNA":
            return pill_rect(v)
        x, y = pos[v]; r = node_size(deg[v]) / 2; w = text_w(lab[v]) + 14
        return (x - r - 8, x + r + w, y - 28, y + 28) if v in col else (x - r - w, x + r + 8, y - 28, y + 28)
    placed = list(col) + order
    def label_rect(v):                                               # 名称文字所占矩形（eRNA 为整个节点）
        if ty[v] == "eRNA":
            return pill_rect(v)[0] + 14, pill_rect(v)[1] - 14, pill_rect(v)[2] + 20, pill_rect(v)[3] - 20
        x, y = pos[v]; r = node_size(deg[v]) / 2; w = text_w(lab[v])
        return (x + r + 6, x + r + 10 + w, y - 16, y + 16) if v in col else (x - r - 10 - w, x - r - 6, y - 16, y + 16)
    def hits(a_, b_, r):                                             # 线段 ab 是否穿过矩形 r（Liang–Barsky 裁剪）
        (x0, y0), (x1, y1) = a_, b_; dx, dy = x1 - x0, y1 - y0; lo, hi = 0.0, 1.0
        for p_, q_ in ((-dx, x0 - r[0]), (dx, r[1] - x0), (-dy, y0 - r[2]), (dy, r[3] - y0)):
            if p_ == 0:
                if q_ < 0: return False
            else:
                s = q_ / p_
                if p_ < 0: lo = max(lo, s)
                else: hi = min(hi, s)
                if lo > hi: return False
        return True
    steps = sorted(((dx, dy) for dx in range(-400, 401, 16) for dy in range(-400, 401, 16)), key=lambda d: d[0] ** 2 + d[1] ** 2)
    for p in sorted(col2, key=lambda p: (-len(shared[p]), lab[p])):
        tx, ty_ = SHRINK * bary[p][0], SHRINK * bary[p][1]
        best_ = None
        for dx, dy in steps:
            x, y = tx + dx, ty_ + dy
            if (x / A1) ** 2 + (y / B1) ** 2 > 0.8:                   # 留在 eRNA 环内
                continue
            pos[p] = (x, y)
            if any(overlap(rect(p), rect(q)) for q in placed):        # 节点与名称不重叠（硬性条件）
                continue
            r = node_size(deg[p]) / 2
            old = [(a_, b_) for a_, b_ in pairs if b_ != p and a_ in pos and b_ in pos]
            mine = [(a_, p) for a_ in shared[p]]
            others = [q for q in placed if q in col or q in col2]
            bad = (10 * sum(seg_dist((x, y), pos[a_], pos[b_]) < r + 10 for a_, b_ in old)      # 节点压在已有连线上
                   + 10 * sum(seg_dist(pos[q], (x, y), pos[a_]) < node_size(deg[q]) / 2 + 8   # 自身连线穿过其他节点
                              for a_, _ in mine for q in others if q != a_)
                   + 3 * sum(hits(pos[a_], pos[b_], label_rect(p)) for a_, b_ in old)          # 已有连线穿过自身名称
                   + 3 * sum(hits(pos[a_], (x, y), label_rect(q))                             # 自身连线穿过其他名称
                             for a_, _ in mine for q in others + order if q != a_))
            score = bad + math.hypot(dx, dy) / 100
            if best_ is None or score < best_[0]:
                best_ = (score, (x, y))
            if bad == 0:
                break
        pos[p] = best_[1]
        if best_[0] >= 1 + 400 * 1.42 / 100:
            print("  compromise placement for", lab[p], round(best_[0], 1))
        placed.append(p)

    look = {}
    for v in bn.node:
        x, y = pos[v]
        if ty[v] == "eRNA":                                          # eRNA：深色圆角矩形，名称写在节点内
            look[v] = dict(w=text_w(lab[v], 30, True) + 34, h=46, rot=0.0, lpos="C,C,c,0.00,0.00")
        elif v in leaf_of:                                           # 外圈：名称沿径向朝外，保持正向可读
            a_ = math.degrees(math.atan2(y, x)); left = x < 0        # 屏幕坐标（y 向下），旋转角顺时针为正
            look[v] = dict(w=node_size(deg[v]), h=node_size(deg[v]), rot=(a_ + 180) % 360 if left else a_,
                           lpos="W,E,r,-8.00,0.00" if left else "E,W,l,8.00,0.00")
        else:                                                        # 圈内：名称水平，中央列在右、其余在左
            look[v] = dict(w=node_size(deg[v]), h=node_size(deg[v]), rot=0.0,
                           lpos="E,W,l,8.00,0.00" if v in col else "W,E,r,-8.00,0.00")

    def extent(v):                                                   # 节点及名称所占矩形（外圈名称按旋转后的端点估计）
        x, y = pos[v]
        if ty[v] == "eRNA":
            w = look[v]["w"]; return x - w / 2, x + w / 2, y - 23, y + 23
        r = look[v]["w"] / 2
        if v in leaf_of:
            a_ = math.atan2(y, x); L_ = r + 8 + text_w(lab[v]); ex, ey = x + L_ * math.cos(a_), y + L_ * math.sin(a_)
            return min(x - r, ex) - 16, max(x + r, ex) + 16, min(y - r, ey) - 16, max(y + r, ey) + 16
        w = text_w(lab[v]) + 8
        return (x - r, x + r + w, y - 18, y + 18) if v in col else (x - r - w, x + r, y - 18, y + 18)
    return pos, look, {v: extent(v) for v in bn.node}, order, core


def build_style():
    if STYLE in req("GET", "/styles"):
        req("DELETE", f"/styles/{STYLE}")
    defaults = {"NODE_SHAPE": "ELLIPSE", "NODE_BORDER_WIDTH": 2.0, "NODE_BORDER_PAINT": "#FFFFFF",
                "NODE_LABEL_FONT_FACE": f"{FONT},plain,{FS}", "NODE_LABEL_FONT_SIZE": FS, "NODE_LABEL_COLOR": "#262626",
                "NODE_LABEL_WIDTH": 2000, "NODE_LABEL_POSITION": "C,C,c,0.00,0.00",
                "EDGE_TARGET_ARROW_SHAPE": "NONE", "EDGE_TRANSPARENCY": 190, "NETWORK_BACKGROUND_PAINT": "#FFFFFF"}
    disc = lambda vp, col, m, kind="String": {"mappingType": "discrete", "mappingColumn": col, "mappingColumnType": kind,
                                             "visualProperty": vp, "map": [{"key": k, "value": v} for k, v in m.items()]}
    pas = lambda vp, col, kind="String": {"mappingType": "passthrough", "mappingColumn": col, "mappingColumnType": kind,
                                         "visualProperty": vp}
    mappings = [disc("NODE_FILL_COLOR", "type", TYPE_COL),
                disc("NODE_SHAPE", "type", {"eRNA": "ROUND_RECTANGLE"}),
                disc("NODE_LABEL_COLOR", "type", {"eRNA": "#FFFFFF"}),
                disc("NODE_LABEL_FONT_FACE", "type", {"eRNA": f"{FONT},bold,30"}),
                disc("NODE_BORDER_WIDTH", "type", {"eRNA": "0"}),
                pas("NODE_LABEL", "label"), pas("NODE_WIDTH", "node_w", "Double"), pas("NODE_HEIGHT", "node_h", "Double"),
                disc("EDGE_STROKE_UNSELECTED_PAINT", "interaction", EDGE_COL),
                {"mappingType": "continuous", "mappingColumn": "abs_rho", "mappingColumnType": "Double", "visualProperty": "EDGE_WIDTH",
                 "points": [{"value": 0.4, "lesser": str(W_LO), "equal": str(W_LO), "greater": str(W_LO)},
                            {"value": 1.0, "lesser": str(W_HI), "equal": str(W_HI), "greater": str(W_HI)}]}]
    req("POST", "/styles", json={"title": STYLE, "defaults": [{"visualProperty": k, "value": v} for k, v in defaults.items()],
                                 "mappings": mappings})
    req("PUT", f"/styles/{STYLE}/dependencies", json=[{"visualPropertyDependency": "nodeSizeLocked", "enabled": False}])


def draw(key, fig):
    bn = pd.read_csv(os.path.join(OUT, f"{key}_backbone_nodes.csv"))
    be = pd.read_csv(os.path.join(OUT, f"{key}_backbone_edges.csv"))
    be["interaction"] = be.sign.str.lower()
    lab = dict(zip(bn.node, bn.label)); ty = dict(zip(bn.node, bn.type))
    pos, look, ext, order, core = layout(bn, be)

    # ---- 导入 Cytoscape ----
    title = f"{fig} revised (polished)"
    sif = os.path.join(OUT, f"{fig}_polished.sif")
    with open(sif, "w", encoding="utf-8") as f:
        for r in be.itertuples():
            f.write(f"{r.source}\t{r.interaction}\t{r.target}\n")
    for n in req("GET", "/networks.names"):
        if n["name"] in (title, os.path.basename(sif)):
            req("DELETE", f"/networks/{n['SUID']}")
    d = cmd("network load file", file=sif)
    suid, view = d["networks"][0], d["views"][0]
    for k in range(5):                                               # 新网络的表格偶尔尚未就绪，稍等后重试
        try:
            time.sleep(2); cmd("network rename", name=title, sourceNetwork=f"SUID:{suid}"); break
        except RuntimeError:
            if k == 4: raise
    req("PUT", f"/networks/{suid}/tables/defaultnode",
        json={"key": "name", "dataKey": "name",
              "data": [{"name": v, "label": lab[v], "type": ty[v], "degree": int(dg),
                        "node_w": float(look[v]["w"]), "node_h": float(look[v]["h"])} for v, dg in zip(bn.node, bn.degree)]})
    req("PUT", f"/networks/{suid}/tables/defaultedge",
        json={"key": "name", "dataKey": "name",
              "data": [{"name": f"{r.source} ({r.interaction}) {r.target}", "rho": float(r.rho), "abs_rho": float(r.abs_rho)}
                       for r in be.itertuples()]})
    req("GET", f"/apply/styles/{STYLE}/{suid}")
    time.sleep(2)                                                    # 等导入时的默认布局完成后再设坐标
    sid = {r["name"]: r["SUID"] for r in req("GET", f"/networks/{suid}/tables/defaultnode/rows")}
    req("PUT", f"/networks/{suid}/views/{view}/nodes",
        json=[{"SUID": sid[v], "view": [{"visualProperty": "NODE_X_LOCATION", "value": pos[v][0]},
                                        {"visualProperty": "NODE_Y_LOCATION", "value": pos[v][1]}]} for v in bn.node])
    req("PUT", f"/networks/{suid}/views/{view}/nodes?bypass=true",   # 节点名称的方向与旋转角逐个设置
        json=[{"SUID": sid[v], "view": [{"visualProperty": "NODE_LABEL_POSITION", "value": look[v]["lpos"]},
                                        {"visualProperty": "NODE_LABEL_ROTATION", "value": look[v]["rot"]}]}
              for v in bn.node if ty[v] != "eRNA"])
    req("PUT", "/networks/currentNetwork", json={"networkSUID": suid})

    # ---- 图例（注释）：放在四角中不压到节点与名称的位置 ----
    X0, X1 = min(e[0] for e in ext.values()), max(e[1] for e in ext.values())
    Y0, Y1 = min(e[2] for e in ext.values()), max(e[3] for e in ext.values())
    V = f"SUID:{view}"
    present = [t for t in ["eRNA"] + TYPE_ORDER if t in set(bn.type)]
    ref_deg = [d_ for d_ in (10, 50, 100) if bn.degree.min() <= d_ <= bn.degree.max()]
    ref_rho = [r_ for r_ in (0.5, 0.9) if be.abs_rho.min() <= r_ <= be.abs_rho.max() + 0.05]
    sizeA = (76 + max(text_w(TYPE_LAB[t], LF) for t in present), RH * (len(present) + 1))
    sizeB = (76 + max(text_w("Degree in full network", TF, True) - 76, text_w("Negative", LF)),
             RH * (3 + 2 + len(ref_rho) + len(ref_deg)) + 2 * 14)
    # 三张图位置统一：节点类型图例在左上角，边与节点大小图例在右下角；
    # 若压到节点或名称，左上角图例向上移、右下角图例向下移，直到不重叠（只增加图高，不改变宽度与字号）
    w, h = sizeA; y = Y0
    while any(overlap((X0, X0 + w, y, y + h), e) for e in ext.values()):
        y -= 10
    placed = [(X0, X0 + w, y, y + h)]
    w, h = sizeB; y = Y1 - h
    while any(overlap((X1 - w, X1, y, y + h), e) for e in ext.values()):
        y += 10
    placed.append((X1 - w, X1, y, y + h))
    print(f"  legend shift: top-left {placed[0][2] - Y0:.0f} px, bottom-right {placed[1][3] - Y1:.0f} px")
    def text(t_, x, cy, size=LF, bold=False):                        # cy 为该行的竖直中心
        cmd("annotation add text", text=t_, x=x, y=cy - size * 0.62, fontSize=size, fontFamily=FONT,
            fontStyle="bold" if bold else "plain", color="#1A1A1A", view=V)
    def shape(kind, cx, cy, w, h, fill):
        cmd("annotation add shape", type=kind, x=cx - w / 2, y=cy - h / 2, width=w, height=h, fillColor=fill,
            borderColor=fill, borderThickness=0, view=V)
    x0 = placed[0][0]; cy = placed[0][2] + RH / 2
    text("Node type", x0, cy, TF, True); cy += RH
    for t in present:
        if t == "eRNA":
            shape("Rounded Rectangle", x0 + 30, cy, 58, 28, TYPE_COL[t])
        else:
            shape("Ellipse", x0 + 30, cy, 26, 26, TYPE_COL[t])
        text(TYPE_LAB[t], x0 + 76, cy); cy += RH
    x0 = placed[1][0]; cy = placed[1][2] + RH / 2
    text("Spearman correlation", x0, cy, TF, True); cy += RH
    for k, lb in (("positive", "Positive"), ("negative", "Negative")):
        shape("Rectangle", x0 + 30, cy, 58, 5, EDGE_COL[k]); text(lb, x0 + 76, cy); cy += RH
    cy += 14; text("Correlation strength", x0, cy, TF, True); cy += RH  # PDF 导出字体不含希腊字母，图注说明边宽与 |ρ| 成正比
    for r_ in ref_rho:
        shape("Rectangle", x0 + 30, cy, 58, edge_width(r_), "#7F7F7F"); text(f"{r_:.1f}", x0 + 76, cy); cy += RH
    cy += 14; text("Degree in full network", x0, cy, TF, True); cy += RH
    for d_ in ref_deg:
        s_ = node_size(d_); shape("Ellipse", x0 + 30, cy, s_, s_, "#A6A6A6"); text(str(d_), x0 + 76, cy); cy += RH
    req("GET", f"/apply/fit/{suid}")
    cmd("view update", network=f"SUID:{suid}")
    time.sleep(5)                                                    # 等界面渲染完成后再导出

    # ---- 导出：Cytoscape 按窗口比例导出整页，用 PyMuPDF 求出图形与文字的实际范围后裁边并缩放 ----
    stem = os.path.join(OUT, f"{fig}_cytoscape_polished")
    for e_ in ("pdf", "svg", "png", "tif"):
        if os.path.exists(f"{stem}.{e_}"): os.remove(f"{stem}.{e_}")
    raw, raw_svg = stem + "_raw.pdf", stem + "_raw.svg"
    cmd("view export pdf", outputFile=raw, view=V, pageSize="Auto", exportTextAsFont=True)
    cmd("view export svg", outputFile=raw_svg, view=V, exportTextAsFont=True)
    d0 = fitz.open(raw); P = d0[0]; H = P.rect.height
    boxes = [fitz.Rect(r_) for k_, r_ in P.get_bboxlog() if fitz.Rect(r_) != P.rect]   # 去掉整页白色背景
    bb = boxes[0]
    for r_ in boxes[1:]:
        bb |= r_
    sizes = [round(s["size"], 2) for b in P.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]]
    s_raw = max(set(sizes), key=sizes.count)                         # 出现最多的字号即 30 px 的节点名称
    d0.close()
    M = 4 * s_raw / 4.3                                              # 四周留白约 1.4 mm
    bx0, by0, bx1, by1 = bb.x0 - M, bb.y0 - M, bb.x1 + M, bb.y1 + M
    k = min(LABEL_PT / s_raw, MAX_W_MM / 25.4 * 72 / (bx1 - bx0))    # 三张图比例尺相同；过宽时缩到 183 mm
    W_PT, H_PT = (bx1 - bx0) * k, (by1 - by0) * k
    pg = PdfReader(raw).pages[0]
    pg.add_transformation(Transformation().translate(-bx0, -(H - by1)).scale(k, k))
    pg.mediabox.lower_left = (0, 0); pg.mediabox.upper_right = (W_PT, H_PT)
    for box_ in ("cropbox", "trimbox", "artbox", "bleedbox"):
        setattr(pg, box_, pg.mediabox)
    w_ = PdfWriter(); w_.add_page(pg); w_.write(stem + ".pdf"); os.remove(raw)
    svg = open(raw_svg, encoding="utf-8", errors="replace").read()   # SVG 与 PDF 坐标一致：改 viewBox 与物理尺寸
    svg = re.sub(r'width="[^"]*"\s+height="[^"]*"\s+viewBox="[^"]*"',
                 f'width="{W_PT / 72 * 25.4:.1f}mm" height="{H_PT / 72 * 25.4:.1f}mm" '
                 f'viewBox="{bx0:.2f} {by0:.2f} {bx1 - bx0:.2f} {by1 - by0:.2f}"', svg, count=1)
    open(stem + ".svg", "w", encoding="utf-8").write(svg); os.remove(raw_svg)
    doc = fitz.open(stem + ".pdf")
    for dpi, e_ in ((600, "tif"), (300, "png")):
        pix = doc[0].get_pixmap(dpi=dpi, alpha=False)
        img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
        img.save(stem + ".tif", dpi=(dpi, dpi), compression="tiff_lzw") if e_ == "tif" else img.save(stem + ".png", dpi=(dpi, dpi))
    sz = sorted({round(s["size"], 1) for b in doc[0].get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]})
    print(f"  {fig}: {len(bn)} nodes, {len(be)} edges; core eRNAs: {sorted(lab[x] for x in core)}; "
          f"{doc[0].rect.width / 72 * 25.4:.0f} x {doc[0].rect.height / 72 * 25.4:.0f} mm; font sizes (pt): {sz}")
    doc.close()


build_style()
for key, fig in NETS:
    print(key)
    draw(key, fig)
cmd("session save as", file=os.path.join(OUT, "Figure 8_S8_S9_cytoscape.cys"))
print("session saved")
