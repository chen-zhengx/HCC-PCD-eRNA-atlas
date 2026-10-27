# =========================================================
# 对应结果：Table S6 与 Supplementary Data（全部 eRNA–基因对在每个数据集中的相关）
# 输入：output/confounding_pairs_by_dataset.csv、output/confounding_summary_5datasets.csv（05 的输出）
#       input/model_risk_eRNAs.csv
# 输出：output/Table_S6_top_gene_per_eRNA.csv、output/Supplementary Data_eRNA-gene correlation pairs.xlsx
# 运行：在模块目录下 python code/06_table_S6.py
# =========================================================
# PCD 相关 eRNA 及其支持性 eRNA–基因相关
# 配对：除 GSE124535 外 5 个数据集均 FDR<1e-5 的 eRNA–PCD 基因对（707/1083/1541 对；196/240/276 个 eRNA）
# Table S6：每个 eRNA 展示一个相关最强的 PCD 基因 = 5 个筛选数据集中合并样本 |rho| 中位数最大者，
#   列出显著相关的基因数，以及校正组织状态的偏相关、仅肿瘤相关的 meta 分析结果（来自 05_tissue_status_confounding.py）；
#   各数据集的合并样本 rho 只放在 Supplementary Data 中
# Supplementary Data（xlsx）：全部配对在每个数据集中的合并、仅肿瘤、偏相关 rho 与 FDR，以及 meta 分析结果
# =========================================================
import os, re, sys, zipfile
from xml.sax.saxutils import escape
import numpy as np, pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
ROOT     = "output"
OUT      = "output"
XLSX     = os.path.join(OUT, "Supplementary Data_eRNA-gene correlation pairs.xlsx")
PCD      = ["Anoikis", "Apoptosis", "Ferroptosis"]
ORDER    = ["GSE148355", "GSE144269", "GSE214846", "GSE124535", "GSE94660", "GSE77509"]   # 与 Table S1 相同
SELECT   = [g for g in ORDER if g != "GSE124535"]                                          # 筛选配对的 5 个数据集

# ---------- 数据 ----------
P = pd.read_csv(os.path.join(ROOT, "confounding_pairs_by_dataset.csv"))
M = pd.read_csv(os.path.join(ROOT, "confounding_summary_5datasets.csv"))
rk = pd.read_csv("input/model_risk_eRNAs.csv")
model = {m.capitalize(): set(rk.eRNA[rk.Model == m.capitalize()])
         for m in ["anoikis", "apoptosis", "ferroptosis"]}
key = ["PCD", "eRNA", "gene"]
W = P.pivot_table(index=key, columns="dataset",
                  values=["rho_pooled", "FDR_pooled", "rho_tumor", "FDR_tumor", "rho_partial", "FDR_partial"])
W.columns = [f"{v}_{g}" for v, g in W.columns]
W = W.reset_index().merge(M[key + ["direction", "meta_rho_partial", "meta_FDR_partial", "retained_partial",
                                   "meta_rho_tumor", "meta_FDR_tumor", "retained_tumor"]], on=key)
# 检查：配对数、筛选规则、eRNA 数
for p, n_pair, n_e in [("Anoikis", 707, 196), ("Apoptosis", 1083, 240), ("Ferroptosis", 1541, 276)]:
    w = W[W.PCD == p]
    assert len(w) == n_pair and w.eRNA.nunique() == n_e, p
    assert (w[[f"FDR_pooled_{g}" for g in SELECT]] < 1e-5).all().all(), p

# ---------- 每个 eRNA 选一个相关最强的基因 ----------
W["median_abs_rho_selection"] = W[[f"rho_pooled_{g}" for g in SELECT]].abs().median(axis=1)
W["abs_rho_GSE144269"] = W["rho_pooled_GSE144269"].abs()
W["n_genes"] = W.groupby(["PCD", "eRNA"]).gene.transform("size")
top = (W.sort_values(["PCD", "eRNA", "median_abs_rho_selection", "abs_rho_GSE144269"], ascending=[True, True, False, False])
         .groupby(["PCD", "eRNA"], sort=False).head(1).copy())
top["model_eRNA"] = [e in model[p] for p, e in zip(top.PCD, top.eRNA)]
def ekey(e):                                     # ENSR 编号在前，坐标形式的 eRNA 按染色体、起点排在后
    if e.startswith("ENSR"): return (0, 0, int(e[4:]))
    chrom, start = e.split(":")[0], int(e.split(":")[1].split("-")[0])
    return (1, 23 if chrom == "X" else 24 if chrom == "Y" else int(chrom), start)
top = top.sort_values(["PCD", "eRNA"], key=lambda s: s.map(ekey) if s.name == "eRNA" else s)
for p in PCD:
    t = top[top.PCD == p]
    print(f"{p}: eRNAs {len(t)} | negative top pairs: {(t.direction == 'negative').sum()} | model eRNAs: {t.model_eRNA.sum()}/{len(model[p])}"
          f" | tissue-adjusted retained: {t.retained_partial.sum()} | tumor-only retained: {t.retained_tumor.sum()}"
          f" | most frequent top genes: {t.gene.value_counts().head(3).to_dict()}")
top.to_csv(os.path.join(OUT, "Table_S6_top_gene_per_eRNA.csv"), index=False)

# ---------- Supplementary Data（xlsx）----------
cols = key + ["direction"]
for g in ORDER:
    cols += [f"rho_pooled_{g}", f"FDR_pooled_{g}", f"rho_partial_{g}", f"FDR_partial_{g}", f"rho_tumor_{g}", f"FDR_tumor_{g}"]
cols += ["meta_rho_partial", "meta_FDR_partial", "retained_partial", "meta_rho_tumor", "meta_FDR_tumor", "retained_tumor"]
allp = W.sort_values(["PCD", "eRNA", "median_abs_rho_selection"], ascending=[True, True, False])[cols]
allp = allp.rename(columns={"gene": "PCD_gene", "direction": "direction_pooled"})
notes = pd.DataFrame({"Item": [
    "Pairs", "Selection rule", "rho_pooled_<dataset>", "FDR_pooled_<dataset>", "rho_partial_<dataset>",
    "rho_tumor_<dataset>", "FDR_partial / FDR_tumor_<dataset>", "meta_rho_partial / meta_rho_tumor",
    "meta_FDR_partial / meta_FDR_tumor", "retained_partial / retained_tumor", "GSE124535"],
    "Description": [
    "All eRNA-gene pairs defining the programmed cell death (PCD)-related eRNAs: 707 anoikis-, 1,083 apoptosis- and 1,541 ferroptosis-related pairs.",
    "Spearman correlation across pooled tumor and non-tumor samples within each dataset, Benjamini-Hochberg (BH)-adjusted within each dataset; pairs with FDR<1e-5 in each of GSE148355, GSE77509, GSE94660, GSE214846 and GSE144269.",
    "Spearman rho across pooled tumor and non-tumor samples.",
    "BH-adjusted P value of the pooled correlation across all candidate pairs in the dataset.",
    "Partial Spearman correlation across all samples adjusted for tissue status (ranks regressed on tissue status).",
    "Spearman rho within tumor samples only.",
    "BH-adjusted within each dataset across the listed pairs.",
    "DerSimonian-Laird random-effects meta-analysis of Fisher z-transformed coefficients across the five selection datasets (variance 1.06/(n-3-k)).",
    "BH-adjusted within each PCD process.",
    "TRUE if meta-analytic FDR<0.05 and the same sign as the pooled correlation.",
    "Not used to select the pairs; shown for comparison."]})
with pd.ExcelWriter(XLSX, engine="openpyxl") as xw:
    notes.to_excel(xw, sheet_name="Notes", index=False)
    for p in PCD:
        allp[allp.PCD == p].drop(columns="PCD").to_excel(xw, sheet_name=f"{p} pairs", index=False)
    for ws in xw.book.worksheets:
        ws.freeze_panes = "B2" if ws.title == "Notes" else "C2"
        for c in ws.columns:
            ws.column_dimensions[c[0].column_letter].width = 110 if (ws.title == "Notes" and c[0].column_letter == "B") else 16
print("saved:", XLSX)

