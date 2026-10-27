# =========================================================
# 对应结果：Table S2（18 个 PCD 基因集）与 Table S3（PCD 过程筛选：超几何检验 + BH 校正）
# 输入：input/pcd_genesets_TableS2.csv、input/common_DE_genes_649.csv
# 输出：output/Table_S2_merged_categories.csv、output/Table_S3_merged_categories.csv
# 运行：在模块目录下 python code/02_merge_categories_hypergeometric.py
# =========================================================
# PCD 基因集：合并重复类别并统一类别名称（基因符号不变），再做超几何检验
#   Entosisname + Entoticcelldeath -> Entosis；NETosis + Netoticcelldeath -> NETosis（合并后同类别内重复的基因只保留一次）
#   其余类别名仅修正书写（Oxeiptosisname、panoptosis、缺少空格的名称）
# 超几何检验：a = 19,416 个编码基因，b = 649 个共有差异基因，c = 该类别条目数，d = c 中的差异基因数
import re, sys, zipfile
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

sys.stdout.reconfigure(encoding="utf-8")
import os
ROOT = "output"; os.makedirs(ROOT, exist_ok=True)
rename = {"Entosisname": "Entosis", "Entoticcelldeath": "Entosis", "Netoticcelldeath": "NETosis",
          "Oxeiptosisname": "Oxeiptosis", "panoptosis": "PANoptosis",
          "Immunogeniccelldeath": "Immunogenic cell death",
          "Lysosome-dependentcelldeath": "Lysosome-dependent cell death"}

# ---------- 1. 读取 PCD 基因集（20 个类别）----------
s2 = pd.read_csv("input/pcd_genesets_TableS2.csv", dtype=str, keep_default_na=False)[["Symbol", "Type"]]
s2["Symbol"] = s2.Symbol.str.strip()
print("20 categories: annotations", len(s2), "| distinct symbols", s2.Symbol.nunique(), "| categories", s2.Type.nunique())
s2["Category"] = s2.Type.replace(rename)
dups = s2[s2.duplicated(["Symbol", "Category"], keep="first")]
print("duplicate gene-category rows created by merging:", dups[["Symbol", "Type", "Category"]].values.tolist())
new = s2.drop_duplicates(["Symbol", "Category"], keep="first")
new = new.assign(_o=new.Category.str.lower()).sort_values("_o", kind="stable")[["Symbol", "Category"]]
print("18 categories: annotations", len(new), "| distinct symbols", new.Symbol.nunique(), "| categories", new.Category.nunique())
multi = new.groupby("Symbol").Category.nunique()
print("symbols in >1 category:", int((multi > 1).sum()), "| extra annotations from multi-category membership:",
      int((multi - 1).sum()))
new.to_csv(ROOT + "/Table_S2_merged_categories.csv", index=False)

# ---------- 2. 超几何检验 ----------
de = set(pd.read_csv("input/common_DE_genes_649.csv").gene)
a, b = 19416, 649
def test(df, col):
    rows = []
    for k, g in df.groupby(col, sort=False):
        c = g.Symbol.nunique(); d = len(set(g.Symbol) & de)
        rows.append([k, a, b, c, d, hypergeom.sf(d - 1, a, c, b)])
    out = pd.DataFrame(rows, columns=["PCD type", "a", "b", "c", "d", "P"])
    out["FDR"] = multipletests(out.P, method="fdr_bh")[1]
    return out.assign(_o=out["PCD type"].str.lower()).sort_values("_o").drop(columns="_o")

orig = test(s2, "Type")
print("\n20 categories:\n", orig.to_string(index=False))
res = test(new, "Category")
res.to_csv(ROOT + "/Table_S3_merged_categories.csv", index=False)
print("\n", res.to_string(index=False))

