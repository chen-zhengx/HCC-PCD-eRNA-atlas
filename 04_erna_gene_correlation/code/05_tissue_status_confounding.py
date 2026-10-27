# =========================================================
# 对应结果：组织状态混杂检查（Table S6 中校正组织状态的偏相关与仅肿瘤相关）
# 输入：output/correlation/<PCD>/<GSE>/input/、output/correlation/<PCD>/<GSE>/output/<GSE>_All_Results.csv.gz（01、02 的输出）
#       input/<GSE>_series_matrix.txt.gz（GEO 组织标签）、input/PRJNA512451_runs.tsv（GSE124535 的 SRR–GSM 对应）
#       input/GSE144269_sample_ids.csv（GSE144269 的 SRR 与组织）、input/model_risk_eRNAs.csv（三个模型的风险 eRNA，= Table S7）
# 输出：output/confounding_pairs_by_dataset.csv、output/confounding_summary_5datasets.csv、output/confounding_summary_6datasets.csv
# 运行：在模块目录下 python code/05_tissue_status_confounding.py
# =========================================================
# 相关性：每个数据集合并肿瘤与癌旁样本计算 Spearman 相关，BH 校正；
#   除 GSE124535 外的 5 个数据集均 FDR<1e-5 的配对为显著配对（失巢凋亡 707、凋亡 1083、铁死亡 1541 对）
# 本脚本对这些配对在每个数据集中计算：仅肿瘤、仅癌旁、校正组织状态的偏 Spearman 相关，并做跨数据集 meta 分析
# =========================================================
import gzip, os, sys
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.stdout.reconfigure(encoding="utf-8")
ROOT = "output"; os.makedirs(ROOT, exist_ok=True)
INP = "input"
B = "output/correlation"
PCD = {"Anoikis": "Anoikis", "Apoptosis": "Apoptosis", "Ferroptosis": "Ferroptosis"}
DS = {g: g for g in ["GSE148355", "GSE77509", "GSE94660", "GSE124535", "GSE214846", "GSE144269"]}
SELECTION = ["GSE148355", "GSE77509", "GSE94660", "GSE214846", "GSE144269"]       # 参与筛选的 5 个数据集

# ---------- 组织标签（True = 肿瘤）----------
def geo_tissue(g):
    L = gzip.open(f"{INP}/{g}_series_matrix.txt.gz", "rt", encoding="utf-8", errors="ignore").read().splitlines()
    gsm = [x.strip('"') for x in next(l for l in L if l.startswith("!Sample_geo_accession")).split("\t")[1:]]
    tis = next([x.strip('"').split(":", 1)[1].strip() for x in l.split("\t")[1:]] for l in L
               if l.startswith("!Sample_characteristics_ch1") and "tissue:" in l)
    return dict(zip(gsm, tis))
TUMOR = {"GSE148355": {"Gastric cancer, liver"},          # GEO 中 HCC 样本被误标为 gastric cancer（分级字段为 G1–G3）
         "GSE77509": {"primary tumor"}, "GSE94660": {"Tumor tissue"},
         "GSE214846": {"hepatocellular carcinoma"}, "GSE124535": {"HCC"}}
labels = {g: {k: v in TUMOR[g] for k, v in geo_tissue(g).items()} for g in ["GSE148355", "GSE77509", "GSE94660", "GSE214846"]}
runs = pd.read_csv(f"{INP}/PRJNA512451_runs.tsv", sep="\t").set_index("run_accession").sample_alias
t124 = geo_tissue("GSE124535"); labels["GSE124535"] = {r: t124[g] in TUMOR["GSE124535"] for r, g in runs.items()}
ids = pd.read_csv(f"{INP}/GSE144269_sample_ids.csv")
labels["GSE144269"] = dict(zip(ids.SRR, ids.tissue == "tumor"))

# ---------- 行向量化的 Spearman 与偏 Spearman ----------
def ranks(M): return np.apply_along_axis(stats.rankdata, 1, M)
def corr_rows(rx, ry, df):
    rx = rx - rx.mean(1, keepdims=True); ry = ry - ry.mean(1, keepdims=True)
    den = np.sqrt((rx ** 2).sum(1) * (ry ** 2).sum(1))
    r = np.divide((rx * ry).sum(1), den, out=np.full(len(den), np.nan), where=den > 0)
    t = r * np.sqrt(df / np.clip(1 - r ** 2, 1e-12, None))
    return r, 2 * stats.t.sf(np.abs(t), df)
def spearman(X, Y): return corr_rows(ranks(X), ranks(Y), X.shape[1] - 2)
def partial_spearman(X, Y, z):
    Z = np.column_stack([np.ones(len(z)), z.astype(float)]); H = Z @ np.linalg.pinv(Z)
    rx, ry = ranks(X), ranks(Y)
    return corr_rows(rx - rx @ H.T, ry - ry @ H.T, X.shape[1] - 3)

rows = []
for p, folder in PCD.items():
    # 显著配对（5 个筛选数据集均 FDR<1e-5）
    allres = {g: pd.read_csv(f"{B}/{folder}/{d}/output/{d}_All_Results.csv.gz").set_index(["eRNA", "Gene"])
              for g, d in DS.items()}
    idx = allres["GSE144269"].index
    F = np.column_stack([allres[g].loc[idx, "FDR"].values for g in SELECTION])
    pairs = idx[(F < 1e-5).all(1)]
    print(f"{p}: significant pairs {len(pairs)} | eRNAs {pairs.get_level_values(0).nunique()} | genes {pairs.get_level_values(1).nunique()}")
    for g, d in DS.items():
        inp = f"{B}/{folder}/{d}/input"
        ef = [f for f in os.listdir(inp) if "eRNA" in f][0]; gf = [f for f in os.listdir(inp) if "eRNA" not in f][0]
        E = pd.read_csv(f"{inp}/{ef}", index_col=0); G = pd.read_csv(f"{inp}/{gf}", index_col=0)
        cols = [c for c in E.columns if c in G.columns]
        z = np.array([labels[g][c] for c in cols])
        X = E.loc[pairs.get_level_values(0), cols].to_numpy(float); Y = G.loc[pairs.get_level_values(1), cols].to_numpy(float)
        rp, _ = spearman(X, Y); rt, pt = spearman(X[:, z], Y[:, z]); rn, pn = spearman(X[:, ~z], Y[:, ~z])
        ra, pa = partial_spearman(X, Y, z)
        df = pd.DataFrame({"PCD": p, "dataset": g, "selection_dataset": g in SELECTION,
                           "eRNA": pairs.get_level_values(0), "gene": pairs.get_level_values(1),
                           "n_tumor": int(z.sum()), "n_adjacent": int((~z).sum()),
                           "rho_pooled": allres[g].loc[pairs, "Correlation"].values, "FDR_pooled": allres[g].loc[pairs, "FDR"].values,
                           "rho_pooled_check": rp, "rho_tumor": rt, "P_tumor": pt, "rho_adjacent": rn, "P_adjacent": pn,
                           "rho_partial": ra, "P_partial": pa})
        for c in ["tumor", "adjacent", "partial"]:
            ok = df[f"P_{c}"].notna()
            df.loc[ok, f"FDR_{c}"] = multipletests(df.loc[ok, f"P_{c}"], method="fdr_bh")[1]
        rows.append(df)
res = pd.concat(rows, ignore_index=True)
print("max |recomputed − pooled rho|:", float((res.rho_pooled_check - res.rho_pooled).abs().max()))
res.drop(columns="rho_pooled_check").to_csv(f"{ROOT}/confounding_pairs_by_dataset.csv", index=False)

# ---------- 跨数据集汇总与 meta 分析 ----------
def meta(r, n, k):
    """DerSimonian–Laird 随机效应合并 Fisher z（Spearman 方差近似 1.06/(n−3−k)）"""
    ok = ~np.isnan(r); r, n = r[ok], n[ok]
    zi = np.arctanh(np.clip(r, -0.999999, 0.999999)); vi = 1.06 / (n - 3 - k); w = 1 / vi
    zf = (w * zi).sum() / w.sum(); Q = (w * (zi - zf) ** 2).sum()
    tau2 = max(0.0, (Q - (len(zi) - 1)) / (w.sum() - (w ** 2).sum() / w.sum())); ws = 1 / (vi + tau2)
    zr = (ws * zi).sum() / ws.sum(); se = np.sqrt(1 / ws.sum())
    return np.tanh(zr), 2 * stats.norm.sf(abs(zr / se))

def summarize(sub, tag):
    out = []
    for (p, e, g), d in sub.groupby(["PCD", "eRNA", "gene"], sort=False):
        nd = len(d); sgn = np.sign(d.rho_pooled.values)
        main = np.sign(np.median(d.rho_pooled))
        ntot = (d.n_tumor + d.n_adjacent).values
        mr_a, mp_a = meta(d.rho_partial.values, ntot, 1)
        mr_t, mp_t = meta(d.rho_tumor.values, d.n_tumor.values, 0)
        out.append({"PCD": p, "eRNA": e, "gene": g, "n_datasets": nd,
                    "direction": "positive" if main > 0 else "negative",
                    "pooled_same_direction_all": bool((sgn == main).all()),
                    "median_rho_pooled": np.median(d.rho_pooled), "median_rho_partial": np.nanmedian(d.rho_partial),
                    "median_rho_tumor": np.nanmedian(d.rho_tumor),
                    "partial_same_direction_n": int((np.sign(d.rho_partial) == main).sum()),
                    "partial_sig_same_direction_n": int(((d.FDR_partial < 0.05) & (np.sign(d.rho_partial) == main)).sum()),
                    "tumor_sig_same_direction_n": int(((d.FDR_tumor < 0.05) & (np.sign(d.rho_tumor) == main)).sum()),
                    "meta_rho_partial": mr_a, "meta_P_partial": mp_a, "meta_rho_tumor": mr_t, "meta_P_tumor": mp_t})
    s = pd.DataFrame(out)
    for c in ["partial", "tumor"]:
        for p in s.PCD.unique():
            m = s.PCD == p; s.loc[m, f"meta_FDR_{c}"] = multipletests(s.loc[m, f"meta_P_{c}"], method="fdr_bh")[1]
        dirn = np.where(s.direction == "positive", 1, -1)
        s[f"retained_{c}"] = (s[f"meta_FDR_{c}"] < 0.05) & (np.sign(s[f"meta_rho_{c}"]) == dirn)
    s.to_csv(f"{ROOT}/confounding_summary_{tag}.csv", index=False)
    return s

risk = {}
rk = pd.read_csv(f"{INP}/model_risk_eRNAs.csv")
for m in ["anoikis", "apoptosis", "ferroptosis"]:
    risk[m.capitalize()] = set(rk.eRNA[rk.Model == m.capitalize()])
for tag, sub in [("5datasets", res[res.selection_dataset]), ("6datasets", res)]:
    s = summarize(sub, tag); nd = sub.dataset.nunique()
    print(f"\n================ {tag} ================")
    for p in PCD:
        q = s[s.PCD == p]; n = len(q); ne = q.eRNA.nunique()
        eret = q.groupby("eRNA").retained_partial.any()
        rk = sorted(risk[p]); rk_ret = [e for e in rk if eret.get(e, False)]
        print(f"{p:11s} pairs {n} | pooled direction consistent across datasets: {q.pooled_same_direction_all.sum()} ({100*q.pooled_same_direction_all.mean():.1f}%)")
        print(f"   median |rho|: pooled {q.median_rho_pooled.abs().median():.2f} | tissue-adjusted {q.median_rho_partial.abs().median():.2f} | tumor-only {q.median_rho_tumor.abs().median():.2f}")
        print(f"   tissue-adjusted: same direction in all {nd} datasets {(q.partial_same_direction_n == nd).sum()} ({100*(q.partial_same_direction_n == nd).mean():.1f}%)"
              f" | FDR<0.05 in >=3 datasets {(q.partial_sig_same_direction_n >= 3).sum()} ({100*(q.partial_sig_same_direction_n >= 3).mean():.1f}%)"
              f" | in all {nd} {(q.partial_sig_same_direction_n == nd).sum()} ({100*(q.partial_sig_same_direction_n == nd).mean():.1f}%)")
        print(f"   meta (random effects) FDR<0.05, same direction: tissue-adjusted {q.retained_partial.sum()} ({100*q.retained_partial.mean():.1f}%)"
              f" | tumor-only {q.retained_tumor.sum()} ({100*q.retained_tumor.mean():.1f}%)")
        print(f"   eRNAs with >=1 retained pair (tissue-adjusted meta): {eret.sum()}/{ne} ({100*eret.mean():.1f}%) | model risk eRNAs retained: {len(rk_ret)}/{len(rk)}"
              + (f" (not retained: {sorted(set(rk) - set(rk_ret))})" if len(rk_ret) < len(rk) else ""))
