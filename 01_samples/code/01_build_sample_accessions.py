# =========================================================
# 对应结果：Table S1（458 个样本的 GSM/SRR 编号、病人配对与组织类型）；../data/samples/ 中的样本表由本脚本的输出生成
# 输入：input/<GSE>_series_matrix.txt.gz、input/ena/<BioProject>.tsv、../data/expression/eRNA_RPM/<GSE>_eRNA_RPM.csv.gz（只读列名）
# 输出：output/Table_S1_sample_accessions.csv
# 运行：在模块目录下 python code/01_build_sample_accessions.py
# =========================================================
# Table S1：逐样本列出分析所用的 458 个样本（229 对）的 GEO/SRA 编号、病人配对与组织类型
# 样本范围：三种 PCD 相关性计算的输入矩阵列名（与 eRNA RPM 矩阵、差异分析所用样本一致）
# 元数据：GEO series matrix（C3/downloads）+ ENA read_run 报告（GSM ↔ SRR、平台、读段布局）
# 配对：按 GEO 中的病人编号确定；同时检查按列位置配对（第 i 列与第 i+n/2 列）是否与 GEO 中的配对一致
# =========================================================
import gzip, os, re, sys
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
ROOT = "input"
OUT  = "output"; os.makedirs(OUT, exist_ok=True)
ENA  = os.path.join(ROOT, "ena")  # ENA run 报告（本地缓存）
CORR = "../data/expression/eRNA_RPM"
DS = {
    "GSE148355": ("GSE148355_eRNA_RPM.csv.gz", "PRJNA624047"), "GSE77509": ("GSE77509_eRNA_RPM.csv.gz", "PRJNA310654"),
    "GSE94660": ("GSE94660_eRNA_RPM.csv.gz", "PRJNA371753"), "GSE124535": ("GSE124535_eRNA_RPM.csv.gz", "PRJNA512451"),
    "GSE214846": ("GSE214846_eRNA_RPM.csv.gz", "PRJNA887382"), "GSE144269": ("GSE144269_eRNA_RPM.csv.gz", "PRJNA603277")}
FIELDS = ("run_accession,sample_alias,sample_title,library_layout,instrument_model,read_count,base_count,fastq_ftp")

# ---------- ENA run 报告（本地缓存，缺失时下载）----------
def ena(project):
    f = os.path.join(ENA, f"{project}.tsv")
    if not os.path.exists(f):
        import urllib.request
        os.makedirs(ENA, exist_ok=True)
        url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={project}"
               f"&result=read_run&fields={FIELDS}&format=tsv")
        urllib.request.urlretrieve(url, f)
    return pd.read_csv(f, sep="\t")

# ---------- GEO series matrix ----------
def geo(g):
    L = gzip.open(f"{ROOT}/{g}_series_matrix.txt.gz", "rt", encoding="utf-8", errors="ignore").read().splitlines()
    cols, seen = {}, {}
    for l in L:
        if not l.startswith("!Sample_"): continue
        parts = [x.strip('"') for x in l.split("\t")]
        k = parts[0]; seen[k] = seen.get(k, 0) + 1
        cols[k if seen[k] == 1 else f"{k}_{seen[k]}"] = parts[1:]
    df = pd.DataFrame(cols).set_index("!Sample_geo_accession")
    ch = {}                                              # characteristics 展开为 key: value
    for c in [c for c in df.columns if c.startswith("!Sample_characteristics_ch1")]:
        for gsm, v in df[c].items():
            if ":" in v:
                key, val = v.split(":", 1); ch.setdefault(key.strip(), {})[gsm] = val.strip()
    return df, pd.DataFrame(ch)

def used_samples(folder):
    return list(pd.read_csv(os.path.join(CORR, folder), nrows=1, index_col=0).columns)

rows = []
for g, (folder, proj) in DS.items():
    used = used_samples(folder)
    runs = ena(proj)
    if used[0].startswith("SRR"):                         # 矩阵列为 SRR：反查 GSM，SRR 即所用 run
        gsm_of = runs.set_index("run_accession").sample_alias
        recs = [(gsm_of[s], [s]) for s in used]
    else:                                                 # 矩阵列为 GSM：列出该 GSM 的全部 run
        by = runs.groupby("sample_alias").run_accession.apply(sorted)
        recs = [(s, by[s]) for s in used]
    meta = runs.set_index("run_accession")
    if g != "GSE144269":
        df, ch = geo(g)
    for order, (gsm, srr) in enumerate(recs):
        r0 = meta.loc[srr[0]]
        if g == "GSE144269":                              # 未下载 series matrix，标题取自 ENA（与 GEO 标题相同）
            title = r0.sample_title
            m = re.match(r"(pat \d+) (tumor|non-tumor)", title)
            patient, tumor, annot = m.group(1), m.group(2) == "tumor", m.group(2)
        else:
            title = df.loc[gsm, "!Sample_title"]
            if g == "GSE148355":
                pid = ch.loc[gsm, "patient id"]; patient = pid.split("_")[0]
                stage = ch.loc[gsm, "tumor stage"] if "tumor stage" in ch else None
                annot = stage if isinstance(stage, str) else ch.loc[gsm, "tissue"]     # 正常肝没有 tumor stage 字段
                tumor = ch.loc[gsm, "tissue"] == "Gastric cancer, liver"
                if tumor:                                  # 核对 GEO 中登记的配对癌旁样本
                    paired = ch.loc[gsm, "paired nontumor"]
            elif g == "GSE77509":
                patient = ch.loc[gsm, "patient id"]; annot = ch.loc[gsm, "tissue"]; tumor = annot == "primary tumor"
            elif g == "GSE94660":
                patient = re.search(r"of (P\d+)$", title).group(1); annot = ch.loc[gsm, "tissue"]; tumor = annot == "Tumor tissue"
            elif g == "GSE124535":
                patient = title[:-1]; annot = ch.loc[gsm, "tissue"]; tumor = annot == "HCC"
            elif g == "GSE214846":
                desc = [df.loc[gsm, c] for c in df.columns if c.startswith("!Sample_description")]
                lab = [d for d in desc if re.match(r"L\d+_sample_[TN]$", d)][0]
                patient = lab.split("_")[0]; annot = ch.loc[gsm, "tissue"]; tumor = annot == "hepatocellular carcinoma"
        rows.append({"dataset": g, "bioproject": proj, "matrix_column": used[order], "column_order": order + 1,
                     "patient_id": patient, "tissue": "Tumor" if tumor else "Non-tumor", "geo_annotation": annot,
                     "sample_title": title, "GSM": gsm, "SRR": ";".join(srr), "n_runs": len(srr),
                     "instrument": r0.instrument_model, "sra_layout": r0.library_layout,
                     "fastq_files_per_run": len(str(r0.fastq_ftp).split(";")),
                     "read_count_M": round(meta.loc[srr, "read_count"].sum() / 1e6, 1),
                     "paired_nontumor_geo": paired if (g == "GSE148355" and tumor) else ""})
S = pd.DataFrame(rows)

# ---------- 核对 ----------
print("samples:", len(S), "| GSM unique:", S.GSM.is_unique, "| SRR unique:", not S.SRR.str.split(";").explode().duplicated().any())
for g, d in S.groupby("dataset", sort=False):
    cnt = d.groupby("patient_id").tissue.agg(lambda x: "".join(sorted(x)))
    ok = (cnt == "Non-tumorTumor").all()
    n = len(d) // 2
    pos = all(d.patient_id.iloc[i] == d.patient_id.iloc[i + n] and d.tissue.iloc[i] != d.tissue.iloc[i + n] for i in range(n))
    first = d.tissue.iloc[:n].unique().tolist()
    print(f"{g}: {len(d)} samples, {d.patient_id.nunique()} patients, each patient one tumor + one non-tumor: {ok} | "
          f"positional pairing (i, i+{n}) matches GEO pairing: {pos} | first half: {first}")
gs = S[(S.dataset == "GSE148355")]
t = gs[gs.tissue == "Tumor"].set_index("patient_id").paired_nontumor_geo
nt = gs[gs.tissue == "Non-tumor"].set_index("patient_id").sample_title
pid = {r.patient_id: r for _, r in gs.iterrows()}
geo148, ch148 = geo("GSE148355")
nt_pid = gs[gs.tissue == "Non-tumor"].set_index("patient_id").GSM.map(ch148["patient id"])
print("GSE148355 tumor 'paired nontumor' == non-tumor patient id:", (t.sort_index() == nt_pid.sort_index()).all())
print("GSE148355 non-tumor tissue types:", gs[gs.tissue == "Non-tumor"].geo_annotation.value_counts().to_dict())
print("GSE148355 tumor annotations:", gs[gs.tissue == "Tumor"].geo_annotation.value_counts().to_dict())
print("SRA layout vs fastq files:", S.groupby("dataset").apply(lambda d: f"{d.sra_layout.unique().tolist()} files/run {d.fastq_files_per_run.unique().tolist()}").to_dict())

S.to_csv(os.path.join(OUT, "Table_S1_sample_accessions.csv"), index=False)
print("saved:", os.path.join(OUT, "Table_S1_sample_accessions.csv"))
