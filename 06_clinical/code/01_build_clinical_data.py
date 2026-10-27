# =========================================================
# 对应结果：Table S4、Table S9 的分析数据（70 例 GSE144269 肿瘤的临床变量 + 三个模型的风险评分和高/低风险组）
# 输入：input/GSE144269_clinical_Candia2020.xlsx
#       input/risk_scores_<module>.csv（05_prognostic_model/code/03_internal_validation.R 输出的 output/<module>/patient_scores_apparent.csv）
# 输出：output/clinical_analysis_data.csv
# 运行：在模块目录下 python code/01_build_clinical_data.py
# =========================================================

# 由临床信息表（input/GSE144269_clinical_Candia2020.xlsx，整理自 Candia et al. 2020）构建 70 例 GSE144269 肿瘤的分析数据
# 并合并 05_prognostic_model 的风险评分与高/低风险组（开发队列中位数切分，即 Figure 3B/S1B/S2B 的高/低风险组）
# 风险评分来自 C1/<module>/patient_scores_apparent.csv（全部 70 例建模的 apparent 线性预测值及其中位数分组）
import os  # 用于创建 output/
import sys
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
os.makedirs("output", exist_ok=True)
u = pd.read_excel("input/GSE144269_clinical_Candia2020.xlsx", na_values=["NA", "#N/A"])
r = u.dropna(subset=["RNASeq_ID"])                                  # 有 RNA-seq 的 70 例
d = pd.DataFrame({
    "SRR": r.iloc[:, 3].str.strip(), "paper_patient_ID": r.Patient.astype(int), "GEO_patient": r.RNASeq_ID.astype(int),
    "age": pd.to_numeric(r.age, errors="coerce"), "male": r["sex(1=male,0=female)"],
    "tnm": pd.to_numeric(r.TNM_staging, errors="coerce"), "stage_III_IV": r.stage,
    "tumor_large": r["tumor_size(1 = '>5cm' ,0 = '<=5cm')"], "multinodular": r.multinodular, "cirrhosis": r.cirrhosis,
    "hbv": r.hbv, "hcv": r.hcv, "hdv": r.hdv, "afp_abn": r.afp, "alb_abn": r.alb, "bil_abn": r.bil, "alt_abn": r.alt,
    "bmi": pd.to_numeric(r.BMI, errors="coerce"), "alcohol": r.alcohol, "smoker": r["smoker(1 = 'ever', 0 = 'never')"],
    "fhx_lc": r["FHx_LC(0 = no, 1 = yes)"], "time_days": r["survival.time"], "status": r["survival.status"].astype(int)})
d["months"] = d.time_days / 30.4375
for m in ["anoikis", "apoptosis", "ferroptosis"]:
    s = pd.read_csv(f"input/risk_scores_{m}.csv")[["Sample", "apparent_lp", "apparent_high"]]
    s.columns = ["SRR", f"app_{m}", f"app_high_{m}"]
    d = d.merge(s, on="SRR", how="left")
assert len(d) == 70 and d.filter(like="app_high").notna().all().all()
d.to_csv("output/clinical_analysis_data.csv", index=False)
print("saved clinical_analysis_data.csv:", d.shape, "| deaths:", d.status.sum())
