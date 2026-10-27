# =========================================================
# 对应结果：Supplementary Table 12 的数据表（Excel，6 个工作表：Figure 8、S8、S9 所示的边；三个完整网络的全部显著边）
# 输入：output/full_network_edge_list.csv、output/network_backbone_edges_shown.csv（07_network_backbone.R 的输出）
# 输出：output/Supplementary Data_full co-expression association networks.xlsx（10 脚本的输入）
# 运行：在模块目录下 python code/09_supp_table_S12_xlsx.py
# =========================================================
import os, sys
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
OUT = "output"
NETS = [("Anoikis", "Figure 8"), ("Apoptosis", "Figure S8"), ("Ferroptosis", "Figure S9")]
REN_LAYER = {"Nearby genes (\u00b11 Mb)": "Candidate target genes (\u00b11 Mb)"}   # 与论文用词一致

full = pd.read_csv(os.path.join(OUT, "full_network_edge_list.csv"))
shown = pd.read_csv(os.path.join(OUT, "network_backbone_edges_shown.csv"))

xlsx = os.path.join(OUT, "Supplementary Data_full co-expression association networks.xlsx")
with pd.ExcelWriter(xlsx) as w:
    for key, fig in NETS:                                            # 图中所示的边：按 eRNA、层、ρ 排序
        s = shown[shown.Network == key].sort_values(["eRNA", "Layer", "Spearman_rho"], kind="stable")
        s = s.drop(columns="Network").replace(REN_LAYER)
        s = s.rename(columns={"Partner": "Correlated node", "Spearman_rho": "Spearman rho"})
        s.to_excel(w, sheet_name=f"{fig} (shown edges)", index=False)
        print(f"{fig} (shown edges): {len(s)} edges")
    for key, fig in NETS:                                            # 完整网络：保持 07 脚本输出的顺序
        x = full[full.Network == key].drop(columns="Network").replace(REN_LAYER)
        x = x.rename(columns={"Spearman_rho": "Spearman rho"})
        x.to_excel(w, sheet_name=f"{key} (full network)", index=False)
        print(f"{key} (full network): {len(x)} edges")
print("saved:", xlsx)
