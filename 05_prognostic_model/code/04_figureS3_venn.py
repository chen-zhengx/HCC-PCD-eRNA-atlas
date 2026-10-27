# =========================================================
# 对应结果：Figure S3（三个预后模型风险 eRNA 的 Venn 图：失巢凋亡 11 个、凋亡 7 个、铁死亡 10 个，共 12 个不同 eRNA）
# 输入：output/<module>/apparent_model_coefficients.csv（由 code/03_internal_validation.R 生成
#       与 Table S7 的 eRNA 相同）
# 输出：output/FigureS3_risk_eRNA_venn.pdf
# 运行：在模块目录下 python code/04_figureS3_venn.py（先运行 03_internal_validation.R 的三个模型）
# =========================================================
from __future__ import annotations

import csv  # 读取模型 eRNA
import os
from pathlib import Path

LOCAL_MPL_CONFIG = Path("output") / ".mplconfig"
LOCAL_MPL_CONFIG.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(LOCAL_MPL_CONFIG))

import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def read_model_eRNAs(module: str) -> set[str]:
    with open(Path("output") / module / "apparent_model_coefficients.csv", newline="", encoding="utf-8") as f:
        return {row["eRNA"] for row in csv.DictReader(f)}


SET_A = read_model_eRNAs("ferroptosis")  # A = Ferroptosis
SET_B = read_model_eRNAs("anoikis")  # B = Anoikis
SET_C = read_model_eRNAs("apoptosis")  # C = Apoptosis


# ============================================================
# Edit values here.
# This script draws a balanced 3-set Venn diagram with all 7 regions.
# ============================================================
CONFIG = {
    "title": "A: Ferroptosis   B: Anoikis   C: Apoptosis",
    "set_labels": ["A", "B", "C"],   # 绘图代码按 A/B/C 定位圆心和标签，保持不变
    "region_values": {
        "A_only": len(SET_A - SET_B - SET_C),
        "B_only": len(SET_B - SET_A - SET_C),
        "C_only": len(SET_C - SET_A - SET_B),
        "AB_only": len((SET_A & SET_B) - SET_C),
        "AC_only": len((SET_A & SET_C) - SET_B),
        "BC_only": len((SET_B & SET_C) - SET_A),
        "ABC": len(SET_A & SET_B & SET_C),
    },
    "output_name": "FigureS3_risk_eRNA_venn.pdf",
}


BACKGROUND_COLOR = "#fbfaf7"
TEXT_COLOR = "#1f1f1f"
EDGE_COLOR = "#333333"
CIRCLE_COLORS = ["#4DBBD5", "#E64B35", "#00A087"]

# Decrease spacing to enlarge the two-set intersection regions.
CIRCLE_RADIUS = 2.25
CIRCLE_CENTERS = {
    "A": (-1.15, 0.55),
    "B": (1.15, 0.55),
    "C": (0.0, -0.95),
}

REGION_POSITIONS = {
    "A_only": (-1.95, 0.55),
    "B_only": (1.95, 0.55),
    "C_only": (0.0, -2.0),
    "AB_only": (0.0, 1.18),
    "AC_only": (-0.86, -0.48),
    "BC_only": (0.86, -0.48),
    "ABC": (0.0, 0.03),
}

LABEL_POSITIONS = {
    "A": (-2.62, 2.52),
    "B": (2.62, 2.52),
    "C": (0.0, -3.45),
}


def validate_config() -> None:
    labels = CONFIG["set_labels"]
    if len(labels) != 3:
        raise ValueError("set_labels must contain exactly 3 labels.")

    required_regions = set(REGION_POSITIONS)
    value_regions = set(CONFIG["region_values"])
    missing = sorted(required_regions - value_regions)
    extra = sorted(value_regions - required_regions)
    if missing or extra:
        raise ValueError(
            f"region_values keys must be {sorted(required_regions)}. "
            f"Missing: {missing}; extra: {extra}"
        )


def draw_region_text(ax: plt.Axes, key: str, label: str, value: object) -> None:
    x, y = REGION_POSITIONS[key]
    ax.text(
        x,
        y,
        f"{label}\n{value}",
        ha="center",
        va="center",
        fontsize=11.5,
        fontweight="bold",
        color=TEXT_COLOR,
        linespacing=1.12,
        bbox={
            "boxstyle": "round,pad=0.24",
            "facecolor": "#fffdf9",
            "edgecolor": "none",
            "alpha": 0.72,
        },
    )


def draw_venn(output_dir: Path) -> Path:
    validate_config()

    labels = CONFIG["set_labels"]
    values = CONFIG["region_values"]

    fig, ax = plt.subplots(figsize=(7.6, 7.2), dpi=220)
    fig.patch.set_facecolor(BACKGROUND_COLOR)
    ax.set_facecolor(BACKGROUND_COLOR)
    ax.set_aspect("equal")
    ax.axis("off")

    for label, color in zip(labels, CIRCLE_COLORS, strict=True):
        circle = Circle(
            CIRCLE_CENTERS[label],
            radius=CIRCLE_RADIUS,
            facecolor=color,
            edgecolor=EDGE_COLOR,
            linewidth=1.35,
            alpha=0.34,
        )
        ax.add_patch(circle)

    region_labels = {
        "A_only": labels[0],
        "B_only": labels[1],
        "C_only": labels[2],
        "AB_only": f"{labels[0]} & {labels[1]}",
        "AC_only": f"{labels[0]} & {labels[2]}",
        "BC_only": f"{labels[1]} & {labels[2]}",
        "ABC": f"{labels[0]} & {labels[1]} & {labels[2]}",
    }
    for key, label in region_labels.items():
        draw_region_text(ax, key, label, values[key])

    for label in labels:
        x, y = LABEL_POSITIONS[label]
        ax.text(
            x,
            y,
            label,
            ha="center",
            va="center",
            fontsize=16,
            fontweight="bold",
            color=TEXT_COLOR,
        )

    ax.set_xlim(-3.9, 3.9)
    ax.set_ylim(-3.75, 3.35)
    ax.set_title(CONFIG["title"], fontsize=16, pad=14, color=TEXT_COLOR)

    output_path = output_dir / str(CONFIG["output_name"])
    fig.savefig(output_path, bbox_inches="tight", facecolor=BACKGROUND_COLOR)
    plt.close(fig)
    return output_path


def main() -> None:
    output_path = draw_venn(Path("output"))
    print(f"Saved: {output_path}")
    print("Region values:", CONFIG["region_values"])


if __name__ == "__main__":
    main()
