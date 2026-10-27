# =========================================================
# 对应结果：Figure 2A（各数据集差异 eRNA 数及 6 个数据集共有的 559 个）、Figure 2B（差异基因，共有 649 个）
# 输入：花瓣数值 = output/common_DE_summary.csv 中各数据集的 DE_total，中心数值 = 最后一个数据集的 Cumulative_intersection
#       （04_common_DE.R 的输出）
# 输出：output/Figure2A_DE_eRNAs_petal_venn.pdf、output/Figure2B_DE_genes_petal_venn.pdf
# 运行：在模块目录下 python code/05_figure2_petal_venn.py
# =========================================================
from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Any

LOCAL_MPL_CONFIG = Path(__file__).resolve().parent / ".mplconfig"
LOCAL_MPL_CONFIG.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(LOCAL_MPL_CONFIG))

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse


# ============================================================
# Fill your data here.
# petal_count controls how many petals are drawn.
# petal_labels and petal_values must match petal_count.
# You can also change the labels, values, rotation, and output file names.
# Use .pdf if you want vector output for publication.
# ============================================================
# 填入 Figure 2 的数值（模板中为 0 占位的三花瓣示例）；花瓣从正上方开始逆时针排列
FIGURE_CONFIGS = [
    {
        "title": "",
        "petal_count": 6,
        "petal_labels": ["GSE148355", "GSE144269", "GSE214846", "GSE124535", "GSE94660", "GSE77509"],
        "petal_values": [13200, 11518, 23705, 5179, 12692, 26645],  # 各数据集差异 eRNA 数
        "center_label": "",
        "center_value": 559,  # 6 个数据集共有的差异 eRNA
        "output_name": "Figure2A_DE_eRNAs_petal_venn.pdf",
        "rotation_offset_deg": 90.0,
    },
    {
        "title": "",
        "petal_count": 6,
        "petal_labels": ["GSE148355", "GSE144269", "GSE214846", "GSE124535", "GSE94660", "GSE77509"],
        "petal_values": [3605, 2594, 3781, 2496, 4640, 2660],  # 各数据集差异基因数
        "center_label": "",
        "center_value": 649,  # 6 个数据集共有的差异基因
        "output_name": "Figure2B_DE_genes_petal_venn.pdf",
        "rotation_offset_deg": 90.0,
    },
]


DEFAULT_ROTATION_OFFSET_DEG = 90.0
BACKGROUND_COLOR = "#fbfaf7"
TEXT_COLOR = "#1f1f1f"
EDGE_COLOR = "#3a3a3a"

# Nature/Science/Cell style inspired palette with lower saturation and cleaner contrast.
PETAL_COLORS = [
    "#4DBBD5",
    "#E64B35",
    "#00A087",
    "#3C5488",
    "#F39B7F",
    "#8491B4",
]


def get_petal_count(config: dict[str, Any]) -> int:
    petal_count = int(config.get("petal_count", len(config["petal_labels"])))
    if petal_count < 2:
        raise ValueError("petal_count must be at least 2.")
    if len(config["petal_labels"]) != petal_count:
        raise ValueError("petal_labels length must match petal_count.")
    if len(config["petal_values"]) != petal_count:
        raise ValueError("petal_values length must match petal_count.")
    return petal_count


def get_geometry(config: dict[str, Any], petal_count: int) -> dict[str, float]:
    if petal_count <= 3:
        defaults = {
            "petal_center_distance": 5.7,
            "petal_width": 14.5,
            "petal_height": 7.4,
            "center_radius": 2.7,
        }
    else:
        defaults = {
            "petal_center_distance": 8.0,
            "petal_width": 16.0,
            "petal_height": max(4.8, 7.4 - 0.15 * petal_count),
            "center_radius": 2.8,
        }

    geometry = {
        key: float(config.get(key, value))
        for key, value in defaults.items()
    }
    geometry["value_radius"] = float(
        config.get(
            "value_radius",
            geometry["petal_center_distance"] + geometry["petal_width"] * 0.28,
        )
    )
    geometry["label_radius"] = float(
        config.get(
            "label_radius",
            geometry["petal_center_distance"] + geometry["petal_width"] * 0.55,
        )
    )
    return geometry


def draw_petal_venn(config: dict[str, Any], output_dir: Path) -> Path:
    petal_count = get_petal_count(config)
    geometry = get_geometry(config, petal_count)

    figure_size = float(config.get("figure_size", 8.0))
    fig, ax = plt.subplots(figsize=(figure_size, figure_size), dpi=200)
    fig.patch.set_facecolor(BACKGROUND_COLOR)
    ax.set_facecolor(BACKGROUND_COLOR)
    ax.set_aspect("equal")
    ax.axis("off")

    angle_step_deg = 360.0 / petal_count
    rotation_offset_deg = float(
        config.get("rotation_offset_deg", DEFAULT_ROTATION_OFFSET_DEG)
    )

    for index, (label, value) in enumerate(
        zip(config["petal_labels"], config["petal_values"], strict=True)
    ):
        angle_deg = index * angle_step_deg + rotation_offset_deg
        angle_rad = math.radians(angle_deg)
        x = geometry["petal_center_distance"] * math.cos(angle_rad)
        y = geometry["petal_center_distance"] * math.sin(angle_rad)
        color = PETAL_COLORS[index % len(PETAL_COLORS)]

        petal = Ellipse(
            (x, y),
            width=geometry["petal_width"],
            height=geometry["petal_height"],
            angle=angle_deg,
            facecolor=color,
            edgecolor=EDGE_COLOR,
            linewidth=1.1,
            alpha=0.42,
        )
        ax.add_patch(petal)

        value_x = geometry["value_radius"] * math.cos(angle_rad)
        value_y = geometry["value_radius"] * math.sin(angle_rad)
        label_x = geometry["label_radius"] * math.cos(angle_rad)
        label_y = geometry["label_radius"] * math.sin(angle_rad)

        ax.text(
            value_x,
            value_y,
            str(value),
            ha="center",
            va="center",
            fontsize=14,
            fontweight="bold",
            color=TEXT_COLOR,
        )
        ax.text(
            label_x,
            label_y,
            str(label),
            ha="center",
            va="center",
            fontsize=13,
            color=TEXT_COLOR,
        )

    center_circle = Circle(
        (0, 0),
        radius=geometry["center_radius"],
        facecolor="#fffdf9",
        edgecolor=EDGE_COLOR,
        linewidth=0.9,
    )
    ax.add_patch(center_circle)
    ax.text(
        0,
        0,
        f'{config["center_label"]}\n{config["center_value"]}',
        ha="center",
        va="center",
        fontsize=13,
        fontweight="bold",
        color=TEXT_COLOR,
    )

    plot_limit = max(
        18.0,
        geometry["label_radius"] + 2.0,
        geometry["petal_center_distance"] + geometry["petal_width"] * 0.7,
    )
    ax.set_xlim(-plot_limit, plot_limit)
    ax.set_ylim(-plot_limit, plot_limit)
    ax.set_title(str(config["title"]), fontsize=16, pad=16, color=TEXT_COLOR)

    output_path = output_dir / str(config["output_name"])
    fig.savefig(output_path, bbox_inches="tight", facecolor=BACKGROUND_COLOR)
    plt.close(fig)
    return output_path


def main() -> None:
    output_dir = Path(__file__).resolve().parent.parent / "output"
    output_dir.mkdir(exist_ok=True)
    saved_files = [draw_petal_venn(config, output_dir) for config in FIGURE_CONFIGS]
    for path in saved_files:
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
