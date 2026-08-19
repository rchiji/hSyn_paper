"""{split}_class_metrics.csv / {split}_confusion_matrix.csv から
クラス別 IoU 図・Dice 図・混同行列ヒートマップを split ごとに描く。

実行: python scripts/exp_val_metrics_plot.py
出力: assets/{split}_class_iou.png, {split}_class_dice.png, {split}_confusion_matrix.png
      (split = train, val, trainval)
"""
import csv
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"   # SVG 内テキストを編集可能なまま保持
plt.rcParams["pdf.fonttype"] = 42       # PDF フォントを TrueType 埋め込みに


def save_all_formats(fig, out_path: Path) -> None:
    """png / svg / pdf の3形式で保存する。out_path は .png のパス。"""
    fig.savefig(out_path, dpi=200)
    fig.savefig(out_path.with_suffix(".svg"))
    fig.savefig(out_path.with_suffix(".pdf"))

ASSETS: Path = Path(__file__).resolve().parent.parent / "assets"
SPLITS: dict[str, str] = {"train": "training set", "val": "validation set", "trainval": "training + validation set"}
METRIC_COLORS: dict[str, str] = {"IoU": "#3B6FB6", "Dice": "#E28934"}


def bar_plot(names: list[str], values: list[float], metric: str, split_label: str, out_path: Path) -> None:
    """クラス別スコアの横棒グラフを 1 メトリクス分描く。"""
    y: np.ndarray = np.arange(len(names))  # (13,)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.barh(y, values, height=0.62, color=METRIC_COLORS[metric])
    for yi, v in zip(y, values):
        ax.text(v + 0.012, yi, f"{v:.3f}", va="center", fontsize=8.5, color="#333")
    ax.set_yticks(y, names)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.1)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel(metric)
    ax.set_title(f"Per-class {metric} ({split_label})", loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#e0e0e0", linewidth=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()
    save_all_formats(fig, out_path)
    plt.close(fig)


def confusion_plot(names: list[str], conf: np.ndarray, split_label: str, out_path: Path) -> None:
    """行正規化した混同行列ヒートマップを描く。conf: (13,13)"""
    norm: np.ndarray = conf / conf.sum(1, keepdims=True)  # (13,13)
    fig, ax = plt.subplots(figsize=(8.5, 7.5))
    im = ax.imshow(norm, cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(13), names, rotation=45, ha="right")
    ax.set_yticks(range(13), names)
    ax.set_xlabel("Predicted class")
    ax.set_ylabel("Ground truth class")
    ax.set_title(f"Confusion matrix (row-normalized, {split_label})", loc="left")
    for i in range(13):
        for j in range(13):
            v: float = norm[i, j]
            if v >= 0.005:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if v > 0.6 else "#333")
    fig.colorbar(im, ax=ax, shrink=0.8, label="Fraction of ground-truth pixels")
    fig.tight_layout()
    save_all_formats(fig, out_path)
    plt.close(fig)


for split, split_label in SPLITS.items():
    names: list[str] = []
    iou: list[float] = []
    dice: list[float] = []
    with open(ASSETS / f"{split}_class_metrics.csv", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            names.append(row["name"])
            iou.append(float(row["IoU"]))
            dice.append(float(row["Dice"]))
    bar_plot(names, iou, "IoU", split_label, ASSETS / f"{split}_class_iou.png")
    bar_plot(names, dice, "Dice", split_label, ASSETS / f"{split}_class_dice.png")
    conf: np.ndarray = np.loadtxt(ASSETS / f"{split}_confusion_matrix.csv", delimiter=",", skiprows=1)  # (13,14)
    confusion_plot(names, conf[:, :13], split_label, ASSETS / f"{split}_confusion_matrix.png")
    print(f"saved: {split}_class_iou.png, {split}_class_dice.png, {split}_confusion_matrix.png")
