#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""07 号 MR 结果森林图"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import csv, os

# 中文字体
for fp in ["C:/Windows/Fonts/msyh.ttc", "C:/Windows/Fonts/simhei.ttf"]:
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=fp).get_name()
        break
plt.rcParams["axes.unicode_minus"] = False

BASE = r"C:/Users/Hero/WorkBuddy/07_rheo_MR"
csv_f = os.path.join(BASE, "04_mr", "mr_summary_all.csv")

rows = []
with open(csv_f, encoding="utf-8") as f:
    rd = csv.DictReader(f)
    for r in rd:
        rows.append(r)

# 按结局分组
outcomes = ["VTE", "DVT", "PE"]
traits = ["HCT", "PLT", "MPV", "RBC", "PCT"]
colors = {"VTE": "#2c7fb8", "DVT": "#d95f0e", "PE": "#31a354"}

fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=True)
for ax, out in zip(axes, outcomes):
    sub = [r for r in rows if r["outcome"] == out]
    sub.sort(key=lambda x: traits.index(x["exposure"]))
    ys = range(len(sub))
    for y, r in zip(ys, sub):
        orv = float(r["OR"]); lo = float(r["OR_low"]); hi = float(r["OR_high"])
        p = float(r["pval"])
        c = "#c0392b" if p < 0.05 else "#7f8c8d"
        ax.errorbar(orv, y, xerr=[[orv - lo], [hi - orv]], fmt="o",
                    color=c, ecolor=c, capsize=3, markersize=6, linewidth=1.5)
        # 标注 OR 和 p
        sig = " *" if p < 0.05 else ""
        ax.text(hi + 0.06, y, f"{orv:.2f}{sig}", va="center", fontsize=8, color=c)
    ax.axvline(1.0, color="#555", linestyle="--", linewidth=1)
    ax.set_yticks(list(ys))
    ax.set_yticklabels([r["exposure"] for r in sub], fontsize=10)
    ax.set_title(f"{out}（静脉血栓）" if out == "VTE" else
                 (f"{out}（深静脉血栓）" if out == "DVT" else f"{out}（肺栓塞）"),
                 fontsize=12, fontweight="bold")
    ax.set_xlabel("OR (95% CI)", fontsize=10)
    ax.set_xlim(0.7, 1.7)
    ax.grid(axis="x", alpha=0.3)

axes[0].set_ylabel("暴露（流变性状）", fontsize=11)
fig.suptitle("流变性状 → 血栓：双样本孟德尔随机化（IVW）", fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])

out_f = os.path.join(BASE, "06_figures", "mr_forest.png")
os.makedirs(os.path.dirname(out_f), exist_ok=True)
fig.savefig(out_f, dpi=150, bbox_inches="tight")
print(f"written -> {out_f}")
