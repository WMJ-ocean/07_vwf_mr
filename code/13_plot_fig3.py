#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""07 号 Fig 3 研究设计示意图（双样本 MR 两层设计）"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

BASE = r"C:/Users/Hero/WorkBuddy/07_rheo_MR"

fig, ax = plt.subplots(figsize=(9, 6.5), dpi=100)
ax.set_xlim(0, 100)
ax.set_ylim(0, 82)
ax.axis("off")

C_BLUE = "#E6F1FB"; C_BLUE_S = "#185FA5"; C_BLUE_T = "#0C447C"
C_PURP = "#EEEDFE"; C_PURP_S = "#534AB7"; C_PURP_T = "#26215C"
C_GRAY = "#F1EFE8"; C_GRAY_S = "#5F5E5A"; C_GRAY_T = "#2C2C2A"
C_GREEN = "#E1F5EE"; C_GREEN_S = "#0F6E56"; C_GREEN_T = "#04342C"

def box(x, y, w, h, fc, ec):
    b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.5",
                       linewidth=1.0, edgecolor=ec, facecolor=fc)
    ax.add_patch(b)

def txt(x, y, s, size=9, weight="normal", color="#222222", ha="left", va="center"):
    ax.text(x, y, s, fontsize=size, fontweight=weight, color=color, ha=ha, va=va)

def arrow(x1, y1, x2, y2, color="#444444"):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=12,
                        linewidth=1.4, color=color)
    ax.add_patch(a)

# ---- Title ----
txt(2, 79, "Figure 3. Two-sample Mendelian randomization study design",
    size=12, weight="bold", color="#111111")

# ---- Layer 1 box (blood-cell traits) ----
box(4, 52, 26, 24, C_BLUE, C_BLUE_S)
txt(6.5, 72, "Layer 1", size=10, weight="bold", color=C_BLUE_T)
txt(6.5, 68.5, "Blood-cell traits", size=10, weight="bold", color=C_BLUE_T)
txt(6.5, 64.5, "HCT  PLT  MPV", size=8.5, color=C_BLUE_S)
txt(6.5, 61.5, "RBC  PCT", size=8.5, color=C_BLUE_S)
txt(6.5, 58.5, "Astle 2016 (Cell)", size=8, color=C_BLUE_S)
txt(6.5, 55.5, "N up to 173,480", size=8, color=C_BLUE_S)

# ---- Layer 2 box (gene cis-eQTL) ----
box(4, 18, 26, 24, C_PURP, C_PURP_S)
txt(6.5, 38, "Layer 2", size=10, weight="bold", color=C_PURP_T)
txt(6.5, 34.5, "Gene cis-eQTL", size=10, weight="bold", color=C_PURP_T)
txt(6.5, 30.5, "24 shear / coag.", size=8.5, color=C_PURP_S)
txt(6.5, 27.5, "VWF  GP1BA  GP6 ...", size=8.5, color=C_PURP_S)
txt(6.5, 24.5, "eQTLGen (Võsa 2021)", size=8, color=C_PURP_S)
txt(6.5, 21.5, "N up to 31,684", size=8, color=C_PURP_S)

# ---- MR analysis box ----
box(42, 38, 24, 22, C_GRAY, C_GRAY_S)
txt(44.5, 56, "Two-sample MR", size=10, weight="bold", color=C_GRAY_T)
txt(44.5, 52, "Layer 1: IVW /", size=8.5, color=C_GRAY_S)
txt(44.5, 49.2, "Egger / weighted", size=8.5, color=C_GRAY_S)
txt(44.5, 46.4, "median", size=8.5, color=C_GRAY_S)
txt(44.5, 43, "Layer 2: Wald ratio", size=8.5, color=C_GRAY_S)
txt(44.5, 40.2, "Steiger direction", size=8.5, color=C_GRAY_S)

# ---- Outcomes box ----
box(82, 52, 15, 24, C_GREEN, C_GREEN_S)
txt(84.5, 72, "Outcomes", size=10, weight="bold", color=C_GREEN_T)
txt(84.5, 68.5, "VTE", size=8.5, color=C_GREEN_S)
txt(84.5, 65.5, "DVT", size=8.5, color=C_GREEN_S)
txt(84.5, 62.5, "PE", size=8.5, color=C_GREEN_S)
txt(84.5, 59.5, "FinnGen R12", size=8, color=C_GREEN_S)

# ---- Arrows ----
arrow(30, 64, 41.5, 55, "#444444")
arrow(30, 30, 41.5, 45, "#444444")
arrow(66, 49, 81.5, 62, "#444444")

# ---- Assumptions box ----
box(4, 2, 93, 16, "#FFFFFF", "#888888")
txt(6.5, 15.3, "Instrumental-variable assumptions", size=9.5, weight="bold", color="#222222")
txt(6.5, 12.0, "(i) Relevance — F > 10", size=8, color="#555555")
txt(6.5, 9.2, "(ii) Independence — strict LD clumping, MR-Egger intercept, Cochran Q", size=8, color="#555555")
txt(6.5, 6.4, "(iii) Exclusion restriction — Steiger filtering, bidirectional MR, pleiotropy-robust estimators", size=8, color="#555555")
txt(6.5, 3.6, "Reporting follows STROBE-MR", size=7.5, color="#888888")

out_f = os.path.join(BASE, "06_figures", "fig3_study_design.png")
fig.savefig(out_f, dpi=600, bbox_inches="tight", facecolor="white")
print(f"written -> {out_f}")
