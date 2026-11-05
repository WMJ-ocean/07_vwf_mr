#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""07 号 B 线结果森林图：基因 cis-eQTL → 血栓 OR（显著结果 + 关键基因）"""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial"]
plt.rcParams["axes.unicode_minus"] = False

r = pd.read_csv(r"C:/Users/Hero/WorkBuddy/07_rheo_MR/04_mr/gene_mr_wald.csv")

# 关键基因（分子对接靶点 + 剪切响应核心），按结局 DVT 优先展示
key_genes = ['VWF', 'GP1BA', 'ITGB1', 'SELP', 'THBD', 'SERPINE1', 'PIEZO1',
             'PROC', 'PROS1', 'F10', 'F5', 'GP6']
outcomes = ['VTE', 'DVT', 'PE']
colors = {'VTE': '#2c7fb8', 'DVT': '#d62728', 'PE': '#2ca02c'}

fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=True)
for ax, oc in zip(axes, outcomes):
    sub = r[r.outcome == oc]
    # 按 OR 排序，只画关键基因
    sub = sub[sub.GeneSymbol.isin(key_genes)].sort_values('OR_perSDU', ascending=False)
    genes = sub.GeneSymbol.tolist()
    ors = sub.OR_perSDU.tolist()
    lows = sub.OR_low.tolist()
    highs = sub.OR_high.tolist()
    ps = sub.wald_p.tolist()
    y = range(len(genes))
    for yi, (o, lo, hi, p) in zip(y, zip(ors, lows, highs, ps)):
        c = colors[oc] if p < 0.05 else '#999999'
        ax.errorbar(o, yi, xerr=[[o-lo], [hi-o]], fmt='o', color=c,
                    capsize=3, markersize=5, elinewidth=1.2)
    ax.axvline(1.0, color='grey', linestyle='--', linewidth=0.8)
    ax.set_yticks(list(y))
    ax.set_yticklabels(genes, fontsize=9)
    ax.set_xlabel('OR per 1-SD gene expression', fontsize=10)
    ax.set_title(f'{oc}', fontsize=12, fontweight='bold')
    ax.set_xlim(0.3, 2.0)
    # 显著标记
    for yi, p in zip(y, ps):
        if p < 0.05:
            ax.annotate(f'p={p:.3g}', xy=(1.9, yi), fontsize=7, color=colors[oc], va='center')
axes[0].set_ylabel('Gene', fontsize=11)
fig.suptitle('cis-eQTL → 血栓结局 MR（eQTLGen，Wald ratio）', fontsize=13, fontweight='bold')
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = r"C:/Users/Hero/WorkBuddy/07_rheo_MR/06_figures/gene_mr_forest.png"
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f'已保存 -> {out}')
