#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
07 号 —— MR 结果汇总
读取 04_mr/*_MR.csv，汇总 IVW 主结果成一个总表。
"""
import os, csv, glob

BASE = r'C:/Users/Hero/WorkBuddy/07_rheo_MR'
MRD = os.path.join(BASE, '04_mr')

rows = []
for f in sorted(glob.glob(os.path.join(MRD, '*_MR.csv'))):
    tag = os.path.basename(f).replace('_MR.csv', '')
    # 跳过 harmonised/heterogeneity/pleiotropy/steiger
    if any(k in tag for k in ['harmonised', 'heterogeneity', 'pleiotropy', 'steiger']):
        continue
    with open(f, encoding='utf-8') as g:
        rd = list(csv.DictReader(g))
    ivw = [r for r in rd if r.get('method') == 'Inverse variance weighted']
    if not ivw:
        continue
    r = ivw[0]
    rows.append({
        'tag': tag,
        'nsnp': r.get('nsnp', ''),
        'b': float(r.get('b', 'nan')),
        'se': float(r.get('se', 'nan')),
        'pval': float(r.get('pval', 'nan')),
        'OR': float(r.get('OR', 'nan')),
        'OR_low': float(r.get('OR_low', 'nan')),
        'OR_high': float(r.get('OR_high', 'nan')),
        'F_stat': float(r.get('F_stat', 'nan')),
    })

if not rows:
    print('无 MR 结果')
else:
    # 排序：按暴露性状，再按结局
    order = {'HCT': 0, 'PLT': 1, 'MPV': 2, 'RBC': 3, 'PCT': 4}
    rows.sort(key=lambda x: (order.get(x['tag'].split('_to_')[0], 9), x['tag']))

    dst = os.path.join(MRD, 'mr_summary_all.csv')
    with open(dst, 'w', encoding='utf-8', newline='') as g:
        w = csv.writer(g)
        w.writerow(['exposure', 'outcome', 'nsnp', 'b', 'se', 'pval', 'OR', 'OR_low', 'OR_high', 'F_stat'])
        for r in rows:
            exp, out = r['tag'].split('_to_')
            w.writerow([exp, out, r['nsnp'], f"{r['b']:.4f}", f"{r['se']:.4f}",
                        f"{r['pval']:.3g}", f"{r['OR']:.3f}", f"{r['OR_low']:.3f}",
                        f"{r['OR_high']:.3f}", f"{r['F_stat']:.1f}"])

    print(f"{'暴露':<6}{'结局':<6}{'nSNP':>6}{'OR':>9}{'95%CI':>22}{'P':>10}{'F':>8}")
    print('-' * 70)
    for r in rows:
        exp, out = r['tag'].split('_to_')
        ci = f"[{r['OR_low']:.2f}-{r['OR_high']:.2f}]"
        sig = '***' if r['pval'] < 0.001 else ('**' if r['pval'] < 0.01 else ('*' if r['pval'] < 0.05 else ''))
        print(f"{exp:<6}{out:<6}{r['nsnp']:>6}{r['OR']:>9.2f}{ci:>22}{r['pval']:>10.3g}{r['F_stat']:>8.1f} {sig}")
    print(f"\nwritten -> {dst}")
