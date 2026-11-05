#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
07 号 流变基因因果推断 MR —— 工具变量提取脚本
从 Astle 2016 narrow_form.tsv.gz 提取 genome-wide significant SNP（p<5e-8），
并做位置窗口 + p 值排序的近似 LD 剪枝（1Mb 窗口，r2 近似由位置去重）。

Astle narrow_form 列（制表符分隔）：
VARIANT ID ID_dbSNP49 CHR BP GENPOS REF ALT ALT_MINOR DIRECTION EFFECT SE P MLOG10P ...
其中 ALT = 效应等位（effect allele），REF = 基线等位（other allele），ALT_FREQ = EAF

输出：03_iv/<trait>_iv_gws.tsv
格式：SNP chr pos effect_allele other_allele beta se pval eaf
"""
import gzip, os, sys, csv

BASE = r'C:/Users/Hero/WorkBuddy/07_rheo_MR'
EXP = os.path.join(BASE, '01_raw_exposure')
OUTD = os.path.join(BASE, '03_iv')
os.makedirs(OUTD, exist_ok=True)

P_THRESH = 5e-8
WINDOW = 1_000_000  # 1 Mb

# 性状 : narrow_form 文件名
TRAITS = {
    'HCT': 'hct_N173039_narrow_form.tsv.gz',
    'PLT': 'plt_N166066_narrow_form.tsv.gz',
    'MPV': 'mpv_N164454_narrow_form.tsv.gz',
    'RBC': 'rbc_N172952_narrow_form.tsv.gz',
    'PCT': 'pct_N164339_narrow_form.tsv.gz',
}

def extract(trait, fn):
    path = os.path.join(EXP, fn)
    if not os.path.exists(path):
        print(f'  [SKIP] {trait}: 文件不存在 {fn}')
        return
    sig = []
    with gzip.open(path, 'rt') as f:
        header = f.readline().rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(header)}
        need = ['ID_dbSNP49', 'CHR', 'BP', 'REF', 'ALT', 'EFFECT', 'SE', 'P', 'ALT_FREQ']
        for c in need:
            if c not in ix:
                print(f'  [ERR] {trait}: 缺少列 {c}，实际列: {header[:12]}')
                return
        n = 0
        for line in f:
            n += 1
            parts = line.rstrip('\n').split('\t')
            try:
                pval = float(parts[ix['P']])
            except (ValueError, IndexError):
                continue
            if pval >= P_THRESH:
                continue
            rsid = parts[ix['ID_dbSNP49']]
            if not rsid or rsid == '.' or not rsid.startswith('rs'):
                continue
            try:
                chr_ = parts[ix['CHR']]
                pos = int(parts[ix['BP']])
                ref = parts[ix['REF']]
                alt = parts[ix['ALT']]
                beta = float(parts[ix['EFFECT']])
                se = float(parts[ix['SE']])
                eaf = float(parts[ix['ALT_FREQ']])
            except (ValueError, IndexError):
                continue
            if se <= 0:
                continue
            sig.append(dict(rsid=rsid, chr=chr_, pos=pos, ea=alt, oa=ref,
                            beta=beta, se=se, pval=pval, eaf=eaf))
        print(f'  {trait}: 扫描 {n} 行，显著(p<5e-8) {len(sig)} 个')

    # 位置窗口剪枝：按 chr 分组，窗口内取最小 p
    from collections import defaultdict
    byc = defaultdict(list)
    for r in sig:
        byc[r['chr']].append(r)
    iv = []
    for c, rs in byc.items():
        rs.sort(key=lambda x: x['pos'])
        cur = []
        for r in rs:
            if cur and r['pos'] - cur[-1]['pos'] > WINDOW:
                iv.append(min(cur, key=lambda x: x['pval']))
                cur = [r]
            else:
                cur.append(r)
        if cur:
            iv.append(min(cur, key=lambda x: x['pval']))
    iv.sort(key=lambda x: x['pval'])

    dst = os.path.join(OUTD, f'{trait}_iv_gws.tsv')
    with open(dst, 'w', encoding='utf-8') as g:
        g.write('SNP\tchr\tpos\teffect_allele\tother_allele\tbeta\tse\tpval\teaf\n')
        for r in iv:
            g.write(f"{r['rsid']}\t{r['chr']}\t{r['pos']}\t{r['ea']}\t{r['oa']}\t"
                    f"{r['beta']:.6g}\t{r['se']:.6g}\t{r['pval']:.6g}\t{r['eaf']:.4g}\n")
    # 统计 F 统计量
    f_stats = [(r['beta']/r['se'])**2 for r in iv]
    mean_f = sum(f_stats)/len(f_stats) if f_stats else 0
    print(f'  -> {trait}: IV = {len(iv)}，平均 F = {mean_f:.1f}，minP = {min(r["pval"] for r in iv):.3g}')
    return len(iv)

if __name__ == '__main__':
    print('=== 07 号 工具变量提取 ===')
    summary = []
    for trait, fn in TRAITS.items():
        n = extract(trait, fn)
        summary.append((trait, n if n else 0))
    print('\n性状\tIV数')
    for t, n in summary:
        print(f'{t}\t{n}')
