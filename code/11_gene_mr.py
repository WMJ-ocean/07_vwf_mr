#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
07 号 B 线（设计 C）：基因 cis-eQTL → 血栓结局 Wald ratio MR

复用 02 号 56_gene_mr.py 的换算逻辑：
  eQTLGen 只给 Zscore 和 NrSamples，换算 beta/se：
    denom = sqrt(2 * f * (1-f) * (n + Z^2))
    beta  = Z / denom
    se    = 1 / denom
  f = AssessedAllele 频率（从 MAF 文件查 AlleleB 频率，判断 AssessedAllele 是 A 还是 B）

输入：
  03_iv/gene_eqtl_iv.tsv  （每基因最强 cis-eQTL）
  02_raw_outcome/{VTE,DVT,PE}_targets.tsv
输出：
  04_mr/gene_mr_wald.csv  （所有基因×结局的 Wald ratio）
"""
import os, gzip, math, csv
import pandas as pd

BASE = r'C:/Users/Hero/WorkBuddy/07_rheo_MR'
EQTL = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/02_raw/eqtlgen/eQTLGen_cis_eQTL_FDR005.txt.gz'
MAF  = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/02_raw/eqtlgen/eQTLGen_SNP_AF_AlleleB.txt.gz'
IV   = os.path.join(BASE, '03_iv', 'gene_eqtl_iv.tsv')
OUTD = os.path.join(BASE, '04_mr')
os.makedirs(OUTD, exist_ok=True)

def main():
    print('=== 07 B线：基因 eQTL → 血栓 Wald ratio MR ===')
    iv = pd.read_csv(IV, sep='\t')
    print(f'读入 {len(iv)} 个基因工具变量')

    # 需要 (SNP, GeneSymbol) -> AssessedAllele/OtherAllele
    need = set(zip(iv.SNP, iv.GeneSymbol))
    al = {}
    with gzip.open(EQTL, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t')
        c = {n: i for i, n in enumerate(h)}
        for line in f:
            fl = line.rstrip('\n').split('\t')
            k = (fl[c['SNP']], fl[c['GeneSymbol']])
            if k in need:
                al[k] = (fl[c['AssessedAllele']], fl[c['OtherAllele']])
    print(f'补到 {len(al)}/{len(need)} 条等位信息')
    iv['AssessedAllele'] = [al.get((s, g), (None, None))[0] for s, g in zip(iv.SNP, iv.GeneSymbol)]
    iv['OtherAllele']    = [al.get((s, g), (None, None))[1] for s, g in zip(iv.SNP, iv.GeneSymbol)]

    # MAF 文件：查 AlleleB 频率
    want = set(iv.SNP)
    rows = []
    with gzip.open(MAF, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t')
        iS, iA, iB, iF = h.index('SNP'), h.index('AlleleA'), h.index('AlleleB'), h.index('AlleleB_all')
        for line in f:
            fl = line.split('\t')
            if fl[iS] in want:
                rows.append((fl[iS], fl[iA], fl[iB], float(fl[iF])))
    maf = pd.DataFrame(rows, columns=['SNP', 'AlleleA', 'AlleleB', 'freqB'])
    iv = iv.merge(maf, on='SNP', how='left')
    print(f'MAF 命中 {len(maf)}/{len(want)} 个 SNP')

    # f = AssessedAllele 频率
    same_B = iv.AssessedAllele == iv.AlleleB
    same_A = iv.AssessedAllele == iv.AlleleA
    iv['f'] = None
    iv.loc[same_B, 'f'] = iv.loc[same_B, 'freqB']
    iv.loc[same_A, 'f'] = 1.0 - iv.loc[same_A, 'freqB']
    iv['f'] = pd.to_numeric(iv.f, errors='coerce')

    # 换算 beta/se
    n = pd.to_numeric(iv.NrSamples)
    z = pd.to_numeric(iv.Zscore)
    den = (2.0 * iv.f * (1.0 - iv.f) * (n + z**2)) ** 0.5
    iv['beta_exp_raw'] = z / den
    iv['se_exp'] = 1.0 / den
    # 一致性自检
    chk = (iv.beta_exp_raw / iv.se_exp - z).abs().max()
    print(f'一致性自检 |beta/se - Z| 最大偏差 = {chk:.2e}  ' + ('OK' if chk < 1e-6 else '!! 异常'))

    # 读结局
    outcomes = {}
    for oc in ['VTE', 'DVT', 'PE']:
        p = os.path.join(BASE, '02_raw_outcome', f'{oc}_gene_targets.tsv')
        df = pd.read_csv(p, sep='\t')
        df = df[['rsid', 'alt', 'ref', 'beta', 'sebeta', 'pval', 'af']].rename(
            columns={'rsid': 'SNP', 'alt': 'out_eff', 'ref': 'out_oth',
                     'beta': 'beta_out', 'sebeta': 'se_out', 'pval': 'p_out', 'af': 'af_out'})
        outcomes[oc] = df
        print(f'结局 {oc}: {len(df)} 个 SNP')

    # Wald ratio
    allres = []
    for oc, od in outcomes.items():
        m = iv.merge(od, on='SNP', how='inner')
        if not len(m):
            print(f'  {oc}: 无 SNP 匹配')
            continue
        # 等位方向对齐：AssessedAllele 应与 FinnGen effect allele 一致，否则翻转表达效应
        flip = (m.AssessedAllele != m.out_eff).astype(int) * 2 - 1
        m['concordance'] = ['same' if v > 0 else 'flipped' for v in flip]
        m['beta_exp'] = m.beta_exp_raw * flip
        m['outcome'] = oc

        bx, sx = m.beta_exp, m.se_exp
        by, sy = m.beta_out, m.se_out
        m['wald_b'] = by / bx
        m['wald_se'] = (sy / bx).abs()
        m['wald_z'] = m.wald_b / m.wald_se
        m['wald_p'] = [math.erfc(abs(v) / math.sqrt(2)) if pd.notna(v) else float('nan') for v in m.wald_z]
        m['OR_perSDU'] = [math.exp(v) if pd.notna(v) else float('nan') for v in m.wald_b]
        m['OR_low'] = [math.exp(v - 1.96*sv) if pd.notna(v) else float('nan') for v, sv in zip(m.wald_b, m.wald_se)]
        m['OR_high'] = [math.exp(v + 1.96*sv) if pd.notna(v) else float('nan') for v, sv in zip(m.wald_b, m.wald_se)]
        allres.append(m)

    R = pd.concat(allres, ignore_index=True)
    print(f'\n生成 {len(R)} 条检验 ({R.GeneSymbol.nunique()} 基因 × {R.outcome.nunique()} 结局)')

    def stars(p):
        if pd.isna(p): return 'na'
        return '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 0.05 else 'ns'
    R['sig'] = R.wald_p.apply(stars)

    cols = ['outcome', 'GeneSymbol', 'SNP', 'AssessedAllele', 'OtherAllele', 'out_eff', 'out_oth',
            'concordance', 'f', 'Zscore', 'NrSamples', 'beta_exp', 'se_exp',
            'beta_out', 'se_out', 'p_out', 'wald_b', 'wald_se', 'wald_p', 'OR_perSDU', 'OR_low', 'OR_high', 'sig']
    out = os.path.join(OUTD, 'gene_mr_wald.csv')
    R[cols].to_csv(out, index=False)
    print(f'\n写入 -> {out}')

    # 结果预览
    print('\n=== 结果预览（按结局分组，p 排序）===')
    for oc in ['VTE', 'DVT', 'PE']:
        sub = R[R.outcome == oc].sort_values('wald_p')
        print(f'\n--- {oc} ---')
        print(f'{"gene":<10}{"SNP":<14}{"Z_exp":>7}{"OR/SDU":>8}{"95%CI":>18}{"P_Wald":>10}  sig')
        for r in sub.itertuples():
            if pd.isna(r.wald_b):
                print(f'{r.GeneSymbol:<10}{r.SNP:<14}{r.Zscore:>7.1f}{"NA":>8}')
                continue
            ci = f'({r.OR_low:.2f}-{r.OR_high:.2f})'
            print(f'{r.GeneSymbol:<10}{r.SNP:<14}{r.Zscore:>7.1f}{r.OR_perSDU:>8.3f}{ci:>18}{r.wald_p:>10.2e}  {r.sig}')

if __name__ == '__main__':
    main()
