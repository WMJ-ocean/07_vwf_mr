#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
07 号 B 线（设计 C）：从 eQTLGen 提取剪切响应/血小板黏附基因的 cis-eQTL 工具变量
一次流式扫描 eQTLGen 主表，按 GeneSymbol 提取目标基因的显著 cis-eQTL。

eQTLGen 主表列（制表符分隔，FDR<0.05 已筛）：
Pvalue SNP SNPChr SNPPos AssessedAllele OtherAllele Zscore Gene GeneSymbol
GeneChr GenePos NrCohorts NrSamples FDR BonferroniP

输出：03_iv/gene_eqtl_raw.tsv（每个基因的全部显著 cis-eQTL）
     03_iv/gene_eqtl_iv.tsv（每个基因选最强 cis-eQTL 作为工具变量，用于 Wald ratio）
"""
import gzip, os, csv

BASE = r'C:/Users/Hero/WorkBuddy/07_rheo_MR'
EQTL = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/02_raw/eqtlgen/eQTLGen_cis_eQTL_FDR005.txt.gz'
MAF  = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/02_raw/eqtlgen/eQTLGen_SNP_AF_AlleleB.txt.gz'
OUTD = os.path.join(BASE, '03_iv')
os.makedirs(OUTD, exist_ok=True)

# 目标基因：剪切响应 + 血小板黏附/整合素 + 凝血（与你分子对接靶点 GP1BA/ITGB1 呼应）
TARGET_GENES = [
    # 血小板黏附/整合素（GP1BA/ITGB1 是分子对接靶点）
    'GP1BA', 'GP1BB', 'GP6', 'ITGB1', 'ITGA2B', 'ITGB3',
    # 内皮剪切响应
    'PIEZO1', 'PECAM1', 'SELE', 'SELP', 'SELPLG', 'THBD', 'VWF',
    # 凝血/纤溶
    'SERPINE1', 'F2', 'F5', 'F7', 'F8', 'F10', 'PROC', 'PROS1', 'FGA', 'FGB', 'FGG',
]

def main():
    print('=== 07 B线：从 eQTLGen 提取目标基因 cis-eQTL ===')
    genes_found = {g: [] for g in TARGET_GENES}

    n = 0
    with gzip.open(EQTL, 'rt') as f:
        header = f.readline().rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(header)}
        need = ['Pvalue', 'SNP', 'SNPChr', 'SNPPos', 'AssessedAllele', 'OtherAllele',
                'Zscore', 'Gene', 'GeneSymbol', 'NrCohorts', 'NrSamples', 'FDR']
        for line in f:
            n += 1
            if n % 5_000_000 == 0:
                print(f'  扫描 {n/1e6:.0f}M 行...', flush=True)
            p = line.rstrip('\n').split('\t')
            gs = p[ix['GeneSymbol']]
            if gs not in genes_found:
                continue
            genes_found[gs].append({
                'Pvalue': p[ix['Pvalue']], 'SNP': p[ix['SNP']],
                'SNPChr': p[ix['SNPChr']], 'SNPPos': p[ix['SNPPos']],
                'AssessedAllele': p[ix['AssessedAllele']], 'OtherAllele': p[ix['OtherAllele']],
                'Zscore': p[ix['Zscore']], 'Gene': p[ix['Gene']], 'GeneSymbol': gs,
                'NrCohorts': p[ix['NrCohorts']], 'NrSamples': p[ix['NrSamples']], 'FDR': p[ix['FDR']],
            })
    print(f'扫描完成，共 {n:,} 行')

    # 写原始记录 + 每基因最强 cis-eQTL
    raw_out = os.path.join(OUTD, 'gene_eqtl_raw.tsv')
    iv_out = os.path.join(OUTD, 'gene_eqtl_iv.tsv')
    fieldnames = ['GeneSymbol', 'Gene', 'SNP', 'SNPChr', 'SNPPos', 'AssessedAllele',
                  'OtherAllele', 'Zscore', 'Pvalue', 'FDR', 'NrCohorts', 'NrSamples']

    with open(raw_out, 'w', newline='') as rf, open(iv_out, 'w', newline='') as ivf:
        rw = csv.DictWriter(rf, fieldnames=fieldnames, delimiter='\t')
        iw = csv.DictWriter(ivf, fieldnames=fieldnames, delimiter='\t')
        rw.writeheader(); iw.writeheader()
        for g in TARGET_GENES:
            rows = genes_found[g]
            if not rows:
                print(f'  {g:12s}: 无 cis-eQTL（FDR<0.05）')
                continue
            for r in sorted(rows, key=lambda x: float(x['Pvalue'])):
                rw.writerow(r)
            # 最强 cis-eQTL（最小 P）
            best = min(rows, key=lambda x: float(x['Pvalue']))
            iw.writerow(best)
            print(f'  {g:12s}: {len(rows):4d} 条 cis-eQTL, 最强 {best["SNP"]} Z={float(best["Zscore"]):+.2f} P={best["Pvalue"]}')

    print(f'\n原始记录 -> {raw_out}')
    print(f'每基因最强 IV -> {iv_out}')

    # 汇总
    n_have = sum(1 for g in TARGET_GENES if genes_found[g])
    print(f'\n有 cis-eQTL 的基因: {n_have}/{len(TARGET_GENES)}')

if __name__ == '__main__':
    main()
