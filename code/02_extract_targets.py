#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
07 号 —— 结局 targets 提取脚本
按 rsID 列表从 FinnGen R12 gz 提取结局统计量。
用法: python 02_extract_targets.py <snp_list.txt> <tag1,tag2,...>
输出: 02_raw_outcome/<tag>_targets.tsv  (chrom pos ref alt rsid pval beta sebeta af)
"""
import gzip, os, sys, time

BASE = r'C:/Users/Hero/WorkBuddy/07_rheo_MR'
RAW = os.path.join(BASE, '02_raw_outcome')
OUTD = RAW
os.makedirs(OUTD, exist_ok=True)

FILES = {
    'VTE': 'finngen_R12_I9_VTE_VTE.gz',
    'DVT': 'finngen_R12_I9_PHLETHROMBDVTLOW_DVT.gz',
    'PE':  'finngen_R12_I9_PULMEMB_PE.gz',
}

def main():
    snpfile = sys.argv[1]
    tags = sys.argv[2].split(',') if len(sys.argv) > 2 else list(FILES)
    want = set()
    for ln in open(snpfile, encoding='utf-8'):
        ln = ln.strip()
        if not ln or ln.startswith('SNP') or ln.startswith('#'):
            continue
        want.add(ln.split()[0])
    print('target SNPs =', len(want), flush=True)

    for tag in tags:
        fn = FILES.get(tag)
        if not fn:
            print('unknown tag', tag); continue
        p = os.path.join(RAW, fn)
        if not os.path.exists(p):
            print(f'  [SKIP] {tag}: 文件不存在 {fn}'); continue
        t0 = time.time()
        found = {}
        n = 0
        with gzip.open(p, 'rt') as f:
            header = f.readline().rstrip('\n').lstrip('#').split('\t')
            ix = {c: i for i, c in enumerate(header)}
            ic = ix.get('chrom') if 'chrom' in ix else ix.get('#chrom')
            for line in f:
                n += 1
                parts = line.rstrip('\n').split('\t')
                rs = parts[ix['rsids']]
                if not rs or rs == '.':
                    continue
                hits = [s for s in rs.split(',') if s in want]
                if not hits:
                    continue
                for s in hits:
                    if s in found:
                        continue
                    found[s] = [parts[ic], parts[ix['pos']], parts[ix['ref']], parts[ix['alt']],
                                s, parts[ix['pval']], parts[ix['beta']], parts[ix['sebeta']],
                                parts[ix['af_alt']]]
        dst = os.path.join(OUTD, f'{tag}_targets.tsv')
        with open(dst, 'w', encoding='utf-8') as g:
            g.write('chrom\tpos\tref\talt\trsid\tpval\tbeta\tsebeta\taf\n')
            for s in sorted(found):
                g.write('\t'.join(found[s]) + '\n')
        print(f'{tag}: scanned {n} lines in {time.time()-t0:.0f}s, found {len(found)} / {len(want)}', flush=True)

if __name__ == '__main__':
    main()
