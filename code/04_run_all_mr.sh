#!/bin/bash
# 07 号 流变基因因果推断 MR —— 主流程
# 暴露(流变性状) x 结局(血栓) 全组合 MR
# 用法: bash 04_run_all_mr.sh

set -u
BASE="C:/Users/Hero/WorkBuddy/07_rheo_MR"
RSCRIPT="/d/R-4.5.2/bin/x64/Rscript.exe"
PY="C:/Users/Hero/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
CODE="$BASE/00_code"
IV="$BASE/03_iv"
OUT="$BASE/02_raw_outcome"
MRD="$BASE/04_mr"
mkdir -p "$MRD"

# 暴露性状
TRAITS=(HCT PLT MPV RBC PCT)
# 结局
OUTCOMES=(VTE DVT PE)

echo "=============================================="
echo "07 号 MR 主分析：流变性状 -> 血栓"
echo "=============================================="

# 合并所有 IV 的 rsID 列表
SNPLIST="$BASE/03_iv/all_iv_snps.txt"
cat "$IV"/*_iv_gws.tsv 2>/dev/null | awk 'NR>1 && $1 ~ /^rs/ {print $1}' | sort -u > "$SNPLIST"
echo "合并 IV SNP 总数: $(wc -l < "$SNPLIST")"

# 提取结局 targets（一次性提取所有 rsID）
for o in "${OUTCOMES[@]}"; do
  "$PY" "$CODE/02_extract_targets.py" "$SNPLIST" "$o"
done

echo ""
echo "========== 逐个组合跑 MR =========="
for t in "${TRAITS[@]}"; do
  ivf="$IV/${t}_iv_gws.tsv"
  if [ ! -f "$ivf" ]; then
    echo "SKIP $t: 无 IV 文件"; continue
  fi
  for o in "${OUTCOMES[@]}"; do
    outf="$OUT/${o}_targets.tsv"
    if [ ! -f "$outf" ]; then
      echo "SKIP ${t}_to_${o}: 无 targets"; continue
    fi
    tag="${t}_to_${o}"
    echo "--- $tag ---"
    "$RSCRIPT" "$CODE/03_mr.R" "$ivf" "$outf" "$tag" 2>&1 | grep -aE "IV =|harmonise|MR results|pleiotropy|Steiger|b |OR|written|MR failed|WARNING" | head -20
  done
done

echo ""
echo "========== 汇总结果 =========="
"$PY" "$CODE/05_summarize.py" 2>&1
