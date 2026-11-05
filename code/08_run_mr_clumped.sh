#!/bin/bash
# 07 号 —— 用严格 LD clump 后的 IV 重跑 MR（复用现有 targets，不重新提取）
# clump 后的 SNP 是原 SNP 的子集，原 targets 文件已包含，无需重新提取
# 用法: bash 08_run_mr_clumped.sh

set -u
BASE="C:/Users/Hero/WorkBuddy/07_rheo_MR"
RSCRIPT="/d/R-4.5.2/bin/x64/Rscript.exe"
CODE="$BASE/00_code"
IV="$BASE/03_iv"
OUT="$BASE/02_raw_outcome"

TRAITS=(HCT PLT MPV RBC PCT)
OUTCOMES=(VTE DVT PE)

echo "=============================================="
echo "07 号 MR（严格 LD clump 后）：流变性状 -> 血栓"
echo "结果 tag 带 _clumped 后缀，输出到 04_mr/"
echo "=============================================="

for t in "${TRAITS[@]}"; do
  ivf="$IV/${t}_iv_clumped.tsv"
  if [ ! -f "$ivf" ]; then
    echo "SKIP $t: 无 clumped IV 文件"; continue
  fi
  for o in "${OUTCOMES[@]}"; do
    outf="$OUT/${o}_targets.tsv"
    if [ ! -f "$outf" ]; then
      echo "SKIP ${t}_to_${o}: 无 targets"; continue
    fi
    tag="${t}_to_${o}_clumped"
    echo "--- $tag ---"
    "$RSCRIPT" "$CODE/03_mr.R" "$ivf" "$outf" "$tag" 2>&1 | grep -aE "IV =|harmonise|MR results|pleiotropy|Steiger|written|MR failed|WARNING|method|Inverse|Weighted|Egger"
    echo ""
  done
done

echo "========== 完成 =========="
echo "clumped MR 结果（tag 带 _clumped）在 04_mr/，与原结果对比即可"
