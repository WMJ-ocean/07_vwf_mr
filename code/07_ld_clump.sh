#!/usr/bin/env bash
# ============================================================
# 07 号 严格 LD clump 加固
# 用 plink 1.9 + 1000G Phase3 EUR 参考面板（HapMap3 交集）
# 对 5 个血细胞性状的 IV 做真实 LD r² clump（r²<0.001、窗口 10Mb）
# 替换 01_extract_iv.py 里的 1Mb 位置近似剪枝
#
# 输出：03_iv/<trait>_iv_clumped.tsv（clump 后独立 SNP 的 IV）
# 用法: bash 07_ld_clump.sh
# ============================================================
set -u

BASE="C:/Users/Hero/WorkBuddy/07_rheo_MR"
TOOLS="$BASE/tools"
PLINK="$TOOLS/plink.exe"
REFDIR="$TOOLS/1000G_EUR_Phase3_plink"   # 解压后的参考面板目录
IV="$BASE/03_iv"

# clump 参数（严格标准）
CLUMP_R2=0.001
CLUMP_KB=10000
CLUMP_P1=5e-8        # 索引 SNP 显著性阈值（与 IV 提取一致）

TRAITS=(HCT PLT MPV RBC PCT)

echo "=============================================="
echo "07 号 严格 LD clump（plink + 1000G EUR）"
echo "r2 < $CLUMP_R2, window = ${CLUMP_KB} kb, p < $CLUMP_P1"
echo "=============================================="

# ---------- 1. 合并参考面板（22 条染色体 -> 单个 bfile）----------
# 已手动合并为 1000G_EUR_merged.{bed,bim,fam}（9,997,231 SNP × 489 人）
MERGED="$REFDIR/../1000G_EUR_merged"
if [ ! -f "$MERGED.bed" ]; then
  echo "[ERR] 未找到合并参考面板 $MERGED.bed，请先合并"
  exit 1
fi
echo "参考面板: $MERGED.bed"

# ---------- 2. 对每个性状做 clump ----------
for t in "${TRAITS[@]}"; do
  ivf="$IV/${t}_iv_gws.tsv"
  if [ ! -f "$ivf" ]; then
    echo "SKIP $t: 无 IV 文件"; continue
  fi

  echo ""
  echo "========== $t =========="
  n_orig=$(($(wc -l < "$ivf") - 1))
  echo "原始 IV 数: $n_orig"

  # 生成 plink clump 输入（SNP CHR BP P，需带 header）
  # IV 列：SNP chr pos effect_allele other_allele beta se pval eaf → pval 是第 8 列
  tmp_p="$BASE/tools/${t}_clump_input.txt"
  awk 'NR>1 {print $1, $2, $3, $8}' "$ivf" > "$tmp_p"
  # 加 header
  sed -i '1i SNP CHR BP P' "$tmp_p" 2>/dev/null || { echo "SNP CHR BP P" > "$tmp_p.tmp"; cat "$tmp_p" >> "$tmp_p.tmp"; mv "$tmp_p.tmp" "$tmp_p"; }

  # 跑 plink clump
  out="$BASE/tools/${t}_clump"
  "$PLINK" --bfile "$MERGED" \
    --clump "$tmp_p" \
    --clump-r2 "$CLUMP_R2" \
    --clump-kb "$CLUMP_KB" \
    --clump-p1 "$CLUMP_P1" \
    --clump-p2 1 \
    --out "$out" 2>&1 | tail -8

  if [ ! -f "$out.clumped" ]; then
    echo "  [WARN] $t: clump 未产出 .clumped 文件"
    continue
  fi

  # 读取 clump 结果里的索引 SNP（.clumped 格式：CHR F SNP BP P ...，rsID 在第 3 列）
  n_clumped=$(awk 'NR>1 && $3 ~ /^rs/ {print $3}' "$out.clumped" | sort -u | wc -l)
  echo "clump 后独立 SNP 数: $n_clumped"

  # 生成 keep 列表 = clump 索引 SNP + 参考面板缺失的 top variant（单独保留）
  awk 'NR>1 && $3 ~ /^rs/ {print $3}' "$out.clumped" | sort -u > "$BASE/tools/${t}_keep_snps.txt"
  # 参考面板缺失的 SNP：原 IV 中不在参考面板的 rsID
  awk 'NR>1 {print $1}' "$ivf" | sort -u > "$BASE/tools/${t}_all_iv.txt"
  comm -23 "$BASE/tools/${t}_all_iv.txt" "$BASE/tools/ref_rsids.txt" > "$BASE/tools/${t}_missing_ref.txt"
  n_missing=$(wc -l < "$BASE/tools/${t}_missing_ref.txt")
  if [ "$n_missing" -gt 0 ]; then
    cat "$BASE/tools/${t}_missing_ref.txt" >> "$BASE/tools/${t}_keep_snps.txt"
    echo "  参考面板缺失的 top variant（单独保留）: $n_missing 个"
  fi
  sort -u "$BASE/tools/${t}_keep_snps.txt" -o "$BASE/tools/${t}_keep_snps.txt"

  # 输出 clumped IV（保留原格式：SNP chr pos effect_allele other_allele beta se pval eaf）
  # 先写 header，再从原 IV 过滤 keep 列表中的 SNP
  head -1 "$ivf" > "$IV/${t}_iv_clumped.tsv"
  awk 'NR==FNR{keep[$1]=1; next} ($1 in keep)' \
    "$BASE/tools/${t}_keep_snps.txt" "$ivf" >> "$IV/${t}_iv_clumped.tsv"

  n_out=$(($(wc -l < "$IV/${t}_iv_clumped.tsv") - 1))
  echo "  -> ${t}_iv_clumped.tsv: $n_out 个独立 IV（原 $n_orig，clump 后独立 $n_clumped + 缺失保留 $n_missing）"

  # 计算 clump 后 F 统计量
  awk -F'\t' 'NR>1 {s+=($6/$7)^2; n++} END {if(n>0) printf "  clump 后平均 F = %.1f (n=%d)\n", s/n, n}' "$IV/${t}_iv_clumped.tsv"
done

echo ""
echo "=============================================="
echo "严格 LD clump 完成。clump 后 IV 在 03_iv/<trait>_iv_clumped.tsv"
echo "下一步：用 clumped IV 重跑 MR（见 08_run_mr_clumped.sh）"
echo "=============================================="
