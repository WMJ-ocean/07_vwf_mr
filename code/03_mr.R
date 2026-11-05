#!/usr/bin/env Rscript
# ============================================================
# 07 号 流变基因因果推断 MR 主分析
# 流变性状(暴露) -> 血栓结局 MR
# 用法: Rscript 03_mr.R <exposure_iv.tsv> <outcome_targets.tsv> <tag>
#   exposure_iv.tsv: SNP chr pos effect_allele other_allele beta se pval eaf
#   outcome_targets.tsv: chrom pos ref alt rsid pval beta sebeta af
#   (FinnGen: alt = 效应等位)
# ============================================================
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("need: exposure_iv.tsv outcome_targets.tsv tag")
exp_f <- args[1]
out_f <- args[2]
tag   <- args[3]

BASE <- "C:/Users/Hero/WorkBuddy/07_rheo_MR"
OUTD <- file.path(BASE, "04_mr")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

# ---------- 读暴露 ----------
e <- read.delim(exp_f, stringsAsFactors = FALSE, check.names = FALSE)
e <- e[!is.na(e$beta) & !is.na(e$se) & e$se > 0, ]
exp_dat <- data.frame(
  SNP           = e$SNP,
  beta.exposure = as.numeric(e$beta),
  se.exposure   = as.numeric(e$se),
  pval.exposure = as.numeric(e$pval),
  effect_allele.exposure = e$effect_allele,
  other_allele.exposure  = e$other_allele,
  eaf.exposure  = as.numeric(e$eaf),
  exposure      = tag,
  id.exposure   = tag,
  samplesize.exposure = 173480,   # Astle 2016 总样本（近似）
  stringsAsFactors = FALSE
)
exp_dat <- exp_dat[!duplicated(exp_dat$SNP), ]

# ---------- 读结局 ----------
o <- read.delim(out_f, stringsAsFactors = FALSE, check.names = FALSE)
out_dat <- data.frame(
  SNP          = o$rsid,
  beta.outcome = as.numeric(o$beta),
  se.outcome   = as.numeric(o$sebeta),
  pval.outcome = as.numeric(o$pval),
  effect_allele.outcome = o$alt,
  other_allele.outcome  = o$ref,
  eaf.outcome  = as.numeric(o$af),
  outcome      = basename(out_f),
  id.outcome   = basename(out_f),
  samplesize.outcome = 440357,   # FinnGen R12 总样本（近似）
  stringsAsFactors = FALSE
)
out_dat <- out_dat[!duplicated(out_dat$SNP), ]

cat(sprintf("exposure IV = %d ; outcome matched rows = %d\n", nrow(exp_dat), nrow(out_dat)))

# ---------- harmonise ----------
dat <- harmonise_data(exp_dat, out_dat, action = 2)
dat <- dat[dat$mr_keep, ]
cat(sprintf("after harmonise (mr_keep) = %d SNPs\n", nrow(dat)))
if (nrow(dat) < 3) {
  cat("WARNING: <3 SNPs after harmonise, MR unreliable\n")
}

# ---------- MR ----------
meth <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
if (nrow(dat) < 4) meth <- c("mr_ivw", "mr_wald_ratio")

res <- tryCatch(mr(dat, method_list = meth), error = function(e) NULL)

if (!is.null(res)) {
  res$OR      <- exp(res$b)
  res$OR_low  <- exp(res$b - 1.96 * res$se)
  res$OR_high <- exp(res$b + 1.96 * res$se)
  res$F_stat  <- mean((dat$beta.exposure / dat$se.exposure)^2, na.rm = TRUE)
  write.csv(res, file.path(OUTD, sprintf("%s_MR.csv", tag)), row.names = FALSE)

  cat("\n===== MR results: %s =====\n" %>% sprintf(tag))
  print(res[, c("method", "nsnp", "b", "se", "pval", "OR", "OR_low", "OR_high")])

  if (nrow(dat) >= 3) {
    het <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
    if (!is.null(het)) write.csv(het, file.path(OUTD, sprintf("%s_heterogeneity.csv", tag)), row.names = FALSE)
    plei <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
    if (!is.null(plei)) {
      write.csv(plei, file.path(OUTD, sprintf("%s_pleiotropy.csv", tag)), row.names = FALSE)
      cat("\n-- pleiotropy (Egger intercept) p --\n")
      print(plei[, c("exposure","outcome","egger_intercept","pval")])
    }
    ss <- tryCatch(directionality_test(dat), error = function(e) NULL)
    if (!is.null(ss)) {
      write.csv(ss, file.path(OUTD, sprintf("%s_steiger.csv", tag)), row.names = FALSE)
      cat("\n-- Steiger directionality --\n")
      print(ss)
    }
  }
  write.csv(dat, file.path(OUTD, sprintf("%s_harmonised.csv", tag)), row.names = FALSE)
  cat("\nwritten ->", OUTD, "\n")
} else {
  cat("MR failed\n")
}
