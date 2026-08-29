#!/usr/bin/env Rscript
# p3_regression.R — structure-preference regression for one comparison
# Usage: Rscript p3_regression.R <sites_structure.tsv> <comparison_name> <direction> <outdir> [covariates.csv]
#   direction: "writer" (loss = log2FC<=-1) or "eraser" (gain = log2FC>=1)
# Models:
#   continuous: diff.log2FC ~ local_pm5 [+ expr_log2FC + expr_level]   (all tested sites)
#   logistic  : changed ~ local_pm5 [+ covariates]                     (changed vs insensitive)
# Cluster-robust (sandwich) SE by gene_id. Seed 42.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4)
tsv <- args[1]; cname <- args[2]; direction <- args[3]; outdir <- args[4]
covfile <- if (length(args) >= 5) args[5] else NA

suppressMessages({library(data.table); library(sandwich); library(lmtest)})
set.seed(42)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

d <- fread(tsv)
d <- d[!is.na(local_pm5) & !is.na(diff.log2FC)]
d[, gene_id := as.factor(gene_id)]

use_cov <- !is.na(covfile) && file.exists(covfile)
if (use_cov) {
  cv <- fread(covfile)
  d <- merge(d, cv, by = "gene_id", all.x = TRUE)
  cat(sprintf("%s: covariates merged for %d/%d sites\n", cname,
              sum(!is.na(d$expr_log2FC)), nrow(d)))
  d <- d[!is.na(expr_log2FC) & !is.na(expr_level)]
}
cat(sprintf("%s: %d sites in analysis\n", cname, nrow(d)))

if (use_cov) {
  fml_cont <- diff.log2FC ~ local_pm5 + expr_log2FC + expr_level
  fml_logi <- changed ~ local_pm5 + expr_log2FC + expr_level
} else {
  fml_cont <- diff.log2FC ~ local_pm5
  fml_logi <- changed ~ local_pm5
}

robust_row <- function(fit, df) {
  V <- vcovCL(fit, cluster = df$gene_id)
  ct <- coeftest(fit, vcov. = V)
  b <- ct["local_pm5", "Estimate"]; se <- ct["local_pm5", "Std. Error"]
  data.frame(term = "local_pm5", beta = b, se_cluster = se,
             z = ct["local_pm5", 3], p = ct["local_pm5", 4],
             ci_lo = b - 1.96 * se, ci_hi = b + 1.96 * se)
}

fit_c <- lm(fml_cont, data = d)
row_c <- robust_row(fit_c, d); row_c$model <- "continuous_log2FC"; row_c$n <- nrow(d)

if (direction == "writer") {
  d[, changed := as.integer(fdr < 0.05 & diff.log2FC <= -1)]
} else {
  d[, changed := as.integer(fdr < 0.05 & diff.log2FC >= 1)]
}
d[, insensitive := as.integer(fdr >= 0.5 & abs(diff.log2FC) < 0.5)]
dl <- d[changed == 1 | insensitive == 1]
cat(sprintf("logistic set: %d changed vs %d insensitive\n",
            sum(dl$changed), sum(dl$insensitive == 1 & dl$changed == 0)))
row_l <- NULL
if (sum(dl$changed) > 50 && sum(dl$changed == 0) > 50) {
  fit_l <- glm(fml_logi, data = dl, family = binomial)
  row_l <- robust_row(fit_l, dl); row_l$model <- "logistic_changed"; row_l$n <- nrow(dl)
}

res <- rbindlist(list(row_c, row_l), fill = TRUE)
res$comparison <- cname
res$covariates <- use_cov
fwrite(res, file.path(outdir, paste0("beta_", cname, ".csv")))
print(res[, c("comparison","model","n","beta","ci_lo","ci_hi","z","p")])

if (direction == "eraser" && nrow(dl) > 100) {
  wt <- wilcox.test(local_pm5 ~ changed, data = dl)
  cat(sprintf("Mann-Whitney p=%.3g | median local_pm5 changed=%.3f insensitive=%.3f\n",
              wt$p.value, dl[changed==1, median(local_pm5)], dl[changed==0, median(local_pm5)]))
}
cat("DONE\n")
