#!/usr/bin/env Rscript
# p3_expression_covariates.R — per-comparison expression covariates from featureCounts
# Usage: Rscript p3_expression_covariates.R <counts.txt> <out.csv> <wt_cols_csv> <mut_cols_csv>
#   col args are comma-separated SRR/CRR prefixes matching bam column names in counts.txt
# Output: gene_id, expr_level (log2 mean CPM WT), expr_log2FC (mut/WT, pseudo 0.5)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
suppressMessages(library(data.table))

cnt <- fread(args[1], skip = 1)
setnames(cnt, 1, "gene_id")
bamcols <- names(cnt)[-(1:6)]
# normalize column keys: last path component without .sorted.bam
keys <- sub("\\.sorted\\.bam$", "", basename(bamcols))
wt_ids  <- strsplit(args[3], ",", fixed = TRUE)[[1]]
mut_ids <- strsplit(args[4], ",", fixed = TRUE)[[1]]
wt_idx  <- which(keys %in% wt_ids); mut_idx <- which(keys %in% mut_ids)
stopifnot(length(wt_idx) == length(wt_ids), length(mut_idx) == length(mut_ids))

m <- as.matrix(cnt[, ..bamcols])
lib <- colSums(m)
cpm <- t(t(m) / lib * 1e6)
wt_cpm  <- rowMeans(cpm[, wt_idx, drop = FALSE])
mut_cpm <- rowMeans(cpm[, mut_idx, drop = FALSE])

out <- data.table(gene_id = cnt$gene_id,
                  expr_level  = log2(wt_cpm + 1),
                  expr_log2FC = log2((mut_cpm + 0.5) / (wt_cpm + 0.5)))
fwrite(out, args[2])
cat("covariates:", nrow(out), "genes ->", args[2], "\n")
