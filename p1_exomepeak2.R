#!/usr/bin/env Rscript
# exomePeak2 peak calling + differential m6A for one dataset
# Usage: Rscript p1_exomepeak2.R <sample_table.csv> <gff> <genome_fa> <ctrl_label> <treat_label> <outdir>
# sample_table.csv columns: sample,group,type,bam
#   group: genotype/condition label (e.g. Col0, alkbh10b)
#   type:  "IP" or "Input"
#   bam:   path to sorted bam
# Output: <outdir>/peaks.bed, <outdir>/diff_peaks.csv, <outdir>/exomePeak2_output/ (full objects)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6)
tab_file <- args[1]; gff <- args[2]; genome <- args[3]
ctrl <- args[4]; treat <- args[5]; outdir <- args[6]

suppressMessages({library(exomePeak2); library(data.table)})
set.seed(42)

tab <- fread(tab_file)
stopifnot(all(c("sample","group","type","bam") %in% names(tab)))
tab[, type := toupper(type)]
cat("Sample table:\n"); print(tab[, .(sample, group, type)])

get_bams <- function(g, t) tab[group == g & type == toupper(t), bam]
bam_ip_c   <- get_bams(ctrl, "IP")
bam_in_c   <- get_bams(ctrl, "INPUT")
bam_ip_t   <- get_bams(treat, "IP")
bam_in_t   <- get_bams(treat, "INPUT")
stopifnot(length(bam_ip_c) > 0, length(bam_in_c) > 0,
          length(bam_ip_t) > 0, length(bam_in_t) > 0)
cat(sprintf("ctrl(%s): %d IP + %d Input | treat(%s): %d IP + %d Input\n",
            ctrl, length(bam_ip_c), length(bam_in_c),
            treat, length(bam_ip_t), length(bam_in_t)))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

res <- exomePeak2(
  bam_ip            = bam_ip_c,
  bam_input         = bam_in_c,
  bam_ip_treated    = bam_ip_t,
  bam_input_treated = bam_in_t,
  gff               = gff,
  genome            = genome,
  strandness        = "unstrand",   # featureCounts -s check: unstranded libraries
  fragment_length   = 100,          # fastp insert peak ~112 bp
  parallel          = 1,            # 11GB RAM WSL: single worker (RSS ~6GB); 2 workers thrash swap
  diff_p_cutoff     = 1,            # keep ALL tested peaks (default 0.01 drops nonsignificant!
                                    # regression needs the full site universe, not just hits)
  save_dir          = outdir,
  experiment_name   = "exomePeak2_output"
)
saveRDS(res, file.path(outdir, "exomePeak2_result.rds"))  # full object, never lose again

# ---- export peak list + differential results (slot names vary by version) ----
cat("result class:", class(res), "\n")
cat("result slots:", paste(slotNames(res), collapse = ", "), "\n")
peaks <- tryCatch(res$peak, error = function(e) NULL)
if (!is.null(peaks) && length(peaks) > 0) {
  rtracklayer::export(peaks, file.path(outdir, "peaks.bed"), format = "BED")
  cat("peaks:", length(peaks), "\n")
}
diff <- tryCatch(res$diffPeak, error = function(e) NULL)
if (is.null(diff)) diff <- tryCatch(res$diff_peak, error = function(e) NULL)
if (!is.null(diff) && length(diff) > 0) {
  df <- as.data.frame(diff)
  fwrite(df, file.path(outdir, "diff_peaks.csv"))
  cat("diff peaks:", nrow(df), "\n")
  sig <- sum(!is.na(df$padj) & df$padj < 0.05 & df$log2FC > 1)
  cat(sprintf("hyper in %s (padj<0.05 & log2FC>1): %d\n", treat, sig))
} else {
  cat("NOTE: diff peaks not found in object slots; inspect exomePeak2_output dir CSVs\n")
  print(list.files(file.path(outdir, "exomePeak2_output"), recursive = TRUE))
}
cat("DONE\n")
