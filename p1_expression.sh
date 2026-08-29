#!/bin/bash
# p1_expression.sh — gene expression counts from Input (= RNA-seq) BAMs via featureCounts
# Usage: p1_expression.sh <BASE_DIR> <GTF> <outdir>
# Counts ALL bams in $BASE/bam (Inputs and IPs; downstream uses Input columns).
set -u
source /home/jxun/miniforge3/etc/profile.d/conda.sh
conda activate m6a

BASE="$1"; GTF="$2"; OUT="$3"
mkdir -p "$OUT"
name=$(basename "$BASE")
featureCounts -T 4 -p --countReadPairs -s 0 -t exon -g gene_id \
  -a "$GTF" -o "$OUT/${name}_gene_counts.txt" "$BASE"/bam/*.sorted.bam 2> "$OUT/${name}_featureCounts.log"
echo "EXPR DONE: $name"
tail -2 "$OUT/${name}_featureCounts.log"
