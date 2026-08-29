#!/bin/bash
# Phase 1: fastp trim + QC + HISAT2 alignment for one dataset (PE fastq.gz in $BASE/fastq)
# Usage: p1_qc_align.sh <BASE_DIR> [n_threads]
#   HISAT2_IDX env var overrides the index (default: Arabidopsis TAIR10)
# Output: $BASE/fastq_trim/<s>_{1,2}.fq.gz, $BASE/qc/{fastp,fastqc}, $BASE/bam/<s>.sorted.bam (+.bai)
# Note: fastp adapter trimming is REQUIRED — MeRIP libraries here have short inserts
# (insert peak ~112 bp < 150 bp read length, ~73% reads adapter-contaminated in PRJCA005164;
# untrimmed HISAT2 rate ~35%, trimmed ~95%).
set -u
source /home/jxun/miniforge3/etc/profile.d/conda.sh
conda activate m6a

BASE="$1"; THREADS="${2:-4}"
IDX=${HISAT2_IDX:-/mnt/j/m6A_plant_data/reference/hisat2_index/TAIR10}
SORT_THREADS=3; SORT_MEM=512M   # 11GB RAM WSL: cap sort buffers (was the OOM source)
mkdir -p "$BASE/fastq_trim" "$BASE/qc/fastp" "$BASE/qc/fastqc" "$BASE/bam" "$BASE/logs"

cd "$BASE/fastq"
for f1 in *_1.fastq.gz; do
  s=${f1%_1.fastq.gz}
  f2="${s}_2.fastq.gz"
  [ -f "$f2" ] || { echo "MISSING MATE: $f2" >> "$BASE/logs/align_failed.log"; continue; }
  t1="$BASE/fastq_trim/${s}_1.fq.gz"; t2="$BASE/fastq_trim/${s}_2.fq.gz"

  if [ ! -s "$t1" ] || [ ! -s "$t2" ]; then
    echo "=== fastp $s ==="
    fastp -w 4 -i "$f1" -I "$f2" -o "$t1" -O "$t2" \
      -j "$BASE/qc/fastp/${s}.json" -h "$BASE/qc/fastp/${s}.html" \
      2> "$BASE/logs/${s}.fastp.log" || { echo "$s fastp" >> "$BASE/logs/align_failed.log"; continue; }
  fi

  if [ ! -s "$BASE/qc/fastqc/${s}_1_fastqc.html" ]; then
    fastqc -t 4 -o "$BASE/qc/fastqc" "$t1" "$t2" || echo "$s fastqc" >> "$BASE/logs/align_failed.log"
  fi

  if [ ! -s "$BASE/bam/${s}.sorted.bam" ]; then
    echo "=== aligning $s ==="
    # Strandness intentionally unset at alignment; infer per dataset via featureCounts -s comparison
    if hisat2 -p "$THREADS" -x "$IDX" \
        -1 "$t1" -2 "$t2" 2> "$BASE/logs/${s}.hisat2.log" \
      | samtools sort -m "$SORT_MEM" -@ "$SORT_THREADS" -o "$BASE/bam/${s}.sorted.bam" -; then
      samtools index "$BASE/bam/${s}.sorted.bam"
    else
      echo "$s hisat2" >> "$BASE/logs/align_failed.log"
    fi
  fi
done

multiqc "$BASE/qc" "$BASE/logs" -o "$BASE/qc" -n multiqc_report 2>/dev/null || true
echo "ALIGN DONE: $BASE"
grep -H "overall alignment rate" "$BASE/logs"/*.hisat2.log 2>/dev/null || true
