#!/bin/bash
# Generic SRA download: prefetch + fasterq-dump --split-files + gzip
# Usage: p1_download_sra.sh <BASE_DIR> <SRR1> [SRR2 ...]
set -u
source /home/jxun/miniforge3/etc/profile.d/conda.sh
conda activate sra

BASE="$1"; shift
mkdir -p "$BASE/sra" "$BASE/fastq" "$BASE/tmp"
cd "$BASE/sra"

for srr in "$@"; do
  if [ -s "$BASE/fastq/${srr}_1.fastq.gz" ] && [ -s "$BASE/fastq/${srr}_2.fastq.gz" ]; then
    echo "[skip] $srr already done"; continue
  fi
  # single-end datasets produce only _1; accept that too
  if [ -s "$BASE/fastq/${srr}_1.fastq.gz" ] && [ ! -f "$BASE/fastq/${srr}.fastq" ]; then
    echo "[skip] $srr already done (single-end)"; continue
  fi
  echo "=== $srr : prefetch ==="
  if ! prefetch "$srr" -O "$BASE/sra" --max-size 100G; then
    echo "$srr prefetch" >> "$BASE/failed.log"; continue
  fi
  echo "=== $srr : fasterq-dump ==="
  if ! fasterq-dump "$srr" -O "$BASE/fastq" -e 8 -m 16384 -t "$BASE/tmp" --split-files; then
    echo "$srr dump" >> "$BASE/failed.log"; continue
  fi
  pigz -p 8 "$BASE/fastq/${srr}"*.fastq 2>/dev/null || gzip "$BASE/fastq/${srr}"*.fastq
  rm -rf "$BASE/sra/$srr" "$BASE/sra/${srr}.sra"
  echo "=== $srr : done ==="
done

echo "DATASET DONE: $BASE"
cat "$BASE/failed.log" 2>/dev/null || echo "no failures"
