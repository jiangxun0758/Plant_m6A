#!/bin/bash
# ENA fastq direct download with aria2c (fast mirror for NCBI-slow networks)
# Usage: p1_download_ena.sh <BASE_DIR> <SRR1> [SRR2 ...]
# ENA URL rule: vol1/fastq/<SRR+first3digits>/[subdir]/<acc>/<acc>_1.fastq.gz
#   subdir: 7-digit accession -> 00<last digit>; 8-digit -> 0<last two digits>
set -u
source /home/jxun/miniforge3/etc/profile.d/conda.sh
conda activate m6a

BASE="$1"; shift
mkdir -p "$BASE/fastq"
cd "$BASE/fastq"

ena_url() {
  local acc=$1
  local prefix=${acc:0:6}          # e.g. SRR294
  local digits=${acc:3}            # e.g. 29462665
  local nd=${#digits}
  local sub=""
  if [ "$nd" -eq 7 ]; then sub="/00${digits: -1}"; fi
  if [ "$nd" -eq 8 ]; then sub="/0${digits: -2}"; fi
  if [ "$nd" -ge 9 ]; then sub="/0${digits: -2}"; fi  # rare; adjust if needed
  echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}${sub}/${acc}"
}

for srr in "$@"; do
  url_base=$(ena_url "$srr")
  for mate in 1 2; do
    out="${srr}_${mate}.fastq.gz"
    if [ -s "$out" ] && gzip -t "$out" 2>/dev/null; then
      echo "[skip] $out ok"; continue
    fi
    echo "=== downloading $out ==="
    if ! aria2c -x 16 -s 16 -k 8M -c --max-tries=5 --retry-wait=10 --summary-interval=60 \
         --console-log-level=warn -o "$out" "${url_base}/${out}"; then
      echo "$out" >> "$BASE/failed.log"; continue
    fi
    if ! gzip -t "$out" 2>/dev/null; then
      echo "$out corrupt" >> "$BASE/failed.log"; rm -f "$out"
    fi
  done
  echo "=== $srr done ==="
done

echo "DATASET DONE: $BASE"
cat "$BASE/failed.log" 2>/dev/null || echo "no failures"
ls -lh "$BASE/fastq" | tail -5
