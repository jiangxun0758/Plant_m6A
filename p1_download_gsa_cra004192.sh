#!/bin/bash
# Download GSA CRA004192 (PRJCA005164, Tang et al. 2021, alkbh10b ABA MeRIP) from China NGDC mirror
# 8 runs PE150: Col0 input/IP x2, alkbh10b input/IP x2
set -u
source /home/jxun/miniforge3/etc/profile.d/conda.sh
conda activate m6a

BASE=/mnt/j/m6A_plant_data/PRJCA005164
mkdir -p "$BASE/fastq"
cd "$BASE/fastq"

for crr in CRR283802 CRR283803 CRR283804 CRR283805 CRR283806 CRR283807 CRR283808 CRR283809; do
  for pair in "f1 1" "r2 2"; do
    set -- $pair; src=$1; mate=$2
    out="${crr}_${mate}.fastq.gz"
    if [ -s "$out" ] && gzip -t "$out" 2>/dev/null; then echo "[skip] $out ok"; continue; fi
    url="https://download.cncb.ac.cn/gsa2/CRA004192/${crr}/${crr}_${src}.fastq.gz"
    echo "=== downloading $out ==="
    if aria2c -x 8 -s 8 -k 4M --max-tries=5 --retry-wait=10 --summary-interval=0 \
         --console-log-level=error -o "$out" "$url" && gzip -t "$out" 2>/dev/null; then
      echo "[done] $out"
    else
      echo "$out" >> "$BASE/failed.log"; rm -f "$out"
    fi
  done
done
echo "GSA CRA004192 DONE"
cat "$BASE/failed.log" 2>/dev/null || echo "no failures"
