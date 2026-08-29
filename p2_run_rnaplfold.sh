#!/bin/bash
# p2_run_rnaplfold.sh — run RNAplfold (-u 5 -W 70 -L 40) on transcript fasta, parallel by chunks
# Usage: p2_run_rnaplfold.sh <workdir_with_transcripts.fa> [n_jobs]
# Output: <workdir>/lunp/<txid>_lunp  (per-sequence unpaired-probability profiles)
set -u
source /home/jxun/miniforge3/etc/profile.d/conda.sh
conda activate colabfold

WD="$1"; NJ=${2:-4}   # 11GB RAM WSL: keep parallel low
mkdir -p "$WD/lunp" "$WD/chunks"

# split fasta into chunks (skip already-computed)
cd "$WD"
[ -f chunks/done.split ] || {
  /home/jxun/miniforge3/envs/colabfold/bin/python - "$WD" << 'EOF'
import sys, os
from Bio import SeqIO
wd = sys.argv[1]
recs = list(SeqIO.parse(os.path.join(wd, "transcripts.fa"), "fasta"))
done = {f[:-5] for f in os.listdir(os.path.join(wd, "lunp")) if f.endswith("_lunp")}
recs = [r for r in recs if r.id not in done]
n = 0
for i in range(0, len(recs), 500):
    with open(os.path.join(wd, "chunks", f"chunk_{i//500:04d}.fa"), "w") as out:
        for r in recs[i:i+500]:
            out.write(f">{r.id}\n{r.seq}\n")
    n += 1
open(os.path.join(wd, "chunks", "done.split"), "w").write(str(n))
print(f"chunks: {n} (remaining {len(recs)} transcripts)")
EOF
}

run_chunk() {
  local ck=$1
  local cdir="$WD/lunp/$(basename $ck .fa)"
  mkdir -p "$cdir"
  ( cd "$cdir" && RNAplfold -u 5 -W 70 -L 40 < "$ck" > /dev/null 2>&1 )
  for f in "$cdir"/*_lunp; do
    [ -f "$f" ] && mv "$f" "$WD/lunp/"
  done
  rm -rf "$cdir" "$ck"
  echo "chunk $(basename $ck) done"
}
export -f run_chunk
export WD

ls "$WD"/chunks/chunk_*.fa 2>/dev/null | xargs -P "$NJ" -I{} bash -c 'run_chunk "$@"' _ {}
echo "RNAPLFOLD DONE: $(ls $WD/lunp/*_lunp 2>/dev/null | wc -l) profiles"
