#!/usr/bin/env python3
"""p2_collect_structure.py — parse RNAplfold lunp profiles, compute site accessibility.

Metrics (aligned with human Paper2 pipeline):
  center_u1 : u1 (channel 1, P(unpaired) over span 1) at the site position
  local_pm5 : mean of u1 over site ±5 nt

Usage: p2_collect_structure.py <workdir>   # workdir contains sites.tsv and lunp/
Output: <workdir>/sites_structure.parquet (+ .tsv for inspection)
"""
import sys, os
import numpy as np
import pandas as pd

wd = sys.argv[1]
sites = pd.read_csv(os.path.join(wd, "sites.tsv"), sep="\t")
prof_dir = os.path.join(wd, "lunp")

cache = {}
def get_u1(txid):
    if txid not in cache:
        f = os.path.join(prof_dir, f"{txid}_lunp")
        if not os.path.exists(f):
            cache[txid] = None
        else:
            try:
                # lunp: exactly 2 header lines; rows = pos, u1..u5; edges may be "NA"
                # (C engine breaks with comment="#" + regex sep -> use skiprows=2)
                df = pd.read_csv(f, sep=r"\s+", header=None, skiprows=2,
                                 na_values=["NA"], engine="c")
                arr = df.to_numpy(dtype=float)
                if arr.shape[0] == 0:
                    raise ValueError("no data rows")
                cache[txid] = arr[:, 1] if arr.shape[1] == 6 else arr[:, 0]
            except Exception:
                cache[txid] = None
    return cache[txid]

center_u1, local_pm5 = [], []
for _, r in sites.iterrows():
    u1 = get_u1(r["tx_id"])
    p = int(r["tx_pos"])  # 1-based
    if u1 is None or p > len(u1):
        center_u1.append(np.nan); local_pm5.append(np.nan); continue
    lo, hi = max(1, p - 5), min(len(u1), p + 5)
    center_u1.append(u1[p - 1])
    local_pm5.append(u1[lo - 1:hi].mean())

sites["center_u1"] = center_u1
sites["local_pm5"] = local_pm5
sites["chr"] = sites["chr"].astype(str)   # chr column mixes int/str across datasets -> unify
sites.to_parquet(os.path.join(wd, "sites_structure.parquet"), index=False)
sites.to_csv(os.path.join(wd, "sites_structure.tsv"), sep="\t", index=False)
ok = sites["local_pm5"].notna().sum()
print(f"annotated: {ok}/{len(sites)} sites")
print(sites["local_pm5"].describe())
