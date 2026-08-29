#!/usr/bin/env python3
"""p2_extract_sites.py — extract m6A sites on longest transcripts for RNAplfold.

Inputs:
  --peaks     exomePeak2 diffPeaks.csv (chr, chromStart, chromEnd, geneID, ...)
  --gtf       annotation GTF (Ensembl Plants / Araport11 gene set)
  --genome    genome fasta
  --outdir    output dir

Outputs:
  sites.tsv         one row per site: site_id, chr, gstart, gend, gene_id, tx_id,
                    tx_pos (1-based position of peak midpoint on spliced transcript),
                    tx_len, plus any extra peak columns (log2FC, fdr)
  transcripts.fa    spliced sequence of each involved longest transcript
"""
import argparse, csv, gzip
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--peaks", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def main():
    args = parse_args()

    # ---- parse GTF: transcripts, exon blocks ----
    tx2gene, tx_strand, tx_chrom = {}, {}, {}
    tx_exons = defaultdict(list)   # tx -> list of (start, end) 1-based inclusive genomic
    with open_maybe_gz(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            attrs = {}
            for kv in f[8].strip().strip(";").split(";"):
                kv = kv.strip()
                if not kv or " " not in kv:
                    continue
                k, v = kv.split(" ", 1)
                attrs[k] = v.strip().strip('"')
            tid, gid = attrs.get("transcript_id"), attrs.get("gene_id")
            if not tid or not gid:
                continue
            tx2gene[tid] = gid
            tx_strand[tid] = f[6]
            tx_chrom[tid] = f[0]
            tx_exons[tid].append((int(f[3]), int(f[4])))

    # longest transcript per gene (by summed exon length)
    gene2txs = defaultdict(list)
    for tid, gid in tx2gene.items():
        gene2txs[gid].append(tid)
    def tx_len(tid):
        return sum(e - s + 1 for s, e in tx_exons[tid])
    gene2best = {g: max(txs, key=tx_len) for g, txs in gene2txs.items()}

    # sort exons in transcript orientation
    for tid in tx_exons:
        tx_exons[tid].sort(reverse=(tx_strand[tid] == "-"))

    def g2tx_pos(tid, gpos):
        """genomic 1-based pos -> transcript 1-based pos, or None if intronic"""
        acc = 0
        for s, e in tx_exons[tid]:
            L = e - s + 1
            if s <= gpos <= e:
                off = gpos - s if tx_strand[tid] == "+" else e - gpos
                return acc + off + 1
            acc += L
        return None

    # ---- read peaks ----
    with open_maybe_gz(args.peaks) as fh:
        rd = csv.DictReader(fh)
        rows = list(rd)
        fieldnames = rd.fieldnames
    # normalize column names
    def col(*names):
        for n in names:
            if n in fieldnames:
                return n
        raise KeyError(f"columns {names} not found in {fieldnames}")
    c_chr, c_s, c_e, c_g = col("chr", "seqnames"), col("chromStart", "start"), col("chromEnd", "end"), col("geneID", "gene_id")

    sites, used_tx = [], set()
    n_skip = 0
    for i, r in enumerate(rows):
        gid = r[c_g]
        tid = gene2best.get(gid)
        if tid is None:
            n_skip += 1; continue
        mid = (int(r[c_s]) + int(r[c_e])) // 2   # peak midpoint as site anchor (exomePeak2 gives no summit)
        tp = g2tx_pos(tid, mid)
        if tp is None:
            n_skip += 1; continue
        used_tx.add(tid)
        sites.append({"site_id": f"S{i}", "chr": r[c_chr], "gstart": r[c_s], "gend": r[c_e],
                      "gene_id": gid, "tx_id": tid, "tx_pos": tp, "tx_len": tx_len(tid),
                      **{k: r[k] for k in fieldnames if k not in (c_chr, c_s, c_e, c_g)}})

    # ---- dump transcript sequences ----
    genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
    n_fa = 0
    import os
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "transcripts.fa"), "w") as out:
        for tid in sorted(used_tx):
            chrom = tx_chrom[tid]
            if chrom not in genome:
                continue
            gseq = genome[chrom].seq
            pieces = [gseq[s - 1:e] for s, e in sorted(tx_exons[tid])]
            seq = sum(pieces, Seq("")) if pieces else None
            if tx_strand[tid] == "-":
                seq = seq.reverse_complement()
            out.write(f">{tid}\n{seq}\n")
            n_fa += 1

    with open(os.path.join(args.outdir, "sites.tsv"), "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=list(sites[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(sites)
    print(f"sites: {len(sites)} (skipped {n_skip}), transcripts: {n_fa}")

if __name__ == "__main__":
    main()
