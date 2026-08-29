# Plant_m6A

Cross-species analysis of RNA structure selectivity of the m6A machinery — **Arabidopsis thaliana** analyses A1 (writers) and A2 (eraser ALKBH10B).

Companion to our human study (writers are structure-blind; FTO, but not ALKBH5, prefers structurally accessible substrates). This repo contains the key analysis code for the plant arm.

## Analyses

| ID | Question | Data |
|---|---|---|
| **A1** | Are plant writers (MTA/MTB/FIP37/VIR/HAKAI) structure-blind? | [GSE174573](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174573) (WT vs fip37-4 / vir-1 / hakai-2 MeRIP-seq, 3 reps + paired inputs) |
| **A2** | Does the plant eraser ALKBH10B show structure preference? | [GSE79523](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79523) (WT vs alkbh10b-1, seedlings); GSA [PRJCA005164](https://bigd.big.ac.cn/gsa) (ABA-treated alkbh10b, stress arm) |

## Pipeline

```
download (SRA/ENA/GSA) → fastp + HISAT2 (TAIR10) → exomePeak2 → RNAplfold → regression
```

| Script | Purpose |
|---|---|
| `p1_download_arabidopsis_all.sh` / `p1_download_sra.sh` / `p1_download_ena.sh` / `p1_download_gsa_cra004192.sh` | Download raw data (GEO/SRA, ENA mirror, NGDC GSA) |
| `p1_qc_align.sh` | fastp adapter trimming (required: short-insert MeRIP libraries) + FastQC + HISAT2 → sorted BAM |
| `p1_exomepeak2.R` | Unified peak calling + differential m6A (exomePeak2), mutant vs WT with paired inputs |
| `p1_expression.sh` | Gene expression counts from Input BAMs (featureCounts) |
| `p2_extract_sites.py` | Map peaks to longest transcripts; extract site windows |
| `p2_run_rnaplfold.sh` | RNAplfold `-u 5 -W 70 -L 40` (parameters identical to the human pipeline) |
| `p2_collect_structure.py` | Site accessibility metrics: `center_u1`, `local_pm5` (mean u1 over site ±5 nt) |
| `p3_expression_covariates.R` | Per-comparison expression covariates (expr level, expr log2FC) |
| `p3_regression.R` | Structure-preference regression: continuous `log2FC(m6A) ~ local_pm5 + covariates` + logistic (changed vs insensitive), gene-level cluster-robust SE |

## Conventions

- Reference: **TAIR10 + Araport11/Ensembl Plants annotation**, unified across all datasets
- Structure metric: RNAplfold u1 channel, site midpoint ±5 nt mean (`local_pm5`), parameter-matched to the human pipeline
- Statistics: continuous log2FC regression primary; logistic sensitivity; cluster-robust SE by gene; seed 42

## Key preliminary findings (2026-08)

- **Plant writers are NOT structure-blind**: sites lost in fip37/vir/hakai hypomorphs are significantly enriched in structured (low-accessibility) regions (continuous β = +0.52 to +1.26, all z > 15) — opposite to human writers (β ≈ 0).
- **ALKBH10B shows no structure preference** under normal conditions (β ≈ 0), paralleling human ALKBH5 rather than FTO; the ABA stress arm shows a weak positive preference (suggestive only).

Preliminary: MeRIP-seq resolution (~200 nt peaks, midpoint-anchored), hypomorphic alleles; in vivo structure (Nuc-SHAPE) validation pending.

## Environment

- WSL2; conda envs `m6a` (fastp/HISAT2/featureCounts/samtools) and `colabfold` (ViennaRNA RNAplfold)
- R with exomePeak2, data.table, sandwich; Python 3 with numpy/pandas/pyarrow
