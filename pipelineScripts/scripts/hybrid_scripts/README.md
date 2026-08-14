# Short vs. long vs. hybrid comparison

A five-sample, three-track comparison of ARG recovery from the same stool samples
sequenced on both Illumina and Oxford Nanopore:

| Track | Assembly | Reads |
|---|---|---|
| `short` | metaWRAP (metaSPAdes + megahit) — existing Duke_short pipeline output | Illumina |
| `long` | metaFlye — existing Nanopore_Resistome_Pipeline output | ONT |
| `hybrid` | OPERA-MS, scaffolding the **short track's own contigs** with ONT reads | both |

## Cohort

`D20248PRE`, `D20248D1`, `D20248D21`, `D21309PRE`, `D21309D10` — two patients,
both baselines included.

These five were chosen because their host-removed ONT inputs (2.1–8.4 GB gzipped)
are the ones that assembled cleanly. There is a sharp cliff above them: the next
smallest ONT input is 73 GB, and every sample at that scale either OOM'd or hit a
wall-clock limit. `D20248D28` (155 GB) is deliberately excluded — it has failed four
separate metaFlye attempts across ~14 days of wall clock.

## Why ARGs are re-called rather than reused

The two source pipelines do **not** call ARGs the same way, so their outputs on disk
are not comparable to each other:

| | `short_scripts/5.5_AA_amr_assembly.sh` | Nanopore `rules/amr_aa.smk` |
|---|---|---|
| ORF calling | `prodigal` (single-genome mode) | `prodigal -p meta` |
| RGI DB load | `rgi load --card_json` | `rgi load --card_json --card_annotation … --local` |

The second row matters: without `--card_annotation … --local`, a packaged protein
reference that is stale relative to `card.json` makes every hit resolve to an absent
model ID and be **silently dropped**, yielding a 0-ARG report with no error. The
short-read calls on disk predate that fix.

`uniform_amr_compare.sbatch` therefore re-calls ARGs across all three assemblies with
one identical stack (`prodigal -p meta` → AMRFinder+ `-p` / RGI `-t protein --local`),
run out of the same pinned conda environment. Differences in the results are then
attributable to the assemblies — i.e. to read type — rather than to tooling. All three
inputs are contig-length filtered at ≥1500 bp before ORF calling.

## Hybrid design: `--contig-file`, not a fresh MEGAHIT

OPERA-MS will assemble the short reads itself with MEGAHIT by default. Instead it is
given the short track's existing metaWRAP contigs via `--contig-file`, which holds the
short-read assembler constant between the `short` and `hybrid` tracks. The
short-vs-hybrid delta is then the contribution of the ONT reads alone, rather than
ONT plus a change of assembler.

## Environment notes (both were hard failures)

**OPERA-MS needs an ncurses 5 shim on RHEL 9.** Its bundled `samtools` is 0.1.19,
linked against `libncurses.so.5`/`libtinfo.so.5`; RHEL 9 ships only the `.so.6`
series. The failure surfaces as `Error in during cov_estimate.pl` — an empty alignment
— rather than as a missing-library error, because samtools' failure to start is
swallowed by a pipeline. `/users/asorgen/PROGRAMS/compat_lib` holds `.so.5 → .so.6`
symlinks and is put on `LD_LIBRARY_PATH`; samtools 0.1.19 references no ncurses
symbols in practice (they are linked only for `tview`), so the version skew is inert.

Substituting a modern samtools is **not** an option: OPERA-MS calls pre-1.0 syntax
(`sort - <prefix>`, `sort -no -`) that samtools 1.x removed.

**Snakemake conda envs currently cannot be rebuilt.** Edits to `envs/*.yaml` on
2026-07-28 changed every environment hash, so Snakemake wants to rebuild all nine —
and conda fails mid-transaction on the project NFS mount (96% full), leaving rolled-back
environments. These scripts sidestep it by invoking the already-built environments by
absolute path, which has the side benefit of pinning identical tool versions across
tracks. The durable fix is `--conda-prefix` on scratch, plus freeing space on
`/projects/afodor_research3`.

## Running

```bash
# 1. hybrid assemblies (OPERA-MS), ~25-40 min/sample
sbatch --array=0-4 opera_ms_compare5.sbatch

# 2. uniform ARG calling, one track at a time
sbatch --array=0-4 uniform_amr_compare.sbatch short
sbatch --array=0-4 uniform_amr_compare.sbatch long
sbatch --array=0-4 uniform_amr_compare.sbatch hybrid

# 3. comparison tables
python3 compare_tracks.py \
    --amr-root /scratch/asorgen/hct_compare/amr \
    --outdir   /scratch/asorgen/hct_compare/summary
```

Hybrid assemblies are archived to
`HCT_Gut_Resistome_Data/unprocessed/Duke/Duke_hybrid/1.1_assembly/<sample>/assembly.fasta`.

## Outputs

`compare_tracks.py` writes `assembly_stats.tsv`, `arg_counts.tsv`,
`gene_matrix_{amrfinder,rgi}.tsv`, `track_overlap.tsv` (pairwise shared/unique counts
and Jaccard, per sample and caller), `track_unique_genes.tsv`, and
`drug_class_by_track.tsv`.
