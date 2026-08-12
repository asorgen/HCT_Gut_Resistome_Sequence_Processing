#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=operams
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=248GB
#SBATCH --time=96:00:00
#SBATCH --array=0-4
#SBATCH --output=/scratch/asorgen/hct_hybrid/logs/operams_%A_%a.log

# ------------------------------------------------------------------------------
# Hybrid (Illumina + ONT) assembly with OPERA-MS for the 5-sample short/long/hybrid
# comparison cohort.
#
# Design note: --contig-file is supplied with the *existing* Duke_short metaWRAP
# (metaSPAdes + megahit) contigs rather than letting OPERA-MS re-assemble the short
# reads with its bundled MEGAHIT. This holds the short-read assembler constant across
# the short and hybrid tracks, so the short-vs-hybrid delta is attributable to the
# ONT reads alone rather than to a change of assembler.
#
# OPERA-MS accepts gzipped short reads but hard-fails on a gzipped long-read file
# (OPERA-MS.pl:257), so the ONT reads are decompressed to scratch first and removed
# on exit.
#
# @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics
# ------------------------------------------------------------------------------

set -euo pipefail

SAMPLES=(D20248PRE D20248D1 D20248D21 D21309PRE D21309D10)
ID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

STUDY=/projects/afodor_research3/asorgen/HCT_Gut_Resistome_Study
SHORT_DIR=${STUDY}/HCT_Gut_Resistome_Data/unprocessed/Duke/Duke_short
ONT_DIR=${STUDY}/Nanopore_Resistome_Pipeline/results
OPERA=/users/asorgen/PROGRAMS/OPERA-MS

WORK=/scratch/asorgen/hct_hybrid
ONT_FLAT=${WORK}/ont_reads/${ID}.fastq
OUTDIR=${WORK}/opera_ms/${ID}

mkdir -p "${WORK}/ont_reads" "${WORK}/logs" "$(dirname "${OUTDIR}")"

R1=${SHORT_DIR}/0.4_host_decontamination/${ID}_1.fastq.gz
R2=${SHORT_DIR}/0.4_host_decontamination/${ID}_2.fastq.gz
CONTIGS=${SHORT_DIR}/1.1_assembly/${ID}/assembly.fasta
ONT_GZ=${ONT_DIR}/${ID}/host_removal/${ID}_host_removed.fastq.gz

for f in "$R1" "$R2" "$CONTIGS" "$ONT_GZ"; do
    [[ -s "$f" ]] || { echo "[ERROR] missing input: $f"; exit 1; }
done

echo "===== ${ID} ====="
echo "[ host      ] $(hostname)"
echo "[ short R1  ] ${R1}"
echo "[ short R2  ] ${R2}"
echo "[ contigs   ] ${CONTIGS}"
echo "[ ONT (gz)  ] ${ONT_GZ}"
echo "[ out       ] ${OUTDIR}"
date
SECONDS=0

# OPERA-MS refuses gzipped long reads -- stage a decompressed copy on scratch.
# Removed on any exit so a failed/cancelled job doesn't leak tens of GB.
cleanup() { rm -f "${ONT_FLAT}"; }
trap cleanup EXIT

if [[ ! -s "${ONT_FLAT}" ]]; then
    echo "[ decompressing ONT reads -> ${ONT_FLAT} ]"
    pigz -dc -p 8 "${ONT_GZ}" > "${ONT_FLAT}" 2>/dev/null || gunzip -c "${ONT_GZ}" > "${ONT_FLAT}"
fi
echo "[ ONT flat size ] $(du -h "${ONT_FLAT}" | cut -f1)"

# OPERA-MS bundles its own tools (bwa, samtools, minimap2, racon, Pilon, MEGAHIT,
# SPAdes, mash) under tools_opera_ms/; its Perl deps live in the user perl5 tree.
export PERL5LIB=/users/asorgen/perl5/lib/perl5
module load R/4.2.2

# The bundled samtools is 0.1.19, linked against libncurses.so.5 / libtinfo.so.5.
# RHEL 9 ships only the .so.6 series, so samtools fails to start -- which manifests
# downstream as "Error in during cov_estimate.pl" with an empty alignment, not as an
# obvious missing-library error. compat_lib/ holds .so.5 -> .so.6 symlinks; samtools
# 0.1.19 references no ncurses symbols in practice (it links them only for `tview`),
# so the version skew is inert here. Swapping in a modern samtools is NOT an option:
# OPERA-MS calls pre-1.0 syntax (`sort - <prefix>`, `sort -no -`) that samtools 1.x
# removed.
export LD_LIBRARY_PATH=/users/asorgen/PROGRAMS/compat_lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

rm -rf "${OUTDIR}"
mkdir -p "${OUTDIR}"

perl ${OPERA}/OPERA-MS.pl \
    --num-processors 32 \
    --short-read1 "${R1}" \
    --short-read2 "${R2}" \
    --long-read "${ONT_FLAT}" \
    --contig-file "${CONTIGS}" \
    --out-dir "${OUTDIR}"

# OPERA-MS emits contigs.polished.fasta when polishing is enabled (default, Pilon),
# otherwise contigs.fasta. Normalize to assembly.fasta so the downstream uniform
# AMR harness can treat all three tracks identically.
if [[ -s "${OUTDIR}/contigs.polished.fasta" ]]; then
    cp "${OUTDIR}/contigs.polished.fasta" "${OUTDIR}/assembly.fasta"
elif [[ -s "${OUTDIR}/contigs.fasta" ]]; then
    cp "${OUTDIR}/contigs.fasta" "${OUTDIR}/assembly.fasta"
else
    echo "[ERROR] no contig output found in ${OUTDIR}"; exit 1
fi

echo "[ final assembly ] $(du -h "${OUTDIR}/assembly.fasta" | cut -f1)"
grep -c '^>' "${OUTDIR}/assembly.fasta" | xargs echo "[ contig count ]"
touch "${OUTDIR}/COMPLETE"

echo "[ elapsed ] $((SECONDS / 3600))h $(((SECONDS % 3600) / 60))m"
date
echo "===== ${ID} DONE ====="
