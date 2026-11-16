#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=12:00:00

# Per-sample "fraction of reads assembled" and MAG counts for one short-read cohort.
#
# 1. Fraction of reads assembled: samtools flagstat on the 3.1 binning BAM
#    (metaWRAP binning maps each sample's clean reads back to its own final
#    assembly with bwa mem; unmapped reads are kept in the BAM). Reported on
#    primary alignments, so supplementary/secondary records don't inflate it.
# 2. MAGs: number of reassembled bins (3.3, >= min_completion / <= max_contam),
#    their median completeness/contamination, and the zero-MAG reason from
#    3.2_refine_bins/NO_BINS/<ID> where present.
#
# Writes to ${PROCESSED_ROOT}/${dataset}_tables/:
#   ${dataset}_assembly_read_mapping.tsv
#   ${dataset}_mag_counts.tsv
#
# Usage (submits itself; run from anywhere):
#   bash assembly_read_fraction_and_mags.sh <Duke_short|UNC_short|Heston_short> [parallel_bams=8]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    CONFIG_NAME=${1:?Usage: assembly_read_fraction_and_mags.sh <config> [parallel_bams]}
    PAR=${2:-8}
    CONFIG_FILE="$(cd "${SCRIPT_DIR}/../../configs" && pwd)/${CONFIG_NAME}-read.config"
    [[ -f "$CONFIG_FILE" ]] || { echo "No such config: $CONFIG_FILE"; exit 1; }
    logDir=$(source "$CONFIG_FILE" >/dev/null 2>&1; echo "${datasetDir}/LOGs")
    sbatch --job-name=asmfrac_${CONFIG_NAME} -o ${logDir}/asm_read_fraction_mags_%j.out \
        --export=ALL,AF_CONFIG=${CONFIG_FILE},AF_PAR=${PAR} "${SCRIPT_DIR}/assembly_read_fraction_and_mags.sh"
    exit $?
fi

# ---------------- Job mode ----------------
set -uo pipefail
source "$AF_CONFIG"
module load samtools/1.19

outDir=${PROCESSED_ROOT}/${dataset}_tables
mkdir -p "$outDir"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/asmfrac_${dataset}_XXXX")
trap 'rm -rf "$tmp"' EXIT

ids=$(ls "${binningDir}/COMPLETE" | sort)
echo "${dataset}: $(echo "$ids" | wc -l) samples with 3.1 COMPLETE"

# --- 1. read mapping (parallel flagstat) ---
export binningDir tmp
echo "$ids" | xargs -P "${AF_PAR}" -I{} bash -c '
    ID={}; bam=${binningDir}/${ID}/work_files/${ID}.bam
    if [[ ! -s $bam ]]; then echo -e "${ID}\tNA\tNA\tNA\tno_bam" > ${tmp}/${ID}.map; exit 0; fi
    samtools flagstat -@ 2 -O tsv "$bam" > ${tmp}/${ID}.flagstat 2>/dev/null
    awk -F"\t" -v id="$ID" "
        \$3==\"primary\"        {tot=\$1}
        \$3==\"primary mapped\" {map=\$1}
        \$3==\"properly paired\"{pp=\$1}
        END{ printf \"%s\t%d\t%d\t%.4f\tok\n\", id, tot, map, (tot>0?map/tot:0) }" \
        ${tmp}/${ID}.flagstat > ${tmp}/${ID}.map
'
mapOut=${outDir}/${dataset}_assembly_read_mapping.tsv
{
    echo -e "SampleID\tprimary_reads\tprimary_mapped\tfrac_reads_assembled\tstatus"
    for ID in $ids; do
        if [[ -s ${tmp}/${ID}.map ]]; then cat ${tmp}/${ID}.map; else echo -e "${ID}\tNA\tNA\tNA\tflagstat_failed"; fi
    done
} > "$mapOut"

# Samples marked NO_BINS by 3.1 never got a BAM; label them rather than "no_bam"
for f in ${refinedbinDir}/NO_BINS/*; do
    [[ -e $f ]] || continue
    ID=$(basename "$f")
    grep -q "3.1_binning" "$f" && sed -i "s/^${ID}\tNA\tNA\tNA\tno_bam$/${ID}\tNA\tNA\tNA\tnot_binned_tiny_assembly/" "$mapOut"
done

# --- 2. MAG counts ---
magOut=${outDir}/${dataset}_mag_counts.tsv
{
    echo -e "SampleID\tn_MAGs\tmedian_completeness\tmedian_contamination\tzero_MAG_reason"
    for ID in $ids; do
        dir=${reassemDir}/${ID}/reassembled_bins
        n=$(ls "$dir"/*.fa 2>/dev/null | wc -l)
        stats=${reassemDir}/${ID}/reassembled_bins.stats
        if [[ $n -gt 0 && -s $stats ]]; then
            med=$(awk -F"\t" 'NR>1{c[++k]=$2; x[k]=$3} END{
                    n=asort(c); asort(x); m=int((n+1)/2);
                    if (n%2) printf "%.2f\t%.3f", c[m], x[m]; else printf "%.2f\t%.3f", (c[m]+c[m+1])/2, (x[m]+x[m+1])/2 }' "$stats")
        else
            med=$'NA\tNA'
        fi
        reason=NA
        if [[ $n -eq 0 ]]; then
            nb=${refinedbinDir}/NO_BINS/${ID}
            if   [[ -e $nb ]] && grep -q "3.1_binning" "$nb"; then reason=assembly_too_small_to_bin
            elif [[ -e $nb ]]; then reason=no_bin_passed_refinement
            else reason=no_bins_unmarked
            fi
        fi
        echo -e "${ID}\t${n}\t${med}\t${reason}"
    done
} > "$magOut"

# --- summary to the job log ---
echo; echo "== ${mapOut}"
awk -F"\t" 'NR>1 && $5=="ok"{v[++n]=$4} NR>1{s[$5]++} END{
    asort(v); printf "frac_reads_assembled (n=%d): p10 %.3f  median %.3f  p90 %.3f\n", n, v[int(n*0.1)+1], v[int((n+1)/2)], v[int(n*0.9)];
    for (k in s) printf "  status %s: %d\n", k, s[k] }' "$mapOut"
echo; echo "== ${magOut}"
awk -F"\t" 'NR>1{v[++n]=$2; t+=$2; if($2==0) z[$5]++} END{
    asort(v); printf "MAGs/sample (n=%d): total %d  median %s  p90 %s  max %s\n", n, t, v[int((n+1)/2)], v[int(n*0.9)], v[n];
    for (k in z) printf "  zero-MAG %s: %d\n", k, z[k] }' "$magOut"
