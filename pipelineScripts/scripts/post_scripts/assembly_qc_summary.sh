#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=assembly_qc_summary
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=02:00:00
# =============================================================================
# assembly_qc_summary.sh — per-sample assembly QC table for one short-read cohort
#
# Usage (from the HCT_Gut_Resistome_Sequence_Processing root):
#   sbatch --output=<log> pipelineScripts/scripts/post_scripts/assembly_qc_summary.sh <config>
#     config   config name (e.g. Heston_short) or path to a *-read.config
#
# For every sample with a 1.2_evaluation COMPLETE flag, writes one row with:
#   - QUAST stats for the metaWRAP final assembly (metaSPAdes + Megahit, contigs >= 1 kb):
#       1.1_assembly/<ID>/QUAST_out/report.tsv
#   - The same stats recomputed on the >= 1.5 kb assembly that feeds 5.5/binning:
#       1.2_evaluation/<ID>_final_assembly.fasta
# Output: ${PROCESSED_ROOT}/${dataset}_tables/${dataset}_assembly_qc.tsv
#
# "Fraction of reads assembled" is NOT computed here (needs a read-mapping step).
# =============================================================================

set -uo pipefail

CONFIG_ARG=${1:?Usage: assembly_qc_summary.sh <config>}
if [[ -f "$CONFIG_ARG" ]]; then CONFIG_FILE=$CONFIG_ARG
else CONFIG_FILE=pipelineScripts/configs/${CONFIG_ARG}-read.config; fi
[[ -f "$CONFIG_FILE" ]] || { echo "No such config: ${CONFIG_FILE}"; exit 1; }

source "$CONFIG_FILE" >/dev/null 2>&1
: "${dataset:?config did not define dataset}" "${datasetDir:?}" "${PROCESSED_ROOT:?}"

asmDir=${datasetDir}/1.1_assembly
evalDir=${datasetDir}/1.2_evaluation
outDir=${PROCESSED_ROOT}/${dataset}_tables
out=${outDir}/${dataset}_assembly_qc.tsv
mkdir -p "$outDir"

echo "dataset=${dataset}  datasetDir=${datasetDir}"
echo "output=${out}"

# Stats from a FASTA: n_contigs, total_bp, largest, N50, L50, GC%
fasta_stats() {
    awk '
        /^>/ { if (len) L[++n]=len; len=0; next }
        { len += length($0); s=toupper($0); gc += gsub(/[GC]/, "", s) }
        END {
            if (len) L[++n]=len
            if (n == 0) { print "0\t0\t0\t0\t0\tNA"; exit }
            tot=0; for (i=1;i<=n;i++) tot+=L[i]
            m=asort(L, S); cum=0; n50=0; l50=0    # asort is ascending (gawk); walk from the top
            for (i=m;i>=1;i--) { cum+=S[i]; l50++; if (cum >= tot/2) { n50=S[i]; break } }
            printf "%d\t%d\t%d\t%d\t%d\t%.2f\n", n, tot, S[m], n50, l50, 100*gc/tot
        }' "$1"
}

# Pull one field from a QUAST report.tsv (2-column: metric<TAB>value)
quast_field() { awk -F'\t' -v k="$2" '$1==k {print $2; f=1} END {if (!f) print "NA"}' "$1"; }

printf "SampleID\tquast_contigs_ge1kb\tquast_total_length\tquast_largest_contig\tquast_N50\tquast_L50\tquast_GC\tfinal1500_contigs\tfinal1500_total_length\tfinal1500_largest_contig\tfinal1500_N50\tfinal1500_L50\tfinal1500_GC\n" > "$out"

n=0; miss_q=0; miss_f=0
for flag in "$evalDir"/COMPLETE/*; do
    ID=$(basename "$flag")
    q=${asmDir}/${ID}/QUAST_out/report.tsv
    f=${evalDir}/${ID}_final_assembly.fasta

    if [[ -s "$q" ]]; then
        qrow=$(printf "%s\t%s\t%s\t%s\t%s\t%s" \
            "$(quast_field "$q" '# contigs')" "$(quast_field "$q" 'Total length')" \
            "$(quast_field "$q" 'Largest contig')" "$(quast_field "$q" 'N50')" \
            "$(quast_field "$q" 'L50')" "$(quast_field "$q" 'GC (%)')")
    else
        qrow=$'NA\tNA\tNA\tNA\tNA\tNA'; miss_q=$((miss_q+1))
    fi

    if [[ -s "$f" ]]; then frow=$(fasta_stats "$f")
    else frow=$'NA\tNA\tNA\tNA\tNA\tNA'; miss_f=$((miss_f+1)); fi

    printf "%s\t%s\t%s\n" "$ID" "$qrow" "$frow" >> "$out"
    n=$((n+1))
done

echo "samples=${n}  missing_QUAST=${miss_q}  missing_final_assembly=${miss_f}"
echo "=== DONE: ${dataset} ==="
