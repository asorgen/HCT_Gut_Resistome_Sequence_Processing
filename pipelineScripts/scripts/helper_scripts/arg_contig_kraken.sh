#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --time=04:00:00
#SBATCH --job-name=arg_contig_kraken

# Kraken2 taxonomy for every ARG-bearing contig (bug-to-bug ARG transfer analysis).
#
# Most ARG hits on mobile genes sit on contigs that never made it into a MAG
# (10-45% binned for the 4.14 short list vs ~69% overall), so MAG membership alone
# gives a host for a minority of them. This classifies the contigs directly.
#
# All ARG contigs are classified, binned ones included: on those, Kraken's call
# can be checked against the MAG's GTDB-Tk taxonomy, which calibrates how far to
# trust Kraken on the unbinned ones.
#
# Inputs : <PROCESSED_ROOT>/<dataset>_tables/<dataset>_arg_hosts.tsv (build_arg_host_tables.py)
#          ${evaluationDir}/<ID>_final_assembly.fasta
# DB     : ${KRAKEN2_DB} from the cohort config (same DB as 2.1_kraken2), default
#          confidence (0), as in 2.1.
# Output : <PROCESSED_ROOT>/<dataset>_tables/<dataset>_arg_contig_kraken.tsv
#          one row per (SampleID, contig): kraken taxid, name, rank, and the
#          NCBI lineage (superkingdom..species) plus the contig's MAG (or NA).
#
# Usage (submits itself; run from anywhere on the login node):
#   bash arg_contig_kraken.sh <log_file> [cohort ...]     default cohorts: Duke_short UNC_short Heston_short

set -uo pipefail
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    log=${1:?Usage: arg_contig_kraken.sh <log_file> [cohort ...]}; shift
    cohorts=${*:-"Duke_short UNC_short Heston_short"}
    sbatch -o "$log" --export=ALL,ARGK_CONFIG_DIR="$(cd "${SCRIPT_DIR}/../../configs" && pwd)",COHORTS="${cohorts}" \
        "${SCRIPT_DIR}/arg_contig_kraken.sh"
    exit $?
fi
CONFIG_DIR=${ARGK_CONFIG_DIR:?}
COHORTS=${COHORTS:?}
PY=/usr/bin/python3
THREADS=${SLURM_CPUS_PER_TASK:-4}

# ARGK_WORK: keep intermediates (fasta, raw kraken output) there; else temp + cleanup
if [[ -n "${ARGK_WORK:-}" ]]; then
    work=$ARGK_WORK
else
    work=${TMPDIR:-/tmp}/argk_${SLURM_JOB_ID}
    trap 'rm -rf "$work"' EXIT
fi
mkdir -p "$work"

# ---------- 1. extract ARG contigs per cohort into one fasta ----------
KDB=""
for c in $COHORTS; do
    read -r ds evalDir procRoot kdb < <(source ${CONFIG_DIR}/${c}-read.config >/dev/null 2>&1; \
        echo "$dataset $evaluationDir $PROCESSED_ROOT $KRAKEN2_DB")
    [[ -z "$KDB" ]] && KDB=$kdb
    [[ "$kdb" == "$KDB" ]] || { echo "Cohorts use different Kraken DBs ($kdb vs $KDB)"; exit 1; }
    hosts=${procRoot}/${ds}_tables/${ds}_arg_hosts.tsv
    [[ -s $hosts ]] || { echo "Missing $hosts - run build_arg_host_tables.py first"; exit 1; }
    echo "${ds}: extracting ARG contigs from $hosts"
    $PY - "$hosts" "$evalDir" "${work}/${ds}.fa" "${work}/${ds}.contigs.tsv" <<'PYEOF'
import csv, sys, os
hosts, evaldir, fa_out, map_out = sys.argv[1:]
want = {}                                   # SampleID -> {contig: MAG}
with open(hosts) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        want.setdefault(r["SampleID"], {})[r["contig"]] = r["MAG"]
n = 0
with open(fa_out, "w") as fo, open(map_out, "w") as mo:
    mo.write("SampleID\tcontig\tMAG\n")
    for sid, contigs in want.items():
        asm = os.path.join(evaldir, f"{sid}_final_assembly.fasta")
        keep = False
        with open(asm) as fh:
            for line in fh:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    keep = name in contigs
                    if keep:
                        fo.write(f">{sid}|{name}\n"); n += 1
                        mo.write(f"{sid}\t{name}\t{contigs[name]}\n")
                elif keep:
                    fo.write(line)
print(f"  {n} contigs written from {len(want)} samples "
      f"({sum(len(v) for v in want.values())} unique contigs requested)")
PYEOF
done

# ---------- 2. kraken2, one run over all cohorts (DB loaded once) ----------
cat ${work}/*.fa > ${work}/all.fa
echo "kraken2 on $(grep -c '>' ${work}/all.fa) contigs, DB ${KDB}"
# conda's activate.d hooks reference unset variables (ADDR2LINE), so relax -u around them
set +u
module load anaconda3/2023.09
source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate metawrap-env
kraken2 --db "$KDB" --threads $THREADS --output ${work}/all.kraken --report ${work}/all.kreport ${work}/all.fa
rc=$?
conda deactivate
set -u
[[ $rc -eq 0 ]] || { echo "kraken2 failed ($rc)"; exit $rc; }

# ---------- 3. lineage from the DB's own NCBI taxonomy; split back per cohort ----------
for c in $COHORTS; do
    read -r ds procRoot < <(source ${CONFIG_DIR}/${c}-read.config >/dev/null 2>&1; echo "$dataset $PROCESSED_ROOT")
    out=${procRoot}/${ds}_tables/${ds}_arg_contig_kraken.tsv
    $PY - "${KDB}/taxonomy" "${work}/all.kraken" "${work}/${ds}.contigs.tsv" "$out" <<'PYEOF'
import sys, csv
taxdir, kraken, contigs, out = sys.argv[1:]
parent, rank, name = {}, {}, {}
with open(f"{taxdir}/nodes.dmp") as fh:
    for line in fh:
        f = line.split("\t|\t")
        parent[f[0]] = f[1]; rank[f[0]] = f[2]
with open(f"{taxdir}/names.dmp") as fh:
    for line in fh:
        f = line.rstrip("\t|\n").split("\t|\t")
        if f[3] == "scientific name":
            name[f[0]] = f[1]
RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
def lineage(t):
    out = {r: "" for r in RANKS}
    seen = 0
    while t in parent and seen < 100:
        if rank.get(t) in out:
            out[rank[t]] = name.get(t, "")
        if t == parent[t]:
            break
        t = parent[t]; seen += 1
    return out
kr = {}
with open(kraken) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        kr[f[1]] = (f[0], f[2], f[3])        # status, taxid, length
cols = ["SampleID", "contig", "MAG", "kraken_status", "kraken_taxid", "kraken_name",
        "kraken_rank", "contig_length"] + [f"ncbi_{r}" for r in RANKS]
n = cl = 0
with open(contigs) as fh, open(out, "w") as fo:
    w = csv.writer(fo, delimiter="\t", lineterminator="\n"); w.writerow(cols)
    for r in csv.DictReader(fh, delimiter="\t"):
        key = f"{r['SampleID']}|{r['contig']}"
        st, tid, ln = kr.get(key, ("NA", "0", ""))
        lin = lineage(tid) if st == "C" else {x: "" for x in RANKS}
        w.writerow([r["SampleID"], r["contig"], r["MAG"], st, tid, name.get(tid, "") if st == "C" else "",
                    rank.get(tid, "") if st == "C" else "", ln] + [lin[x] for x in RANKS])
        n += 1; cl += st == "C"
print(f"  {out}: {n} contigs, {cl} classified ({100*cl/max(n,1):.1f}%)")
PYEOF
done
echo "done"
