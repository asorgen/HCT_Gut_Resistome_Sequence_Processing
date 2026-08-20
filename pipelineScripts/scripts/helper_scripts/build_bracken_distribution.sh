#!/bin/bash
#SBATCH --job-name=build_bracken_distrib
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180GB
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/asorgen/bracken_build/build_bracken_distrib_%j.out
#SBATCH --mail-type=END,FAIL

# Build an additional Bracken k-mer distribution for the Kraken2 DB.
#
#   sbatch build_bracken_distribution.sh [READ_LEN]      (default 125)
#
# Why: Bracken's kmer_distrib is built for a specific read length, and the DB
# only shipped 100/150/200/250. Duke and UNC are 150-151bp so they match, but the
# Heston cohort (PRJNA890666) is largely ~126bp NextSeq — 260 of its samples sit
# in the 125-149bp band, where 2.1_kraken2.sh must otherwise fall back to the
# 100mer distribution. A 125mer build gives those samples a matched distribution.
#
# Cost: the expensive step is classifying the whole 152 GB library against the
# 76 GB hash to produce database.kraken. That is cached — once it exists, further
# read lengths are cheap, so building e.g. 75mers later costs minutes not hours.
#
# Writes: everything large goes to scratch via a symlink farm; the real DB is
# only touched at the end, and only to add the ~3 MB kmer_distrib file.

# Deliberately no `set -u`: conda's activation scripts reference unset variables
# and would abort the job before kraken2 ever runs.
set -o pipefail

# SLURM copies the submitted script to a spool dir, so ${BASH_SOURCE[0]} does NOT
# point back into the repo — resolve the pipeline root explicitly instead.
# Order: PIPELINE_ROOT env override -> the dir sbatch was invoked from -> default.
DEFAULT_ROOT=/projects/afodor_research3/asorgen/HCT_Gut_Resistome_Study/HCT_Gut_Resistome_Sequence_Processing

if [[ -n "${PIPELINE_ROOT:-}" && -f "${PIPELINE_ROOT}/pipelineScripts/configs/private.config" ]]; then
    REPO_ROOT=${PIPELINE_ROOT}
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/pipelineScripts/configs/private.config" ]]; then
    REPO_ROOT=${SLURM_SUBMIT_DIR}
elif [[ -f "${DEFAULT_ROOT}/pipelineScripts/configs/private.config" ]]; then
    REPO_ROOT=${DEFAULT_ROOT}
else
    echo "Could not locate pipelineScripts/configs/private.config."
    echo "Submit from the sequence-processing root, or run:"
    echo "  sbatch --export=ALL,PIPELINE_ROOT=/path/to/HCT_Gut_Resistome_Sequence_Processing $0"
    exit 1
fi

source "${REPO_ROOT}/pipelineScripts/configs/private.config"

READ_LEN=${1:-125}
KMER_LEN=35          # matches the DB build: opts.k2d reports k=35, l=31

ROOT=${HPC_PROJECTS}/HCT_Gut_Resistome_Study/HCT_Gut_Resistome_Sequence_Processing
REAL_DB=${ROOT}/databases/KRAKEN2-TESSA-DB
WORK=${HPC_SCRATCH}/bracken_build
WORK_DB=${WORK}/KRAKEN2-TESSA-DB
BRACKEN=${HPC_HOME}/PROGRAMS/Bracken-2.7

mkdir -p "${WORK_DB}"

echo "[$(date)] Bracken distribution build"
echo "  read length : ${READ_LEN}"
echo "  kmer length : ${KMER_LEN}"
echo "  real DB     : ${REAL_DB}"
echo "  work dir    : ${WORK_DB}"
echo "  host        : $(hostname), ${SLURM_CPUS_PER_TASK} cpus"

if [[ -s "${REAL_DB}/database${READ_LEN}mers.kmer_distrib" ]]; then
    echo "[$(date)] ${REAL_DB}/database${READ_LEN}mers.kmer_distrib already exists — nothing to do."
    exit 0
fi

# Symlink farm: kraken2 reads the index and library through the links, while every
# large output lands on scratch. Avoids copying ~228 GB and avoids writing multi-GB
# intermediates onto the 94%-full project filesystem.
for item in hash.k2d taxo.k2d opts.k2d seqid2taxid.map taxonomy library; do
    if [[ ! -e "${REAL_DB}/${item}" ]]; then
        echo "[$(date)] FAILED — missing ${REAL_DB}/${item}"
        exit 1
    fi
    ln -sfn "${REAL_DB}/${item}" "${WORK_DB}/${item}"
done

# bracken-build calls `python`, which resolves to Python 2.7 inside metawrap-env.
# Shim it to python3 for the duration of this job.
SHIM=${WORK}/bin
mkdir -p "${SHIM}"
printf '#!/bin/bash\nexec /usr/bin/python3 "$@"\n' > "${SHIM}/python"
chmod +x "${SHIM}/python"

module load anaconda3/2023.09
source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate metawrap-env          # same kraken2 (2.0.9-beta) the pipeline classifies with
export PATH=${SHIM}:${PATH}

echo "[$(date)] kraken2: $(command -v kraken2)"
echo "[$(date)] python : $(command -v python) -> $(python --version 2>&1)"

if [[ -s "${WORK_DB}/database.kraken" ]]; then
    echo "[$(date)] database.kraken already present ($(du -h "${WORK_DB}/database.kraken" | cut -f1)) — reusing"
else
    echo "[$(date)] database.kraken not found; bracken-build will classify the library (long step)"
fi

SECONDS=0
"${BRACKEN}/bracken-build" \
    -d "${WORK_DB}" \
    -t "${SLURM_CPUS_PER_TASK}" \
    -k "${KMER_LEN}" \
    -l "${READ_LEN}"
STATUS=$?

conda deactivate

OUT=${WORK_DB}/database${READ_LEN}mers.kmer_distrib
if [[ $STATUS -ne 0 || ! -s "$OUT" ]]; then
    echo "[$(date)] FAILED — bracken-build exited ${STATUS}; ${OUT} not produced"
    echo "  Intermediates kept in ${WORK_DB} for inspection."
    exit 1
fi

echo "[$(date)] Built ${OUT} ($(du -h "$OUT" | cut -f1)) in ${SECONDS}s"

# Only now touch the real DB, and only to add the new distribution.
cp "$OUT" "${REAL_DB}/database${READ_LEN}mers.kmer_distrib"
echo "[$(date)] SUCCESS — installed:"
ls -lh "${REAL_DB}"/database*mers.kmer_distrib

echo
echo "database.kraken is cached at ${WORK_DB}/database.kraken."
echo "Keep it to make future read lengths cheap; delete it to reclaim scratch."
