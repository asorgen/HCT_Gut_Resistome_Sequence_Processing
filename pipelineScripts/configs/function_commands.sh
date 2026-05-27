# Source private config for HPC paths (if not already sourced)
if [[ -z "${HPC_PROJECTS}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/private.config"
fi

# Function to easily run pipelines
run_pipeline() {
    if [[ "$1" == "-h" || "$1" == "--help" ]]; then
        echo "Usage: run_pipeline <cohort> <read> [stop_at]"
        echo ""
        echo "Launches a pipeline wrapper in the background via nohup and logs output."
        echo ""
        echo "Arguments:"
        echo "  cohort    Cohort name matching the wrapper filename (e.g. Duke, UNC)"
        echo "  read      Read type matching the wrapper filename (e.g. short, long, hybrid)"
        echo "  stop_at   Optional pipeline tag passed to the wrapper to halt after that module"
        echo "            (e.g. trim, kraken, bins)"
        echo ""
        echo "Wrapper called: pipelineScripts/pipeline_wrappers/<cohort>_<read>_pipeline.sh"
        echo "Log written to: HCT_Gut_Resistome_Data/unprocessed/<cohort>/<cohort>_<read>/LOGs/"
        return 0
    fi
    cohort=$1
    read=$2
    dataset=${cohort}_${read}
    export version=$(date +"%Y.%m.%d")
    root=${HPC_PROJECTS}/HCT_Gut_Resistome_Study
    pipeline_root=${root}/HCT_Gut_Resistome_Sequence_Processing
    data_root=${root}/HCT_Gut_Resistome_Data
    cd $pipeline_root
    
    local log_file="${data_root}/unprocessed/${cohort}/${dataset}/LOGs/${dataset}_pipeline_$version.out"
    local script="${pipeline_root}/pipelineScripts/pipeline_wrappers/${dataset}_pipeline.sh"
    
    if [ -n "$3" ]; then
        echo "nohup sh $script $3 >> $log_file 2>&1 &"
        nohup sh "$script" "$3" >> "$log_file" 2>&1 &
    else
        echo "nohup sh $script >> $log_file 2>&1 &"
        nohup sh "$script" >> "$log_file" 2>&1 &
    fi
    
    cd -
}

# Check the number of completed jobs
check_jobs() {
    if [[ "$1" == "-h" || "$1" == "--help" ]]; then
        echo "Usage: check_jobs <cohort> <read> [module]"
        echo ""
        echo "Prints the number of samples that have completed processing."
        echo ""
        echo "Arguments:"
        echo "  cohort    Cohort name (e.g. Duke, UNC)"
        echo "  read      Read type (e.g. short, long, hybrid)"
        echo "  module    Optional module prefix to check (e.g. 5.6_AA_amr_bins)."
        echo "            Counts COMPLETE flags inside matching module directories."
        echo "            Omit to count top-level COMPLETE flags (full pipeline completions)."
        return 0
    fi
    local cohort=$1
    local read=$2
    local dataset=${cohort}_${read}
    local root=${HPC_PROJECTS}/HCT_Gut_Resistome_Study
    local pipeline_root=${root}/HCT_Gut_Resistome_Sequence_Processing
    local data_root=${root}/HCT_Gut_Resistome_Data

    if [ -n "$3" ]; then
        local module=$3
        local completed_dir=${data_root}/unprocessed/${cohort}/${dataset}/${module}_*/COMPLETE
        local completed_samples=$(ls ${completed_dir} | wc -l)
    else
        local completed_dir=${data_root}/unprocessed/${cohort}/${dataset}/COMPLETE
        local completed_samples=$(ls ${completed_dir} | wc -l)
        local completed_samples=$(( completed_samples - 1 ))
    fi
    
    echo $completed_samples
}

# Run post-pipeline summarization
run_summary() {
    if [[ "$1" == "-h" || "$1" == "--help" ]]; then
        echo "Usage: run_summary <cohort> <read> [modules]"
        echo ""
        echo "Runs pipeline_summary.sh to aggregate per-sample results into count tables."
        echo ""
        echo "Arguments:"
        echo "  cohort    Cohort name (e.g. Duke, UNC)"
        echo "  read      Read type (e.g. short, long, hybrid)"
        echo "  modules   Optional comma-separated list of modules to summarize"
        echo "            (e.g. rgi_bwt,asm_gene_profiling)."
        echo "            Omit to summarize all modules."
        echo ""
        echo "Output tables written to: <PROCESSED_ROOT>/<cohort>_<read>_tables/"
        return 0
    fi
    local cohort=$1
    local read=$2
    local dataset=${cohort}_${read}
    local root=${HPC_PROJECTS}/HCT_Gut_Resistome_Study
    local pipeline_root=${root}/HCT_Gut_Resistome_Sequence_Processing
    local data_root=${root}/HCT_Gut_Resistome_Data
    cd $pipeline_root

    local script="${pipeline_root}/pipelineScripts/scripts/post_scripts/pipeline_summary.sh"

    if [ -n "$3" ]; then
        $script -p $dataset -m $3
    else
        $script -p $dataset
    fi
    
    cd -
}