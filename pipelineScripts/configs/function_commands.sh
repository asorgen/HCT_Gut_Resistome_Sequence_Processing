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
        echo "Short-read cohorts run through the generic wrapper:"
        echo "  pipelineScripts/pipeline_wrappers/run_short_pipeline.sh <cohort>_<read>"
        echo "Long/hybrid still use their own wrappers: <cohort>_<read>_pipeline.sh"
        echo "Log written to the config's own datasetDir/LOGs/"
        return 0
    fi
    cohort=$1
    read=$2
    dataset=${cohort}_${read}
    export version=$(date +"%Y.%m.%d")
    root=${HPC_PROJECTS}/HCT_Gut_Resistome_Study
    pipeline_root=${root}/HCT_Gut_Resistome_Sequence_Processing
    cd $pipeline_root

    local config_file="${pipeline_root}/pipelineScripts/configs/${dataset}-read.config"
    if [[ ! -f "$config_file" ]]; then
        echo "No such config: ${config_file}"
        cd - >/dev/null
        return 1
    fi

    # Take the log location from the config rather than reconstructing it — DATA_ROOT
    # is no longer always under HCT_Gut_Resistome_Data (Heston lives beside its FASTQs
    # on /projects/afodor_research). Read it in a subshell so sourcing the config
    # doesn't leak its variables into the caller's interactive shell.
    local datasetDir
    datasetDir=$(source "$config_file" >/dev/null 2>&1; echo "$datasetDir")
    if [[ -z "$datasetDir" ]]; then
        echo "Config ${config_file} did not define datasetDir"
        cd - >/dev/null
        return 1
    fi

    # Build the command as an array — the short-read form carries the config name as
    # a second word, which a quoted "$script" would collapse into one filename.
    local -a cmd
    if [[ "$read" == "short" ]]; then
        cmd=("${pipeline_root}/pipelineScripts/pipeline_wrappers/run_short_pipeline.sh" "$dataset")
    else
        cmd=("${pipeline_root}/pipelineScripts/pipeline_wrappers/${dataset}_pipeline.sh")
    fi
    [ -n "$3" ] && cmd+=("$3")

    mkdir -p "${datasetDir}/LOGs"
    local log_file="${datasetDir}/LOGs/${dataset}_pipeline_$version.out"

    echo "nohup sh ${cmd[*]} >> $log_file 2>&1 &"
    nohup sh "${cmd[@]}" >> "$log_file" 2>&1 &
    
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