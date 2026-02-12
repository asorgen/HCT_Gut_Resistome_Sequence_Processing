# Function to easily run pipelines
run_pipeline() {
    cohort=$1
    read=$2
    dataset=${cohort}_${read}
    export version=$(date +"%Y.%m.%d")
    root=/projects/afodor_research3/asorgen/HCT_Gut_Resistome_Study
    pipeline_root=${root}/HCT_Gut_Resistome_Pipeline/sequence_processing
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
