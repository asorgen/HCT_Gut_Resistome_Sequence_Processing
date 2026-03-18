# Resolve path to print_header.py relative to this file
    PRINT_HEADER="$(cd "$(dirname "${BASH_SOURCE[0]}")/../scripts/helper_scripts" && pwd)/print_header.py"

# Set function for output comments
    H1 () { "$PRINT_HEADER" "$1" "H1"; }
    H2 () { "$PRINT_HEADER" "$1" "H2"; }
    H3 () { "$PRINT_HEADER" "$1" "H3"; }
    comment () { "$PRINT_HEADER" "$1" "#"; echo; }
    print () { "$PRINT_HEADER" "$1" "#"; }
    error () { echo $1; exit 1; }
    pFunc () { echo $1; echo; }

# Find out if the job is already queued
    job_lookup() { squeue -u ${HPC_USER} --format='%.18i   %.9P   %.40j   %.1T   %.12M   %.10l   %.6D   %R' | awk -v id="$jobID" 'match($3,id) {print $1}'; }

# Set up module directory, variables, and flags
    module_setup() {
        script=$1
        
        moduleDir="${script%.sh}"

        if [ ! -z "$2" ]; then moduleDir=$moduleDir/$2; fi

        if [[ ! -d ${moduleDir} ]]; then mkdir -p $moduleDir; fi
        export moduleDir
        
        logDir=${moduleDir}/logs
        if [[ ! -d ${logDir} ]]; then mkdir -p $logDir; fi

        Complete_tag=(${moduleDir}/COMPLETE/${ID})

        jobID="${script%%_*}"
        jobID=${jobID}_${ID}_${dataset}
        if [ ! -z "$2" ]; then jobID=${jobID}_$2; fi

    }

# Find out if this is the first module submitted from this SampleID
    first_ID() {
        if [[ $module == 0 ]]; then 
            H2 "${ID}"
            ((module++))
            ((count++))
        fi    
    }

# Run a pipeline module step with all standard boilerplate in one call.
# Usage: run_module_step <script> <header> <dep_job_id> <hpc_opts> <tag> <result_var>
#   script      — module script filename (e.g., "0.3_sequence_trim.sh")
#   header      — display name printed by H3 (e.g., "Sequence trim")
#   dep_job_id  — job ID of the preceding step, or "COMPLETE" if no dependency
#                 (always pass the stripped ID using ${PREV_JOB##* })
#   hpc_opts    — slurm resource string (e.g., $trim_opts from config)
#   tag         — short pipeline tag for single-step reruns (e.g., "trim")
#   result_var  — name of the variable to store Current_Job in (e.g., "TRIM_JOB")
#
# Returns 1 if $pipeline_stop_tag matches this step's tag (signals caller to `continue`).
# Returns 0 otherwise.
    run_module_step() {
        local script="$1"
        local header="$2"
        local dep_id="$3"
        local opts="$4"
        local tag="$5"
        local result_var="$6"

        module_setup "$script"
        header3="$header"
        DEPENDENT_JOB=("$dep_id")
        hpc_opts="$opts"
        pipeline_tag="$tag"

        run_module

        printf -v "$result_var" '%s' "$Current_Job"
        CLEAN_UP_DEP+=("${Current_Job##* }")

        [[ "$pipeline_stop_tag" == "$tag" ]] && return 1
        return 0
    }

# Run the module script
    run_module() {
        
        # Set flag for if the completion file exists
        all_exist=true
        for file in "${Complete_tag[@]}"; do
          if [[ ! -e "$file" ]]; then
            # echo $file
            all_exist=false
            break
          fi
        done

        if $all_exist; then
            Current_Job=COMPLETE
        else
            queued=$(job_lookup $jobID)
            # If the job is not queued
            if [ -z "$queued" ]; then
                ready_to_run=true
                dependencies=""
                for job in "${DEPENDENT_JOB[@]}"; do
                    if [[ $job != "COMPLETE" ]]; then # if dependent job isn't COMPLETE
                        ready_to_run=false # set flag to false
                        if [[ -z "$dependencies" ]]; then # if dependencies is empty
                            dependencies="${job}" # set dependencies to $job
                        else
                            dependencies+=",${job}" # add this $job to dependencies
                        fi
                    fi
                done
                # echo $dependencies

                # If the dependent job is complete
                if $ready_to_run; then
                    first_ID
                    H3 "$header3"
                    Current_Job=$(sbatch \
                        --chdir=${datasetDir} \
                        $hpc_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ${moduleScriptPath}/${script})
                    echo $Current_Job
                    echo "[ Log file ] -> ${dataset}/${logDir}/${ID}.${Current_Job##* }.log"
                else
                    H3 "$header3"
                    # echo -e "sbatch --dependency=afterok:${dependencies} \
                    #     $hpc_opts \
                    #     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ${moduleScriptPath}/${script}"
                    Current_Job=$(sbatch --dependency=afterok:${dependencies} \
                        --chdir=${datasetDir} \
                        $hpc_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ${moduleScriptPath}/${script})
                    echo $Current_Job
                    echo "[ Log file ] -> ${dataset}/${logDir}/${ID}.${Current_Job##* }.log"
               fi 
            else
                Current_Job=$queued
           fi # if [ -z "$queued" ]
        fi     
    }
