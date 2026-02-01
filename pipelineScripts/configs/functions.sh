# Set function for output comments
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; }
    error () { echo $1; exit 1; }

# Find out if the job is already queued
    job_lookup() { squeue -u asorgen --format='%.18i   %.9P   %.40j   %.1T   %.12M   %.10l   %.6D   %R' | awk -v id="$jobID" 'match($3,id) {print $1}'; }

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
                echo $dependencies

                # If the dependent job is complete
                if $ready_to_run; then
                    first_ID
                    H3 "$header3"
                    Current_Job=$(sbatch \
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
