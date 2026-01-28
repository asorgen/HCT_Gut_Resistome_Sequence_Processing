# Logical test to determine if files exist
    test_for_output() {
        local File_list=("$@")  # grab all args as one big array
        all_complete=true
        for file in "${File_list[@]}"; do 
            if [[ ! -e "$file" ]]; then 
                all_complete=false
                break
            fi
        done
        echo "$all_complete"    
    }


# Test if the function ran properly and generated output
    step_completion() {
        File_list=("$@")  # capture all arguments as an array
        
        all_complete=true
        for file in "${File_list[@]}"; do 
            if [[ ! -e "$file" ]]; then 
                all_complete=false
                break
            fi
        done

        if $all_complete; then 
            comment "SUCCESS: $STEP Complete"
        else
            error "[ $STEP ERROR! ] - Exiting..."
        fi
        Complete_tag+=("${File_list[@]}")
    }



# Test if the function ran properly and generated intermediate files
    substep_completion() {
        File_list=("$@")  # capture all arguments as an array
        
        all_complete=true
        for file in "${File_list[@]}"; do 
            if [[ ! -e "$file" ]]; then 
                all_complete=false
                break
            fi
        done

        if $all_complete; then 
            comment "SUCCESS: $STEP Complete"
        else
            error "[ $STEP ERROR! ] - Exiting..."
        fi
        Intermediate_files+=("${File_list[@]}")
    }


# Generate completion tag and remove intermediate files
    module_completion() {
        output_exists=$(test_for_output "${Complete_tag[@]}")
        if $output_exists; then 
            touch ${moduleDir}/COMPLETE/$ID

            # Remove intermediate files
            for int_file in "${Intermediate_files[@]}"; do
                rm $int_file
            done
        fi
    }


