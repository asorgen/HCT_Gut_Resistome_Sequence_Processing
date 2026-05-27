#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. Filters final assembled contigs to >=1000 bp
    # 2. Searches contigs against the ChocoPhlAn pangenome database (waafle_search)
    # 3. Calls ORFs on contigs (waafle_genecaller)
    # 4. Scores contigs for lateral gene transfer events (waafle_orgscorer)
    #
    # Output: per-sample LGT-classified contig tables:
    #   {ID}.lgt.tsv         — putative LGT (two-clade) contigs
    #   {ID}.no_lgt.tsv      — single-clade contigs
    #   {ID}.unclassified.tsv
    #
    # Prerequisite: build the ChocoPhlAn BLAST database once before running any samples:
    #   conda activate waafle-env
    #   zcat /projects/datasets/humann/3.8/chocophlan/*.ffn.gz > /tmp/chocophlan_concat.fna
    #   makeblastdb -in /tmp/chocophlan_concat.fna -dbtype nucl -out ${CHOCOPHLAN_BLAST_DB}
    #   rm /tmp/chocophlan_concat.fna
    #
    # Created on 2026-05-26
    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics
    # Version: 1

# Config files -------------------------------------------------------------------------------------------------------------
    source $pipelineConfig
    source $config_file
    source $bashrc
    source $bash_profile
    source $module_functions
    source $print_functions

# Print script information to log ------------------------------------------------------------------------------------------
    H1 "Description: 5.7_waafle.sh"
        echo -e "Detects lateral gene transfer (LGT) events from assembled contigs using WAAFLE v1.1."
        echo -e "1. Filter contigs >=1000 bp"
        echo -e "2. waafle_search  — BLAST contigs vs ChocoPhlAn"
        echo -e "3. waafle_genecaller — call ORFs on contigs"
        echo -e "4. waafle_orgscorer  — classify contigs as LGT / no-LGT / unclassified"

    H1 "Job Context"
        print "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
        print "Running on host: $(hostname)"
        Total_Gb=$(( SLURM_MEM_PER_NODE / 1024 ))
        JobTime=$(squeue -h -j $SLURM_JOBID -o "%l")
        print "----- Resources Requested -----"
        print "Nodes:           $SLURM_NNODES"
        print "Cores / node:    $SLURM_CPUS_PER_TASK"
        print "Total memory:    ${Total_Gb} GB"
        print "Wall-clock time: $JobTime"
        print "-------------------------------"

    H1 "Variables"
        comment "SampleID: ${ID}"

        H2 "Input"
            inputFile=${evaluationDir}/${ID}_final_assembly.fasta
            echo -e "$inputFile"
            [[ ! -f "$inputFile" ]] && { echo "ERROR: Input file not found: $inputFile"; exit 1; }

        H2 "Output"
            if [[ ! -d ${moduleDir}/${ID} ]]; then mkdir -p ${moduleDir}/${ID}; fi
            filteredContigs=${moduleDir}/${ID}/${ID}_contigs_1kb.fasta
            blastOut=${moduleDir}/${ID}/${ID}
            gffOut=${moduleDir}/${ID}/${ID}.gff
            lgtOut=${moduleDir}/${ID}/${ID}.lgt.tsv
            noLgtOut=${moduleDir}/${ID}/${ID}.no_lgt.tsv
            unclassOut=${moduleDir}/${ID}/${ID}.unclassified.tsv
            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir -p ${moduleDir}/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=("$lgtOut")
    Intermediate_files=("$filteredContigs" "$blastOut" "$gffOut")

# Load environments --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    source $miniforge_init
    conda activate $WAAFLE_ENV
    module load blast/2.15.0+
    module load seqtk

# Run functions ------------------------------------------------------------------------------------------------------------

func="Filter contigs >=1000 bp"
    H1 "$func"
    if [[ -s "$filteredContigs" ]]; then
        comment "Output already found. Skipping..."
    else
        start=$SECONDS
        seqtk seq -L 1000 "$inputFile" > "$filteredContigs"
        [[ $? -ne 0 ]] && error "seqtk failed for $ID"
        N_CONTIGS=$(grep -c "^>" "$filteredContigs" || echo 0)
        comment "Contigs >=1000 bp: $N_CONTIGS"
        end=$SECONDS; comment "Elapsed: $(elapsed_time $((end - start)))"
    fi
    step_completion "$filteredContigs"


func="WAAFLE search (blastn vs ChocoPhlAn)"
    H1 "$func"
    if [[ -s "$blastOut" ]]; then
        comment "Output already found. Skipping..."
    else
        start=$SECONDS
        waafle_search \
            "$filteredContigs" \
            "$CHOCOPHLAN_BLAST_DB" \
            --out "$blastOut" \
            --threads $SLURM_CPUS_PER_TASK
        [[ $? -ne 0 ]] && error "waafle_search failed for $ID"
        HITS=$(wc -l < "$blastOut")
        comment "BLAST hits: $HITS"
        end=$SECONDS; comment "Elapsed: $(elapsed_time $((end - start)))"
    fi
    step_completion "$blastOut"


func="WAAFLE gene caller (ORF prediction)"
    H1 "$func"
    if [[ -s "$gffOut" ]]; then
        comment "Output already found. Skipping..."
    else
        start=$SECONDS
        waafle_genecaller "$filteredContigs" --out "$gffOut"
        [[ $? -ne 0 ]] && error "waafle_genecaller failed for $ID"
        N_GENES=$(wc -l < "$gffOut")
        comment "Predicted ORFs: $N_GENES"
        end=$SECONDS; comment "Elapsed: $(elapsed_time $((end - start)))"
    fi
    step_completion "$gffOut"


func="WAAFLE org scorer (LGT classification)"
    H1 "$func"
    if [[ -s "$lgtOut" ]]; then
        comment "Output already found. Skipping..."
    else
        start=$SECONDS
        waafle_orgscorer \
            "$filteredContigs" \
            "$blastOut" \
            "$gffOut" \
            --basename "${moduleDir}/${ID}/${ID}"
        [[ $? -ne 0 ]] && error "waafle_orgscorer failed for $ID"

        N_LGT=$(awk 'NR>1' "$lgtOut" | wc -l)
        N_NO_LGT=$(awk 'NR>1' "$noLgtOut" | wc -l)
        N_UNCLASS=$(awk 'NR>1' "$unclassOut" | wc -l)
        comment "LGT contigs:          $N_LGT"
        comment "Single-clade contigs: $N_NO_LGT"
        comment "Unclassified:         $N_UNCLASS"
        end=$SECONDS; comment "Elapsed: $(elapsed_time $((end - start)))"
    fi
    step_completion "$lgtOut"


    conda deactivate
    module unload seqtk

# Completion ---------------------------------------------------------------------------------------------------------------
    output_exists=$(test_for_output "${Complete_tag[@]}")
    if $output_exists; then
        touch ${moduleDir}/COMPLETE/$ID
        for int_file in "${Intermediate_files[@]}"; do
            echo -e "rm $int_file"; rm -f "$int_file"
        done
    else
        error "Not all outputs were created."
    fi

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
