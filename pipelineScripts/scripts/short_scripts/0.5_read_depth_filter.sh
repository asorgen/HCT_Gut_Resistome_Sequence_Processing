#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. Counts the read pairs remaining after deduplication, trimming, and host decontamination (module 0.4).
    # 2. Records the count for every sample, whether it passes or not.
    # 3. Marks samples below MIN_CLEAN_READS as EXCLUDED so downstream modules skip them.

    # Why this module exists:
    #   Samples that lose most of their reads in QC have too little data for a
    #   meaningful resistome call. 5.4_rgi_bwt treats "no ARG hits" as a fatal
    #   error, so a depleted sample fails, never writes its COMPLETE flag, and is
    #   resubmitted on every pass of the throttled driver - which also prevents
    #   the driver from ever reaching its exit condition. Filtering on depth up
    #   front removes the cause rather than the symptom, and keeps the cohort to
    #   samples with enough reads to plausibly detect anything.

    # This module never fails a sample by exiting non-zero. Excluded samples are a
    # legitimate outcome, not an error, so it always writes its COMPLETE flag and
    # signals exclusion via a marker file. Failing here would strand every
    # downstream job in DependencyNeverSatisfied.

    # Created on Aug 27, 2026

    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

# Config files -------------------------------------------------------------------------------------------------------------
    source $pipelineConfig
    source $config_file
    source $bashrc
    source $bash_profile
    source $module_functions
    source $print_functions

# Print script information to log ------------------------------------------------------------------------------------------
    H1 "Usage"
        echo -e "This script does the following:"
        echo -e "1. Counts read pairs surviving QC (module 0.4)."
        echo -e "2. Records the count for every sample."
        echo -e "3. Marks samples below MIN_CLEAN_READS as EXCLUDED so downstream modules skip them."

    H1 "Job Context"
        comment "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
        comment "Running on host: `hostname`"

    H1 "Variables"
        comment "SampleID: ${ID}"

        : "${MIN_CLEAN_READS:=500000}"

        H2 "Input"
            R1=${clean_readDir}/${ID}_1.fastq.gz; echo -e "${R1}"
            R2=${clean_readDir}/${ID}_2.fastq.gz; echo -e "${R2}"

        H2 "Output"
            out=${moduleDir}
            mkdir -p ${out}/COMPLETE ${out}/EXCLUDED ${out}/counts
            echo -e "${out}/counts/${ID}.count"

        comment "Minimum clean read pairs required: ${MIN_CLEAN_READS}"

    H2 "[ Start ]"
    /bin/date
    SECONDS=0

# Run functions ------------------------------------------------------------------------------------------------------------

STEP="Read depth filter"
    H1 "$STEP"

    if [[ ! -s "$R1" || ! -s "$R2" ]]; then
        error "Missing decontaminated reads for ${ID} - run module 0.4 first"
    fi

    # Count from R1 only; reads are paired so R1 count == read pairs.
    n_pairs=$(zcat "$R1" | awk 'END {print int(NR/4)}')

    if [[ -z "$n_pairs" ]]; then
        error "Could not count reads in ${R1}"
    fi

    printf '%s\t%s\n' "${ID}" "${n_pairs}" > ${out}/counts/${ID}.count
    comment "${ID}: ${n_pairs} read pairs"

    if [[ ${n_pairs} -lt ${MIN_CLEAN_READS} ]]; then
        touch ${out}/EXCLUDED/${ID}
        comment "EXCLUDED - ${n_pairs} < ${MIN_CLEAN_READS} read pairs"
        comment "Downstream modules will skip ${ID}. This is a recorded outcome, not a failure."
    else
        # Clear any stale marker so raising the threshold and re-running is idempotent.
        rm -f ${out}/EXCLUDED/${ID}
        comment "PASS - ${n_pairs} >= ${MIN_CLEAN_READS} read pairs"
    fi

H1 "Completion"

    # Always complete: exclusion is a result, not an error.
    touch ${out}/COMPLETE/${ID}

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
/bin/date
