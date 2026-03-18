#!/bin/bash

declare -A name_map
name_map=(
    [D21309D22]="D21309d22_dorado_0.4.1_dna_r10.4.1_e8.2_400bps_hac_v4.2.0.fastq.gz"
    [D21309PRE]="21309PRE_PASS.fastq.gz"
    [D20248D21]="D20248D21_08_11_23_Promethion.pass.fastq.gz"
    [D21309D10]="21309d10_PASS.fastq.gz"
    [D20248D1]="D20248D1_10_31_23_dorado_0.4.1_dna_r10.4.1_e8.2_400bps_sup_v4.2.0.fastq.gz"
    [D20248PRE]="20248PRE_PASS.fastq.gz"
    [D21309D15]="01_08_24_D21309D15_Promethion_dorado_0.5.1_dna_r10.4.1_e8.2_400bps_hac_v4.3.0.fastq.gz"
    [D20248D54]="12_15_23_D20248D57_Promethion.fastq.gz"
    [D20248D14]="12_15_13_D20248D14_Promethion.fastq.gz"
    [D20248D8]="D20248d8_11_3_23_dorado_0.4.1_dna_r10.4.1_e8.2_400bps_hac_v4.2.0.fastq.gz"
    [D20248D28]="D20248D28_10_31_23_dorado_0.4.1_dna_r10.4.1_e8.2_400bps_hac_v4.2.0.fastq.gz"
    [D21309D29]="12_01_23_D21309D29_Promethion.fastq.gz"
    [D21309D0]="01_08_24_D21309D0_Promethion_dorado_0.5.1_dna_r10.4.1_e8.2_400bps_hac_v4.3.0.fastq.gz"
    [D21309D58]="12_01_23_D21309D58_Promethion.fastq.gz"
)

ont_file="${name_map[$ID]}"

if [[ -n "$ont_file" ]]; then
    if [[ ! -f "${moduleDir}/${ID}.fastq.gz" ]]; then
        # echo -e "Command: cp ${seqPath}/${ont_file} ${moduleDir}/${ID}.fastq.gz"
        cp "${seqPath}/${ont_file}" "${moduleDir}/${ID}.fastq.gz"
    fi
fi