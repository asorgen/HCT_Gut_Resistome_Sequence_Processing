#!/bin/bash

ont_names=(D21309d22_dorado_0.4.1_dna_r10.4.1_e8.2_400bps_hac_v4.2.0.fastq.gz 21309PRE_PASS.fastq.gz D20248D21_08_11_23_Promethion.pass.fastq.gz 21309d10_PASS.fastq.gz D20248D1_10_31_23_dorado_0.4.1_dna_r10.4.1_e8.2_400bps_sup_v4.2.0.fastq.gz 20248PRE_PASS.fastq.gz 01_08_24_D21309D15_Promethion_dorado_0.5.1_dna_r10.4.1_e8.2_400bps_hac_v4.3.0.fastq.gz 12_15_23_D20248D57_Promethion.fastq.gz 12_15_13_D20248D14_Promethion.fastq.gz D20248d8_11_3_23_dorado_0.4.1_dna_r10.4.1_e8.2_400bps_hac_v4.2.0.fastq.gz D20248D28_10_31_23_dorado_0.4.1_dna_r10.4.1_e8.2_400bps_hac_v4.2.0.fastq.gz 12_01_23_D21309D29_Promethion.fastq.gz 01_08_24_D21309D0_Promethion_dorado_0.5.1_dna_r10.4.1_e8.2_400bps_hac_v4.3.0.fastq.gz 12_01_23_D21309D58_Promethion.fastq.gz)

new_names=(D21309D22 D21309PRE D20248D21 D21309D10 D20248D1 D20248PRE D21309D15 D20248D54 D20248D14 D20248D8 D20248D28 D21309D29 D21309D0 D21309D58)

for i in "${!new_names[@]}"; do
  if [[ "${new_names[$i]}" == "$ID" ]]; then
    if [[ ! -f "${moduleDir}/${new_names[$i]}.fastq.gz" ]]; then
	    # echo -e "Command: cp ${seqPath}/${ont_names[$i]} ${moduleDir}/${new_names[$i]}.fastq.gz"
	    cp ${seqPath}/${ont_names[$i]} ${moduleDir}/${new_names[$i]}.fastq.gz
    fi
  fi
done