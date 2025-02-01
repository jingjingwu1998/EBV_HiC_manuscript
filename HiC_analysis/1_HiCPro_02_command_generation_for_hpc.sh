#!/bin/bash

# This is the script to create HiC-Pro commands for HPC systems

# initial working directories
# export PROJECT_DIR="/share/lab_teng/xiangliu/hic_analysis/hicpro_mapping_2"
# export RAW_DATA="/share/lab_teng/data/sequencing/hic_bo/raw_data"
# export MAP_DATA="/share/lab_teng/xiangliu/hic_analysis/hicpro_mapping_2"
export PROJECT_DIR="/share/lab_teng/trainee/JingjingWu/EBV/hicpro_mapping_2"
export RAW_DATA="/share/lab_teng/trainee/JingjingWu/EBV/raw_data"
# export MAP_DATA="/share/lab_teng/trainee/JingjingWu/EBV/hicpro_mapping_2"
export MAP_DATA="/share/lab_teng/trainee/JingjingWu/EBV/hicpro_mapping_2/hicpro_mapping_2"


# module load gcc/5.4.0
module load gcc/11.2.0

cd $PROJECT_DIR

# Step 1: HiC-Pro mapping command
# /share/lab_teng/data/software/HiC-Pro_3.1.0/bin/HiC-Pro -c 0_config-hicpro.txt -i $RAW_DATA -o hicpro_mapping_2 -s mapping -s quality_checks -p

/share/lab_teng/trainee/JingjingWu/EBV/HiC_performance/HiC-Pro_install_path/HiC-Pro_3.1.0/bin/HiC-Pro -c 0_config-hicpro.txt -i $RAW_DATA -o hicpro_mapping_2 -s mapping -s quality_checks -p

#/share/lab_teng/trainee/JingjingWu/EBV/HiC_performance/HiC-Pro_install_path/HiC-Pro_3.1.0/bin/HiC-Pro -c 0_config-hicpro.txt -i /share/lab_teng/trainee/JingjingWu/EBV/raw_data -o hicpro_mapping_2 -s mapping -s quality_checks -p

# Step 2: HiC-Pro proc_hic command
/share/lab_teng/trainee/JingjingWu/EBV/HiC_performance/HiC-Pro_install_path/HiC-Pro_3.1.0/bin/HiC-Pro -c 0_config-hicpro.txt -i ${MAP_DATA}/bowtie_results/bwt2 -o proc_hic_out -s proc_hic -s quality_checks -p

# Step 3: HiC-Pro merge_persample command 
/share/lab_teng/trainee/JingjingWu/EBV/HiC_performance/HiC-Pro_install_path/HiC-Pro_3.1.0/bin/HiC-Pro -i ${MAP_DATA}/proc_hic_out/hic_results/data -o ${MAP_DATA}/merge_validpairs -c ${MAP_DATA}/0_config-hicpro.txt -s merge_persample -p

# Step 4: HiC-Pro build_contact_maps and ice_norm commands 
/share/lab_teng/trainee/JingjingWu/EBV/HiC_performance/HiC-Pro_install_path/HiC-Pro_3.1.0/bin/HiC-Pro -i ${MAP_DATA}/merge_validpairs/hic_results/data -o ${MAP_DATA}/contact_maps -c ${MAP_DATA}/0_config-hicpro.txt -s build_contact_maps -s ice_norm -p
