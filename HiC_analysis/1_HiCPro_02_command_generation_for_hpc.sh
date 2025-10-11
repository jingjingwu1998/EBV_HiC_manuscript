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

# 1.0 for this command you can submit it locally using ./
# 1.1 after submission, you will get /share/lab_teng/trainee/JingjingWu/EBV/hicpro_mapping_2/hicpro_mapping_2
# you will run HiCPro_step1_hic_bo.sh after step 1.2
# 1.2 remember to edit input files in hicpro_mapping_2 file as following format, and ends up with _1 and _2 is important
# Day28/Day28_CKDL220009104-1a_HN2CMDSX3_L3_1.fq.gz
# Day28/Day28_CKDL220009104-1a_HN2CMDSX3_L3_2.fq.gz
# run HiCPro_step2_hic_bo.sh and check log files

/share/lab_teng/trainee/JingjingWu/EBV/HiC_performance/HiC-Pro_install_path/HiC-Pro_3.1.0/bin/HiC-Pro -c 0_config-hicpro.txt -i $RAW_DATA -o hicpro_mapping_2 -s mapping -s quality_checks -p

#/share/lab_teng/trainee/JingjingWu/EBV/HiC_performance/HiC-Pro_install_path/HiC-Pro_3.1.0/bin/HiC-Pro -c 0_config-hicpro.txt -i /share/lab_teng/trainee/JingjingWu/EBV/raw_data -o hicpro_mapping_2 -s mapping -s quality_checks -p

# Step 2: HiC-Pro proc_hic command
# update the map data directory by export MAP_DATA="/share/lab_teng/trainee/JingjingWu/EBV/hicpro_mapping_2/hicpro_mapping_2, before running this command
/share/lab_teng/trainee/JingjingWu/EBV/HiC_performance/HiC-Pro_install_path/HiC-Pro_3.1.0/bin/HiC-Pro -c 0_config-hicpro.txt -i ${MAP_DATA}/bowtie_results/bwt2 -o proc_hic_out -s proc_hic -s quality_checks -p
# after running this command, please go to file proc_hic_out for next steps
# to run HiCPro_step1_hic_bo.sh, please adjust the file format in inputfiles_hic_bo.txt into "rawdata/cMCL1/cMCL1_read_1_GRCh37.bwt2merged"

# Step 3: HiC-Pro merge_persample command 
# change this directory before running export MAP_DATA="/share/lab_teng/trainee/JingjingWu/EBV/hicpro_mapping_2"
/share/lab_teng/trainee/JingjingWu/EBV/HiC_performance/HiC-Pro_install_path/HiC-Pro_3.1.0/bin/HiC-Pro -i ${MAP_DATA}/proc_hic_out/hic_results/data -o ${MAP_DATA}/merge_validpairs -c ${MAP_DATA}/0_config-hicpro.txt -s merge_persample -p

# Step 4: HiC-Pro build_contact_maps and ice_norm commands 
# before run the command please upgrade memory to 128G
/share/lab_teng/trainee/JingjingWu/EBV/HiC_performance/HiC-Pro_install_path/HiC-Pro_3.1.0/bin/HiC-Pro -i ${MAP_DATA}/merge_validpairs/hic_results/data -o ${MAP_DATA}/contact_maps -c ${MAP_DATA}/0_config-hicpro.txt -s build_contact_maps -s ice_norm -p
