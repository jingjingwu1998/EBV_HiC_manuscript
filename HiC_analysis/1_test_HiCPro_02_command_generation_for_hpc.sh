#!/bin/bash

# Define the project and raw data paths
export PROJECT_DIR="/share/lab_teng/trainee/JingjingWu/EBV/test_hicpro_mapping"  # Updated PROJECT_DIR for the test run
export RAW_DATA="/share/lab_teng/trainee/JingjingWu/EBV/test_dataset"
export MAP_DATA="/share/lab_teng/trainee/JingjingWu/EBV/test_dataset/test_output"  # Updated output directory to avoid overwriting

# Change to the project directory
cd $PROJECT_DIR

# Step 1: HiC-Pro mapping command (for Bowtie2 alignment and quality checks)
# Ensure output is directed to the new directory specified by $MAP_DATA

/share/lab_teng/trainee/JingjingWu/EBV/HiC_performance/HiC-Pro_install_path/HiC-Pro_3.1.0/bin/HiC-Pro \
    -c 0_config-hicpro.txt \
    -i $RAW_DATA \
    -o $MAP_DATA \  # Use the new output directory to avoid overwriting
    -s mapping \
    -s quality_checks \
    -p
