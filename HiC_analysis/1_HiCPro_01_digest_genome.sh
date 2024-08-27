#!/bin/bash

#PBS -N 01_hicpro_digest 
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=4,mem=32gb
#PBS -m bea
#PBS -M Xiang.Liu@moffitt.org

# This script is used to generate genome fragment based on arima cut site (^GATC ^GANTC) for Hicpro usage

export PROJECT_DIR="/share/lab_teng/data/software/HiC-Pro_3.1.0/bin/utils"

cd $PROJECT_DIR

python digest_genome.py -r ^GATC G^ANTC -o /share/lab_teng/data/software/HiC-Pro_3.1.0/annotation/arima_resfrag_hg19.bed ../../../../genomes/GRCh37/GRCh37.primary_assembly.genome.fa
