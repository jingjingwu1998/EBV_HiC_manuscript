#!/bin/bash
#PBS -l nodes=1:ppn=1,mem=72gb,walltime=120:00:00
#PBS -M xiang.liu@moffitt.org
#PBS -m ae
#PBS -j eo
#PBS -N mapping_quality_checks
#PBS -q batch
#PBS -V

module load gcc/5.4.0

cd $PBS_O_WORKDIR
make --file /share/lab_teng/data/software//HiC-Pro_3.1.0/scripts/Makefile CONFIG_FILE=/share/lab_teng/xiangliu/hic_analysis/hicpro_mapping_2/config-hicpro.txt CONFIG_SYS=/share/lab_teng/data/software//HiC-Pro_3.1.0/config-system.txt quality_checks 2>&1
