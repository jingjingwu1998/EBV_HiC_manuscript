#!/bin/bash
#PBS -l nodes=1:ppn=12,mem=72gb,walltime=120:00:00
#PBS -M xiang.liu@moffitt.org
#PBS -m ae
#PBS -j eo
#PBS -N proc_hic
#PBS -q batch
#PBS -V
#PBS -t 1-2
cd $PBS_O_WORKDIR
FASTQFILE=$PBS_O_WORKDIR/inputfiles_hic_bo.txt; export FASTQFILE
make --file /share/lab_teng/data/software//HiC-Pro_3.1.0/scripts/Makefile CONFIG_FILE=/share/lab_teng/xiangliu/hic_analysis/hicpro_mapping_2/config-hicpro.txt CONFIG_SYS=/share/lab_teng/data/software//HiC-Pro_3.1.0/config-system.txt proc_hic  2>&1
