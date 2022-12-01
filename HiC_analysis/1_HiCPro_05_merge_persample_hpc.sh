#!/bin/bash
#PBS -l nodes=1:ppn=1,mem=72gb,walltime=120:00:00
#PBS -M xiang.liu@moffitt.org
#PBS -m ae
#PBS -j eo
#PBS -N merge_persample
#PBS -q batch
#PBS -V

cd $PBS_O_WORKDIR
make --file /share/lab_teng/data/software//HiC-Pro_3.1.0/scripts/Makefile CONFIG_FILE=/share/lab_teng/xiangliu/hic_analysis/hicpro_mapping_2/config-hicpro.txt CONFIG_SYS=/share/lab_teng/data/software//HiC-Pro_3.1.0/config-system.txt merge_persample 2>&1
