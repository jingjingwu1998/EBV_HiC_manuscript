#!/bin/bash
# by Xiang Liu

# This script is used to convert *.pairs files to *.hic files for Juicebox visulization

# cd /share/lab_teng/xiangliu/hic_analysis/hicpro_mapping_2/contact_maps/hicpro2juice_out
# java -Xmx32g -jar /share/lab_teng/xiangliu/hic_analysis/juicer/scripts/juicer_tools.jar pre -r 5000,10000,25000,50000,100000,500000,1000000 -d -t Day0_tmp -n -v \
#    Day0_5000_iced.pairs Day0_5000_iced.hic hg19

# java -Xmx32g -jar /share/lab_teng/xiangliu/hic_analysis/juicer/scripts/juicer_tools.jar pre -r 5000,10000,25000,50000,100000,500000,1000000 -d -t Day28_tmp -n -v \
#    Day28_5000_iced.pairs Day28_5000_iced.hic hg19



# by Jingjing Wu

#PBS -l nodes=1:ppn=12,mem=72gb,walltime=12:00:00
#PBS -M Jingjing.Wu@moffitt.org
#PBS -N juicer_tool_GCBC_5000_iced.hic

java -Xmx32g -jar /share/lab_teng/trainee/JingjingWu/EBV/juicer/scripts/juicer_tools_linux_0.8_1.5.3.jar pre -r 5000,10000,25000,50000,100000,500000,1000000 -d -t tmp_GCBC -n -v \/share/lab_teng/trainee/JingjingWu/EBV/hicpro2juice_out/GCBC_5000_iced.pairs /share/lab_teng/trainee/JingjingWu/EBV/hicpro2juice_out/GCBC_5000_iced.hic /share/lab_teng/trainee/JingjingWu/EBV/HiC_performance/chrom_hg19.sizes
