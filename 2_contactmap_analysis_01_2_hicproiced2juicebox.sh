#!/bin/bash

# This script is used to convert *.pairs files to *.hic files for Juicebox visulization

cd /share/lab_teng/xiangliu/hic_analysis/hicpro_mapping_2/contact_maps/hicpro2juice_out
java -Xmx32g -jar /share/lab_teng/xiangliu/hic_analysis/juicer/scripts/juicer_tools.jar pre -r 5000,10000,25000,50000,100000,500000,1000000 -d -t Day0_tmp -n -v \
    Day0_5000_iced.pairs Day0_5000_iced.hic hg19

java -Xmx32g -jar /share/lab_teng/xiangliu/hic_analysis/juicer/scripts/juicer_tools.jar pre -r 5000,10000,25000,50000,100000,500000,1000000 -d -t Day28_tmp -n -v \
    Day28_5000_iced.pairs Day28_5000_iced.hic hg19
