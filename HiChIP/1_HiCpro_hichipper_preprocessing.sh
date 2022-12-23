#!/bin/bash

# 1. Run HiC PRo
HiC-Pro -i raw_fastq/ -o /HiC-Pro_out -c config-hicpro.txt

# 2. Run HiChipper
hichipper --out EBNA3A_H3K27ac_hichipper_merged  --make-ucsc  yaml/merged.yaml 