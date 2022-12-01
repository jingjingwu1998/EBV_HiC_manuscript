#!/bin/bash

# This script is used to calculate 3D genome structures with bedpe files by miniMDS tools

`for a in `cat ../../data/bedpe_list`; do pythonw minimds.py -l 1 -o ../structure_out/"$a"_str.tsv ../../data/bedpe_100k/"$a".bed; done`
