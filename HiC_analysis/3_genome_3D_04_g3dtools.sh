#!/bin/bash

# This script is used to convert genome 3D structure bed files to g3d format for washU genome browser.

`for a in `cat chr_list`; do python g3dtools load g3d_input/"$a"_structure.bed -o washU_input/$a; done`
