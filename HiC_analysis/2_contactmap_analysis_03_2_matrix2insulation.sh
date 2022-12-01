#!/bin/bash

# this command is used to calculate insulation vector from matrix file using matrix2insulation.pl tool


`perl matrix2insulation.pl -i sysmat.chr1.matrix -is 500000 -ids 200000 -im mean -bmoe 3 -nt 0.1 -v`
