#!/bin/bash

polyftsdir=/home/tquah/PolyFTS_ALL/PolyFTS/bin/Release

IDIR=`pwd`

source activate mp

python ~/toolbox_github/newsubmit/runner_general.py -PolyFTS ${polyftsdir} -stat DGC -p DIS LAM -nt 1 1 -sp $IDIR/SEEDS/ -ds 0.1 -cl AB-Bottlebrush -dt 0.00100000 -nref 1 -ss 0.001 -fs 1.0 0.1 -L0 5.0 5.0 -npw 512 512 -ftol 1e-7 -stol 1e-8 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -nbbmin 199 -nbbmax 199 -dnbb 1 -nscmin 40 40 -nscmax 40 40 -ndsc 1 1 -chimin 0.00220273 -chimax 0.00180223 -dchi -0.00000610 -dirs nbb nsc f chi -dirnl None 0 0 0 -runpoly False -subpoly False -chain False

