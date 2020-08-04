#!/bin/bash
IDIR=`pwd`

source activate mp

python ~/toolbox_github/newsubmit/runner_general.py -stat DGC -p GYR -nt 1 -sp $IDIR/SEEDS/ -ds 0.1 -cl AB-Bottlebrush -dt 0.001 -nref 1 -ss 0.001 -fs 1.0 0.5 -L0 5.0 -npw 32 -ftol 1e-5 -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -nbbmin 99 -nbbmax 99 -dnbb 1.0 -nscmin 20 20 -nscmax 10 30 -ndsc -2 2 -chimin 0.0289 -chimax 0.0289 -dchi -0.0005 -fmin 0.33 0.67 -fmax 0.33 0.67 -df 0.1 0.1 -dirs chi nbb f nsc -dirnl 0 None 0 0 -runpoly false  -subpoly False -chain false

