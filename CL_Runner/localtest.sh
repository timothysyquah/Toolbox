#!/bin/bash



IDIR=`pwd`

source activate mp

python /home/tquah/toolbox_github/CL_Runner/runner_general.py -stat DGC -p LAM -nt 1 1 -sp $IDIR/SEEDS/ \
                                                            -ds 0.1 -cl AB-Bottlebrush -dt 0.0001 -nref 1 \
                                                            -ss 0.001 -fs 0.1 1.0 -L0 5.0 -npw 512 \
                                                            -ftol 1e-6 -stol 1e-5 -bsp 1 2 -sas 1 2 \
                                                            -sac 1.0 1.0 -space 1 1 -kuhn 1 \
                                                            -gwidth 0.2 -cdensity 1 \
                                                            -nbbmin 99 -nbbmax 500 -dnbb 10.0 \
                                                            -nscmin 20 20 -nscmax 20 20 -ndsc 1 1 \
                                                            -fmin 0.5 0.5 -fmax 0.5 0.5 -df 0.01 0.01\
                                                             -chimin 0.1 -chimax 0.1 -dchi -0.005 \
                                                             -dirs  f chi nsc nbb -dirnl 0 0 0 None\
                                                              -runpoly True -subpoly False -chain True