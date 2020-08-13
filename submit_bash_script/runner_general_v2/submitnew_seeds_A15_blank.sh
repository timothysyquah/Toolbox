#!/bin/bash

IDIR=`pwd`

source activate mp

fA=0.33000
fB=0.67000
phase=A15
L0=5
npw=32

~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/SEEDS/ -ds 0.1 -cl AB-Bottlebrush -dt 0.01 \
    -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
    -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 0.5 \
    -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
    -nscmin 0 0 -nscmax 0 0 -ndsc 2 -2 \
    -chimin 30 -chimax 30 -dchi -0.0005 \
    -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
    -nref_list 100 \
    -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 0 \
    -runpoly False -subpoly False -chain False



~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/SEEDS/  -ds 0.1 -cl AB-Bottlebrush -dt 0.01 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 0.5 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 0 0 -nscmax 0 0 -ndsc 2 -2 \
            -chimin 3 -chimax 3 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 10 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 0 \
            -runpoly False -subpoly False -chain False



~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/SEEDS/  -ds 0.1 -cl AB-Bottlebrush -dt 0.01 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 0.5 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 0 0 -nscmax 0 0 -ndsc 2 -2 \
            -chimin .3 -chimax .3 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 0 \
            -runpoly False -subpoly False -chain False



~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/SEEDS/  -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 1 2 -nscmax 1 2 -ndsc 2 -2 \
            -chimin 0.15 -chimax 0.15 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 0 \
            -runpoly False -subpoly False -chain False


~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/SEEDS/  -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 2 4 -nscmax 2 4 -ndsc 2 -2 \
            -chimin 0.1 -chimax 0.1 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 0 \
            -runpoly False -subpoly False -chain False 


~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/SEEDS/  -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 4 8 -nscmax 4 8 -ndsc 2 -2 \
            -chimin 0.06 -chimax 0.06 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 0 \
            -runpoly False -subpoly False -chain False 


~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/SEEDS/  -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 8 16 -nscmax 8 16 -ndsc 2 -2 \
            -chimin 0.033 -chimax 0.033 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 0 \
            -runpoly False -subpoly False -chain False 


~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/SEEDS/  -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 12 28 -nscmax 12 28 -ndsc 2 -2 \
            -chimin 0.023 -chimax 0.029 -dchi 0.001 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is nbb nsc f nref chi -drst nref/nref NscA-NscB/nsc fA/f chiAB/chi -dirnl None 0-1 0 0\
            -runpoly False -subpoly False -chain False 


