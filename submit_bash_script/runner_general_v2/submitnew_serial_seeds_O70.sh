#!/bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00
#PBS -V
#PBS -j oe
#PBS -N BB-AB_chain
######################################


cd $PBS_O_WORKDIR
outdir=${PBS_O_WORKDIR}
rundir=${outdir}
username=`whoami`
cat $PBS_NODEFILE > nodes


IDIR=`pwd`

source activate mp

fA=0.36000
fB=0.64000
phase=O701
L0=4-8-13
npw=32-32-64

~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/SEEDS/ -spt Main -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 0 0 -nscmax 0 0 -ndsc 2 -2 \
            -chimin 30 -chimax 30 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 100 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly true  -subpoly False -chain true 

PATH=nref100.0/NscA_0.0_NscB_0.0/fA${fA}/O701Phase


~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/${PATH}/ -spt Job -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 0 0 -nscmax 0 0 -ndsc 2 -2 \
            -chimin 3 -chimax 3 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 10 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly true  -subpoly False -chain true 

PATH=nref10.0/NscA_0.0_NscB_0.0/fA${fA}/O701Phase


~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/${PATH}/ -spt Job -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 0 0 -nscmax 0 0 -ndsc 2 -2 \
            -chimin .3 -chimax .3 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly true  -subpoly False -chain true 

PATH=nref1.0/NscA_0.0_NscB_0.0/fA${fA}/O701Phase


~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/${PATH}/ -spt Job -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 1 1 -nscmax 1 1 -ndsc 2 -2 \
            -chimin 0.15 -chimax 0.15 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly true  -subpoly False -chain true 

PATH=nref1.0/NscA_1.0_NscB_1.0/fA${fA}/O701Phase

~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/${PATH}/ -spt Job -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 2 2 -nscmax 2 2 -ndsc 2 -2 \
            -chimin 0.1 -chimax 0.1 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly true  -subpoly False -chain true 

PATH=nref1.0/NscA_2.0_NscB_2.0/fA${fA}/O701Phase

~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/${PATH}/ -spt Job -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 5 5 -nscmax 5 5 -ndsc 2 -2 \
            -chimin 0.05 -chimax 0.05 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly true  -subpoly False -chain true 

PATH=nref1.0/NscA_10.0_NscB_10.0/fA${fA}/O701Phase

~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/${PATH}/ -spt Job -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 10 10 -nscmax 10 10 -ndsc 2 -2 \
            -chimin 0.027 -chimax 0.027 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly true  -subpoly False -chain true 
PATH=nref1.0/NscA_10.0_NscB_10.0/fA${fA}/O701Phase

~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/${PATH}/ -spt Job -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 15 15 -nscmax 15 15 -ndsc 2 -2 \
            -chimin 0.0187 -chimax 0.0187 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly true  -subpoly False -chain true 

PATH=nref1.0/NscA_15.0_NscB_15.0/fA${fA}/O701Phase

~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt 1 -sp ${IDIR}/${PATH}/ -spt Job -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 20 20 -nscmax 20 20 -ndsc 2 -2 \
            -chimin 0.0143 -chimax 0.0143 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly true  -subpoly False -chain true 

# Copy back results
if [ "$rundir" != "$outdir" ]; then
  mv * ${outdir}
fi

# Force good exit code here - e.g., for job dependency
exit 0

