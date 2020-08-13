#!/bin/bash
#PBS -q gpuq
#PBS -l nodes=1:ppn=6
#PBS -l walltime=60:00:00
#PBS -V
#PBS -j oe
#PBS -N bottlebrush-AB
######################################

cd $PBS_O_WORKDIR
outdir=${PBS_O_WORKDIR}
rundir=${outdir}
username=`whoami`

############# TO USE LOCAL SCRATCH FOR INTERMEDIATE IO, UNCOMMENT THE FOLLOWING
#if [ ! -d /scratch_local/${username} ]; then
#  rundir=/scratch_local/${username}/${PBS_JOBID}
#  mkdir -p $rundir
#  cp ${PSB_O_WORKDIR}/* $rundir
#  cd $rundir
#fi
#####################################################

# Fetch the device ID for the GPU that has been assigned to the job
GPUDEV=`cat $PBS_GPUFILE | awk '{print $1}'`

IDIR=`pwd`

source activate mp

fA=0.30000
fB=0.70000
phase=O701
L0=4-8-13
npw=32-32-64
~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt -1 -sp ${IDIR}/SEEDS/ -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 0 0 -nscmax 0 0 -ndsc 2 -2 \
            -chimin 30 -chimax 30 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 100 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly False  -subpoly False -chain False -gd ${GPUDEV}

PATH=nref100.0/NscA_0.0_NscB_0.0/fA${fA}/O701Phase


~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt -1 -sp ${IDIR}/${PATH}/ -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 0 0 -nscmax 0 0 -ndsc 2 -2 \
            -chimin 3 -chimax 3 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 10 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly False  -subpoly False -chain False -gd ${GPUDEV}

PATH=nref10.0/NscA_0.0_NscB_0.0/fA${fA}/O701Phase


~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt -1 -sp ${IDIR}/${PATH}/ -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 0 0 -nscmax 0 0 -ndsc 2 -2 \
            -chimin .3 -chimax .3 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly False  -subpoly False -chain False -gd ${GPUDEV}

PATH=nref1.0/NscA_0.0_NscB_0.0/fA${fA}/O701Phase


~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt -1 -sp ${IDIR}/${PATH}/ -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 1 1 -nscmax 1 1 -ndsc 2 -2 \
            -chimin 0.15 -chimax 0.15 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly False  -subpoly False -chain False -gd ${GPUDEV}

PATH=nref1.0/NscA_1.0_NscB_1.0/fA${fA}/O701Phase

~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt -1 -sp ${IDIR}/${PATH}/ -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 2 2 -nscmax 2 2 -ndsc 2 -2 \
            -chimin 0.07 -chimax 0.07 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly False  -subpoly False -chain False -gd ${GPUDEV}

PATH=nref1.0/NscA_2.0_NscB_2.0/fA${fA}/O701Phase

~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt -1 -sp ${IDIR}/${PATH}/ -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 5 5 -nscmax 5 5 -ndsc 2 -2 \
            -chimin 0.028 -chimax 0.028 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly False  -subpoly False -chain False -gd ${GPUDEV}

PATH=nref1.0/NscA_10.0_NscB_10.0/fA${fA}/O701Phase

~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt -1 -sp ${IDIR}/${PATH}/ -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 10 10 -nscmax 10 10 -ndsc 2 -2 \
            -chimin 0.015 -chimax 0.015 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly False  -subpoly False -chain False -gd ${GPUDEV}
PATH=nref1.0/NscA_10.0_NscB_10.0/fA${fA}/O701Phase

~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p ${phase} -nt -1 -sp ${IDIR}/${PATH}/ -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 ${L0} -npw ${npw} -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -b 1.0 1.0 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 10 10 -nscmax 10 10 -ndsc 2 -2 \
            -chimin 0.015 -chimax 0.015 -dchi -0.0005 \
            -fmin ${fA} ${fB} -fmax ${fA} ${fB} -df 0.1 -0.1 \
            -nref_list 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly False  -subpoly False -chain False -gd ${GPUDEV}

# Copy back results
if [ "$rundir" != "$outdir" ]; then
  mv * ${outdir}
fi

# Force good exit code here - e.g., for job dependency
exit 0

