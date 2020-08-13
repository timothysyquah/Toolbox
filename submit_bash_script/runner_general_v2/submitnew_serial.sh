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

~/toolbox/newsubmit/runner_general_v2.py \
    -stat DGC -p GYR -nt 1 -sp $IDIR/SEEDS/ -ds 0.1 -cl AB-Bottlebrush -dt 0.001 \
            -ss 0.001 -fs 1.0 0.5 -L0 12.6 -npw 32 -ftol 1e-5 \
            -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 \
            -nbbmin 99 -nbbmax 99 -dnbb 1.0 \
            -nscmin 20 20 -nscmax 20 20 -ndsc 2 -2 \
            -chimin 0.0289 -chimax 0.0289 -dchi -0.0005 \
            -fmin 0.2 0.8 -fmax 0.8 0.2 -df 0.1 -0.1 \
            -nref_list 100 10 1 \
            -is chi nbb nsc f nref -drst nref/nref NscA-NscB/nsc fA/f -dirnl None 0-1 1 \
            -runpoly False  -subpoly False -chain False




if [ "$rundir" != "$outdir" ]; then
  mv * ${outdir}
fi

# Force good exit code here - e.g., for job dependency
exit 0

~                                                                                                                                                                                                           
