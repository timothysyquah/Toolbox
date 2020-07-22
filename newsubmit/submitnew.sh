#!/bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -V
#PBS -j oe
#PBS -N bottlebrush-ABC
######################################
polyftsdir=/home/tquah/PolyFTS_ALL/PolyFTS/bin/Release

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

cat $PBS_NODEFILE > nodes

# Run the job



IDIR=`pwd`

source activate mp

python ./runner_general.py -PolyFTS ${polyftsdir} -stat DGC -p DIS LAM -nt 1 1 -sp $IDIR/SEEDS/ -ds 0.1 -cl AB-Bottlebrush -dt 0.001 -nref 1 -ss 0.001 -fs 1.0 1.0 -L0 5.0 5.0 -npw 512 512 -ftol 1e-5 -stol 1e-4 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -nbbmin 99 -nbbmax 99 -dnbb 1.0 -nscmin 20 20 -nscmax 20 20 -ndsc 1 1 -chimin 0.1 -chimax 0.005 -dchi -0.005 -dirs nbb nsc f chi -dirnl None 0 0 0 -runpoly False -subpoly False -chain False



# Copy back results
if [ "$rundir" != "$outdir" ]; then
  mv * ${outdir}
fi

# Force good exit code here - e.g., for job dependency
exit 0
