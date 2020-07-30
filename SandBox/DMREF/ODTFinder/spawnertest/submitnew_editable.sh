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

python ~/toolbox/newsubmit/runner_general.py -stat DGC -p DIS LAM -nt 1 1 -sp $IDIR/SEEDS/ -ds 0.1 -cl AB-Bottlebrush -dt __dt__ -nref 1 -ss 0.001 -fs 1.0 0.1 -L0 5.0 5.0 -npw 512 512 -ftol 1e-7 -stol 1e-8 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -nbbmin __nbbmin__ -nbbmax __nbbmax__ -dnbb __dnbb__ -nscmin __nscmin__ __nscmin__ -nscmax __nscmax__ __nscmax__ -ndsc __dnsc__ __dnsc__ -chimin __chimin__ -chimax __chimax__ -dchi __dchi__ -dirs nbb nsc f chi -dirnl None 0 0 0 -runpoly True -subpoly False -chain True



# Copy back results
if [ "$rundir" != "$outdir" ]; then
  mv * ${outdir}
fi

# Force good exit code here - e.g., for job dependency
exit 0
