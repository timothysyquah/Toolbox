#!/bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00
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

python /home/tquah/toolbox_github/CL_Runner/runner_general.py -stat __chain_stats__ -p LAM -nt 1 -sp $IDIR/SEEDS/ \
                                                            -ds 0.1 -cl AB-Bottlebrush -dt __dt__ -nref 1 \
                                                            -ss 0.001 -fs 0.1 1.0 -L0 5.0 -npw __npw__ \
                                                            -ftol 1e-5 -stol 1e-4 -bsp 1 2 -sas 1 2 \
                                                            -sac 1.0 1.0 -space 1 1 -kuhn 1 \
                                                            -gwidth __gwidth__ -cdensity 1 --invzeta  __invzeta__\
                                                            -nbbmin __nbbstart__ -nbbmax 500 -dnbb 10.0 \
                                                            -nscmin __nsc__ __nsc__ -nscmax __nsc__ __nsc__ -ndsc 1 1 \
                                                            -fmin 0.5 0.5 -fmax 0.5 0.5 -df 0.01 0.01\
                                                             -chimin __chi__ -chimax __chi__ -dchi -0.005 \
                                                             -dirs  f chi nsc nbb -dirnl 0 0 0 None\
                                                              -runpoly True -subpoly False -chain True


# Copy back results
if [ "$rundir" != "$outdir" ]; then
  mv * ${outdir}
fi

# Force good exit code here - e.g., for job dependency
exit 0
