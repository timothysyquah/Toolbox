#!/bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -V
#PBS -j oe
#PBS -N bottlebrush-ABC
######################################
inputfile=LAM.in
outputfile=LAM.out
polyftsdir=/home/tquah/PolyFTS_ALL/PolyFTS/bin/Release
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

cat $PBS_NODEFILE > nodes

# Run the job
${polyftsdir}/PolyFTS.x ${inputfile} > ${outdir}/${outputfile}
# Copy back results
if [ "$rundir" != "$outdir" ]; then
  mv * ${outdir}
fi

# Force good exit code here - e.g., for job dependency
exit 0
