#!/bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=__NTHREADS__
#PBS -l walltime=012:00:00
#PBS -V
#PBS -j oe
#PBS -N bottlebrush
######################################
inputfile=__INFILE__
outputfile=__OUTFILE__
polyftsdir=__PFTPATH__
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

# Generate the nodes file and compute the number of
# nodes and processors per node that we have requested
cat $PBS_NODEFILE > nodes
# How many cores total do we have?
NCORES=`cat $PBS_NODEFILE | egrep -v '^#'\|'^$' | wc -l | awk '{print $1}'`
NNODES=`cat $PBS_NODEFILE | sort | uniq -c | wc -l | awk '{print $1}'`
#NTDS=`cat $PBS_NODEFILE | sort | uniq -c | head -n 1 | awk '{print $1}'`
NTDS=${PBS_NUM_PPN}

# Substitute the number of OpenMP threads into the input file
grep OpenMP_nthreads ${inputfile} > /dev/null
if [ $? -ne 0 ]; then
  echo "OpenMP_nthreads line not found in $inputfile"
  exit 1
fi
sed -i "s/\(OpenMP_nthreads\).*/\1 = ${NTDS}/g" ${inputfile}


if [ ${NNODES} -gt 1 ]; then
  mpirun -np $NNODES -machinefile nodes ${polyftsdir}/PolyFTSPLL.x ${inputfile} > ${outdir}/${outputfile}
else
  ${polyftsdir}/PolyFTSPLL.x ${inputfile} > ${outdir}/${outputfile}
fi
# Copy back results to outdir
if [ "$rundir" != "$outdir" ]; then
  mv * ${outdir}
fi

# Force good exit code here - e.g., for job dependency
exit 0
