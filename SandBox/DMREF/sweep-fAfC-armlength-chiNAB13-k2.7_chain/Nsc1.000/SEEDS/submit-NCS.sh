#!/bin/bash

#SBATCH --job-name=NCS-pentablock
#SBATCH --time=02:00:00
#SBATCH --partition=batch
#SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --mail-type=NONE
##SBATCH --mail-user=lequieu@mrl.ucsb.edu

##PBS -N __JOBNAME__

######################################
inputfile=__INFILE__
outputfile=__OUTFILE__
#polyftsdir=/home/lequieu/PolyFTS/bin/Release
polyftsdir=/home/lequieu/PolyFTS-develop/bin/Release
######################################

cd $SLURM_SUBMIT_DIR
outdir=${SLURM_SUBMIT_DIR}
rundir=${outdir}
username=`whoami`


# Run the job
${polyftsdir}/PolyFTS.x ${inputfile} > ${outdir}/${outputfile}

# Force good exit code here - e.g., for job dependency
exit 0
