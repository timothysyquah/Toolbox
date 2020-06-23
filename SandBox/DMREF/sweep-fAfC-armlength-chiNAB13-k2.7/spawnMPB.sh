#!/bin/bash

idir=`pwd`

epsA=7.8
epsB=1.0
epsC=1.0

PhaseList=(AGYR32 ADIA32 ACYL)
for d in Nsc0.040/fAfC*; do
    for phase in ${PhaseList[@]}; do
        wdir=$d/${phase}Phase
        if [ ! -e $wdir ]; then
            echo "$wdir does not exist. Skipping!"; continue
        fi
        mpb_output_file="mpb-${epsA}-${epsB}-${epsC}.out"
        if [ -e $wdir/${mpb_output_file} ]; then
          echo "$wdir/${mpb_output_file} exists...skipping"
          continue
        fi
        if [ -e $wdir/STATUS ]; then
          status=`cat $wdir/STATUS`
          if [ $status -ne 2 ]; then
            echo "${wdir} not converged (STATUS=$status)...skipping"
            continue
          fi
        fi

        echo "Running $wdir"

        #mkdir -p $wdir
        cd $wdir

        # generate h5 file
        ~/tools/photonics/PolyFTS_to_epsilon_h5.py --epsilons $epsA $epsB $epsC --input density.bin

        mpb_input_file="mpb_${phase}.py"
        cp $idir/SEEDS/${mpb_input_file} .

        cmd="sed -e 's/__INFILE__/$mpb_input_file/g' 
              -e 's/__OUTFILE__/$mpb_output_file/g' 
                < $idir/SEEDS/submit_mpb.sh > submit_mpb.sh"
        #echo $cmd
        eval $cmd

        qsub submit_mpb.sh
        cd $idir
    done
done
