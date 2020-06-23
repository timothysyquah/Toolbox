#!/bin/bash

PhaseList=(AGYR32 ADIA32 CsCl ACYL)
idir=`pwd`
for phase in ${PhaseList[@]}; do
  #for i in `find . -name "${phase}Phase"`; do
  for i in Nsc0.100/fAfC*/${phase}Phase/; do
    s=`cat $i/STATUS`
    if [ $s -eq 1 ]; then 
      cd $i
      echo $i
      sed -i '/DT/s/=.*/= 0.01/g' ${phase}.in
      sed -i '/walltime/s/=.*/=12:00:00/g' submit.sh
      qsub submit.sh
      cd $idir
    fi
  done
done
