#!/bin/bash
idir=`pwd`
for fAfC in `seq 0.10 0.05 0.40`; do
  #wdir="../delta0.00/fAfC${fAfC}/ADIA1Phase/"
  wdir="../../Nsc0.040/fAfC${fAfC}/ADIAC64Phase/"
  echo $wdir
  cd $wdir
 
  label=`echo "scale=0;$fAfC*100 /1" | bc`
  echo $label
  vtkname="density_${label}.vtk"
  ~/tools/plot/PolyFTS_to_VTK.py density.bin

  cp density.vtk ${idir}/${vtkname}

  cd $idir
  
done
