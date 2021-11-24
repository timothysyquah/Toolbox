#!/bin/bash





IDIR=`pwd`


for i in L*;
do
  cd $i
  python ~/PolyFTS_ALL/PolyFTS/tools/plot/PolyFTS_to_VTK.py DensityOperator.dat
  python ~/PolyFTS_ALL/PolyFTS/tools/plot/PolyFTS_to_VTK.py 
  cd ../


done

