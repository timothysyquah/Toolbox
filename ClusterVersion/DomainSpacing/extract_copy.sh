#!/bin/bash





IDIR=`pwd`

rm -r density

mkdir density


for i in L*;
do
  cd $i
  cp DensityOperator.vtk ../density/${i}_densityoperator.vtk
  cp density.vtk ../density/${i}_density.vtk
  cd ../


done

