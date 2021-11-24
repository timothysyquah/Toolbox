#!/bin/bash

output_file=operators_result.dat
#output_full_file=operators_result_full.dat
#for maindir in ./CL/f*/chi*/nsc*/nbb*/LAM*; do

for maindir in ./CL_c_2.0/f*/chi*/nsc*/nbb*/LAM*; do
      	cd ${maindir}
  rm ${output_file}
  touch ${output_file}
  rm ${output_full_file}
#  touch ${output_full_file}
  echo "L Hamiltonian.Real Hamiltonian.Real_Error StressXX.Real StressXX.Real_Error StressYY.Real StressYY.Real_Error StressZZ.Real StressZZ.Real_Error StressXY.Real StressXY.Real_Error StressXZ.Real StressXZ.Real_Error StressYZ.Real StressYZ.Real_Error ChemicalPotential.Real ChemicalPotential.Real_Error" > ${output_file}
  for dir in L*; do

    cd $dir

    s=`cat STATUS`
    if [ $s -eq 0 ] || [ $s -eq 2 ]; then 
      if [ -f "operators.dat" ]; then
        echo ${maindir} ${dir} ${s}
        Lvalue=$(echo "$dir" | sed 's/[^0-9]*//g')
#        var=$(python3 /home/tquah/toolbox_github/Analysis/stats.py -f operators.dat -o Hamiltonian.Real StressXX.Real StressYY.Real StressZZ.Real StressXY.Real StressXZ.Real StressYZ.Real -a -q) 
#        var="${var%"${var##*[![:space:]]}"}"   
#        echo "${Lvalue} ${var}" >> ../${output_file}
#        python3 /home/tquah/toolbox/ClusterVersion/DomainSpacing/PreprocessADT.py 
        var=$(python3 /home/tquah/toolbox/Analysis/stats.py -o Hamiltonian.Real StressXX.Real StressYY.Real StressZZ.Real StressXY.Real StressXZ.Real StressYZ.Real ChemicalPotential.Real -a -q)
        echo "${Lvalue} ${var}" >> ../${output_file}
      fi 
    fi 

    cd ../
    done
  cd ../../../../../../
done
