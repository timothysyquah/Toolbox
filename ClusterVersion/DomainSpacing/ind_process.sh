#!/bin/bash

output_file=averaged_operators.dat
rm ${output_file}
touch ${output_file}
echo "Hamiltonian.Real Hamiltonian.Real_Error StressXX.Real StressXX.Real_Error StressYY.Real StressYY.Real_Error StressZZ.Real StressZZ.Real_Error StressXY.Real StressXY.Real_Error StressXZ.Real StressXZ.Real_Error StressYZ.Real StressYZ.Real_Error ChemicalPotential.Real ChemicalPotential.Real_Error" > ${output_file}
echo ${maindir} ${dir} ${s}
python3 /home/tquah/toolbox/ClusterVersion/DomainSpacing/PreprocessADT.py 
var=$(python3 /home/tquah/toolbox/Analysis/stats.py -f RW_operators.dat -o Hamiltonian.Real StressXX.Real StressYY.Real StressZZ.Real StressXY.Real StressXZ.Real StressYZ.Real ChemicalPotential.Real -a -q -dt)
echo "${var}" >> ${output_file}

