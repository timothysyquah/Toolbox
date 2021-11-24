#!/bin/bash

output_file=averaged_operators.dat
python3 /home/tquah/toolbox/ClusterVersion/DomainSpacing/PreprocessADT.py 
var=$(python3 /home/tquah/toolbox/Analysis/stats.py -f RW_operators.dat -o Hamiltonian.Real StressXX.Real StressYY.Real StressZZ.Real StressXY.Real StressXZ.Real StressYZ.Real ChemicalPotential.Real -a -q -dt)

echo "Hamiltonian.Real Hamiltonian.Real_Error StressXX.Real StressXX.Real_Error StressYY.Real StressYY.Real_Error StressZZ.Real StressZZ.Real_Error StressXY.Real StressXY.Real_Error StressXZ.Real StressXZ.Real_Error StressYZ.Real StressYZ.Real_Error ChemicalPotential.Real ChemicalPotential.Real_Error" > ${output_file}
echo "${var}" >> ${output_file}
