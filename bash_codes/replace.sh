#!/bin/bash

echo $PWD
original='spacegroup.dat'
targets=($(grep -r -l $original *))
replace='spacegroup.in'
length=${#targets[@]}
for ((i = 1; i != length; i++)); do
#	echo "target $i: '${targets[i]}'"
#	echo $replace
	echo "${targets[i]}"
#	echo "$i"
#	echo "s|$original|$replace|g "${targets[i]}"		
done
 


