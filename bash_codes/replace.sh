#!/bin/bash

echo $PWD
original='~/PolyFTS-update-symmetrizer/bin/Release'
targets=($(grep -r -l $original *))
replace='/home/lequieu/PolyFTS-update-symmetrizer/bin/Release'
length=${#targets[@]}
for ((i = 0; i != length; i++)); do
#	echo "target $i: '${targets[i]}'"
#	echo $replace
#	echo "${targets[i]}"
	echo "$i"
	echo "s|$original|$replace|g "${targets[i]}"		
done
 


