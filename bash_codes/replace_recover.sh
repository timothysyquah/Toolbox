#!/bin/bash

echo $PWD
original='~/'
targets=($(grep -r -l $original *))
replace='~/PolyFTS_ALL/'
length=${#targets[@]}
for ((i = 0; i != length; i++)); do
#	echo "target $i: '${targets[i]}'"
#	echo $replace
#	echo "${targets[i]}"
	echo "$i"
	sed -i "s|$original|$replace|g" "${targets[i]}"		
done
 


