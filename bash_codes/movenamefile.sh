#!/bin/bash

path=$PWD

directname="bandexport"
directoryname=$path/$directname
mkdir $directoryname

targetname="band*"

targets=$(find . -name $targetname)

length=${#targets[@]}

for ((i=0; i!=length; i++)); do

	echo "${targets[i]}"

	echo 
done
