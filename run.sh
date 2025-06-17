#!/bin/bash

input=$1
output=$2

echo "Running analysis"
echo "Input edm4eic data is at [ ${input} ] "
echo "Output at [ ${output} ]"

root -b -q diffractive_vm_analysis.cxx+\(\"${input}\",\"${output}\"\)
