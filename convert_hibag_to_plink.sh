#!/bin/bash

out_file_prefix=$1
MAF_filter=$2

#Shift to list of files
shift
shift

#Create PLINK ped file
cat /home/analyst/convert_hibag_to_plink.R | R --vanilla --args ped_${out_file_prefix} "$@"

#Convert to PLINK binary file, applying the MAF filter
plink --file ped_${out_file_prefix} --maf $MAF_filter --make-bed --out $out_file_prefix
