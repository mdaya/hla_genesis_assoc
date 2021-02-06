#!/bin/bash

HLA_LA_output_dir=$1
min_Q1=$2
min_avg_cov=$3
min_MAF=$4

incl_genes=$5   #Comma delimitted list of genes to include
out_file_prefix=$6

#Create PLINK ped file
cat /home/analyst/convert_hla_la_to_plink.R | R --vanilla --args $HLA_LA_output_dir $min_Q1 $min_avg_cov $incl_genes $out_file_prefix

#Convert to PLINK binary file, applying the MAF filter
plink --file ped_${out_file_prefix} --maf $min_MAF --make-bed --out $out_file_prefix
