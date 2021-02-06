#!/bin/bash

null_rdata_file=$1
plink_bed_file=$2
plink_bim_file=$3
plink_fam_file=$4
out_file_name=$5

cat /home/analyst/run_genesis_assoc.R | R --vanilla --args $null_rdata_file $plink_bed_file $plink_bim_file $plink_fam_file $out_file_name
