#!/bin/bash
# -*- coding: utf-8 -*-

##get the present directory
WD=$(pwd)

# get the parent directory
#https://stackoverflow.com/questions/3790101/bash-script-regex-to-get-directory-path-up-nth-levels
PARENT=$(dirname $PWD)

#BASDIR="/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6"
BASDIR="/home/sknu003/uoa00029_runs"
PRIDIR="/Analyse_MONIS5_6_2024Oct/Analysis_MONIS5_to_6_v01/output11_get_species_from_priority_table"
PRJDIR="/Analyse_MONIS5_6_v03_2024Oct/Analysis_MONIS5_to_6_v01"
#PRIDIR="/priority_list_prepared_from_worms"
OUTDIR1="/Analyse_MONIS5_6_v03_2024Oct/split_priority_list_prepared_from_worms"
OUTDIR2="/Analyse_MONIS5_6_v03_2024Oct/all_aphiaIDs_and_locations"
INFL="/priority_spc.csv"

# see how to How to concatenate string variables in Bash
#https://stackoverflow.com/questions/4181703/how-to-concatenate-string-variables-in-bash
# this will be the full path to the input file
FPINF="${BASDIR}${PRIDIR}${INFL}"

#make a path to the input directory
IPD="${BASDIR}${OUTDIR1}"

#make a path to the output directory
OPD="${BASDIR}${OUTDIR2}"

echo "${OPD}"
# remove 
rm -rf ${OPD}
# make a new version of the output directory
mkdir ${OPD}
for dir in ${IPD}/*/     # list directories in the form "/tmp/dirname/"
do
	cd ${WD}
    #dir=${dir%*/}      # remove the trailing "/"
    #echo "${dir##*/}"    # print everything after the final "/"
    echo ${dir}
    cd ${dir}
    cp ApSe_* ${OPD}/.
    cp location_* ${OPD}/.
done
cd ${WD}

rm -rf all_aphiaIDs_and_locations.tar.gz
# tar <options> <archive-file> <files-or-directory>
tar -czvf  all_aphiaIDs_and_locations.tar.gz all_aphiaIDs_and_locations/*



