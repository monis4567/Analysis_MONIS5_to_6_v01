#!/bin/bash
# -*- coding: utf-8 -*-

##get the present directory
WD=$(pwd)

# get the parent directory
#https://stackoverflow.com/questions/3790101/bash-script-regex-to-get-directory-path-up-nth-levels
PARENT=$(dirname $PWD)

BASDIR="/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6"
BASDIR="/home/sknu003/uoa00029_runs"
PRIDIR="/Analyse_MONIS5_6_2024Oct/Analysis_MONIS5_to_6_v01/output11_get_species_from_priority_table"
#PRIDIR="/priority_list_prepared_from_worms"
OUTIDIR="Analyse_MONIS5_6_v03_2024Oct/split_priority_list_prepared_from_worms"
INFL="/priority_spc.csv"

# see how to How to concatenate string variables in Bash
#https://stackoverflow.com/questions/4181703/how-to-concatenate-string-variables-in-bash
# this will be the full path to the input file
FPINF="${BASDIR}${PRIDIR}${INFL}"

#make a path to the output directory
OPD="${BASDIR}${OUTIDIR}"

echo "$OPD"
# remove 
rm -rf $OPD
mkdir $OPD

# Omit the first line
# https://stackoverflow.com/questions/7318497/omitting-the-first-line-from-any-linux-command-output
#echo "without first line"
#head -13 "${FPINF}" | cut -d';' -f1-3 | tail -n +2
#echo "only first line"
#head -13 "${FPINF}" | cut -d';' -f1-3 | head -1

# get the header line from the input fil and place in a variable
HDLN=$(head -13 "${FPINF}" | head -1)
# see the header line
echo "${HDLN}"
# use the split command to split the 'priority_spc.csv' file
# https://linuxhandbook.com/split-command/
# there are almost 20 k Aphia IDs in the 'priority_spc.csv' file, that I need to look up
# over 72 hours, I could get the Rscript doing a query for 5500 Aphia IDs. That is nearly 77 Aphia IDs per hour
# Which also equals 1833 Aphia IDs over 24 hours.
# So perhaps if I try splitting the 'priority_spc.csv' file into files with 1500 Aphia IDs each?
# cd "${OPD}"
# https://unix.stackexchange.com/questions/659018/use-output-of-cat-with-split-command-and-specified-output-directory
# Use - as the input filename. e.g. # cat file.csv | tail -n +2 | split -l 500 - /mnt/outdir
# or use '/dev/stdin' the '"${OPD}"/splt_priority_' specifies the output directory with PREFIX added
# the '-d' option dictates the split command to use digits in the filenames instead of letters
# the '-l 3' option dictates the split command to put 3 lines per file created
# the '-a 3' option dictates the split command to have 3 digits in the filenames created 
# for each file created, e.g. ..001.., ..002.. and so on
head -13 "${FPINF}" | tail -n +2 | split -l 3 -d -a 3 --additional-suffix=.txt /dev/stdin "${OPD}"/splt_priority_ 
#get a list of files in the directory
#ls -lh "${OPD}" 

# grep for 'txt' among the files in the output directory, and get the 1st line in this list
TXTF1=$(ls "${OPD}" | grep 'txt' | head -1)

#echo "$TXTF1"
# loop over files
for filename in "${OPD}"/*.txt; do
	#echo "this is the filenames"
    #echo "$filename"
    # insert a text at the beginning of a file?
	# https://stackoverflow.com/questions/9533679/how-to-insert-a-text-at-the-beginning-of-a-file
	sed -i "1s/^/${HDLN}\n/" "${filename}"
    # modify the filename to get the number on the file
    NofF=$(echo "${filename}" | sed  "s:${OPD}/splt_priority_::g" | sed 's/\.txt//g')
    #
    #get only the filename
    FlNm=$(echo "${filename}" | sed  "s:${OPD}/::g")
    echo $FlNm
    # Use the input file name to make a filename that can be inserted in the Rcode12 for the output part
    # that results from the aphia search
    ApSe_filename=$(echo "${filename}" | 
    	sed -e "s:${OPD}/splt_priority_:ApSe_:g" | 
    	sed 's/\.txt/\.csv/g')
    # append new directory and filename
    ApSe_filename=$(echo ""${OPD}"/aphia_retrieve_"$NofF"/"$ApSe_filename"")
    FlNm=$(echo ""${OPD}"/aphia_retrieve_"$NofF"/"$FlNm"")
    # make a directory to move the splitted priority list into
    mkdir "${OPD}"/aphia_retrieve_"$NofF"
    # move the splitted priority list into
    mv $filename "${OPD}"/aphia_retrieve_"$NofF"
    # copy the 'worms_safe.R' file into the directory where the splitted priority list is placed
    cp "$BASDIR"/Analysis_MONIS5_to_6_v01/Rcode_scripts/worms_safe.R "${OPD}"/aphia_retrieve_"$NofF"/. 
    # get the Rcode12 that is supposed to query the worms database with the aphiaIDs
    # write the Rcode12-file out and replace in it
    
    cat "$BASDIR"/Analysis_MONIS5_to_6_v01/Rcode_scripts/Rcode12_limit_priority_spclist_w_areaofdistr_v01.R |
    # by using sed, and append a number (that equals to splitted file) to the end of the file
    # notice that sed msut escape the quotation marks and the punctuation mark
    
    sed -e "s:flNm <- \"priority_spc\.csv\":flNm <- \"$FlNm\":g" |
    sed -e "s:flNm<-\"limited_priority_spc\.csv\":flNm <- \"$ApSe_filename\":g" > "${OPD}"/aphia_retrieve_"$NofF"/Rcode12_limit_priority_spclist_w_areaofdistr_v01.R
    
    #
    # cat "${OPD}"/aphia_retrieve_"$NofF"/Rcode12_limit_priority_spclist_w_areaofdistr_v01_"$NofF".R | grep ApSe

    # also get the sbatch submission file 
   	cat "$BASDIR"/Analysis_MONIS5_to_6_v01/bash_code/bash_remote_server_01_sbatch_rcode_get_worms_aphia_01.sh |
   	# and replace in this file using sed
	sed -e "s:#SBATCH -J 01:#SBATCH -J $NofF:g" |
	sed -e "s:#SBATCH -o stdout_01\.txt:#SBATCH -o stdout_$NofF.txt:g" |
	sed -e "s:#SBATCH -e stderr_01\.txt:#SBATCH -o stdeer_$NofF.txt:g" > "${OPD}"/aphia_retrieve_"$NofF"/bash_remote_server_01_sbatch_rcode_get_worms_aphia.sh

	# modify the files to ensure they are executable when the are to be started as jobs
	chmod 755 "${OPD}"/aphia_retrieve_"$NofF"/Rcode12_limit_priority_spclist_w_areaofdistr_v01.R
	chmod 755 "${OPD}"/aphia_retrieve_"$NofF"/bash_remote_server_01_sbatch_rcode_get_worms_aphia.sh
done

#ls -lh "${OPD}"/aphia_retrieve_001

