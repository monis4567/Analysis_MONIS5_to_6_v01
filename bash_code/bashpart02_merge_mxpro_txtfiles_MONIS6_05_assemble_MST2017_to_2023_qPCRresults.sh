#!/bin/bash
# -*- coding: utf-8 -*-

##get the present directory
WD=$(pwd)

# get the parnet directory
#https://stackoverflow.com/questions/3790101/bash-script-regex-to-get-directory-path-up-nth-levels
PARENT=$(dirname $PWD)

# #INDIR1="MONIS3_qpcrruns_20180515_w_txt_rep_data"
# INDIR1="input03a_txtreport_files_from_MxPro"
# INDIR2="input03b_txtreport_csvfiles_from_MxPro"
DATADIR="data" 
INDIR_PREFIX="qpcr_csv_files_for_"
PDDIR=$(echo "$PARENT"/"$DATADIR")
echo $PDDIR
cd "$PARENT"/"$DATADIR"
# INDIR1="input03c_txtreport_from_MxPro_MST2022_samples"
# INDIR2="input03b_txtreport_csvfiles_from_MxPro_forMST2022_samples"

# OUTDIR1="output02_merged_txtfiles_from_mxpro_for_MONIS6_MST2017_to_2022_samples"
OUTDIR1="output02_merged_txtfiles_from"
#OUTFILE1="outfile02_merged_mxpro_csvfls_MONIS_MST2017_to_2022_samples.csv"

#OUTDIR2="output03_merged_qPCRdata_txt_reports_from_MONIS_MST2017_to_2022_samples"
TMPDIR="tmp_dir"

#remove the old versions of the in- and output directory
rm -rf output02_merged_txtfiles_from*
rm -rf output02_merged_txtfiles_for_MONIS_*
#rm -rf ${OUTDIR2}
rm -rf ${INDIR2}
rm -rf ${TMPDIR}
#make new versions of the in- and output directory
#mkdir ${OUTDIR1}
#mkdir ${OUTDIR2}
#mkdir ${INDIR2}
mkdir ${TMPDIR}

# make a list that holds the names of the directories with qpcr text report files
LSSMPL=$(ls | grep "${INDIR_PREFIX}")
#make the list of samples an array you can iterate over
declare -a SMPLARRAY=($LSSMPL)

#iterate over directories with text reports
for D in ${SMPLARRAY[@]}
do
	echo "${D}"
	#done
	# use 'sed ' to remove the variable name and get the year for the sample
	YS=$(echo "${D}" | sed s:"$INDIR_PREFIX"::g)
	# make a new name for the output directory
	OUTDIR3=$(echo "output02_merged_txtfiles_for_MONIS_"$YS"")
	echo ${OUTDIR3}

	#make a new directory based on the name
	mkdir ${OUTDIR3}
	echo "${PARENT}"/"${DATADIR}"/"${D}"
	#copy the input files to the new directory
	cp "${PARENT}"/"${DATADIR}"/"${D}"/*MONIS* "${PARENT}"/"${DATADIR}"/"${OUTDIR3}"/.
	# change directory to the present output directory
	cd "${PARENT}"/"${DATADIR}"/"${OUTDIR3}"/
	# Make an if then test, to check if the file name comprises 'biorad' and 'csv' 
	#https://unix.stackexchange.com/questions/417691/how-to-check-if-file-exist-based-on-file-partial-name
	# if it does it needs to be modified in a different way 
	for FILE in *biorad*csv
	do
	    if [ -f "$FILE" ]
	    	then
	    		#echo "$FILE"
	        #printf 'at least one file exists (%s)\n' "$FILE"
	        # replace in the file name to remove '.csv' ending, and to remove the underscore between 'MONIS5' and 'biorad'
	        NWFILE=$(echo ${FILE} | sed -e 's:\.csv::g' | sed -e 's/MONIS5_biorad/MONIS5biorad/g' | sed -e 's/_rundate/_qpcrrundate/g' )
	        	#echo "$NWFILE"
			cat ${FILE} | \
				# see this weblink for how to ignore foreign characters:	#https://stackoverflow.com/questions/19242275/re-error-illegal-byte-sequence-on-mac-os-x/19770395#19770395
				#LC_ALL=C sed '/multiplex/d' | \ # delete lines with multiplex occuring
				#alter all occurences of 'not_used' to 'notused'
				LC_ALL=C sed 's/not_used/notused/g' | \
				#see how to remove tabs here: https://stackoverflow.com/questions/5398395/how-can-i-insert-a-tab-character-with-sed-on-os-x?noredirect=1&lq=1
				sed -E $'s/\t/;/g' | \
				sed -E $'s/,/;/g' | \
				#use sed to replace spaces with nothing
				sed 's/ //g' | \
				#sed -E sed '1,/Melt step/d' | \
				#print all lines after a match including the line matched 
				sed -n  '/Target/,$p' | \
				#delete line with 'ROX'
				awk '!/ROX/' | \
				# replace in the 3rd column using sed
				sed 's/;;/;NaN;/g' | \
				# replace in 4th column using awk see this question
				# this question : https://stackoverflow.com/questions/42004482/shell-script-replace-a-specified-column-with-sed
				awk 'BEGIN{FS=OFS=";"} {gsub("-","",$4)}1' | \
				# use awk and field separator 'FS=";"' to define delimeter between colunms, and then put back together columns in a different
				# order that matches the format of the MxPro text reports. Notice that no output field separator (OFS="") is defined
				# since the output needs a different arrangement of columns 
				# the BioRad file has both replicate and well type information in the 4th column 'Content' , which is why this column
				# is repeated a couple of times in awk command below

				awk 'BEGIN {FS=";"; OFS=""} {print $1,";", $4,"_",$5,"_3uL_Qty or ID",";",$4,";",$4,";Treshold;",$6,";",$7,";","RSq (dRn);Slope (dRn)"}' | \
				# the next 'sed' command now needs to edit the multiple occurences of the 'Content' columns
				# use sed to replace in the header row that now has column headers that differ from the MxPro file format
				sed -E $'s/Well;Content_Sample_3uL_Qty or ID;Content;Content;Treshold;Cq;StartingQuantity\(SQ\);RSq \(dRn\)/Well;Well Name;Replicate;Well Type;Threshold \(dRn\);Ct \(dRn\);Quantity \(copies\);RSqdRn/g' > "${NWFILE}".txt 

	#echo "H07;Unkn24_MST0527_3uL_Qty or ID;Unkn24;Unkn24;" | awk 'BEGIN{FS=OFS=";"} {gensub(/Un(k*)2/, "Upside \\1", "g", $2;}'

				# remove the original input file to prevent it being merged in later on
				rm ${FILE}
				# # copy the new file to a temporary directory to allow for inspection of the resulting file
				#cp "${NWFILE}".txt "$WD"/"$TMPDIR"/"${NWFILE}".txt

	    fi
	done

	#loop over all MONIS qpcr txt report-files and modify the filename-endings, and move the resulting files
	for FILENAME in *MONIS*txt
	do
		NEWFILENAME=$(echo "${FILENAME}" | sed 's/ /_/g' | sed 's/_-_Text_Report_Data//g' | sed 's/\#/no/g'| \
		#sed 's/\.csv/_02\.csv/g' | \
		sed 's/\.txt/\.csv/g')
		mv "${FILENAME}" "${NEWFILENAME}"
	done





#OUTFILE1_MIX_EFT="MONIS3_AssID21_Corcas_efteraar_qpcrrundate20180515_02.csv"
#OUTFILE1_MIX_FOR="MONIS3_AssID21_Corcas_foraar_qpcrrundate20180515_02.csv"

##modify the files in a loop, the loop can be used if there are more incorrect qPCR run files
#MIX_FOR_EFT_FILES="MONIS3_AssID21_Corcas_efteraar_qpcrrundate20180515_retest_07080910.csv"
#for FILE in ${MIX_FOR_EFT_FILES}
#do
#	head -1 ${FILE} > ${OUTFILE1_MIX_EFT}
#	grep efteraar ${FILE} | sed 's/efteraar//g' >> ${OUTFILE1_MIX_EFT}
#	head -1 ${FILE} > ${OUTFILE1_MIX_FOR}
#	grep foraar ${FILE} | sed 's/foraar//g' >> ${OUTFILE1_MIX_FOR}
#done


#INFILE1_MIX_FOR="MONIS3_AssID21_Corcas_foraar_qpcrrundate20180515.csv"
#cat everything, apart from the first line into the outfile
#cat ${INFILE1_MIX_FOR} | sed '1d' >> ${OUTFILE1_MIX_FOR}

#move the original mixed qPCR run file
#mv MONIS3_AssID21_Corcas_foraar_qpcrrundate20180515.csv MONIS3_AssID21_Corcas_foraar_qpcrrundate20180515_original.csv

#move the two outfiles , and give them the required filenames
#mv ${OUTFILE1_MIX_FOR} MONIS3_AssID21_Corcas_foraar_qpcrrundate20180515.csv
#mv ${OUTFILE1_MIX_EFT} MONIS3_AssID21_Corcas_efteraar_qpcrrundate20180515.csv


	# get only the first part of the filename
	CSV_FILENAMES=$(ls *.csv | sed s/\.csv//g)

	#loop over all csv-files and modify the contents, and end with making new csv-files
	for ENDING in ${CSV_FILENAMES}
	do
	# for how to convert line endings in DOS-file
	#see this website: https://stackoverflow.com/questions/2613800/how-to-convert-dos-windows-newline-crlf-to-unix-newline-lf-in-a-bash-script
		#echo "inputfile format is:"
		#file ${ENDING}.csv #check the file format of the input file 
		tr -d '\015' <${ENDING}.csv >${ENDING}.01.csv # replace DOS CRLF-end-of-lines
		#echo "outputfile format is:"
		#file ${ENDING}.01.csv #check the file format of the output file
	# use sed command on every line except the first line, but instead add ";MONISprojectno_AssayID_season_qpcrrundate" to first line
		sed '1 s,$,;qpcrno_MONISprojectno_AssayID_speciesabbr_plateno_qpcrrundate,; 1! s,$,'';'${ENDING}',' ${ENDING}.01.csv > "${ENDING}".02.csv
		#use sed to remove spaces in values
		sed 's/ //g' "${ENDING}".02.csv  | \

		# The qpcr setups I prepared for 2023 samples had an underscore between MST-SX extraction number and the subnumber for the extraction
		# This was no good as I need to be able to split by the underscore as a delimiter later on. I should have used the dash instead.
		# But I then went for a solution where I replace the underscore after the STEX number with a dash
		# following this question : https://stackoverflow.com/questions/45018156/using-sed-for-search-and-replace-multi-digits#45018256
		# I can get a dash instead of the underscore  between the MST-SX extraction number and the subnumber for extraction
		# I decided to use : [0-9]\+ which will match one or more digits.
		sed 's/\(MST[0-9]\+STEX[0-9]\+\)_/\1-/g' | \
		# The qpcr setups I prepared for 2022 samples did not have the STEX prefix added. I decided I needed to search for
		# the MST-SX extraction number and the subnumber for the extraction and then
		# replace in the sample names that have the extraction numbers added with underscores,
		#  as the underscore is used later on to split the sample name
		# to workaround this I found this webiste
		# https://unix.stackexchange.com/questions/247068/how-to-replace-values-in-a-string-using-sed-but-keep-the-string-intact
		# that explains how to substitute and retain
		# to identify the sterivex extractions I here give them the 'STEX' prefix 
		sed 's/_\([0-9]\{2\}\)_\([0-9]\{2\}\)/STEX\1-\2/g' | \
		#> "${ENDING}".03.csv

		# # # replace the MST number that has been typed incorrect in the  MXPro setup in qpcr970 to qpcr1015 for MST2022 samples
		# but MST2021-064 collected in 2021 appears to be correct in qpcr0926 to qpcr0968
	 	sed 's/MST-2021-064STEX27-17/MST-2022-064STEX27-13/g'  > "${ENDING}".03.csv

		# substitute in the intermediate file to replace tabs with semi colons
		sed -E 's/\t/;/g' ${ENDING}.03.csv > ${ENDING}.04.csv	
		
		#only print the last column in the file
		#awk -F ";" '{print $NF}' ${ENDING}03.csv
		#replace in a specified column, see: https://stackoverflow.com/questions/42004482/shell-script-replace-a-specified-column-with-sed
		# for some reason this replaces semicolons w spaces in the other columns - I do not know why
		awk -F";" '{gsub("_",";",$NF)}1' ${ENDING}.04.csv | \

		#but with sed the spaces can be replaced back to semicolons
		sed 's/ /;/g' | \
		
		#and the remaining underscores can be replaced
		sed 's/_/;/g' | \
		sed 's/Well Name/WellName/g' | \
		#and the ;WellName; can be split into new column names and added semicolons
		sed 's/;WellName;/;replno;smpltp;templvol;QtyID;/g' | \
		#sed 's/WellName/replno;smpltp;templvol;QtyID/g' | \
		#use square brackets to allow sed to replace parentheses
		sed 's/[(]//g' | sed 's/[)]//g' | \
		#the csv-file had decimal numbers with commas instead of points, replace commas w points
		sed 's/,/./g' | \
		#use sed to replace double separators with single separators
		sed 's/;;/;/g' | \
		#delete line with 'NotinUse'
		awk '!/NotinUse/' | \
		#delete line with 'notused'
		awk '!/notused/' | \
		#delete line with '---;---;Unknown;0.1'
		awk '!/---;---;Unknown;0.1/' | \
		#delete line with 'retest' - to remove the rerun tests on Corcas 
		awk '!/retest/' | \
		#delete line with 'original' - to remove the rerun tests on Corcas 
		awk '!/original/' | \
			#delete the not needed NTCs
		#awk '!/0;Qty;or;ID;NTC/' | \
		#awk '!/3uL;QtyorID;Unknown/' | \
			
		#delete line with 'Reference'
		awk '!/Reference/' > ${ENDING}.06.csv 
	done	

	#uncomment part below if some of the individual input files have different columns
	##in a loop check if the input .csv-files have the word 'Dye' included
	#CSV04_FILENAMES=$(ls *04.csv | sed s/04.csv//g)
	#
	#for ENDING in ${CSV04_FILENAMES}
	#do
	#	## check if the input file has the word 'Dye' included
	#if grep -Fq Dye ${ENDING}04.csv; then
	#	while IFS= read -r line
	#	do
	#		## Assuming the fifth column always holds the 'Dye' column
	#		## With cut fields 1 to 4 and from 6 an onwards are retained : see this website https://www.cyberciti.biz/faq/unix-linux-bsd-appleosx-skip-fields-command/
	#	     cut -d ';' -f1-4,6- <<<"$line"
	#	    ### same stuff with awk ###
	#	    ### awk '{print substr($0, index($0,$3))}' <<< "$line" ###
	#	done < "${ENDING}"04.csv > ${ENDING}05.csv
	#else NEWFILENAME=$(echo "${ENDING}05.csv")
	#	cp "${ENDING}"04.csv "${NEWFILENAME}"
	#fi
	#done
	#
	#CSV05_FILENAMES=$(ls *05.csv | sed s/05.csv//g)
	##in a loop check if the input .csv-files have the word 'Replicate' included
	#for ENDING in ${CSV05_FILENAMES}
	#do
	#	## check if the input file has the word 'Replicate' included
	#if grep -Fq Replicate ${ENDING}05.csv; then
	#	while IFS= read -r line
	#	do
	#		## Assuming the sixth column always holds the 'Replicate' column
	#		## With cut fields 1 to 5 and from 7 an onwards are retained : see this website https://www.cyberciti.biz/faq/unix-linux-bsd-appleosx-skip-fields-command/
	#	     cut -d ';' -f1-5,7- <<<"$line"
	#	    ### same stuff with awk ###
	#	    ### awk '{print substr($0, index($0,$3))}' <<< "$line" ###
	#	done < "${ENDING}"05.csv > ${ENDING}06.csv
	#else NEWFILENAME=$(echo "${ENDING}06.csv")
	#	cp "${ENDING}"05.csv "${NEWFILENAME}"
	#fi
	#done


	#see note about quotes around variables: https://stackoverflow.com/questions/2462385/getting-an-ambiguous-redirect-error
	# especially for writing to a file in a path that incl. spaces
	for FILE in *.06.csv
	do
		#write the first line of every csv-file into a temporary file
		head -1 ${FILE} >> tmp01.txt
		#echo ${FILE} >> "${PARENT}"/"$DATADIR"/"${OUTDIR1}"/tmp01.txt
	done



	OUTFILE1=$(echo "outfile02_merged_qpcr_csvfls_from_MONIS_MST_"$YS".txt")
	#get the unique lines from the tmp01.txt file, 
	#if all csv-files are set up in the same way, this should return only a single line
	#this line can put into the outputfile and serve as a header with column names
	cat tmp01.txt | uniq > "${OUTFILE1}"

	#see this website on how to use sed to get all lines apart from the first line: https://unix.stackexchange.com/questions/55755/print-file-content-without-the-first-and-last-lines/55757
	for FILE in *.06.csv
	do
		sed '1d' ${FILE} >> "${OUTFILE1}"
	done



	#head -10 "${PARENT}"/"$DATADIR"/"${OUTDIR1}"/"${OUTFILE1}"
	#tail -10 "${PARENT}"/"$DATADIR"/"${OUTDIR1}"/"${OUTFILE1}"

	# see this website about : Echo newline in Bash prints literal \n
	# https://stackoverflow.com/questions/8467424/echo-newline-in-bash-prints-literal-n
	# -e flag did it for me, which "enables interpretation of backslash escapes"

	echo -e " \n make sure there only is one unique line for headers \n"

	#see the content of the tmp01.txt file, to check all input files have the same header
	#using the uniq command in the end , will make sure it only returns the unique lines 
	cat tmp01.txt | uniq

	echo -e " \n make sure there only is one unique species name per species \n"
	#Print the third column
	# and sort the output, and get only uniq values
	#awk -F";" '{print $3}' "${PARENT}"/"$DATADIR"/"${OUTDIR1}"/"${OUTFILE1}" | sort | uniq
	awk -F";" '{print $16}' "${OUTFILE1}" | sort | uniq

	echo -e " \n make sure there only is one unique sample type name per sample type \n"
	#Print the fourth column
	# and sort the output, and get only uniq values
	#awk -F";" '{print $2}' "${PARENT}"/"$DATADIR"/"${OUTDIR1}"/"${OUTFILE1}" | sort | uniq
	awk -F";" '{print $3}' "${OUTFILE1}" | sort | uniq
	# 
	printf "${OUTFILE1}"
	# 
	#cat "${PARENT}"/"$DATADIR"/"${OUTDIR1}"/"${OUTFILE1}" | head -30
	#delete all the temporary files
	rm *.01.csv
	rm *.02.csv
	rm *.03.csv
	rm *.04.csv
	#rm *.05.csv
	rm *.06.csv
	rm tmp01.txt
# change back to the directory that holds the directories with input files
cd "$PARENT"/"$DATADIR"
done