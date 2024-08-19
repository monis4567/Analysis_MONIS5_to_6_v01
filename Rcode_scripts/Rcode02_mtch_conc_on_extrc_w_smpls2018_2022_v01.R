#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

########################################################################################
# R-code for assembly of  "MST samples "

########################################################################################
getwd()
########################################################################################
#remove everything in the working environment, without a warning!!
#rm(list=ls())
########################################################################################
#install packages
#get readxl package
# if(!require(readxl)){
#   install.packages("readxl")
# }  

# if(!require(readr)){
#   install.packages("readr")
# }  
#get ggplot package
# if(!require(ggplot2)){
#   install.packages("ggplot2")
# }  
#get pdp package
# if(!require(pdp)){
#   install.packages("pdp")
# }  
#get chron package
# if(!require(chron)){
#   install.packages("chron")
# }  
library(chron)
library(cowplot)
library(dplyr)
library(ggforce)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(googleway)
library(maps)
library(openxlsx)
library(pdp)
library(readr)
library(readxl)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sf)

#name working directory
#wd00 ="/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6/Analysis_MONIS5_to_6"
#wd00 = "/Users/steenknudsen/Documents/Documents/NIVA_Ansaettelse_2020/NOVANA_proever_2018_2019"
# define dir w output files
wd03 <- "output03_match_MONIS_conc_on_extraction_with_water_sample"
#wd03 <- "output"
#set working directory
#setwd(wd00)
wd00 <- getwd()
wd01 <- paste0(wd00,"/data")
#make complete path to output dir
wd00_wd03 <- paste(wd00,"/",wd03,sep="")
# delete previous versions of the outpur directory
unlink(wd00_wd03, recursive=TRUE)
# # create a new version of the directory you just deleted
dir.create(wd00_wd03)
#check the wd

#_______________________________________________________________________
# section 01 -  start -  read in qpcr text report files
#_______________________________________________________________________
# define prefix for  w input files
prefix_ind_01 = "output02_merged_txtfiles_for_MONIS"
#list all directories in pwd - 
ls.dirs01 <- list.dirs(path = wd01, full.names = TRUE, recursive = TRUE)
# get a list of directories that has the prefix for input files
ind_01 <- ls.dirs01[grep(prefix_ind_01,ls.dirs01)]
#make an empty list to 
ls.txtreps <- list()
# make a count that can increase for each iteration
i <- 1
# iterate over directories that holds the input files
for (dir in ind_01)
  {
  #print(dir)
  #}
  # define file prefix to grep for
  prefix_indf_02 <- "outfile02_merged_qpcr_csvfls_from_MONIS"
  # get a list of all the files in the directory
  ls.inf01 <- list.files(path = dir,full.names = T, recursive = T)
  # among these files, then grep for the file with the prefix above
  infile <- ls.inf01[grep(prefix_indf_02, ls.inf01   )]
  # grep only the text files   
  infile <- ls.inf01[grep("\\.txt", infile)]
  # then read this file, if more files in the directory match this prefix
  # then this part will fail
  # see https://www.statology.org/r-duplicate-row-names-not-allowed/
  # # if the read csv otherwise throws an error about column order
  #
  tibl_txtrep <- read.csv2(file=infile,
                           sep=";")
  #colnames(tibl_txtrep)
  # add the file read in to the list of text reports
  ls.txtreps[[i]] <- tibl_txtrep
  # increase the number that keeps count of iterations
  i <- i+1
# end the iteration over input directories
}
# Use 'plyr::ldply' instead of 'datatable::rbindlist'
# as 'datatable::rbindlist' was hysterical about columns
# see this question: https://stackoverflow.com/questions/55706560/make-rbindlist-skip-ignore-or-change-class-attribute-of-the-column
#bind the rows in each list in to one data frame
df_txtr01 <- plyr::ldply(ls.txtreps, data.frame)

# ##___ERROR CHECK
# df_txtr01.3 <- df_txtr01[grepl("MST-2021-064",df_txtr01$smpltp),]
# unique(df_txtr01.3$smpltp)
# ##___ERROR CHECK

# use gsub in columns to standardize replicate numbers
df_txtr01$replno<- gsub("[A-Za-z]","",df_txtr01$replno)
# use gsub to substitute the WellType
df_txtr01$WellType <- gsub("[0-9]","",df_txtr01$WellType)
df_txtr01$WellType <- gsub("^Unkn$","Unknown",df_txtr01$WellType)
df_txtr01$WellType <- gsub("^Std$","Standard",df_txtr01$WellType)
df_txtr01$WellType <- gsub("^NegCtrl$","NTC",df_txtr01$WellType)
# use gsub in columns to standardize replicate numbers
df_txtr01$Replicate<- gsub("[A-Za-z]","",df_txtr01$Replicate)
#pad with zeros to 4 characters
#see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  qpNo <- stringr::str_pad(gsub("qpcr","",df_txtr01$qpcrno), 4, pad = "0")
# paste 'qpcr' back
df_txtr01$qpcrno <- paste0("qpcr",qpNo)
# use gsub in columns to standardize AssID names
df_txtr01$AssayID <- gsub("assid","AssID",df_txtr01$AssayID)
# remove the extra zeroes added by typing error for the 2022 runs
df_txtr01$qpcrrundate <- gsub("qpcrrundate202200","qpcrrundate20220",df_txtr01$qpcrrundate)
# identify grep elements and place in an object
MSn<- df_txtr01$smpltp[grep("^MST[0-9]{3}$",df_txtr01$smpltp)]
MSn <- gsub("MST","",MSn)
#pad with zeros to 4 characters
#see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
MSn <- stringr::str_pad(gsub("qpcr","",MSn), 4, pad = "0")
# replace the grep elements
df_txtr01$smpltp[grep("^MST[0-9]{3}$",df_txtr01$smpltp)] <- paste0("MST",MSn,sep="")


#
#write.csv(df_txtr01,file=paste0(wd00,"/tmp_dir/read_qpcr_tables.csv"),sep = ";")
#_______________________________________________________________________
# section 01 -  end -  read in qpcr text report files
#_______________________________________________________________________

#_______________________________________________________________________
# section 02 -  start -  read in excel files with concentrations on extractions from filters
#_______________________________________________________________________

# https://stackoverflow.com/questions/46305724/merging-of-multiple-excel-files-in-r
# https://medium.com/@niharika.goel/merge-multiple-csv-excel-files-in-a-folder-using-r-e385d962a90a
# https://stackoverflow.com/questions/54474404/how-to-detect-time-when-reading-from-an-excel-sheet-using-r
# get list of excel files
lst_xlsfl_in_wd01 <- list.files(wd01)[grepl("xls",list.files(wd01))]
# get the xls files that has "MST_ekstraktioner_fra_filterproever_" 
# in the filename
fl_ekstr <- lst_xlsfl_in_wd01[grepl("MST_ekstraktioner_fra_filterproever_",
                                    lst_xlsfl_in_wd01)]
# make an empty list to add read in tables to
ls.conc.ex <- list()
# make a number that can increase for each iteration
j <-1
#fi <- fl_ekstr[1]
# iterate over files with concentrations on the extractions from filters
for (fi in fl_ekstr)
  {
  print(fi)
  fi <- paste0(wd01,"/",fi)
  #read in excel file to a tibble
  tb_co01 <- read_xls(fi)
  #tb_co01 <- openxlsx::read.xlsx(fi)
  #tb_co01$Sub_Nmb_for_extraction
  # make the tibble a data frame
  df_co01 <- as.data.frame(tb_co01)
  # make all columns characters
  #library(dplyr)
  df_co01 <- df_co01 %>%
    dplyr::mutate(across(everything(), as.character))
  #View(df_co01)
    # add the data frame to the list of tables with concentrations
  ls.conc.ex[[j]]<- df_co01
  # increase the count by 1
  j<- j+1
# end iteration over files
}  
# Use 'plyr::ldply' instead of 'datatable::rbindlist'
# as 'datatable::rbindlist' was hysterical about columns
# see this question: https://stackoverflow.com/questions/55706560/make-rbindlist-skip-ignore-or-change-class-attribute-of-the-column
#bind the rows in each list in to one data frame
df_cex01 <- plyr::ldply(ls.conc.ex, data.frame)

df_cex01 <- df_cex01[order(df_cex01$MST.nummer),]
#View(df_cex01)
# identify numbers that do not have an 'MST' in front
MST_Numb.miss.MST<- df_cex01$MST.nummer[grep("^[0-9]+{2}$",df_cex01$MST.nummer)]
#pad with zeros to 4 characters
#see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
MST_Numb.miss.MST <- stringr::str_pad(MST_Numb.miss.MST, 4, pad = "0")
# use paste to add the MST lettercode in front
df_cex01$MST.nummer[grep("^[0-9]+{2}$",df_cex01$MST.nummer)] <- paste("MST-",
                                                MST_Numb.miss.MST,sep = "")
# use gsub to ensure that all MST numbers start with 'MST-' 
df_cex01$MST.nummer <- gsub("^MST_","MST",df_cex01$MST.nummer)
df_cex01$MST.nummer <- gsub("^MST-","MST",df_cex01$MST.nummer)
# identify MSt samples with too few leading zeroes
MSTwfz <- df_cex01$MST.nummer[grep("^MST[0-9]{3}$",df_cex01$MST.nummer)]
# and then replace with a new MST number
df_cex01$MST.nummer[grep("^MST[0-9]{3}$",
                         df_cex01$MST.nummer)] <- gsub("MST",
                                                       "MST0",MSTwfz)
# make acolumn with unique proeve numbers
df_cex01$U_Pr_Nr <- gsub("-","",df_cex01$MST.nummer)
# # identify duplicated MST numbers
# dupl_MSTnos <- df_cex01$MST.nummer[duplicated(df_cex01$MST.nummer)]
# # check of the vector is in the rows with MST numbers
# #https://stackoverflow.com/questions/11612235/select-rows-from-a-data-frame-based-on-values-in-a-vector
# # and only retain rows where the MST numbers in the vector is also present
# df_co01.dpl <- df_co01[df_co01$`MST nummer` %in% dupl_MSTnos,]

# substitute to remove the ' mark in front of dates
df_cex01$Dato_for_indsml <- gsub("'","",df_cex01$Dato_for_indsml)
# identify non-NA in collection-data column in data frame
# with filter extractions
dcfil<- df_cex01$Dato_for_indsml[!is.na(df_cex01$Dato_for_indsml)]
# make a data frame with numbers for months
df_mntabb <- data.frame(cbind(month.abb,seq(1,12,1)))
colnames(df_mntabb) <- c("month.abb", "month.no")
# get the month part of the date
m.dcfil <- gsub("^(.*)-(.*)-(.*)$","\\2",dcfil)
ml.dcfil <- m.dcfil[grepl("[A-Za-z]",m.dcfil)]
mn.dcfil <- df_mntabb$month.no[match(ml.dcfil,df_mntabb$month.abb)]
#pad with zeros to 2 characters
#see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
mn.dcfil <- stringr::str_pad(mn.dcfil, 2, pad = "0")
# get only letter code months
m.dcfil[grepl("[A-Za-z]",m.dcfil)] <- mn.dcfil
# paste years and modified months and dates back together
yyfl <- gsub("^(.*)-(.*)-(.*)$","\\1",dcfil)
ddfl <- gsub("^(.*)-(.*)-(.*)$","\\3",dcfil)
dcfil <- paste(yyfl,m.dcfil,ddfl,sep="-")
# replace the non NA in the dates column
df_cex01$Dato_for_indsml[!is.na(df_cex01$Dato_for_indsml)] <- dcfil

# paste together the MST number and the date for collection
df_cex01$MSTno_di <- paste(df_cex01$U_Pr_Nr,df_cex01$Dato_for_indsml,sep="::")
# re order the data frame
df_cex01 <- df_cex01[order(df_cex01$MSTno_di),]

df_cex01_MSTno_di.2021028 <- df_cex01$MSTno_di[grepl("2021028",df_cex01$MSTno_di)]
unique(df_cex01_MSTno_di.2021028)
#_______________________________________________________________________
# section 02 -  end -  read in excel files with concentrations on extractions from filters
#_______________________________________________________________________

#_______________________________________________________________________
# section 03 -  start -  read in excel files with sample locations
#_______________________________________________________________________

# define infile prefix
inf_prefix2 <- "sample_locations_MST_"
# get list of 'excel'csv' files
lst_xlsfl_in_wd01 <- list.files(wd01)[grepl("csv",list.files(wd01))]
# get the xls files that has "inf_prefix2" 
# in the filename
fls <- lst_xlsfl_in_wd01[grepl(inf_prefix2,
                                    lst_xlsfl_in_wd01)]
# make an empty list to add read in tables to
ls.smplMST <- list()
# make a number that can increase for each iteration
j <-1
fi <- fls[1]
# iterate over files with concentrations on the extractions from filters
for (fi in fls)
{
  
  inf2 <- paste0(wd01,"/",fi)
  #print(inf2)}
  #read in 'csv' file to a tibble
  df_sMST <- readr::read_delim(inf2,delim=",")
  # make it a data frame
  df_sMST <- as.data.frame(df_sMST)
  # check if there is 'lokalitet_vanda' column in the data frame
  vandaclmprsent <- T %in% grepl("vanda",colnames(df_sMST))
  # if there is NOT, then 
  if (vandaclmprsent==F)
  {
    # add an empty column with this name
    df_sMST$lokalitet_vanda <- NA
  # end the if check for the 'lokalitet_vanda' column
  }
  # substitute in column names to standardize column names
  colnames(df_sMST)[grepl("Vandprvenummer_U",colnames(df_sMST))] <- "U_Pr_Nr"
  colnames(df_sMST)[grepl("Vandprvenummer_u",colnames(df_sMST))] <- "U_Pr_Nr"
  colnames(df_sMST)[grepl("Ansvarlig_",colnames(df_sMST))] <- "Navn_Inds"
  colnames(df_sMST)[grepl("Skib_",colnames(df_sMST))] <- "skib"
  colnames(df_sMST)[grepl("Institu",colnames(df_sMST))] <- "Inst"
  colnames(df_sMST)[grepl("^.ato_dato",colnames(df_sMST))] <- "Dato_inds"
  colnames(df_sMST)[grepl("^.okalitet_na",colnames(df_sMST))] <- "Lok_omr01"
  colnames(df_sMST)[grepl("^.osition_l.*pos$",colnames(df_sMST))] <- "lok_pos_lat"
  colnames(df_sMST)[grepl("^.osition_b.*pos$",colnames(df_sMST))] <- "lok_pos_lon"

  colnames(df_sMST)[grepl("^Max_vandd",colnames(df_sMST))] <- "max_Dyb"
  colnames(df_sMST)[grepl("^Dybde_hvor_v",colnames(df_sMST))] <- "Dyb_Str_m"
  colnames(df_sMST)[grepl("^Volumen_vand_f",colnames(df_sMST))] <- "Vwf_mL"
  colnames(df_sMST)[grepl("^Er_vands",colnames(df_sMST))] <- "vandsoejllagd"
  colnames(df_sMST)[grepl("^Dybde_for_la",colnames(df_sMST))] <- "Dyb_lagd01"
  colnames(df_sMST)[grepl("^.ato_hvor_vandprve_er_o",colnames(df_sMST))] <- "dato_80C"
  colnames(df_sMST)[grepl("^.emperatur_p",colnames(df_sMST))] <- "Temp_Inds"
  colnames(df_sMST)[grepl("^.emrkninger",colnames(df_sMST))] <- "Bemaerkn"
  
  # define columns to keep
  ctkee <- c("U_Pr_Nr",
  "Navn_Inds",
  "skib",
  "Inst",
  "Dato_inds",
  "Lok_omr01",
  "lokalitet_vanda",
  "lok_pos_lat",
  "lok_pos_lon",
  "max_Dyb",
  "Dyb_Str_m",
  
  "Vwf_mL",
  "vandsoejllagd",
  "Dyb_lagd01",
  "dato_80C",
  "Temp_Inds",
  "Bemaerkn")
  df_sMST <- df_sMST[ctkee]
  
  # add the data frame with samples to a list
  ls.smplMST[[j]] <- df_sMST
  
  #increase the count
  j <- j+1
  # end the iteration over files that are to become data frames

}  

# Use 'plyr::ldply' instead of 'datatable::rbindlist'
# as 'datatable::rbindlist' was hysterical about columns
# see this question: https://stackoverflow.com/questions/55706560/make-rbindlist-skip-ignore-or-change-class-attribute-of-the-column
#bind the rows in each list in to one data frame
require(plyr)
df_smplMST <- plyr::ldply(ls.smplMST, data.frame)
# get the library for working with excel files
library(openxlsx)
# replce MST sample numbers that only have 3 number characters
MSTnosw3c <- df_smplMST$U_Pr_Nr[grep("MST[0-9]{3}$",df_smplMST$U_Pr_Nr)]
MSTnosw3c <- gsub("MST","MST0",MSTnosw3c)
df_smplMST$U_Pr_Nr[grep("MST[0-9]{3}$",df_smplMST$U_Pr_Nr)] <- MSTnosw3c

# modify the samples that do not have a collection date
# for better or worse , use the date for transferring the sample to -80C
# then it at least has date for seasons and collection

df_smplMST$Dato_inds[is.na(df_smplMST$Dato_inds)] <- gsub(":","-",df_smplMST$dato_80C[is.na(df_smplMST$Dato_inds)])
# for now appears to be only sample 'MST2022040' that is missing the collection
# date. This part might need modification if more samples are without sampling date
df_smplMST[grepl("MST2022040", df_smplMST$U_Pr_Nr),]

#_______________________________________________________________________________

# make the collection coordinates numeric
df_smplMST$lok_pos_lat <- as.numeric(df_smplMST$lok_pos_lat)
df_smplMST$lok_pos_lon <- as.numeric(df_smplMST$lok_pos_lon)

#_______________________________________________________________________________
# section 02 - start - swap wrong longitude and wrong latitude around
#_______________________________________________________________________________

# Latitude is N-S positions
shifted_lat<- df_smplMST$lok_pos_lat[!(50<df_smplMST$lok_pos_lat & 
                                         df_smplMST$lok_pos_lat<60)]
# get row index numbers for the latitudes, that do not match Denmark
rwnshlat <- which(!(50<df_smplMST$lok_pos_lat & 
                      df_smplMST$lok_pos_lat<60))
# Longitude is E-W positions
shifted_lon <- df_smplMST$lok_pos_lon[!(4<df_smplMST$lok_pos_lon & 
                                          df_smplMST$lok_pos_lon<18)]
# get row index numbers for the longitudes, that do not match Denmark
rwnshlon <- which(!(4<df_smplMST$lok_pos_lon & 
                      df_smplMST$lok_pos_lon<18))
# find the row index numbers for the longitudes that are odd for both
# latitude and longitude
oddll <- setdiff(rwnshlon,rwnshlat)
# use these index row numbers to check the contents of odd  
df_smplMST.oddll <- df_smplMST[c(oddll),]
# compare two vectors and find the common index row numbers
#https://stackoverflow.com/questions/3695677/how-to-find-common-elements-from-multiple-vectors
cidx <- Reduce(intersect, list(rwnshlon,rwnshlat))
# for the index numbers that are common in having a wrong latitude outside 50-60
# and the longitude that are outside the 4-18, the index numbers in common,
# can be used to swap around the latitude and the longitude
wlat <- df_smplMST[cidx,]$lok_pos_lat
wlon <- df_smplMST[cidx,]$lok_pos_lon
# put the wrong lon into the latitude, and the wrong lat in to the longitiude
df_smplMST[cidx,]$lok_pos_lat <- wlon
df_smplMST[cidx,]$lok_pos_lon <- wlat

#check the positions in Sweden.
# Longitude is E-W positions
Swed_lon_lat <- df_smplMST[(13<df_smplMST$lok_pos_lon & 
                                          56<df_smplMST$lok_pos_lat),]
# get index row numbers for the Swedish wrong  positions
wswps <- which(13<df_smplMST$lok_pos_lon & 
              56<df_smplMST$lok_pos_lat)
# get the locations names for the odd swedish samples
lkwswps <- df_smplMST$Lok_omr01[wswps]
# substitute to get the vanda location number
lvwswps<- gsub("(^[0-9]{+})_(.*)","\\1",lkwswps)
# match the vanda location numbers to other samples in the data frame
# to get longitude and latitude
swdlat <- df_smplMST$lok_pos_lat[match(lvwswps,df_smplMST$lokalitet_vanda)]
swdlon <- df_smplMST$lok_pos_lon[match(lvwswps,df_smplMST$lokalitet_vanda)]
#use the lat and long for the corresponding vanda number locations
df_smplMST$lok_pos_lat[wswps] <- swdlat
df_smplMST$lok_pos_lon[wswps] <- swdlon

# define boundaries for the wrong latitudes to look for in Jylland
uplatJW <- 56.6
lolatJW <- 55.5
uplonJW <- 9.5
lolonJW <- 8.8
  
#check the positions in Jylland
# Longitude is E-W positions
Jyll_wrong01 <- df_smplMST[(lolatJW<df_smplMST$lok_pos_lat & 
                              uplatJW>df_smplMST$lok_pos_lat &
                              lolonJW<df_smplMST$lok_pos_lon &
                              uplonJW>df_smplMST$lok_pos_lon),]

# define boundaries for the wrong latitudes to look for in Jylland
uplatJW2 <- 56.6
lolatJW2 <- 56.0
uplonJW2 <- 9.8
lolonJW2 <- 8.8


Jyll_wrong02 <- df_smplMST[(lolatJW2<df_smplMST$lok_pos_lat & 
                              uplatJW2>df_smplMST$lok_pos_lat &
                              lolonJW2<df_smplMST$lok_pos_lon &
                              uplonJW2>df_smplMST$lok_pos_lon),]
# bind together wrong locations in Jylland
Jyll_wrong03 <- rbind(Jyll_wrong02,Jyll_wrong01)

indxnJW <- match(Jyll_wrong03$U_Pr_Nr,df_smplMST$U_Pr_Nr)
JWlv <- df_smplMST[indxnJW,"lokalitet_vanda"]
JWlk <- df_smplMST[indxnJW,"Lok_omr01"]
JWlk <- gsub("(^.*)_.*_.*","\\1",JWlk)
# make a subsetted version of the MST sample data frame, but without the 
# wrong Jylland positions
df_smplMST.nJW <- df_smplMST[!(rownames(df_smplMST)  %in% indxnJW),]
# use the vanda numbers for the wrong Jylland inland positions, to 
# get other vanda numbers and their latitudes and longitudes
JW_latrpl <- df_smplMST.nJW$lok_pos_lat[match(Jyll_wrong03$lokalitet_vanda,df_smplMST.nJW$lokalitet_vanda)]
JW_lonrpl <- df_smplMST.nJW$lok_pos_lon[match(Jyll_wrong03$lokalitet_vanda,df_smplMST.nJW$lokalitet_vanda)]
#use the lat and long for the corresponding vanda number locations
df_smplMST$lok_pos_lat[indxnJW] <- JW_latrpl
df_smplMST$lok_pos_lon[indxnJW] <- JW_lonrpl

# check for more longitude positions, that could be latitudes
idxfwlon <- which(df_smplMST$lok_pos_lon>50)
# make a subsetted version of the MST sample data frame, but without the 
# wrong longitude
df_smplMST.nlon <- df_smplMST[!(rownames(df_smplMST)  %in% idxfwlon),]
# use the vanda numbers for the wrong longitude positions, to 
# get other vanda numbers and their latitudes and longitudes

W_latrpl <- df_smplMST.nlon$lok_pos_lat[match(df_smplMST[idxfwlon,"lokalitet_vanda"],df_smplMST.nJW$lokalitet_vanda)]
W_lonrpl <- df_smplMST.nlon$lok_pos_lon[match(df_smplMST[idxfwlon,"lokalitet_vanda"],df_smplMST.nJW$lokalitet_vanda)]
#use the lat and long for the corresponding vanda number locations
df_smplMST$lok_pos_lat[indxnJW] <- W_latrpl
df_smplMST$lok_pos_lon[indxnJW] <- W_lonrpl

#_______________________________________________________________________________
# section 02 - end - swap wrong longitude and wrong latitude around
#_______________________________________________________________________________
# First (before making maps)  make sure all the required packages are loaded
# https://uchicagoconsulting.wordpress.com/tag/r-ggplot2-maps-visualization/
#install packages needed
# if(!require(maps)){
#   install.packages("maps")
# }
library(maps)
# if(!require(ggplot2)){
#   install.packages("ggplot2")
# }
library(ggplot2)
library(maps)
# # #https://www.r-spatial.org/r/2018/10/25/ggplot2-sf-2.html
# To get rgdal and googleway to work,
#first run these in a terminal:

# $ sudo apt install netcdf-*
# $   sudo apt install libnetcdf-dev
# $ sudo apt install libjq-dev
# $ sudo apt install gdal-bin libgdal-dev libproj-dev
# $ sudo apt install libudunits2-dev
# if(!require(cowplot)){
#   install.packages("cowplot")
#   library(cowplot)
# }
library(cowplot)
# if(!require(googleway)){
#   install.packages("googleway")
#   library(googleway)
# }
library(googleway)
# if(!require(ggrepel)){
#   install.packages("ggrepel")
#   }
library(ggrepel)
# if(!require(ggspatial)){
#   install.packages("ggspatial")
# }
library(ggspatial)
# if(!require(libwgeom)){
#   install.packages("libwgeom")
#   library(libwgeom)
# }
# if(!require(sf)){
#   install.packages("sf")
#   library(sf)
# }
library(sf)
# if(!require(rnaturalearth)){
#   install.packages("rnaturalearth")
# }
library(rnaturalearth)
# if(!require(rnaturalearthdata)){
#   install.packages("rnaturalearthdata")
#   library(rnaturalearthdata)
# }

library(rnaturalearthdata)
# if(!require(ggforce)){
#   install.packages("ggforce")
# }
#get 'rnaturalearthhires' installed
# if(!require(rnaturalearthhires)){
#   #install.packages("rnaturalearthhires")
#   install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
#   library(rnaturalearthhires)
# }
# # 
library(rnaturalearthhires)
library(rnaturalearthdata)
library(rnaturalearth)
library("ggplot2")
library("sf")
library(ggforce)


theme_set(theme_bw())

#install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
# # 
library("rnaturalearth")
library("rnaturalearthdata")
library("rnaturalearthhires")
# # Get a map, use a high number for 'scale' for a coarse resolution
# use a low number for scale for a high resolution
# if the map 'world' does not exist, then download it
world <- ne_countries(scale = 10, returnclass = "sf")
library(ggplot2)

# subset data frame to only comprise locations that have 'lokalitet_vanda'
# number
df_smplMST.lkv <- df_smplMST[!is.na(df_smplMST$lokalitet_vanda),]
#
df_smplMST.lkv$dec_lon<- df_smplMST.lkv$lok_pos_lon
df_smplMST.lkv$dec_lat<- df_smplMST.lkv$lok_pos_lat
# add a level for jittering points
jitlvl <- 0.03
#make plot

#make plot
p05 <- ggplot(data = world) +
  geom_sf(color = "black", fill = "azure3", lwd=0.4) +
  #https://ggplot2.tidyverse.org/reference/position_jitter.html
  #https://stackoverflow.com/questions/15706281/controlling-the-order-of-points-in-ggplot2
  # use 'geom_jitter' instead of 'geom_point' 
  geom_point(data = Swed_lon_lat,
             aes(x =lok_pos_lon, 
                 y = lok_pos_lat, 
                 label=lokalitet_vanda),
             color=alpha(c("black"),c(0.8)),
             size=2.6,
             fill=alpha(c("yellow"),c(0.8)),
             shape=22,
             width = jitlvl, #0.07, jitter width 
             height = jitlvl) + #, #0.07, # jitter height
  geom_point(data = Jyll_wrong03,
             aes(x =lok_pos_lon, 
                 y = lok_pos_lat, 
                 label=lokalitet_vanda),
             color=alpha(c("black"),c(0.8)),
             size=2.6,
             fill=alpha(c("cyan"),c(0.8)),
             shape=22) +
  geom_text_repel(data=df_smplMST.lkv,
                  aes(x =lok_pos_lon, 
                      y = lok_pos_lat,  
                      label=lokalitet_vanda),
                  color="red") +
  geom_hline(yintercept=uplatJW, linetype="dashed", color="red") +
  geom_hline(yintercept=lolatJW, linetype="dashed", color="red") +
  geom_vline(xintercept=uplonJW, linetype="dashed", color="red") +
  geom_vline(xintercept=lolonJW, linetype="dashed", color="red") +
  
  geom_hline(yintercept=uplatJW2, linetype="dashed", color="blue") +
  geom_hline(yintercept=lolatJW2, linetype="dashed", color="blue") +
  geom_vline(xintercept=uplonJW2, linetype="dashed", color="blue") +
  geom_vline(xintercept=lolonJW2, linetype="dashed", color="blue") +
  geom_text_repel(data=Swed_lon_lat,
                  aes(x =lok_pos_lon, 
                      y = lok_pos_lat,  
                      label=Lok_omr01),
                  color="orchid4") +
  #define limits of the plot 
  ggplot2::coord_sf(xlim = c(6.6, 17.2),
                    ylim = c(54.2, 58.4), 
                    expand = FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#change axis labels

p05t <- p05 + xlab("longitude") + ylab("latitude")
#change the header for the legend on the side, 
#this must be done for both 'fill', 'color' and 'shape', to avoid 
#getting separate legends
p05t <- p05t + labs(color='')
p05t <- p05t + labs(fill='')
p05t <- p05t + labs(shape='')
p05t <- p05t + labs(size='')
# https://github.com/tidyverse/ggplot2/issues/3492
#adjust tick marks on axis
p05t <- p05t + scale_y_continuous(breaks=seq(54.0, 58.4,0.5))
p05t <- p05t + scale_x_continuous(breaks=seq(6.0, 17.2,1.0))
# see the plot
p05t

#_______________________________________________________________________________
# subset the data frame to only comprise the locations that have a 
# 'lokalitet_vanda' number
df_smplMST.lkv <- df_smplMST[!is.na(df_smplMST$lokalitet_vanda),]
# Many of the harbour location names for the collection localities
# have different names assigned when the sample was collected.
# use 'dplyr::distinct' to only have the distinct rows for 
# column 'Lok_omr01' - to be able to see where they are located

df_Lk <-  df_smplMST.lkv %>% 
  dplyr::distinct(Lok_omr01, .keep_all = TRUE)
# add 'lokalitet_vanda' for the 'Lok_omr01' that can be matched 
lkmsVN <- df_smplMST$Lok_omr01[is.na(df_smplMST$lokalitet_vanda)]
print("no 1")
length(lkmsVN)
lkmsVN <- df_Lk$lokalitet_vanda[match(lkmsVN,df_Lk$Lok_omr01)]
df_smplMST$lokalitet_vanda[is.na(df_smplMST$lokalitet_vanda)] <- lkmsVN
# add 'lokalitet_vanda' numbers for the rows where the 'Lok_omr01'
# holds the 'lokalitet_vanda' number
df_smplMST$lokalitet_vanda[grepl("^[0-9]{8}$",df_smplMST$Lok_omr01)] <- df_smplMST$Lok_omr01[grepl("^[0-9]{8}$",df_smplMST$Lok_omr01)]
# check again for the locations missing a 'lokalitet_vanda' number
lkmsVN <- df_smplMST$Lok_omr01[is.na(df_smplMST$lokalitet_vanda)]
# substitute in the 'Lok_omr01' where the 1st part of the name
# has the 'lokalitet_vanda' number
lkmsVN <- gsub("([0-9]{8})_.*","\\1",lkmsVN)
# if it only has the 'lokalitet_vanda' number, then use this,
# otherwise return NA
lkmsVN <- ifelse(grepl("^[0-9]{8}$",lkmsVN),lkmsVN,NA)
# add back the 'lokalitet_vanda' number
df_smplMST$lokalitet_vanda[is.na(df_smplMST$lokalitet_vanda)] <- lkmsVN
# check which 'Lok_omr01' still are missing a 'lokalitet_vanda' number
lkmsVN <- df_smplMST$Lok_omr01[is.na(df_smplMST$lokalitet_vanda)]
# substitute to get the missing 'lokalitet_vanda' number
lkmsVN <- gsub(".*_([0-9]{8})$","\\1",lkmsVN)
lkmsVN <- gsub(".*_([0-9]{8})_.*$","\\1",lkmsVN)
lkmsVN <- ifelse(grepl("^[0-9]{8}$",lkmsVN),lkmsVN,NA)
# add back the 'lokalitet_vanda' number
df_smplMST$lokalitet_vanda[is.na(df_smplMST$lokalitet_vanda)] <- lkmsVN




df_smplMST$lokalitet_vanda[grepl("arrebaeksminde",df_smplMST$Lok_omr01)] <- "95220322"
df_smplMST$lokalitet_vanda[grepl("DMU444",df_smplMST$Lok_omr01)] <- "99000052"
df_smplMST$lokalitet_vanda[grepl("DMU_444",df_smplMST$Lok_omr01)] <- "99000052"
df_smplMST$lokalitet_vanda[grepl("Hovedeskov",df_smplMST$Lok_omr01)] <- "99210007"
df_smplMST$lokalitet_vanda[grepl("Hovedskov",df_smplMST$Lok_omr01)] <- "99210007"
df_smplMST$lokalitet_vanda[grepl("Ringkoebing",df_smplMST$Lok_omr01)] <- "91320001"
df_smplMST$lokalitet_vanda[grepl("RKB1",df_smplMST$Lok_omr01)] <- "91320001"
df_smplMST$lokalitet_vanda[grepl("NOR7715",df_smplMST$Lok_omr01)] <- "92200001"
df_smplMST$lokalitet_vanda[grepl("Laes",df_smplMST$Lok_omr01)] <- "93010001"
df_smplMST$lokalitet_vanda[grepl("NOR5503",df_smplMST$Lok_omr01)] <- "93610032"
df_smplMST$lokalitet_vanda[grepl("attegataabnedel",df_smplMST$Lok_omr01)] <- "93000001"
df_smplMST$lokalitet_vanda[grepl("FYN6900017",df_smplMST$Lok_omr01)] <- "94230001"
df_smplMST$lokalitet_vanda[grepl("Lister",df_smplMST$Lok_omr01)] <- "91650035"
df_smplMST$lokalitet_vanda[grepl("RIB1510007",df_smplMST$Lok_omr01)] <- "91500007"
df_smplMST$lokalitet_vanda[grepl("Soenderborg",df_smplMST$Lok_omr01)] <- "95700008"
df_smplMST$lokalitet_vanda[grepl("VejleFjord",df_smplMST$Lok_omr01)] <- "95130028"
df_smplMST$lokalitet_vanda[grepl("VIB3708",df_smplMST$Lok_omr01)] <- "93730002"
df_smplMST$lokalitet_vanda[grepl("VIB3727",df_smplMST$Lok_omr01)] <- "93740007"
df_smplMST$lokalitet_vanda[grepl("VIB_3727",df_smplMST$Lok_omr01)] <- "93740007"
df_smplMST$lokalitet_vanda[grepl("ARH170006",df_smplMST$Lok_omr01)] <- "94400007"
df_smplMST$lokalitet_vanda[grepl("HB01",df_smplMST$Lok_omr01)] <- "93500119"
df_smplMST$lokalitet_vanda[grepl("evring",df_smplMST$Lok_omr01)] <- "93500119"
df_smplMST$lokalitet_vanda[grepl("NOR409",df_smplMST$Lok_omr01)] <- "93000005"
df_smplMST$lokalitet_vanda[grepl("KBH431",df_smplMST$Lok_omr01)] <- "97200002"
df_smplMST$lokalitet_vanda[grepl("Franks",df_smplMST$Lok_omr01)] <- "99150002"
df_smplMST$lokalitet_vanda[grepl("Loegstoer",df_smplMST$Lok_omr01)] <- "93730002"

df_smplMST$lokalitet_vanda[grepl("ragskov",df_smplMST$Lok_omr01)] <- "93920064"
df_smplMST$lokalitet_vanda[grepl("Mariager",df_smplMST$Lok_omr01)] <- "93610032"
df_smplMST$lokalitet_vanda[grepl("0101142",df_smplMST$Lok_omr01)] <- "96220363"
df_smplMST$lokalitet_vanda[grepl("Roskilde",df_smplMST$Lok_omr01)] <- "93220011"
df_smplMST$lokalitet_vanda[grepl("Skive",df_smplMST$Lok_omr01)] <- "93740007"
df_smplMST$lokalitet_vanda[grepl("alborg",df_smplMST$Lok_omr01)] <- "93000005"
df_smplMST$lokalitet_vanda[grepl("Bornholm",df_smplMST$Lok_omr01)] <- "99000049"
df_smplMST$lokalitet_vanda[grepl("Fornaes",df_smplMST$Lok_omr01)] <- "93420310"
df_smplMST$lokalitet_vanda[grepl("ROS",df_smplMST$Lok_omr01)] <- "93220011"
df_smplMST$lokalitet_vanda[grepl("FYN6300043",df_smplMST$Lok_omr01)] <- "95600002"
df_smplMST$lokalitet_vanda[grepl("VEJ0006870",df_smplMST$Lok_omr01)] <- "94300001"
df_smplMST$lokalitet_vanda[grepl("FFSB0011",df_smplMST$Lok_omr01)] <- "95700008" 
df_smplMST$lokalitet_vanda[grepl("20925",df_smplMST$Lok_omr01)]  <- "93000001"
df_smplMST$lokalitet_vanda[grepl("6300043",df_smplMST$Lok_omr01)] <- "95600002"
df_smplMST$lokalitet_vanda[grepl("DSoe",df_smplMST$Lok_omr01)]  <- "96510125"
df_smplMST$lokalitet_vanda[grepl("oeresund",df_smplMST$Lok_omr01)] <- "97200002"
df_smplMST$lokalitet_vanda[grepl("006870",df_smplMST$Lok_omr01)] <-   "94300001"
df_smplMST$lokalitet_vanda[grepl("Gniben",df_smplMST$Lok_omr01)] <-   "94300001"
df_smplMST$lokalitet_vanda[grepl("Sundet",df_smplMST$Lok_omr01)]  <-  "97200002"
df_smplMST$lokalitet_vanda[grepl("Odense",df_smplMST$Lok_omr01)] <-   "94230001"
df_smplMST$lokalitet_vanda[grepl("Soenderhooest",df_smplMST$Lok_omr01)] <- "91500007"
df_smplMST$lokalitet_vanda[grepl("Flensborg",df_smplMST$Lok_omr01)] <- "95700008"
df_smplMST$lokalitet_vanda[grepl("Boderne",df_smplMST$Lok_omr01)]   <- "99150002"
df_smplMST$lokalitet_vanda[grepl("D0027",df_smplMST$Lok_omr01)]  <- "91650035"
df_smplMST$lokalitet_vanda[grepl("ydfynske",df_smplMST$Lok_omr01)] <- "96510125"


ulkv <- unique(df_smplMST$Lok_omr01[is.na(df_smplMST$lokalitet_vanda)])
print(ulkv)

lst_lov <- list()
j<- 1
for (Ml in ulkv)
{  print(Ml)
  lov<- df_smplMST[grepl(  Ml,df_smplMST$Lok_omr01),]
  df_lov <- lov[c("U_Pr_Nr",
                  "Dato_inds",
                  "Lok_omr01",
                  "lok_pos_lon",
                  "lok_pos_lat")]
  lst_lov[[j]] <- df_lov
  j <- j+1
}
# Use 'plyr::ldply' instead of 'datatable::rbindlist'
# as 'datatable::rbindlist' was hysterical about columns
# see this question: https://stackoverflow.com/questions/55706560/make-rbindlist-skip-ignore-or-change-class-attribute-of-the-column
#bind the rows in each list in to one data frame
df_lov01 <- plyr::ldply(lst_lov, data.frame)
#
write.table(df_lov01,file=paste0(wd00_wd03,
                                 "/table01b_unmatched_vanda_lokaliteter.csv"),
                                 sep=";",
            row.names = F)
# read in table with the vanda numbers that are not appearing in the 
# 2023 samples
df_missing_vlok <- readxl::read_excel(paste0(wd01,"/table01c_unmatched_vanda_lokaliteter.xlsx"))
df_missing_vlok <- as.data.frame(df_missing_vlok)
# only retain unique vanda numbers per location
df_missing_vlok <-  df_missing_vlok %>% 
  dplyr::distinct(Lok_omr01, VandaNr, .keep_all = TRUE)
# get index row numbers for missing vanda numbers
indxfm <- which(is.na(df_smplMST$lokalitet_vanda))
# add the missing vanda number for the locations index row missing this number
df_smplMST[indxfm,"lokalitet_vanda"] <- df_missing_vlok$VandaNr[match(df_smplMST[indxfm,"Lok_omr01"],df_missing_vlok$Lok_omr01)]


# Many of the harbour location names for the collection localities
# have different names assigned when the sample was collected.
# But at least their lat long coordinates should roughly match
# use 'dplyr::distinct' to only have the distinct rows for 
# column 'Harbour' - to be able to see where they are located

df_smplMST.dstc.1 <-  df_smplMST %>% 
  dplyr::distinct(Lok_omr01, .keep_all = TRUE)

df_smplMST.dstc.1 <-  df_smplMST %>% 
  dplyr::distinct(lok_pos_lat, lok_pos_lon, .keep_all = TRUE)
#df_smplMST.dstc.1 <- df_smplMST[is.na(df_smplMST$lokalitet_vanda),]
#nrow(df_smplMST.dstc.1)
df_smplMST.dstc.2 <- df_smplMST[is.na(df_smplMST$lokalitet_vanda),]
# df_smplMST.dstc.2 <-  df_smplMST %>% 
#   dplyr::distinct(lokalitet_vanda, .keep_all = TRUE)

df_smplMST.dstc.2 <-  df_smplMST %>% 
  dplyr::distinct(Lok_omr01, .keep_all = TRUE)

df_smplMST.dstc.2 <- df_smplMST[is.na(df_smplMST$lokalitet_vanda),]

df_smplMST.dstc3<-  rbind(df_smplMST.dstc.1,df_smplMST.dstc.2)
nrow(df_smplMST.dstc3)
#View(df_smplMST.dstc3)
df_smplMST.4 <- df_smplMST
df_smplMST.4$lokalitet_vanda[which(duplicated(df_smplMST.4$lokalitet_vanda))] <- ""


#View(df_smplMST.dstc3)
# add a level for jittering points
jitlvl <- 0.03
#make plot
#make plot
p05 <- ggplot(data = world) +
  geom_sf(color = "black", fill = "azure3", lwd=0.4) +
  #https://ggplot2.tidyverse.org/reference/position_jitter.html
  #https://stackoverflow.com/questions/15706281/controlling-the-order-of-points-in-ggplot2
  # use 'geom_jitter' instead of 'geom_point' 
  geom_point(data = df_smplMST.dstc.1 ,
             aes(x =lok_pos_lon, y = lok_pos_lat),
             color=alpha(c("black"),c(0.8)),
             size=2.6,
             fill=alpha(c("yellow"),c(0.8)),
             shape=22,
             width = jitlvl, #0.07, jitter width 
             height = jitlvl) + #, #0.07, # jitter height
  geom_text_repel(data=df_smplMST.dstc.1 ,
                  aes(x=lok_pos_lon, 
                      y = lok_pos_lat, 
                      label=lokalitet_vanda),
                  color="brown", label.size=2.1,
                  max.overlaps=30, size=1.2) +
  geom_text_repel(data=df_smplMST.dstc.2,
                  aes(x=lok_pos_lon, 
                      y = lok_pos_lat, 
                      label=Lok_omr01 ),
                  color="darkgreen", label.size=1.1,
                  max.overlaps=10, size=1.6) +
  #define limits of the plot 
  ggplot2::coord_sf(xlim = c(6.6, 17.2),
                    ylim = c(54.2, 58.4), 
                    expand = FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#change axis labels

p05t <- p05 + xlab("longitude") + ylab("latitude")
#change the header for the legend on the side, 
#this must be done for both 'fill', 'color' and 'shape', to avoid 
#getting separate legends
p05t <- p05t + labs(color='')
p05t <- p05t + labs(fill='')
p05t <- p05t + labs(shape='')
p05t <- p05t + labs(size='')
# https://github.com/tidyverse/ggplot2/issues/3492
#adjust tick marks on axis
p05t <- p05t + scale_y_continuous(breaks=seq(54.0, 58.4,0.5))
p05t <- p05t + scale_x_continuous(breaks=seq(6.0, 17.2,1.0))
# see the plot
p05t
#set variable to define if figures are to be saved
bSaveFigures<-T
if(bSaveFigures==T){
  ggsave(plot = p05t, 
         filename = paste0(wd00_wd03,"/Fig03_map_of_vanda_lok_no.png"),
         #width=210,height=297,
         #width=297,height=210,
         width=297,height=1.2*210,
         #width=1.4*297,height=210,
         units="mm",dpi=300)
}

#_______________________________________________________________________________
# replace MST number because of comment
comment2021_027 <- df_smplMST$Bemaerkn[grepl("MST2021027",df_smplMST$U_Pr_Nr)]
df_smplMST$U_Pr_Nr[grepl("ErangivetMST2021_029",df_smplMST$Bemaerkn)] <- "MST2021029"
# replace MST number because of comment
comment2021_008 <- df_smplMST$Bemaerkn[grepl("MST2021008",df_smplMST$U_Pr_Nr)]
df_smplMST$U_Pr_Nr[grepl("ErangivetMST2021_009",df_smplMST$Bemaerkn)] <- "MST2021009"

# # make vectors with location names that are missing coordinates
# # and another location name that comes from the same area
loctrpl<- c("lok.wc","lok.mc",
            "97230016_Nivaa","Nivaaoeresund",
            "93220011_322Roskildefjord","ROS60_322RoskildeFjord",
            "FYN6900017","94230001_FYN6900017")
# make a data frame of locations to replace
df_lokNmstrpl <- as.data.frame(matrix(loctrpl,ncol=2, byrow = T))
colnames(df_lokNmstrpl) <- df_lokNmstrpl[1,]
df_lokNmstrpl <- df_lokNmstrpl[-1,]
# iterate over locations missing coordinates
for (lok.mc in df_lokNmstrpl$lok.mc)
{
  #match to get the wrong coordinate location names
  lok.wc <- df_lokNmstrpl$lok.wc[match(lok.mc,df_lokNmstrpl$lok.mc)]
  replc.lat <- df_smplMST$lok_pos_lat[(lok.mc==df_smplMST$Lok_omr01)]
  replc.lon <- df_smplMST$lok_pos_lon[(lok.mc==df_smplMST$Lok_omr01)]
  replc.lat <- unique(replc.lat)[1]
  replc.lon <- unique(replc.lon)[1]
  df_smplMST$lok_pos_lat[(lok.wc==df_smplMST$Lok_omr01)] <- replc.lat
  df_smplMST$lok_pos_lon[(lok.wc==df_smplMST$Lok_omr01)] <- replc.lon
}
# replace dates that are in incorrect format
odddat <- df_smplMST$Dato_inds[grep("/",df_smplMST$Dato_inds)]
odddat <- gsub("^([0-9]+)/([0-9]+)/([0-9]+)","\\3-\\2-\\1",odddat)
df_smplMST$Dato_inds[grep("/",df_smplMST$Dato_inds)] <- odddat



# remove duplicated rows
df_smplMST <- df_smplMST[!duplicated(df_smplMST),]
# add MST letter code to numbers missing MST in front
MSTlettmiss <- df_smplMST$U_Pr_Nr[grep("^[0-9]",as.character(df_smplMST$U_Pr_Nr))]
df_smplMST$U_Pr_Nr[grep("^[0-9]",as.character(df_smplMST$U_Pr_Nr))] <- paste0("MST",MSTlettmiss)
# 
df_smplMST.miss.date <- df_smplMST[is.na(df_smplMST$Dato_inds),]
#
df_smplMST.unk.MSTno <-  df_smplMST[(df_smplMST$U_Pr_Nr=="MST"),]
df_smplMST.miss.lat <- df_smplMST[!(df_smplMST$lok_pos_lat>51 &df_smplMST$lok_pos_lat<60),]
df_smplMST.miss.lon <- df_smplMST[!(df_smplMST$lok_pos_lon>3 &df_smplMST$lok_pos_lon<19),]
df_smplMST$Temp_Inds <- as.numeric(df_smplMST$Temp_Inds)
df_smplMST.notNA <-  df_smplMST[!is.na(df_smplMST$Temp_Inds),]
df_smplMST.miss.indstemp <- df_smplMST.notNA[(df_smplMST.notNA$Temp_Inds>27),]
df_smplMST.NA.indstemp <- df_smplMST[is.na(df_smplMST$Temp_Inds),]
df_smplMST.missinfo <- rbind(df_smplMST.miss.date,
                             df_smplMST.unk.MSTno,
                              df_smplMST.miss.lat,
                              df_smplMST.miss.lon,
                              df_smplMST.miss.indstemp,
                              df_smplMST.NA.indstemp)
# remove duplicated rows
df_smplMST.missinfo <- df_smplMST.missinfo[!duplicated(df_smplMST.missinfo),]


# find the MST numbers that appear more than once
# i.e. the MST sample numbers that have been used more than once in the field
dupl_MSTnos <- df_smplMST$U_Pr_Nr[duplicated(df_smplMST$U_Pr_Nr)]
# check of the vector (the vector with duplicated MST numbers) is in the rows with MST numbers
#https://stackoverflow.com/questions/11612235/select-rows-from-a-data-frame-based-on-values-in-a-vector
# and only retain rows where the MST numbers in the vector is also present
df_smplMST.dpl <- df_smplMST[df_smplMST$U_Pr_Nr %in% dupl_MSTnos,]
# sort the data fram using 'order' to get the rows arranged by the MST no column 
df_smplMST.dpl <- df_smplMST.dpl[order(df_smplMST.dpl$U_Pr_Nr),]

# paste together MST sample number and collection date
df_smplMST.dpl$MSTn.di <- paste(df_smplMST.dpl$U_Pr_Nr,
                                df_smplMST.dpl$Dato_inds,
                                sep="::")
# evaluate if the MSTnumber plus the sampling date is present in the
# data frame with concentrations for extractions
df_cex01$wdt <- (df_cex01$MSTno_di  %in%  df_smplMST.dpl$MSTn.di)
# limit to MStnumbers and sample dates that can be matched
df_cex01T <- df_cex01[(df_cex01$wdt==T),]
df_cex01F <- df_cex01[(df_cex01$wdt==F),]
#match(df_smplMST.dplT$MSTn.di,df_cex01$MSTno_di)


# copy the column to have columns with the same name
df_smplMST.dpl$MSTno_di <- df_smplMST.dpl$MSTn.di
# use dplyr::left_join to match between data frames 
dfg01 <- dplyr::left_join(df_cex01T,
                          df_smplMST.dpl %>% dplyr::select(
                          everything()
                        ),
                        by = "MSTno_di")


# make a common column for the coming left_join
df_smplMST$U_Pr_Nr <- gsub("MSR","",df_smplMST$U_Pr_Nr)
df_smplMST$MST.nummer <- df_smplMST$U_Pr_Nr
# use dplyr::left_join to match between data frames 
# for the MST numbers that could not be matched with both MSTnumber and date

dfg02 <- dplyr::left_join(df_cex01F,
                          df_smplMST ,
                          by = "MST.nummer")

# dfg02 <- dplyr::left_join(df_cex01F,
#                           df_smplMST %>% dplyr::select(
#                             everything()
#                           ),
#                           by = "MST.nummer")

# grepl for the columns that do not have '.y' in the column name
# as these are duplicates, and are not needed
cnmtk1 <- colnames(dfg01)[!grepl("\\.y",colnames(dfg01))]
cnmtk2 <- colnames(dfg02)[!grepl("\\.y",colnames(dfg02))]
# only keep the columns listed in the line above 
dfg01 <- dfg01[cnmtk1]
dfg02 <- dfg02[cnmtk2]
# use gsub to remove the '.x' added ad the end, as this is not needed
colnames(dfg01) <- gsub("(.*)\\.x$","\\1",colnames(dfg01))
colnames(dfg02) <- gsub("(.*)\\.x$","\\1",colnames(dfg02))

# Use 'plyr::ldply' instead of 'datatable::rbindlist'
# as 'datatable::rbindlist' was hysterical about columns
# see this question: https://stackoverflow.com/questions/55706560/make-rbindlist-skip-ignore-or-change-class-attribute-of-the-column
#bind the rows in each list in to one data frame
require(plyr)
dfg03 <- plyr::ldply(list(dfg01,
                          dfg02), data.frame)
# re oder the data frame
dfg03 <- dfg03[order(dfg03$MST.nummer),]

# make the column characters to be able to replace with the dat_for_when_the_sample_was_collected_
# from the other data frame
dfg03$Dato_inds<- as.character(dfg03$Dato_inds)
dfg03$MSTno.sbExtr <- paste(dfg03$U_Pr_Nr,dfg03$Sub_Nmb_for_extraction,sep=".") 

# find the MST numbers that appear more than once
# i.e. the MST sample numbers that have been used more than once in the field
dupl_dfg03 <- dfg03$U_Pr_Nr[duplicated(dfg03$U_Pr_Nr)]
# check of the vector (the vector with duplicated MST numbers) is in the rows with MST numbers
#https://stackoverflow.com/questions/11612235/select-rows-from-a-data-frame-based-on-values-in-a-vector
# and only retain rows where the MST numbers in the vector is also present
dfg03.dpl <- dfg03[dfg03$U_Pr_Nr %in% dupl_dfg03,]
# sort the data fram using 'order' to get the rows arranged by the MST no column 
dfg03.dpl <- dfg03.dpl[order(dfg03.dpl$U_Pr_Nr),]

# exclude the samples where the "Sub_Nmb_for_extraction" is not zero
# and the "U_Pr_Nr" is not duplicated, and the 'Noter.på.plastic.posen'
# is not empty
dfg03 <- dfg03[!(dfg03$Sub_Nmb_for_extraction!="0" &
         duplicated(dfg03$U_Pr_Nr) &
           is.na(dfg03$Noter.på.plastic.posen)),]

# get the vanda numbers associated with the data frame with all samples
# among the files in the working directory, grep for the file
# that is from 2023, and that is 'Pøvetagningsskema'
# as this file contains the vanda numbers for collection locations
lstfwd01 <- list.files(wd01)
lstfwd01 <- lstfwd01[grepl("2023",lstfwd01)]
lstfwd01 <- lstfwd01[grepl("tagningsskema",lstfwd01)]
Prtgnsk2023 <- paste0(wd01,"/", lstfwd01)
# read in the excel xlsx file, 
df_vx01 <- readxl::read_xlsx(Prtgnsk2023)
# and make it a data frame instead of a tibble
df_vx01 <- as.data.frame(df_vx01)
# find the row number that has 'Obs' in it
obsrw <- which(grepl("Obs",df_vx01[,1]))
# and then get the data frame without this row
df_vx01 <- df_vx01[-obsrw,]
# find the row number that has 'Format' in it
formrw <- which(grepl("Format",df_vx01[,1]))
# and then get the data frame without this row
df_vx01 <- df_vx01[-formrw,]
# substitute in the column names
colnames(df_vx01) <- gsub("\r\n","",colnames(df_vx01))
# and replace the specific column name with 'Vandpr.venummer' 
colnames(df_vx01)[grepl("Vandpr.venummer",colnames(df_vx01))] <- "U_Pr_Nr"
colnames(df_vx01)[grepl("^.ato.*ato_inds",colnames(df_vx01))] <- "Dato_inds"
# substitute in the MST number
df_vx01$MSTno <- gsub(" ","",df_vx01$U_Pr_Nr)
df_vx01$MSTno <- gsub("-","",df_vx01$MSTno)
# get dates for the sampling dates read in when the excel file was read in
df_vx01$Dinds <- as.Date(as.numeric(df_vx01$Dato_inds), origin = "1899-12-30")
# paste the MST number and the collection date together
df_vx01$MSTno_Dinds<- paste0(df_vx01$MSTno,"::",df_vx01$Dinds)

# #View(dfg01)
# 
# min(unique(dfg01$lok_pos_lat))
# max(unique(dfg01$lok_pos_lat))
# min(unique(dfg01$lok_pos_lon))
# max(unique(dfg01$lok_pos_lon))
# 
# min(unique(dfg02$lok_pos_lat))
# max(unique(dfg02$lok_pos_lat))
# min(unique(dfg02$lok_pos_lon))
# max(unique(dfg02$lok_pos_lon))
# 
# min(unique(dfg03$lok_pos_lat))
# max(unique(dfg03$lok_pos_lat))
# min(unique(dfg03$lok_pos_lon))
# max(unique(dfg03$lok_pos_lon))
# 
# 
#_______________________________________________________________________
# section 03 - end -  read in excel files with sample locations
#_______________________________________________________________________
dfg02[grepl("MST2022040", dfg02$U_Pr_Nr),]
dfg03[grepl("MST2022040", dfg03$U_Pr_Nr),]

dfg03[grepl("MST2021028", dfg03$U_Pr_Nr),]
# convert collection dates
dfg03$Dato_inds <- anytime::anydate(dfg03$Dato_inds)
# identify dates that are NA
dtfsmpl <- dfg03$Dato_inds[!is.na(dfg03$Dato_inds)]
# The 'as.Date' function cannot hande NAs
dtfsmpl <- as.Date(dtfsmpl, format = "%YYYY-%mm-%dd")
# add back to data frame for rows in column that are not NAs
dfg03$Dato_inds[!is.na(dfg03$Dato_inds)] <-  dtfsmpl
# reorder (sort) the data frame by two columns 
dfg03 <- dfg03[with(dfg03, 
             order(dfg03$Dato_inds, 
                   MSTno.sbExtr)), ]
# Identify MST samles that are missing a collection date
smcld<- dfg03$U_Pr_Nr[is.na(dfg03$Dato_inds)]
# identify missing dates
missdt <- df_smplMST$MST.nummer[match(smcld,df_smplMST$MST.nummer)]
missdt <- missdt[!is.na(missdt)]
# return the dates from the 'df_smplMST' data frame to the 'dfg03'
# data frame, by matching those isolated missing dates you 
# identified above
dtrplc <- as.character(df_smplMST$Dato_inds[match(missdt,df_smplMST$MST.nummer)])
dtrplc[is.na(dtrplc)] <- "0000-00-00"
dfg03$Dato_inds[match(missdt,dfg03$U_Pr_Nr)] <- dtrplc
# paste together the MST number and the date 
# for collection i.e. : 'Dato_inds'
df_smplMST$MST.n_Dinds <- paste0(df_smplMST$MST.nummer,"_",df_smplMST$Dato_inds)
# find dates and locations for samples '2021-010'
dt_F_2021_010 <- df_smplMST$Dato_inds[grepl("MST2021010",df_smplMST$MST.n_Dinds)]
lok_F_2021_010 <- df_smplMST$Lok_omr01[df_smplMST$Dato_inds %in% dt_F_2021_010]
#subset to only comprise one sample per collection date



df_smplMST02 <- df_smplMST[duplicated(df_smplMST$Dato_inds), ]

# find odd placed sample. Sample 'MST-2021-074' was placed at 
# an odd place
ddfmsmpl <- df_smplMST$Dato_inds[grepl("MST2021074",df_smplMST$MST.n_Dinds)]
locfmsmpl <- df_smplMST$Lok_omr01[df_smplMST$U_Pr_Nr=="MST2021074"]

smpl_notextracted <- df_smplMST[df_smplMST$Dato_inds==ddfmsmpl & 
                                   df_smplMST$Lok_omr01==locfmsmpl,]

# copy column and make it as characters
df_smplMST$Dato_inds <- as.character(df_smplMST$Dato_inds)
# find missing dates, add no date
df_smplMST$Dato_inds[is.na(df_smplMST$Dato_inds)] <- "0000-00-00"
# get only unique per sample location and per date sampled 
# -  to identify one representative of a sample set 
df_smplMST02 <- dplyr::distinct(dfg03, 
                               Dato_inds, 
                               Lok_omr01,
                               .keep_all= TRUE)
# check that they have different number of rows
nrow(df_smplMST02)!=nrow(dfg03)
#reorder the data frame based on extraction number
df_smplMST02 <- df_smplMST02[order(df_smplMST02$Sub_Nmb_for_extraction),]
#define columns to keep
ke<- c("MSTn.di", "Dato_inds","Sub_Nmb_for_extraction")
df_smplMST03 <- df_smplMST02[ke]
#
#https://stackoverflow.com/questions/11996135/create-a-sequential-number-counter-for-rows-within-each-group-of-a-dataframe
# make an consecutive sequential count of sample ID Nos
library(dplyr)
# add numbering to get a count that can be compared 
# to strip-of-eight tubes
df_smplMST03 <- df_smplMST03 %>% dplyr::mutate(id = row_number())
nrMST03 <- nrow(df_smplMST03)
nrMST03 <- ceiling(nrMST03)
df_smplMST03$nofstrp <- rep(seq(1,8,1),
                            (nrMST03/8))[seq(1,nrow(df_smplMST03),1)]


#df_smplMST[grepl("MST2022040", df_smplMST$MST.n_Dinds),]

#_______________________________________________________________________
# section 04 - start - write out the data frames
#_______________________________________________________________________
# # write out the combined text-reports from the qpcr runs
# write.table(df_txtr01,file=paste0(wd00_wd03,
#           "/table01_combined_qpcr_textreports.csv"),
#           sep = ";")
# write out the combined concentrations from filter extractions
write.table(df_cex01,file=paste0(wd00_wd03,
              "/table02_combined_extractions_from_filters.csv"),
          sep = ";")
# write out the combined concentrations from filter extractions
write.table(df_smplMST,file=paste0(wd00_wd03,
            "/table03_combined_MST_samples_locations.csv"),
          sep = ";")

# copy the data frame
df_smplMST05 <- df_smplMST
##
##
##
# ERROR HERE
df_smplMST05$yearSmpl <- gsub("(^[0-9]{+})-.?.*","\\1",df_smplMST05$Dato_inds)
df_smplMST05$montSmpl <- gsub("(^[0-9]{+})-([0-9]{+})-.?.*","\\2",df_smplMST05$Dato_inds)
# make the sampling year numeric
df_smplMST05$yearSmpl <- as.numeric(df_smplMST05$yearSmpl)
df_smplMST05$montSmpl <- as.numeric(df_smplMST05$montSmpl)
# evaluate which season was sampled
df_smplMST05$seasSmpl <- ifelse(df_smplMST05$montSmpl>=7,"1stseason","2ndseason")
# https://stackoverflow.com/questions/1660124/how-to-sum-a-variable-by-group
  # summarise specific variables, not all
# df_smplMST05.01 <- df_smplMST05 %>% 
#           dplyr::group_by(lokalitet_vanda,
#                           Dato_inds) %>% 
#   dplyr::summarise(across("Dato_inds",
#                           list(sumDtVa = sum)))
##
##
##

#_______________________________________________________________________
# section 04 - end - write out the data frames
#_______________________________________________________________________

#________________________________________________________________________
# section 05 - start - combine the data frames qPCR and MST sample
# data frames
#________________________________________________________________________
#-------------------------------------------------------------------------
# section 05.1 - standardize MST numbers across data frame that holds 
# sample location data and the data frame that holds information on the 
# extraction of the sample

# get the 'U_Pr_Nr' that only has 3 digits and add a 0
ms3n <- df_smplMST$U_Pr_Nr[grepl("^MST[0-9]{3}$",df_smplMST$U_Pr_Nr)]
ms3n <-gsub("MST","MST0",ms3n)
df_smplMST$U_Pr_Nr[grepl("^MST[0-9]{3}$",df_smplMST$U_Pr_Nr)] <- ms3n
# # get the 'MST.nummer' that only has 3 digits and add a 0
ms3n <- df_smplMST$MST.nummer[grepl("^MST[0-9]{3}$",df_smplMST$MST.nummer)]
ms3n <-gsub("MST","MST0",ms3n)
df_smplMST$MST.nummer[grepl("^MST[0-9]{3}$",df_smplMST$MST.nummer)] <- ms3n
# make subsetted data frames that holds either the duplicated
# rows or the non-duplicated rows
MSTno.dpl <- df_smplMST$U_Pr_Nr[duplicated(df_smplMST$U_Pr_Nr)]
MSTcex01.dpl <- df_cex01$U_Pr_Nr[duplicated(df_cex01$U_Pr_Nr)]
df_smplMST.dpl <- df_smplMST[which(df_smplMST$U_Pr_Nr %in% MSTno.dpl),]
df_smplMST.ndpl <- df_smplMST[which(!df_smplMST$U_Pr_Nr %in% MSTno.dpl),]
# use the row numbers for the duplicated rows to subset the data frames
df_cex01.dpl <- df_cex01[which(df_cex01$U_Pr_Nr %in% MSTcex01.dpl),]
df_cex01.ndpl <- df_cex01[which(!df_cex01$U_Pr_Nr %in% MSTcex01.dpl),]

# re order the data frames by the 'U_Pr_Nr'
df_smplMST.ndpl <- df_smplMST.ndpl[order(df_smplMST.ndpl$U_Pr_Nr),]
df_smplMST.dpl <- df_smplMST.dpl[order(df_smplMST.dpl$U_Pr_Nr),]
# find duplicated and non duplicated elements for the MST sample numbers
df_cex01.ndpl <- df_cex01.ndpl[order(df_cex01.ndpl$U_Pr_Nr),]
df_cex01.dpl <- df_cex01.dpl[order(df_cex01.dpl$U_Pr_Nr),]
#-------------------------------------------------------------------------
# section 05.2 - combine data frames that holds the extracted samples plus the
# concentration measured of DNA and the sub numbering added for extraction

# combine by 'U_Pr_Nr' for the non-duplicated MST numbers
df_cex02 <- dplyr::left_join(df_cex01.ndpl ,
                         df_smplMST.ndpl %>% dplyr::select(
                           everything()),
                         by = "U_Pr_Nr")
# substitute to have a common column between data frames
# that have the MST number and the sampling date combined
df_smplMST.dpl$MSTno_di <- gsub("_","::",df_smplMST.dpl$MST.n_Dinds)
# exclude MST samples that have not been extracted
df_cex01.dpl <- df_cex01.dpl[!is.na(df_cex01.dpl$Dato.for.ekstraktion),]
# join data frames for the duplicated MST numbers
df_cex03 <- dplyr::left_join(df_cex01.dpl ,
                              df_smplMST.dpl %>% dplyr::select(
                                everything()),
                              by = "MSTno_di")
#
df_cex03[grepl("2021028",df_cex03$MSTno_di),]
# copy the columns into a new column without the .x appended
# on the column name
df_cex03$MST.n_Dinds <- df_cex03$MST.n_Dinds.y
df_cex03$U_Pr_Nr <- df_cex03$U_Pr_Nr.x
df_cex03$MST.nummer <-  df_cex03$MST.nummer.x
# get rid of the doubled columns
ctkee <- colnames(df_cex03)[!grepl("\\.x",colnames(df_cex03))]
ctkee <- ctkee[!grepl("\\.y",ctkee)]
df_cex03 <- df_cex03[ctkee]
# make an inner_join to combine data frames by 'MST.nummer' 
df_cex03.1 <- dplyr::inner_join(df_cex01.ndpl ,
                                df_smplMST.dpl %>% dplyr::select(
                                  everything()),
                                by = "MST.nummer")
# get columns without the '.y' appended
df_cex03.1$U_Pr_Nr <- df_cex03.1$U_Pr_Nr.y
df_cex03.1$MSTno_di <- df_cex03.1$MSTno_di.y
# get rid of the doubled columns
ctkee <- colnames(df_cex03.1)[!grepl("\\.x",colnames(df_cex03.1))]
ctkee <- ctkee[!grepl("\\.y",ctkee)]
df_cex03.1 <- df_cex03.1[ctkee]
# make an inner_join using the "MST.nummer"
df_cex03.2 <- dplyr::inner_join(df_cex01.dpl ,
                                df_smplMST.ndpl %>% dplyr::select(
                                  everything()),
                                by = "MST.nummer")
# get columns without the '.y' and '.x' appended 
df_cex03.2$MST.n_Dinds <- df_cex03.2$MST.n_Dinds.y
df_cex03.2$U_Pr_Nr <- df_cex03.2$U_Pr_Nr.x
# get rid of the doubled columns
ctkee <- colnames(df_cex03.2)[!grepl("\\.x",colnames(df_cex03.2))]
ctkee <- ctkee[!grepl("\\.y",ctkee)]
df_cex03.2 <- df_cex03.2[ctkee]
# use inner_join to match up  data frames
df_cex03.3 <- dplyr::inner_join(df_cex01.ndpl ,
                                df_smplMST.ndpl %>% dplyr::select(
                                  everything()),
                                by = "MST.nummer")
# copy the 'U_Pr_Nr' toe nex column without the appended '.x'
df_cex03.3$U_Pr_Nr <- df_cex03.3$U_Pr_Nr.x
# get rid of the doubled columns
ctkee <- colnames(df_cex03.3)[!grepl("\\.x",colnames(df_cex03.3))]
ctkee <- ctkee[!grepl("\\.y",ctkee)]
df_cex03.3 <- df_cex03.3[ctkee]

# reorder the data frames
df_cex01.dpl <- df_cex01.dpl[order(df_cex01.dpl$Dato.for.ekstraktion),]
df_smplMST.dpl <- df_smplMST.dpl[order(df_smplMST.dpl$U_Pr_Nr),]
# Substitute to get a  'MST.n_Dinds' with and underscore
df_cex01.dpl$MST.n_Dinds <- gsub("::","_",df_cex01.dpl$MSTno_di)
# join data frames for the duplicated MST numbers
df_cex04 <- dplyr::inner_join(df_cex01.dpl ,
                             df_smplMST %>% dplyr::select(
                               everything()),
                             by = "MST.n_Dinds")
#
df_cex04[grepl("2021028",df_cex04$MSTno_di),]

# make a common column to use for inner_join
df_cex01$MST.n_Dinds<-gsub("::","_",df_cex01$MSTno_di)
# join data frames for the duplicated MST numbers
df_cex04.1 <- dplyr::inner_join(df_cex01 ,
                              df_smplMST %>% dplyr::select(
                                everything()),
                              by = "MST.n_Dinds")

# make new columns that are copies of already existing columns, but without
# the .x nad .y appended
df_cex04.1$U_Pr_Nr <- df_cex04.1$U_Pr_Nr.y
df_cex04.1$MST.nummer <- df_cex04.1$MST.nummer.y
# get rid of the doubled columns
ctkee <- colnames(df_cex04.1)[!grepl("\\.x",colnames(df_cex04.1))]
ctkee <- ctkee[!grepl("\\.y",ctkee)]
df_cex04.1 <- df_cex04.1[ctkee]

# remove duplicated rows
# see this question: https://stackoverflow.com/questions/60861391/remove-duplicate-rows-based-on-multiple-columns-using-dplyr-tidyverse
df_cex04 <- df_cex04 %>%
  dplyr::distinct(MST.n_Dinds, .keep_all = TRUE)
# bind rows together
df_cex05 <- bind_rows(df_cex02, 
                      df_cex03, 
                      df_cex03.1,
                      df_cex03.2,
                      df_cex03.3,
                      df_cex04,
                      df_cex04.1)

#
# replace if 'MST.nummer' is missing and if 'U_Pr_Nr' is missing
df_cex05$MST.nummer[is.na(df_cex05$MST.nummer)] <- df_cex05$MST.nummer.y[is.na(df_cex05$MST.nummer)]
df_cex05$MST.nummer[is.na(df_cex05$MST.nummer)] <- df_cex05$MST.nummer.x[is.na(df_cex05$MST.nummer)]
df_cex05$U_Pr_Nr.x[is.na(df_cex05$U_Pr_Nr.x)] <- df_cex05$U_Pr_Nr.y[is.na(df_cex05$U_Pr_Nr.x)]
df_cex05$U_Pr_Nr.x[is.na(df_cex05$U_Pr_Nr.x)] <- df_cex05$U_Pr_Nr[is.na(df_cex05$U_Pr_Nr.x)]
df_cex05$U_Pr_Nr[is.na(df_cex05$U_Pr_Nr)] <- df_cex05$U_Pr_Nr.x[is.na(df_cex05$U_Pr_Nr)]
df_cex05$U_Pr_Nr[is.na(df_cex05$U_Pr_Nr)]  <- df_cex05$MST.nummer[is.na(df_cex05$U_Pr_Nr)]
df_cex05$U_Pr_Nr[is.na(df_cex05$U_Pr_Nr)] <-  df_cex05$MST.nummer.x[is.na(df_cex05$U_Pr_Nr)]
df_cex05$U_Pr_Nr[is.na(df_cex05$U_Pr_Nr)] <-  df_cex05$MST.nummer.y[is.na(df_cex05$U_Pr_Nr)]
df_cex05$MST.nummer[is.na(df_cex05$MST.nummer)] <- df_cex05$U_Pr_Nr[is.na(df_cex05$MST.nummer)]

# remove the NA rows for 'MST.nummer.x'
df_cex05 <- df_cex05[!is.na(df_cex05$MST.nummer.x) ,]


# select columns to keep, if they do NOT have '.x' and '.y'
ctkee <- colnames(df_cex05)[!grepl("\\.x",colnames(df_cex05))]
ctkee <- ctkee[!grepl("\\.y",ctkee)]
# only keep columns without '.y' and without '.x'
df_cex05 <- df_cex05[ctkee]
# reorder the data frame by the 'U_Pr_Nr'
df_cex05 <- df_cex05[order(df_cex05$MST.n_Dinds),]
# Sort data by number of NA's in each line
# https://stackoverflow.com/questions/24093446/sort-data-by-number-of-nas-in-each-line
df_cex05 <- df_cex05[order(rowSums(is.na(df_cex05))), ]

# remove duplicated rows
# see this question: https://stackoverflow.com/questions/60861391/remove-duplicate-rows-based-on-multiple-columns-using-dplyr-tidyverse
df_cex05 <- df_cex05 %>%
  dplyr::filter(!duplicated(.))

# # add the missing subnumber for extraction to the 'MST.nummer
# # numbers that do not have the subnumber for extraction
df_cex05$MST_STEXno <- paste0(df_cex05$MST.nummer,
                              "STEX",
                              df_cex05$Sub_Nmb_for_extraction)
# substitute "-" in 'MST_STEXno'
df_cex05$MST_STEXno <- gsub("-","",df_cex05$MST_STEXno)
# remove duplicated rows
# see this question: https://stackoverflow.com/questions/60861391/remove-duplicate-rows-based-on-multiple-columns-using-dplyr-tidyverse
df_cex05 <- df_cex05 %>%
  dplyr::filter(!duplicated(.))
#
#-------------------------------------------------------------------------
# section 05.3 - standardize MST numbers across data frame that holds 
# the combined sample location data and  information on the 
# extraction of the sample, and standardize in the data table that holds all
# the combined data from the qPCR runs

# global substitute  in the well sample number
smplMSTnos <- gsub("MST-([0-9]{4})-","MST\\1",df_txtr01$smpltp)
smplMSTnoMST.nummer.ysmplMSTnos <- gsub("MST([0-9]{4})-([0-9]{3})","MST\\1\\2",smplMSTnos)
smplMSTnos <- gsub("MST-([0-9]{+3})","MST\\1",smplMSTnos)
df_txtr01$smplNo <- gsub("STEX([0-9]{+2})-","STEX\\1_",smplMSTnos)

# substitute the '-' in the MST number 
df_cex05$smplNo <- gsub("-","",df_cex05$U_Pr_Nr)
df_cex05$smplNoSTEX <- paste0(df_cex05$smplNo,"STEX",df_cex05$Sub_Nmb_for_extraction)

# differentiate between duplicated and non-duplicated sample numbers 
dpl.c05 <- df_cex05$smplNo[duplicated(df_cex05$smplNo)]
df_cex05.dpl <- df_cex05[which(df_cex05$smplNo %in%  dpl.c05 ),]
df_cex05.ndpl <- df_cex05[!duplicated(df_cex05$smplNo),]

# copy the 'smplNo' column but wihtout the 
df_txtr01$smplNo02 <- gsub("-","",df_txtr01$smplNo)
df_txtr01$smplNoSTEX <- df_txtr01$smplNo02
# find duplicates in the assmebled txtreports
dpl.c05.m <- dpl.c05[(dpl.c05 %in% df_txtr01$smplNo02)]
df_txtr01.f <- df_txtr01[!(df_txtr01$smplNo02 %in% dpl.c05.m),]
df_txtr01.m <- df_txtr01[(df_txtr01$smplNo02 %in% dpl.c05.m),]
# It appears the non-duplicated  (in df_txtr01.f) are from qpcr s 
# performed from 2023 and onwards
# and the the nonduplicated  (in df_txtr01.m) are from qpcr s 
# performed from 2022 and prior
unique(df_txtr01.f$qpcrrundate[grepl("2021028",df_txtr01.f$smplNo02)])
unique(df_txtr01.f$smpltp[grepl("2021028",df_txtr01.f$smplNo02)])
unique(df_txtr01.m$qpcrrundate[grepl("2021028",df_txtr01.m$smplNo02)])
unique(df_txtr01.m$smpltp[grepl("2021028",df_txtr01.m$smplNo02)])
#unique(df_txtr01.m$smplNo02)
# identify the 'smplNoSTEX numbers already found
alrfMST <- unique(df_txtr01.f$smplNoSTEX)
# use this to subset the extractions table
# to make it only inlcude the 'not found'
df_cex05.dpl.nf <- df_cex05.dpl[!(df_cex05.dpl$smplNoSTEX %in% alrfMST),]
# match to get the STEX MST numbers for sample numbers missing the
# STEX information
nfSTEXnos <- df_cex05.dpl.nf$smplNoSTEX[match(df_txtr01.m$smplNo02,df_cex05.dpl.nf$smplNo)]
# use these not found STEX number to replace the STEX numbers
df_txtr01.m$smplNoSTEX <- nfSTEXnos 
# and also replace the sample- type number to get at STEX number format
df_txtr01.m$smpltp <- gsub("_","-",nfSTEXnos)
unique(df_txtr01.m$smplNoSTEX)

df_cex05[grepl("MST2021085",df_cex05$smplNo02),]
# bind the rows together again into one data frame
df_txtr01 <- rbind(df_txtr01.m,df_txtr01.f)
# Here I have to add a vector containing the 'df_txtr01$smpltp'
# that end up being unmatched
unmsmplnos <- c("MST2021-031", "MST2021-032", "MST2021-040", "MST2021-044", 
  "MST2021-049", "MST2021-051", "MST2021-052", "MST2021-062", "MST2021-071", 
  "MST2021-075", "MST2021-081", "MST2021-088", "MST2021-089", "MST2021-099", 
  "MST2021-008", "MST2021-039", "MST2021-043", "MST2021-055", "MST2021-058", 
  "MST2021-064", "MST2021-070", "MST2021-076", "MST2021-079", "MST2021-082", 
  "MST2021-083", "MST2021-085")
# differentiate between matched and unmatched
df_txtr01.unm <- df_txtr01[(df_txtr01$smpltp %in% unmsmplnos),]
df_txtr01.mtc <- df_txtr01[!(df_txtr01$smpltp %in% unmsmplnos),]
# identify unique smplNoSTEX in the unmatched data frame
unmsmplnos2 <- unique(df_txtr01.unm$smplNoSTEX)
unm.smplNoSTEX <- df_cex05$smplNoSTEX[match(df_txtr01.unm$smplNoSTEX,df_cex05$smplNo)]
df_txtr01.unm$smplNoSTEX <- unm.smplNoSTEX 
df_cex05.unm <- df_cex05[match(unmsmplnos2,df_cex05$smplNo),]
df_cex05.unm <- df_cex05.unm %>% dplyr::arrange(smplNoSTEX)
# bind back the matched and the unmatched together to a data frame
df_txtr01 <- rbind(df_txtr01.unm,df_txtr01.mtc)
# Check if the 'MST2021028' sample still has qpcrs from 2022 and 2023
df_txtr01.MST2021028 <- df_txtr01[grepl("2021028",df_txtr01$smplNoSTEX),]
unique(df_txtr01.MST2021028$qpcrrundate)
unique(df_txtr01.MST2021028$smplNoSTEX)


# limit the data frame with all qpcr runs to only comprise
# the MST samples that do NOT have an STEX number
df_txtr02 <- df_txtr01[!grepl("STEX",df_txtr01$smplNo),]
# and also limit to the numbers that DO have STEX numbers
df_txtr03 <- df_txtr01[grepl("STEX",df_txtr01$smplNo),]
# copy column to a column with a different name
df_txtr03$smplNoSTEX <- df_txtr03$smplNo
# check sample 'smplNoSTEX' with '028'
df_cex05.dpl$smplNoSTEX[grepl("028",df_cex05.dpl$smplNo02)]


# write out the combined text-reports from the qpcr runs
write.table(df_txtr01,file=paste0(wd00_wd03,
                                  "/table01_combined_qpcr_textreports.csv"),
            sep = ";")
#-------------------------------------------------------------------------
# section 05.4 - combine the combined data frame that 
# holds the extracted samples combined with the
# concentration measured of DNA and the sub numbering added for extraction, 
# and combine this with the data frame that holds all the qPCR run data

#df_cex05.ndpl$smplNo[grepl("2021",df_cex05.ndpl$smplNo)]
#df_cex05.ndpl$smplNo
# join data frames by the sample number 
df_txtr04 <- dplyr::left_join(df_txtr02 ,
                              df_cex05.ndpl %>% dplyr::select(
                                everything()),
                              by = "smplNo")

# join data frames by the sample number and STEX number
df_txtr05 <- dplyr::left_join(df_txtr03 ,
                              df_cex05.ndpl %>% dplyr::select(
                                 everything()),
                               by = "smplNoSTEX")
# replace if NA in 'df_txtr04'
df_txtr04$MST_STEXno[is.na(df_txtr04$MST_STEXno)] <- df_txtr04$smplNoSTEX.x[is.na(df_txtr04$MST_STEXno)]
# replace if NA in 'df_txtr05'
df_txtr05$MST_STEXno[is.na(df_txtr05$MST_STEXno)] <- df_txtr05$smplNoSTEX[is.na(df_txtr05$MST_STEXno)]

unique(df_txtr04$MST_STEXno)[grepl("028",unique(df_txtr04$MST_STEXno))]

#

# bind the data frames together
df_txtr06 <- bind_rows(df_txtr04,df_txtr05)
# replace if 'smplNo' is NA
df_txtr06$smplNo[is.na(df_txtr06$smplNo)] <- df_txtr06$smplNo.x[is.na(df_txtr06$smplNo)]
df_txtr06$smplNo[is.na(df_txtr06$smplNo)] <- df_txtr06$smplNo.y[is.na(df_txtr06$smplNo)]
# select columns to keep, if they do NOT have '.x' and '.y'
ctkee <- colnames(df_txtr06)[!grepl("\\.x",colnames(df_txtr06))]
ctkee <- ctkee[!grepl("\\.y",ctkee)]
df_txtr06 <- df_txtr06[ctkee]
# get rows that have NA in the 'U_Pr_Nr'
# and have NA in the  'MSTno_di', and that do not have
# "std" in the 'smplNo' column
df_txtr07 <- df_txtr06[is.na(df_txtr06$U_Pr_Nr)	& 
                      is.na(df_txtr06$MSTno_di) &
                      !grepl("std",df_txtr06$smplNo) &
                      !grepl("NTC",df_txtr06$smpltp) &
                      !grepl("NEC",df_txtr06$smpltp)  ,]

# get rid of the doubled columns
ctkee <- colnames(df_txtr07)[!grepl("\\.x",colnames(df_txtr07))]
df_txtr07 <- df_txtr07[ctkee]
colnames(df_txtr07) <- gsub("\\.y","",colnames(df_txtr07))

df_txtr07_2021028 <- df_txtr07$MST_STEXno[grepl("2021028",df_txtr07$MST_STEXno)]
df_txtr06_2021028 <- df_txtr06$MST_STEXno[grepl("2021028",df_txtr06$MST_STEXno)]
unique(df_txtr07_2021028)
unique(df_txtr06_2021028)
# Remove columns from dataframe where ALL values are NA
#https://stackoverflow.com/questions/2643939/remove-columns-from-dataframe-where-all-values-are-na
df_txtr07 <- df_txtr07[,colSums(is.na(df_txtr07))<nrow(df_txtr07)]

# if 'MST_STEXno' is NA then replace with 'smplNoSTEX'
df_txtr07$MST_STEXno[is.na(df_txtr07$MST_STEXno)] <- df_txtr07$smplNoSTEX[is.na(df_txtr07$MST_STEXno)]
# find the MST numbers with only 3 digits after the MST text
MSTnmw3dig <- df_txtr07$MST_STEXno[grepl("^MST[0-9]{3}STEX[0-9]{2}_[0-9]{2}",df_txtr07$MST_STEXno)]
# modify the MST numbers with only 3 digits after the MST text 
MSTnmw3dig <- gsub("^MST([0-9]{3})(STEX[0-9]{2}_[0-9]{2})","MST0\\1\\2",MSTnmw3dig)
df_txtr07$MST_STEXno[grepl("^MST[0-9]{3}STEX[0-9]{2}_[0-9]{2}",df_txtr07$MST_STEXno)] <- MSTnmw3dig

# join data frames by the sample number and STEX number
df_txtr08 <- dplyr::inner_join(df_txtr07 ,
                              df_cex05.dpl %>% dplyr::select(
                                everything()),
                              by = "smplNoSTEX")


df_txtr08_2021028 <- df_txtr08$MST_STEXno[grepl("2021028",df_txtr08$MST_STEXno)]
unique(df_txtr08_2021028)

# replace if 'smplNo' is NA
df_txtr08$smplNo[is.na(df_txtr08$smplNo)] <- df_txtr08$smplNo.x[is.na(df_txtr08$smplNo)]
df_txtr08$smplNo[is.na(df_txtr08$smplNo)] <- df_txtr08$smplNo.y[is.na(df_txtr08$smplNo)]
# select columns to keep, if they do NOT have '.x' and '.y'
ctkee <- colnames(df_txtr08)[!grepl("\\.x",colnames(df_txtr08))]
ctkee <- ctkee[!grepl("\\.y",ctkee)]
df_txtr08 <- df_txtr08[ctkee]
# join data frames by the sample number and smplNo number
df_txtr09 <- dplyr::inner_join(df_txtr07 ,
                               df_cex05.dpl %>% dplyr::select(
                                 everything()),
                               #relationship = "many-to-many",
                               by = "smplNoSTEX")

df_txtr09_2021028 <- df_txtr09$MST_STEXno[grepl("2021028",df_txtr09$MST_STEXno)]
unique(df_txtr09_2021028)

# replace if 'smplNo' is NA
df_txtr09$smplNo[is.na(df_txtr09$smplNo)] <- df_txtr09$smplNo.x[is.na(df_txtr09$smplNo)]
df_txtr09$smplNo[is.na(df_txtr09$smplNo)] <- df_txtr09$smplNo.y[is.na(df_txtr09$smplNo)]
# select columns to keep, if they do NOT have '.x' and '.y'
ctkee <- colnames(df_txtr09)[!grepl("\\.x",colnames(df_txtr09))]
ctkee <- ctkee[!grepl("\\.y",ctkee)]
df_txtr09 <- df_txtr09[ctkee]
# Exclude samples that have not been
# extracted  -  they have no  date in the "Dato.for.ekstraktion"
df_cex05.ndpl <- df_cex05.ndpl[!is.na(df_cex05.ndpl$Dato.for.ekstraktion),]
# join data frames by the sample number and smplNo number
# for this last 'inner_join' there were MST sample numbers
# that still had not been matched , to work around this I added the 
#relationship = "many-to-many",
df_txtr10 <- dplyr::inner_join(df_txtr07 ,
                               df_cex05.dpl %>% dplyr::select(
                                 everything()),
                               #relationship = "many-to-many",
                               by = "smplNoSTEX")

df_txtr10_2021028 <- df_txtr10$MST_STEXno[grepl("2021028",df_txtr10$MST_STEXno)]
unique(df_txtr10_2021028)

# keep only unique MST_STEXno numbers
df_cex05.1 <- df_cex05 %>%
  dplyr::distinct(MST_STEXno, .keep_all = TRUE)

# match with innner_join
df_txtr11 <- dplyr::inner_join(df_txtr07 ,
                               df_cex05.1 %>% dplyr::select(
                                 everything()),
                               #relationship = "many-to-many",
                               by = "MST_STEXno")
# copy to rename columns with 'y.' appended
df_txtr11$smplNoSTEX  <-  df_txtr11$smplNoSTEX.y
df_txtr11$smplNo <- df_txtr11$smplNo.y
# select columns to keep, if they do NOT have '.x' and '.y'
ctkee <- colnames(df_txtr11)[!grepl("\\.x",colnames(df_txtr11))]
ctkee <- ctkee[!grepl("\\.y",ctkee)]
df_txtr11 <- df_txtr11[ctkee]


df_txtr11_2021028 <- df_txtr11$MST_STEXno[grepl("2021028",df_txtr11$MST_STEXno)]
unique(df_txtr11_2021028)

# copy a column so that there is a common column for joining
df_cex05$smplNo02<- df_cex05$smplNo
df_cex05.dpl$smplNo02<- df_cex05.dpl$smplNo
df_cex05.ndpl$smplNo02<- df_cex05.ndpl$smplNo
# use inner_join to match between data frames
df_txtr13 <- dplyr::inner_join(df_txtr07 ,
                               df_cex05.ndpl %>% dplyr::select(
                                 everything()),
                               #relationship = "many-to-many",
                               by = "smplNo02")

# copy columns so that they do not have the '.y' appended 
df_txtr13$smplNoSTEX <- df_txtr13$smplNoSTEX.y
df_txtr13$smplNo <- df_txtr13$smplNo.y
df_txtr13$MST_STEXno <- df_txtr13$MST_STEXno.y
# select columns to keep, if they do NOT have '.x' and '.y'
ctkee <- colnames(df_txtr13)[!grepl("\\.x",colnames(df_txtr13))]
ctkee <- ctkee[!grepl("\\.y",ctkee)]
df_txtr13 <- df_txtr13[ctkee]
# match with innner_join
df_txtr11 <- dplyr::inner_join(df_txtr07 ,
                               df_cex05.1 %>% dplyr::select(
                                 everything()),
                               #relationship = "many-to-many",
                               by = "MST_STEXno")
# copy to rename columns with 'y.' appended
df_txtr11$smplNoSTEX  <-  df_txtr11$smplNoSTEX.y
df_txtr11$smplNo <- df_txtr11$smplNo.y
# select columns to keep, if they do NOT have '.x' and '.y'
ctkee <- colnames(df_txtr11)[!grepl("\\.x",colnames(df_txtr11))]
ctkee <- ctkee[!grepl("\\.y",ctkee)]
df_txtr11 <- df_txtr11[ctkee]
# replace if 'smplNo' is NA
df_txtr10$smplNoSTEX[is.na(df_txtr10$smplNoSTEX)] <- df_txtr10$smplNoSTEX.x[is.na(df_txtr10$smplNo)]
df_txtr10$smplNoSTEX[is.na(df_txtr10$smplNoSTEX)] <- df_txtr10$smplNoSTEX.y[is.na(df_txtr10$smplNo)]
# make a 'MST_STEXno' colunm
df_txtr10$MST_STEXno <- df_txtr10$smplNoSTEX
# get 'smplNoSTEX' if 'smplNo.x' is NA
df_txtr10$smplNo.x[is.na(df_txtr10$smplNo.x)] <- df_txtr10$smplNo.y[is.na(df_txtr10$smplNo.x)]
df_txtr10$smplNo.x[is.na(df_txtr10$smplNo.x)] <- df_txtr10$smplNoSTEX[is.na(df_txtr10$smplNo.x)]
# make a new 'smplNo' column that does not have the x appended in the
# column name
df_txtr10$smplNo <- df_txtr10$smplNo.x
# select columns to keep, if they do NOT have '.x' and '.y'
ctkee <- colnames(df_txtr10)[!grepl("\\.x",colnames(df_txtr10))]
ctkee <- ctkee[!grepl("\\.y",ctkee)]
df_txtr10 <- df_txtr10[ctkee]

# df_cex05.1[grepl("MST2023015",df_cex05.1$MST.nummer),]
#df_txtr10[grepl("STEX32-06",df_txtr10$smpltp),]

# bind all rows together into one data frame
df_ctxtr <- bind_rows(df_txtr06,
                       df_txtr07,
                       df_txtr08,
                       df_txtr09,
                       df_txtr10,
                       df_txtr11)#,
                       #df_txtr13)


# remove duplicated rows
# see this question: https://stackoverflow.com/questions/60861391/remove-duplicate-rows-based-on-multiple-columns-using-dplyr-tidyverse
df_ctxtr <- df_ctxtr %>%
  dplyr::filter(!duplicated(.))

# Sort data by number of NA's in each line
# https://stackoverflow.com/questions/24093446/sort-data-by-number-of-nas-in-each-line
df_ctxtr <- df_ctxtr[order(rowSums(is.na(df_ctxtr))), ]

df_ctxtr <- df_ctxtr %>% dplyr::arrange(qpcrno, Well)


df_ctxtr.dpl <- df_ctxtr %>%
  dplyr::filter(duplicated(qpcrno, Well))

#View(df_ctxtr.dpl)
df_ctxtr.NEC <- df_ctxtr[grepl("NEC",df_ctxtr$smpltp),]
#View(df_ctxtr.NEC)

df_ctxtr_2021028 <- df_ctxtr[grepl("2021028",df_ctxtr$MST_STEXno),]
unique(df_ctxtr_2021028$smplNoSTEX)

#View(df_ctxtr_2021028)

# there can only be a unique "Well" per "qpcrno", so
# # using dplyr::distinct on  "Well" per "qpcrno" that occurs more than once
# # are removed
df_ctxtr <- df_ctxtr %>%
  dplyr::distinct(qpcrno, Well, .keep_all = TRUE)


# Sort (order) data frame rows by multiple columns
# https://stackoverflow.com/questions/1296646/sort-order-data-frame-rows-by-multiple-columns
df_ctxtr <- df_ctxtr[with(df_ctxtr, order(qpcrno, 
                                             smpltp,
                                             replno,
                                             Well)), ]


# ##___ERROR CHECK
# define columns to keep to be able to check resulting data frame
ctkee <- c("replno",	
           "smpltp",
           "qpcrno",
           "qpcrrundate",
           "lok_pos_lat",	"lok_pos_lon",
           "MST.n_Dinds","MST.nummer","MST_STEXno",
           "smplNoSTEX"
)
# susbet data frame
df_ctxtr.1 <- df_ctxtr[ctkee]

#
df_ctxtr.1$MST_STEXno[grepl("2021028",df_ctxtr.1$MST_STEXno)]

#
# exclude std1 rows
df_ctxtr.1 <- df_ctxtr.1[!grepl("std",df_ctxtr.1$smpltp),]
# exclude std1 rows
df_ctxtr.1 <- df_ctxtr.1[!grepl("NEC",df_ctxtr.1$smpltp),]
df_ctxtr.1 <- df_ctxtr.1[!grepl("NEK",df_ctxtr.1$smpltp),]
df_ctxtr.1 <- df_ctxtr.1[!grepl("NTC",df_ctxtr.1$smpltp),]
df_ctxtr.1 <- df_ctxtr.1[!grepl("^MSTEX",df_ctxtr.1$smpltp),]
df_ctxtr.1 <- df_ctxtr.1[!grepl("^NaN",df_ctxtr.1$smpltp),]
#
df_ctxtr.1 <- df_ctxtr.1[is.na(df_ctxtr.1$lok_pos_lat),]
notfnd <- unique(df_ctxtr.1$smpltp)
notfnd <- notfnd[order(notfnd)]
notfnd
dput(notfnd)
#
#df_ctxtr[grepl("1026",df_ctxtr$qpcrno),]


# dput(notfnd)
# length(unique(notfnd))
#View(df_ctxtr.1)
# ##___ERROR CHECK

# write out the combined concentrations from filter extractions
write.table(df_ctxtr,file=paste0(wd00_wd03,
                    "/table04_combined_MST_samples_and_qpcr_data.csv"),
          sep = ";", row.names = F)

df_ctxtr.MST2021028 <- df_ctxtr[grepl("MST2021028",df_ctxtr$MST.n_Dinds),]
#View(df_ctxtr.MST2021028)
unique(df_ctxtr.MST2021028$Dato.for.ekstraktion)
unique(df_ctxtr.MST2021028$Lok_omr01)
#________________________________________________________________________
# section 05 - end - combine the data frames qPCR and MST sample
# data frames
#________________________________________________________________________
