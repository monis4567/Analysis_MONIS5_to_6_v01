#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
# R-code provided for the project:
# “MONIS6”
# Authors: Steen Wilhelm Knudsen.
# Change the working directory to a path on your own computer , and run
# the individual parts below to reproduce the diagrams presented in the paper
#
# All input data required needs to be available as csv-files in the same directory 
# as this R-code use for working directory.
#
# Occassionally the code will have difficulties producing the correct diagrams,
# if the packages and libraries are not installed.
# Make sure the packages are installed, and libraries are loaded, if the R-code
# fails in producing the diagrams.
#
#________________IMPORTANT!!_________________________________________________#
# (1)
#You have to change the path to the working directory before running this code
#
# (2)
# The 4 data input files required:
#
# must be located in the same working directory - as specified in the code below
#
#This code is able to run in:
#
#R studio: Version 0.98.994 – © 2009-2013 RStudio, Inc.
#R version:
# > version
# _                           
#____________________________________________________________________________#
# #remove everything in the working environment, without a warning!!
# rm(list=ls())
#see this
#website
#on how to only install required packages
# #https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
# if (!require("pacman")) install.packages("pacman")
# pacman::p_load(
#   scales, 
#   fields, 
#   gplots,
#   plyr)#,
library(plyr)
library(scales)
library(gplots)
library(fields)
#ReporteRs)
## install the package 'scales', which will allow you to make points on your plot more transparent
# #install.packages("scales")
# if(!require(scales)){
#   install.packages("scales")
#   library(scales)
# }
library(scales)
# #install.packages("fields")
# if(!require(fields)){
#   install.packages("fields")
#   library(fields)
# }
library(fields)
# ## install the package 'gplots', to be able to translate colors to hex - function: col2hex
# #install.packages("gplots")
# if(!require(gplots)){
#   install.packages("gplots")
#   library(gplots)
# }
library(gplots)
## install the package 'glad', to be able to color using the function 'myPalette'
#install.packages("glad")
#library(glad)
require(graphics)
#get package to do count number of observations that have the same value at earlier records:
# see this website: https://stackoverflow.com/questions/11957205/how-can-i-derive-a-variable-in-r-showing-the-number-of-observations-that-have-th
#install.packages("plyr")
# if(!require(plyr)){
#   install.packages("plyr")
#   library(plyr)
# }
library(plyr)
# I had difficulties getting 'xlsx' installed
# and found I needed "rJava"
# but "rJava" required this in a terminal first :
# $ sudo apt-get install openjdk-8-jdk
# $ sudo apt-get install default-jre
# $ sudo apt-get install default-jdk
# $ sudo R CMD javareconf
# Solution was here: https://stackoverflow.com/questions/3311940/r-rjava-package-install-failing
# install.packages("rJava")
# Once "rJava" was installed
#I was able to install "xlsx" 
# library("rJava")
# if(!require(xlsx)){
#   install.packages("xlsx")
#   library(xlsx)
# }
library(xlsx)

#get package to make maps - see this website: http://www.molecularecologist.com/2012/09/making-maps-with-r/
#install.packages("mapdata")
library(mapdata)
#get package to make maps - see this website: http://www.molecularecologist.com/2012/09/making-maps-with-r/
#install.packages("maps")
library(maps)
# #get package for shapefiles see this website: http://www.molecularecologist.com/2012/09/making-maps-with-r/
# install.packages(maptools)
# library(maptools)  #for shapefiles
# #get package for adding pies on the map
#install.packages("mapplots")
library(mapplots)
# devtools::install_github("davidgohel/ReporteRs")
# devtools::install_github("davidgohel/officer")
# if(!require(officer)){
#   install.packages("officer")
#   library(officer)
# }
library(officer)

# if(!require(splitstackshape)){
#   install.packages("splitstackshape")
#   library(splitstackshape)
# }
library(splitstackshape)
#install.packages("tableHTML")
# #https://cran.r-project.org/web/packages/tableHTML/vignettes/tableHTML.html
# if(!require(tableHTML)){
#   install.packages("tableHTML")
#   library(tableHTML)
# }
library(tableHTML)
# install package if required
# if(!require(envDocument)){
#   install.packages("envDocument")
#   
# }
library(envDocument)
# define working diretory

#wd00 ="/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6/Analysis_MONIS5_to_6"
wd00 <- getwd()
#wd00 <- "/Users/steenknudsen/Documents/Documents/NIVA_Ansaettelse_2020/NOVANA_proever_2018_2019/"
# define out directory
wd05 <- "output05_color_tables_from_Rcode_for_MST2017_2022samples"
#paste dirs together
wd00_wd05 <- paste0(wd00,"/",wd05)
#Delete any previous versions of the output directory
unlink(wd00_wd05, recursive=TRUE)
#Create a directory to put resulting output files in
dir.create(wd00_wd05)
# # set working directory
# setwd(wd00)
# #getwd()
wd01 <- paste0(wd00,"/data")
#define filename for list of assays used
infl01 <- "lst_assays_MO5_2021apr.csv"
#paste together path and file name
pth_fl01 <- paste(wd01,"/",infl01,sep="")
#read excel with species names
df_spcas01 <-as.data.frame(read.csv(pth_fl01, header = TRUE, 
                                   sep = ",",quote = "\"",
                                   dec = ".", fill = TRUE, 
                                   comment.char = "", 
                                   stringsAsFactors = FALSE))
scpnmames <- df_spcas01
# the list of assays incl an old species name
#grep for this name and make a separate data frame
df_Psever <- scpnmames[grepl("Pse.*serruculata",scpnmames$Lat_Species), ]
#then replace in this df
df_Psever$species <- gsub("serruculata","verruculosa",df_Psever$species)
df_Psever$Lat_Species <- gsub("serruculata","verruculosa",df_Psever$Lat_Species)
#bind the subsetted df back to the main data frame
scpnmames02 <- rbind(scpnmames,df_Psever)
#identify assay numbers for species
assIDnoParcam <- unique(scpnmames02$Assay_ID[grepl("Paralit",scpnmames02$Lat_Species)])
assIDnoHomame <- unique(scpnmames02$Assay_ID[grepl("Homarus",scpnmames02$Lat_Species)])
assIDnoNeomel <- unique(scpnmames02$Assay_ID[grepl("Neogobius",scpnmames02$Lat_Species)])
#exclude these assays from the list of assays
# the assays to excl are those that performed suboptimal
scpnmames02 <- scpnmames02[!grepl("17",scpnmames02$Assay_ID),]
scpnmames02 <- scpnmames02[!grepl("15",scpnmames02$Assay_ID),]
scpnmames02 <- scpnmames02[!grepl("09A",scpnmames02$Assay_ID),]
#replace for the old species name for 'Crassostrea gigas'
scpnmames02$Lat_Species <- gsub("Magallana","Crassostrea",scpnmames02$Lat_Species)
scpnmames02$Genus <- gsub("Magallana","Crassostrea",scpnmames02$Genus)

#replace the first read in data frame
scpnmames <- scpnmames02
df_spcas01 <- scpnmames
# # set the working directory 
# setwd(wd00)
# #getwd()
# define the input file directory
wd03 <- "output03_match_MONIS_conc_on_extraction_with_water_sample"
wd00_03 <- paste0(wd00,"/",wd03)
#make a list of the files in the directory - get the first element
# that matches the 'csv' bit -  make sure there only is one file
# that matches this criteria
lst.infl <- list.files(wd00_03)[grep("csv",list.files(wd00_03))]
infl03 <- lst.infl[grep("table03",lst.infl)]
infl03 <- paste0(wd00_03,"/",infl03)
#read in collection positions
df_MSTsc01 <-as.data.frame(read.csv(infl03,
                          header = TRUE, sep = ";", 
                          quote = "\"", dec = ".", 
                          fill = TRUE, comment.char = "", 
                          stringsAsFactors = FALSE))
# split text - see: https://stevencarlislewalker.wordpress.com/2013/02/13/remove-or-replace-everything-before-or-after-a-specified-character-in-r-strings/
# and concatenate text - see: https://stackoverflow.com/questions/7201341/how-can-2-strings-be-concatenated 
# to get 6 letter abbr of latin speciesnames
ls.abbr.spcnm <-  paste(
  substr(sub('\\ .*', '', df_spcas01$Lat_Species), 1, 3),
  substr(sub('.*\\ ', '', df_spcas01$Lat_Species), 1, 3),
  sep="."
)
#add back on to latin name dataframe
df_spcas01$abbr.nm <- ls.abbr.spcnm

#define input file  with merged qPCR results
wd04 <- "output04_stdcrv_plots_and_tables_from_Rcode_for_MST2017_2023_samples"
#define input file directory
wd00_04 <- paste0(wd00,"/",wd04)
#make a list of files in the 'wd00_04' directory 
lst.infl <- list.files(wd00_04)[grep("csv",list.files(wd00_04))]
infl04 <- lst.infl[grep("table05_1",lst.infl)]
infl04 <- paste0(wd00_04,"/",infl04)
#read csv with all merged mxpro results
smpls01 <- read.csv(infl04, header = TRUE, sep = ";", 
                    quote = "\"", dec = ".", fill = TRUE, 
                    comment.char = "", stringsAsFactors = FALSE)


# set working directory
wd00_05 <- paste0(wd00,"/",wd05)
#setwd (wd00_05)
#getwd()
#make a list of the files in the directory - get the first element
# that matches the 'csv' bit -  make sure there only is one file
# that matches this criteria
#define input file  with merged qPCR results
# check only part of dataframe that holds results for 'Mya arenaria'
mya <- smpls01[smpls01$speciesabbr=="Myaare",]
mya <- mya[c("smpltp","CtdRn")]


# define path for out out directory
wd00_wd05 <- paste(wd00,"/",wd05,sep="")
#Delete any previous versions of the output directory
unlink(wd00_wd05, recursive=TRUE)
#Create a directory to put resulting output files in
dir.create(wd00_wd05)
# set working directory - this will be the output directory
#setwd(wd00_wd05)
#getwd()
#copy the dataframe
smpls02 <- smpls01
#remove blanks
#NOTE!! This will remove all NTC's with "No Ct"
#smpls02<-na.omit(smpls02)
#remove "No Ct"
#smpls02<-smpls02[!grepl("NoCt", smpls02$Quantitycopies),]
#smpls02<-smpls02[!grepl("NaN", smpls02$Quantitycopies),]

#limit to only FAM colur observations
#smpls02 <- smpls02[(smpls02$Fluor=="FAM"),]
smpls02$Sample <- smpls02$smpltp
smpls02$Sample2 <- smpls02$smpltp
#grep not starting with letter,after MST number has been stripped
#then pad with 4 zeroes
smpls02$Sample[!grepl("^[A-Za-z]",gsub("MST","",smpls02$Sample))] <- paste("MST",stringr::str_pad(gsub("MST","",smpls02$Sample[!grepl("^[A-Za-z]",gsub("MST","",smpls02$Sample))]), 4, pad = "0"),sep="")
# replace MST and pad with zero
mstno01 <- gsub("MST","",smpls02$Sample2[!grepl("^[A-Za-z]",gsub("MST","",smpls02$Sample2))])
mstno01 <- gsub("-","",mstno01)
mstno01 <- paste("MST",stringr::str_pad((mstno01), 7, pad = "0"),sep="")
# ad back to column where the sample Name holds MST in the name
smpls02$Sample2[!grepl("^[A-Za-z]",gsub("MST","",smpls02$Sample2))] <- mstno01
#change x into numeric variable
#smpls02$CtdRn <- smpls02$Cq
#smpls02$Quantitycopies <- smpls02$StartingQuantitySQ
smpls02$CtdRn=as.numeric(as.character(smpls02$CtdRn))
smpls02$Quantitycopies=as.numeric(as.character(smpls02$Quantitycopies))

smpls02$Quantitycopies[is.na(smpls02$Quantitycopies)] <- 0
smpls02$CtdRn[is.na(smpls02$smpls02$CtdRn)] <- 0
#remove points '.' with gsub from text string
df_spcas01$six_lett_spec_abbrv <- gsub("\\.","",df_spcas01$abbr.nm)
#match between dataframes to add latin species names and DK common names
smpls02$Lat_Species <- df_spcas01$Lat_Species[match(smpls02$speciesabbr, df_spcas01$six_lett_spec_abbrv)]



#paste a new column based on variables separated by point
smpls02$Well.ssn.no1 <- paste(gsub("^(.*)-(.*)$","\\1",smpls02$Content), smpls02$ssn, "1",  sep=".")
#paste a new column based on variables separated by point
smpls02$Lat_Species.season <- paste(smpls02$Lat_Species, smpls02$ssn,  sep=".")

#get month year, by using gsub to substitute, while also retaining specific parts
#convert to numeric, and get matching month name
smpls02$MSTsmpl.ye<- as.numeric(gsub("^(.*)-(.*)-(.*)$","\\1",smpls02$Dato_inds))

# paste together the sample number and the month for collection
smpls02$MSTno.mnt <-  paste(smpls02$Sample2,smpls02$mnt,sep=".")
smpls02$MSTno.ye.mnt <-  paste(smpls02$Sample2,smpls02$yea,smpls02$mnt,sep=".")

#paste a new column based on variables separated by point
smpls02$welltp <- smpls02$WellType
smpls02$local.Welltype <- paste(smpls02$locnm, smpls02$welltp,  sep=".")

#get the unique smpl names for df_MSTsc01 and WellTypes
unHaWT <- unique(smpls02$local.Welltype)
# make a transparent color
transp_col <- rgb(0, 0, 0, 0)
#transp_col <- as.character("#FFFFFF")
HaWTnoNA <- addNA(unHaWT)
col.01<-as.numeric(as.factor(unHaWT))
#make a small dataframe w df_MSTsc01 and standards and numbers assigned, 
#use the col2hex in gplot pacakge to convert the 'red' color name to hex-color
col.02 <- col2hex(palette(rainbow(length(col.01))))
harbourcols <- cbind(unHaWT,col.01, col.02)

#replace the colour for the standard dilution sample type with the transparent colour
col.03<-replace(col.02, col.01==15, transp_col)
col.04 <- cbind(harbourcols,col.03)
colforharb <- as.data.frame(col.04)
#match to main data frame and add as new color
smpls02$col.06 <- colforharb$col.03[match(smpls02$local.Welltype, colforharb$unHaWT)]
#insert the transparent color for all matches with "NA.Standard"
smpls02$col.06[smpls02$local.Welltype=="NA.Standard"] <- transp_col
smpls02$col.06[smpls02$local.Welltype=="NA.Std"] <- transp_col
# check samples for Mya arenaria - I wanted to know if all triplactes are still included
mya <- smpls02[grepl("Mya",smpls02$Lat_Species.season),]
mya <- mya[c("Sample","Quantitycopies")]
#get the unique species names
latspecnm <- unique(smpls02$Lat_Species)
#remove any NAs
latspecnm <- unlist(as.list(na.omit(latspecnm)))
#match the assay number to the data frame with species
AIfps <- df_spcas01$Assay_ID[match(latspecnm, df_spcas01$Lat_Species)]
#pad with zeros to two characters
#see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
AIfps <-stringr::str_pad(AIfps, 2, pad = "0")
#make a new data frame with assay Id No and species
nlspnm <- data.frame(AIfps,latspecnm)
#reorder by the column 'AssayIDNo'
nlspnm<- nlspnm[order(nlspnm$AIfps),]
#make a list of numbers for the unique species
no.latspc <- seq(1:length(latspecnm))
#add a new column with no to use for appendix numbering
nlspnm <- cbind(nlspnm, no.latspc) 
#use the new order of latin species names for producing plots
latspecnm <- unique(nlspnm$latspecnm)
#define the two seasons
dk.seasons <- c("foraar","efteraar")
eng.seasons <- c("spring","autumn")
ab.eng.seasons <- c("spri","fall")
sblt <- letters[1:2]
tr.seasons <-  data.frame(dk.seasons, eng.seasons,ab.eng.seasons,sblt)
#tr.seasons
amp <- smpls02

###############################################################################################
# start -calculate copies per L of filtered water
###############################################################################################
#smpls02$MSTno_di

#copy the data frame
smpls02.1 <- smpls02
#set NA blanks to zero
smpls02.1$CtdRn[is.na(smpls02.1$CtdRn)] <- 0
smpls02.1$Quantitycopies[is.na(smpls02.1$Quantitycopies)] <- 0
#make sure numbers are numbers
smpls02.1$Quantitycopies <- as.numeric(as.character(smpls02.1$Quantitycopies))
smpls02.1$volfilt_mL <- as.numeric(as.character(smpls02.1$volfilt_mL))
#add column with copies per Liter of filtered water
#Ae = (Cqpcr /Fe) /Vwf. 
#’Ae’ number of  eDNA-copies per volumen filtered water, 
#’Cqpcr’ number of copies detected in the qPCR-well, #smpls02.1$meanQuantitycopies 
#’Fe’ the ratio of the eluted extrated filtrate used in a qPCR-well #5/350
# For the MONIS5 project I eluted in 350 uL and used 3 uL template per well
# For the MONIS6 project I eluted in 400 uL and used 3 uL template per well

#’Vwf’ is volumen of seawater filtered. #smpls02.1$volfilt_mL
#these volumes will change depending on how the project
# is carried out. The values used here are specific for the
# MONIS projects
tmpl.volused <- as.numeric(gsub("uL","",smpls02.1$templvol))
#volume in uL used in qpcr reaction
uLpq <- tmpl.volused
#volume in uL used for elution
elvol <- 400
elvol<- smpls02.1$Elueringvolumen.AE.buffer.uL
#per mL
smpls02.1$copies_per_mLwater <- (smpls02.1$Quantitycopies/(uLpq/elvol))/smpls02.1$volfilt_mL
#per Liter
smpls02.1$copies_per_Lwater <- smpls02.1$copies_per_mLwater*1000
#replace nas with zeros
smpls02.1$copies_per_Lwater[is.na(smpls02.1$copies_per_Lwater)]<-0
#add one to be able to do logarithmic scales
smpls02.1$copies_per_Lwater_plone<- smpls02.1$copies_per_Lwater+1
#take log10 to all copies
smpls02.1$log.10_copies_L <- log10(smpls02.1$copies_per_Lwater_plone)

#substitute and paste to get new columns to subset by later on
#smpls02.1$smpltp <- gsub("^(.*)-.*$","\\1",smpls02.1$Content)
smpls02$Content <- smpls02$smpltp
smpls02.1$Content <- smpls02.1$smpltp
#smpls02.1$WellType #<- smpls02.1$smpltp
# make column with latin species name, the qpcr number and the plate number
smpls02.1$gen_specnm.qpcrno.pltn <- paste(smpls02.1$Lat_Species,smpls02.1$qpcrno,smpls02.1$plateno,sep=".")

#subset data frame to only matched values in column
smpls02.2 <- smpls02.1[which(smpls02.1$WellType=="Standard"),]
#remove samples that did not amplify - i.e. returned NaN for Cq
#smpls02.2 <- smpls02.2[!is.na(smpls02.2$Cq),] # use this line for the Biorad Machine
#smpls02.2 <- smpls02.2[!is.na(smpls02.2$CtdRn),] # use this line for the MxProMachine
#smpls02.2 <- smpls02.2[!smpls02.2$CtdRn==0,]
smpls02.2 <- smpls02.2[!smpls02.2$Quantitycopies==0,]
#get the LOD for each qPCR run
#use the function aggregate to get the minimum value for a group
lodtable1 <- aggregate(smpls02.2[, "Quantitycopies"], list(smpls02.2$gen_specnm.qpcrno.pltn, smpls02.2$welltp), min)
#subset this table by group
#lodtable2 <- lodtable1[ which(lodtable1$Group.2=="Std"), ] # use this line for the Biorad Machine
lodtable2 <- lodtable1[ which(grepl("Standard",lodtable1$Group.2)), ] # use this line for the MxPro Machine
#lodtable2[grepl("Mya",lodtable2$Group.1),]
#rename the column names - last column holds 'Limit Of Detection' per
# qPCR run
colnames(lodtable2) <- c("gen_specnm.qpcrno.pltn","WellT","LOD")
#subset the data frame to only incl standard dilutions
# this time from the data frame that includes the 'failed' 
# standard - i.e the standards that di not amplify
fg <- smpls02.1[ which(smpls02.1$smpltp=="Std"), ]
fg <- smpls02.1[ which(smpls02.1$WellType=="Standard"), ]
#count the number of replicate used per plate
#see this webpage: https://www.miskatonic.org/2012/09/24/counting-and-aggregating-r/
#and this webpage: https://stackoverflow.com/questions/9809166/count-number-of-rows-within-each-group
lodtb3 <- fg %>% dplyr::count(gen_specnm.qpcrno.pltn,smpltp)
#lodtb3[grepl("Neogo",lodtb3$gen_specnm.qpcrno.pltn),]
# paste together  species name and plate number and welltype content
lodtable2$nrepl <- lodtb3$n[match(lodtable2$gen_specnm.qpcrno.pltn, lodtb3$gen_specnm.qpcrno.pltn)]
#__________________________________________________________________
#Now identify LOQ for each qPCR run
#limit the dataframe to only well type that equals standard
oc <- smpls02.2[(smpls02.2$WellType=='Std'),] #use with BioRad
oc <- smpls02.2[(smpls02.2$WellType=='Standard'),] #use with MxPro
#exclude samples in the 'Std' Standard curve that did no amplify
#oc <- oc[!is.na(as.numeric(as.character(oc$Cq))),] #use with BioRad
oc <- oc[!is.na(as.numeric(as.character(oc$CtdRn))),] #use with MxPro
oc <- oc[oc$CtdRn!=0,]
#add a new column that merges two columns for species and qPCR plate no
oc$Content.spc.pltn <- paste(oc$Content, oc$gen_specnm.qpcrno.pltn,  sep=".")

#count the occurences of dilution steps - i.e. the number of succesful replicates
#see this webpage: https://www.miskatonic.org/2012/09/24/counting-and-aggregating-r/
#and this webpage: https://stackoverflow.com/questions/9809166/count-number-of-rows-within-each-group
od <- oc %>% dplyr::count(Content, Content.spc.pltn)
#turn this into a dataframe
oe<-as.data.frame(od)
#The 'n' column now holds the count of wells at this specific
# dilution level that successfully amplified
# lowest number below the number of total replicates will be
# the level just below the LOQ
#add a new column that merges two columns
#match the dilution step to the number of occurences -i.e. match between the two dataframes
no.occ <- oe$n[match(oc$Content.spc.pltn,oe$Content.spc.pltn)]
#add this column with counted occurences to the limited dataframe
og <- cbind.data.frame(oc,no.occ)
# check out the results for one species
og$Content.spc.pltn
#get the number of replicates used
og$nrepl <- lodtable2$nrepl[match(og$gen_specnm.qpcrno.pltn,lodtable2$gen_specnm.qpcrno.pltn)]
#exlude all observations where 
#less than '3'number of replicates' amplified
oh<-og[(og$no.occ>=og$nrepl),]
oh$Quantitycopies <- as.numeric(oh$Quantitycopies)

#get the lowest dilution step that succesfully amplified on all 3 repliactes
#use aggregate to get the minimum for each
loqtable1 <- aggregate(oh[, "Quantitycopies"], list(oh$gen_specnm.qpcrno.pltn), min)
#change the column names
colnames(loqtable1) <- c("spc.plt","LOQ")
#copy the LOD table and add the corresponding LOQ values
loq.lod.table <- lodtable2
loq.lod.table$LOQ <- loqtable1$LOQ[match(lodtable2$gen_specnm.qpcrno.pltn,loqtable1$spc.plt)]
#append limit of quantification back to main data frame
smpls02.1$LOQ <- loq.lod.table$LOQ[match(smpls02.1$gen_specnm.qpcrno.pltn,loq.lod.table$gen_specnm.qpcrno.pltn) ]
#append limit of detection back to main data frame
smpls02.1$LOD <- loq.lod.table$LOD[match(smpls02.1$gen_specnm.qpcrno.pltn,loq.lod.table$gen_specnm.qpcrno.pltn) ]
#add an empty column with just NAs
smpls02.1[,"eDNA_eval"] <- NA
#replace in the empty column, the order is important, as you otherwise will end up with the last evaluations
smpls02.1$eDNA_eval[smpls02.1$Quantitycopies>=smpls02.1$LOQ] <- "aboveLOQ"
smpls02.1$eDNA_eval[smpls02.1$Quantitycopies<=smpls02.1$LOQ] <- "AbLOD_BeLOQ"
smpls02.1$eDNA_eval[smpls02.1$Quantitycopies<=smpls02.1$LOD & !smpls02.1$Quantitycopies==0] <- "belowLOD"
smpls02.1$eDNA_eval[smpls02.1$Quantitycopies==0] <- "NoCt"
# check if all lines have been assigned an evalution
smpls02.1[is.na(smpls02.1$eDNA_eval),]
#smpls02.1$eDNA_eval
######################################################################################################
# Get mean for each set of 3 technical qPCR replicates per species per season per harbour
######################################################################################################
# copy the data frame into a new object
smpls04 <- smpls02.1
# use dplyr
library(dplyr)
# to select columns to work on, and to group by, and to count up for 
# see this examples:https://stackoverflow.com/questions/22767893/count-number-of-rows-by-group-using-dplyr
#Notice that all these columns are not required for the the dplyr::summarise part
# Simply just selecting the 'qpcr.rplno' should be sufficient to summarise
# over the qpcr wells, but I need all the additional information for all the
# later steps in the next pieces of code
smpls04.1 <- smpls04 %>% 
  dplyr::select(qpcrno,qpcr.rplno, lokalitet_vanda,volfilt_mL,declat,declon,plateno,Lat_Species,mnt,smpltp,smplNoSTEX ,welltp,Genus,species, yea,ssn, MST.n_Dinds,eDNA_eval) %>% 
  dplyr::group_by(qpcrno,qpcr.rplno, lokalitet_vanda,volfilt_mL,declat,declon,plateno,Lat_Species,mnt,smpltp,smplNoSTEX ,welltp,Genus,species, yea,ssn, MST.n_Dinds,eDNA_eval) %>% 
  dplyr::summarise(n.eval = n())

# use tidyr::pivot_wider to get the eDNA_eval category columns
# notice that I use the values_fill = 0 option, to get 0 instead of NAs
smpls04.1 <- smpls04.1 %>% tidyr::pivot_wider(names_from = "eDNA_eval",
                                values_from = "n.eval",
                                values_fill = 0)
# check the rows that are missing a 'MST.n_Dinds'
smpls04.2 <- smpls04.1[is.na(smpls04.1$MST.n_Dinds),]
# just to check that the 'smplNoSTEX' that does not have additional information
# are the extraction blanks
unique(smpls04.2$smplNoSTEX)
# paste together to get numbered values for categories 
smpls04.1$freq_repl_eval <- paste( smpls04.1$NoCt,
                                   smpls04.1$belowLOD,
                                   smpls04.1$AbLOD_BeLOQ,
                                   smpls04.1$aboveLOQ,
                                   sep="/")
#add quotation marks in the end and the beginning to make the 
# color search in the html table creation work
smpls04.1$freq_repl_eval <- paste0("'",smpls04.1$freq_repl_eval,"'")


# ______________________________________________________________________________

smpls08 <- smpls04.1
#then paste together columns to be able to match correctly back to 
# the aggregated data frame
smpls08$gen_specnm.qpcrno.pltn.MSTno.mnt <- paste(smpls08$Lat_Species,
                                                  smpls08$qpcrno,
                                                  smpls08$pltno,
                                                  smpls08$MSTno.mnt,sep=".")


# write a copy of the table
write.csv(smpls08,paste0(wd00_wd05,"/table07_out_smpls08.csv"))

####################################################################################
# Start Appendix C
####################################################################################

##########################################################################################
#
# Make a table that looks somewhat similar to Table 5 presented by :
#  Li, J, Hatton‐Ellis, TW, Lawson Handley, L‐J, et al. Ground‐truthing of a fish‐based environmental DNA metabarcoding method for assessing the quality of lakes. J Appl Ecol. 2019; 56: 1232– 1244. https://doi.org/10.1111/1365-2664.13352 
#
# Make a table for each year sampled
##########################################################################################
#copy the data frame
df_MO5_10 <- smpls08
# paste to get genus and species name and sampling year
df_MO5_10$gen_specnm.year_inds <- paste(df_MO5_10$Lat_Species,df_MO5_10$yea,sep="_")
# get MST number only by substituting and keeping first element before '.'
# #_
# df_MO5_10$MSTno <- gsub("^(.*)\\.(.*)\\.(.*)$","\\1",df_MO5_10$MSTno.mnt)
# df_MO5_10$month_inds2 <- gsub("^(.*)\\.(.*)\\.(.*)$","\\3",df_MO5_10$MSTno.mnt)
# df_MO5_10$smpltp <- df_MO5_10$Content
# #_
df_MO5_10$season_cat <- df_MO5_10$ssn
df_MO5_10$yer <- df_MO5_10$yea
df_MO5_10$month_inds2 <- df_MO5_10$mnt
#View(df_MO5_10$MST.n_Dinds)
#define the columns to keep 
keeps <- c(#"smpltp",
           "gen_specnm.year_inds",
           "freq_repl_eval",
           #"eDNA_eval",
           "yer",
           "plateno",
           "month_inds2",
           #"MSTno",
           "MST.n_Dinds",
           "season_cat")
#keep only selected columns
df_MO5_11 <- df_MO5_10[keeps]
#keep unique rows only
df_MO5_12 <- df_MO5_11 %>% dplyr::distinct(MST.n_Dinds, 
                    gen_specnm.year_inds, 
                    freq_repl_eval, .keep_all = TRUE)
# Sort by vector name [smpltp] then [gen_specnm.year_inds] 
# https://chartio.com/resources/tutorials/how-to-sort-a-data-frame-by-multiple-columns-in-r/
df_MO5_11 <- df_MO5_11[
  with(df_MO5_11, order(MST.n_Dinds, gen_specnm.year_inds)),]
# Sort by vector name [smpltp] then [gen_specnm.year_inds] 
# https://chartio.com/resources/tutorials/how-to-sort-a-data-frame-by-multiple-columns-in-r/
df_MO5_12 <- df_MO5_12[
  with(df_MO5_12, order(MST.n_Dinds, gen_specnm.year_inds)),]
# make a new column that fuses MST sample number together with month
df_MO5_12$smpltp.month <- paste(df_MO5_12$MST.n_Dinds,".",
                                df_MO5_12$month_inds2,sep="")
#split column by delimiter, and turn in to data frame # https://www.rdocumentation.org/packages/splitstackshape/versions/1.4.8/topics/cSplit
df_MO5_13 <- as.data.frame(splitstackshape::cSplit(df_MO5_12,
              "gen_specnm.year_inds", sep = "."))
# copy the column with the year sampled
df_MO5_13$yrs_smpl <- df_MO5_13$yer
#Rename specific column # see :  https://stackoverflow.com/questions/7531868/how-to-rename-a-single-column-in-a-data-frame
#remove rows where there is no sampling year
df_MO5_13 <- df_MO5_13[!is.na(df_MO5_13$yrs_smpl),]
df_MO5_13 <- df_MO5_13[(df_MO5_13$yrs_smpl!=0),]
#remove rows where there is no sampling season
df_MO5_13 <- df_MO5_13[!is.na(df_MO5_13$season_cat),]
# remove rows that are not measurements on water samples

#get year
year.dv <- unique(df_MO5_13$yrs_smpl)
#use the year listed in a vector previously
yrs <- year.dv[order(year.dv)]
# use just a single year to start with for testing the loop
#get season category
cat.of.seasons <- unique(df_MO5_13$season_cat)
cat.of.seasons <- cat.of.seasons[!is.na(cat.of.seasons)]

#yrs <- 2017
#cat.of.seasons <- "fall" 
#loop over years sampled -  to produce individual tables per year sampled
for (yr_smpl in yrs){
  print(yr_smpl)
  #}
  #subset based on variable values - only retain rows where the column that match the criteria 
  sbs.MO13y <- df_MO5_13[ which(df_MO5_13$yrs_smpl==yr_smpl), ]
  #to try out the loop assign only one category
  #categories.of.seasons <- "season_1"
  #loop over the seasons
  for (season in cat.of.seasons){
    print(season)
    #}
    #subset based on variable values - only retain rows where the column that match the criteria 
    sbs.MO14ym <- sbs.MO13y[ which(sbs.MO13y$season_cat==season), ]
    
    #head(sbs.MO14ym,4)
    #define the columns to keep 
    keeps <- c(#"gen_specnm",
               "gen_specnm.year_inds_1",
               "freq_repl_eval",
               "MST.n_Dinds",
               "smpltp.month",
               "season_cat")
    #keep only selected columns
    sbs.MO15ym <- sbs.MO14ym[keeps]
    sbs.MO15ym$mnth <- gsub("^(.*)\\.(.*)$","\\2",sbs.MO15ym$smpltp.month)
    #sbs.MO15ym$MSTno.mnth <- paste(sbs.MO15ym$MSTno,sbs.MO15ym$mnth,sep=".")
    
    sbs.MO15ym <- subset (sbs.MO15ym, select = -c(season_cat,
                                                  smpltp.month,
                                                  #MST.n_Dinds,
                                                  mnth))
    # copy column
    sbs.MO15ym$MSTno <- sbs.MO15ym$smpltp.month
    sbs.MO15ym$MSTno <- sbs.MO15ym$MST.n_Dinds
    # and remove columns no longer needed
    sbs.MO15ym$smpltp.month <- NULL
    sbs.MO15ym$mnth <- NULL
    sbs.MO15ym$MST.n_Dinds <- NULL
    #View(sbs.MO15ym)
    #reshape the data frame to have smpls for columns
    sbs.MO16ym <- reshape(sbs.MO15ym, idvar = "gen_specnm.year_inds_1", 
                          timevar = "MSTno", 
                          direction = "wide")
    # change the column header
    colnames(sbs.MO16ym)[1] <- "genus_species_aar"
    #Replace characters in column names gsub : 
    # https://stackoverflow.com/questions/39670918/replace-characters-in-column-names-gsub
    names(sbs.MO16ym) <- gsub(x = names(sbs.MO16ym), 
                              pattern = "freq_repl_eval\\.", 
                              replacement = "")  
    #count the number of columns
    nc.MO16 <- ncol(sbs.MO16ym)
    #use match to match the season with a data frame and get the name for the season
    spcfc_seaon_name <- tr.seasons$eng.seasons[match(season, tr.seasons$ab.eng.seasons)]
    spcfc_seaon_name <- as.character(spcfc_seaon_name)
    tr.seasons$nfss <- seq(1,nrow(tr.seasons))
    #use match to match the season with a data frame and get the category number for the season
    spcfc_seaon_no <- tr.seasons$nfss[match(season, tr.seasons$ab.eng.seasons)]
    spcfc_seaon_no <- as.numeric(spcfc_seaon_no)
    
    if (dim(sbs.MO16ym)[1] == 0) {
      print(paste("data frame for",spcfc_seaon_name,yr_smpl,"is empty", sep=" "))
      sbs.MO16ym <- as.data.frame(rbind(c("MST_smpl01","MST_smpl02"),c("no_data","sampled")))
    } else {
    
    #replace values in entire data frame
    sbs.MO16ym[sbs.MO16ym=="NoCt"]<-"NoCq"
    sbs.MO16ym[sbs.MO16ym=="belowLOD"]<-"bLOD" 
    sbs.MO16ym[sbs.MO16ym=="AbLOD_BeLOQ"]<-"aLODbLOQ"
    sbs.MO16ym[sbs.MO16ym=="1aboveLOQ"]<-"1aLOQ"
    sbs.MO16ym[sbs.MO16ym=="3aboveLOQ"]<-"3aLOQ"
    #get the unique years sampled
    yrs <- yr_smpl
    # Remove columns from dataframe where ALL values are NA # https://stackoverflow.com/questions/2643939/remove-columns-from-dataframe-where-all-values-are-na
    sbs.MO16ym <- sbs.MO16ym[,colSums(is.na(sbs.MO16ym))<nrow(sbs.MO16ym)]
    # delete the first column from the data frame
    #sbs.MO_df[,1] <- NULL
    #order data frame alphabetically by column w species
    sbs.MO16ym <- sbs.MO16ym[order(sbs.MO16ym$genus_species_aar, decreasing = F), ] 
    
    sbs.MO <- sbs.MO16ym
    # #https://cran.r-project.org/web/packages/tableHTML/vignettes/tableHTML.html
    # if(!require(tableHTML)){
    #   install.packages("tableHTML")
    #   library(tableHTML)
    # }
    library(tableHTML)
    
    #try the tableHTML with no border
    tableHTML <- sbs.MO %>% 
      tableHTML(border = 0,rownames=F) 
   #get the number of columns
    nc.MO <- ncol(sbs.MO)
    # get unique values in the data frame from a selected range of columns # see : https://stackoverflow.com/questions/40003028/extracting-unique-values-from-data-frame-using-r
    unqval.MO16 <- unique(unlist((sbs.MO16ym)[,2:nc.MO]))
    #remove NA from vector - this will exclde the non-evaluated categories
    unqval.MO16 <- unqval.MO16[!is.na(unqval.MO16)]
    
    # With the grep function in R the different elements in 'freq_repl_eval' 
    # can be categorized for 
    # identification later on in the preparation of the html table that is to show
    # the colored categories for eDNA levels detected.
    # The coding for the elements in 'freq_repl_eval' are:
    # ' no of replicates with no Ct/no of replicates below LOD /no of replicates above LOD but below LOQ / no of replicates above LOQ ' 
    # the color coding for these elements in 'freq_repl_eval' are
    # ' white/yellow /orange / red or black '
    # red if a minimum if 1 replicate is above LOQ (disregarding if any lower levels are detected)
    # black if all replicates are above LOQ (this will automatically equal all lower categories being zero)
    # To try out the grep function, I have here below tried grepping for different elements in the list
    #grep for elements that begin with '0
    grep("^'0", unqval.MO16, value=T)
    #grep for elements that end with 0'
    grep("0'$", unqval.MO16, value=T)
    # grep for elements starting '0/0/0/ and then 1 or any higher number
    b_cat <- grep("'0/0/0/[1-9]+", unqval.MO16, value=T) # this will equal the black category with all replicates amplyfiying
    # then grep for all elements with not zero in the last category
    br_cat <- grep("'[0-9]+/[0-9]+/[0-9]+/[1-9]+", unqval.MO16, value=T) # this will equal both the red and balck category
    # find the difference between these two vectors - i.e. subtract the black category from the fused red-black category
    r_cat <- setdiff(br_cat,b_cat)
    # grep for all elements with not zero in the first category
    w_cat <- grep("'[1-9]+/0/0/0'", unqval.MO16, value=T) # this will equal all replicates not amplifying - i.e. this will equal the white category
    
    w_cat <- c(c(w_cat),"NA")
    # find the difference between two vectors
    wyo_cat <- setdiff(unqval.MO16,br_cat)
    # find the difference between two vectors
    yo_cat <- setdiff(wyo_cat,w_cat)
    # grep for all elements with not zero in the last and second last category
    y_cat <- grep("'[0-9]+/[1-9]+/0/0'", unqval.MO16, value=T)
    # grep for all elements with not zero in thelast category
    o_cat <- grep("'[0-9]+/[0-9]+/[1-9]+/0'", unqval.MO16, value=T)
    #these categories are used below to identify the color coding in the html tables
    
    sbs.MO_df <- sbs.MO16ym
    # #https://cran.r-project.org/web/packages/tableHTML/vignettes/tableHTML.html
    # if(!require(tableHTML)){
    #   install.packages("tableHTML")
    #   
    # }
    library(tableHTML)
    #try the tableHTML with no border
    tableHTML <- sbs.MO_df %>% 
      tableHTML(border = 0,rownames = F) 
    #count the number of columns in the dataframe
    l.s.MO <- length(sbs.MO_df)
    #get unique cell values in dataframe : see : http://r.789695.n4.nabble.com/Retrieve-distinct-values-within-a-whole-data-frame-td1460205.html
    #apart from the first column
    unique(unlist(sbs.MO_df[2:l.s.MO]))
    #make lists of the words in the cells to color using the 'add_css_conditional_column' function
    #class(eDNA.lvl01)
    eDNA.lvl01 <- w_cat #c("/0/0/0'") #white
    eDNA.lvl02 <- y_cat #c("'0/") #yellow
    eDNA.lvl03 <- o_cat #c("aLODbLOQ") #orange
    eDNA.lvl04 <- r_cat #c("/1'") #red
    eDNA.lvl05 <- b_cat #c("'0/0/0/") #black
    #place the data frame in a tableHTML object
    tableHTML <- sbs.MO_df %>% 
      tableHTML(rownames = F)
    # for eDNA.lvl01 <- c("NoCq") #white
    words <- eDNA.lvl01 #<- c("NoCq") #white
    col.f.cell <- "white"
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color"),
                                              c(col.f.cell)))
    }
    # for eDNA.lvl02 <- c("beLOD") #yellow
    words <- eDNA.lvl02
    col.f.cell <- "yellow"
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color"),
                                              c(col.f.cell)))
    }
    # for eDNA.lvl03 <- c("abLODbeLOQ") #orange
    words <- eDNA.lvl03
    col.f.cell <- "orange"
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color"),
                                              c(col.f.cell)))
    }
    # for eDNA.lvl04 <- c("1abLOQ") #red
    words <- eDNA.lvl04
    col.f.cell <- "red"
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color"),
                                              c(col.f.cell)))
    }
    # for eDNA.lvl05 <- c("3abLOQ") #black
    words <- eDNA.lvl05
    col.f.cell <- "black" # use this color for the cell
    col.f.font <- "white" #use this color for the font
    for (word in words) {
      tableHTML <- tableHTML %>% 
        add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                   conditional = "contains",
                                   value = word,
                                   css = list(c("background-color","color"),
                                              c(col.f.cell,col.f.font)))
    }
    # copy html table to a new object
    t.HTML17 <- tableHTML
      
    df_t.HT17 <- sbs.MO_df
    #pad with zeros to two characters
    #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
    no.e4 <-stringr::str_pad(spcfc_seaon_no, 2, pad = "0")
    #and to export in a file
    write_tableHTML(t.HTML17, file = paste(wd00_05,"/table08_App_C_",
                                           yr_smpl,
                                           "_",
                                           no.e4,
                                           "_",
                                           spcfc_seaon_name,
                                           "_table_eDNA_evalu.html",
                                           sep=""))
    # write out a csv file
    write.csv(df_t.HT17,file=paste(wd00_05,"/table09_App_C_",
                                   yr_smpl,
                                   "_",
                                   no.e4,
                                   "_",
                                   spcfc_seaon_name,
                                   "_table_eDNA_evalu.csv",
                                   sep=""))
    # write to an excel file
    # http://www.sthda.com/english/wiki/writing-data-from-r-to-excel-files-xls-xlsx
    write.xlsx(df_t.HT17, 
               file=paste(wd00_wd05,"/table10_App_C_",
                          yr_smpl,
                          "_",
                          no.e4,
                          "_",
                          spcfc_seaon_name,
                          "_table_eDNA_evalu.xlsx",
                          sep=""), sheetName = "Sheet1", 
               col.names = TRUE, row.names = TRUE, append = FALSE)
     
    # end if test ti check if data frame is empty   
    }
    #end loop over seasons
  }
  #end loop over years
}

####################################################################################
# End Appendix C
####################################################################################

