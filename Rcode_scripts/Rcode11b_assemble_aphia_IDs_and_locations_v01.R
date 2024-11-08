#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# path for externmal dir for data
finpdir <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6/data_prepared_on_remote_server/all_aphiaIDs_and_locations"

#wd00 <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6/Results_from_ABI7500_for_MONIS6_species_specific_assays"
wd00 <- getwd()
setwd (wd00)
#define directory with output flies
wd01 <- "output11b_assembled_aphia_ID_and_location_files"
# define full path for input directory
outdir01 <- paste(wd00,wd01, sep="/")
# remove previous versions of the 'outdir01'
unlink(outdir01, force = T, recursive = T)
# create the 'outdir01' in a new version
dir.create(outdir01)

#list all files in external wd - all the xls-files for which you want to 
# prepare plots from 
ls.fl01 <- list.files(finpdir,full.names = T)

ApSe.fls <- ls.fl01[grepl("\\/ApSe_[0-9]{3}\\.csv",ls.fl01)]
location.fls <- ls.fl01[grepl("\\/location_for_ApSe_[0-9]{3}\\.csv",ls.fl01)]

#
library("dplyr") # Load dplyr package
library("plyr")  # Load plyr package
library("readr")  # Load readr package
# read all files in the list
# make an empty list
lst_ApSe <- list()
# make a sequence of numbers from 1 to the length of the list
sqAf <- seq(1,length(ApSe.fls),1)
# iterate over the sequence
for (f in sqAf){
  #f<- 1 
  # read in the csv-file
  df_Af <- read_delim(ApSe.fls[f],delim=";") %>% # Store all files in list
  bind_rows                       # Combine data sets into one data set
  # add the data set to the list
  lst_ApSe[[f]] <- df_Af
  }
# combine all data sets into one data set
df_ApSe <- do.call(rbind.fill, lst_ApSe)  
# make sure that the AphiaID is numeric
df_ApSe$AphiaID <- as.numeric(df_ApSe$AphiaID)
# remove rows with NA in the AphiaID column
df_ApSe <- df_ApSe[!is.na(df_ApSe$AphiaID),]
# remove duplicated rows for AphiaID  
df_ApSe <- df_ApSe[!duplicated(df_ApSe$AphiaID),]
#
nrow(df_ApSe)

# make an empty list
lst_locs <- list()
# make a sequence of numbers from 1 to the length of the list
sqlf <- seq(1,length(location.fls),1)
# iterate over the sequence
for (f in sqlf){
  #f<- 1 
  # read in the csv-file
  df_lf <- read_delim(location.fls[f],delim=";") %>% # Store all files in list
    bind_rows                       # Combine data sets into one data set
  # add the data set to the list
  lst_locs[[f]] <- df_lf
}
# combine all data sets into one data set
df_locat <- do.call(rbind.fill, lst_locs)  
# make sure that the AphiaID is numeric
df_locat$AphiaID <- as.numeric(df_locat$AphiaID)
# remove rows with NA in the AphiaID column
df_locat <- df_locat[!is.na(df_locat$AphiaID),]
# remove duplicated rows for AphiaID  
df_locat <- df_locat[!duplicated(df_locat$AphiaID),]
# order the data frames
df_ApSe <- df_ApSe %>% dplyr::arrange( kingdom, phylum, 
                            class, order, family, genus,scientificname)
df_locat <- df_locat %>% dplyr::arrange(Kingdom,Phylum, 
                                       Class, Order, Family, Genus,ScientificName)
# remove rows with status not equal to "accepted"
df_ApSe<-df_ApSe[(df_ApSe$status=="accepted"),]
# write out the data sets
write_delim(df_ApSe, paste(outdir01,"/df_ApSe.csv",sep=""),delim=";")
write_delim(df_locat, paste(outdir01,"/df_locat.csv",sep=""),delim=";")

#View(df_locat)
#View(df_ApSe)

nrow(df_locat)
nrow(df_ApSe)

