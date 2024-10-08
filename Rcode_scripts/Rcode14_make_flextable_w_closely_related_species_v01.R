#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# To be able to install 'flextable', it is required to have  'gdtools'
# and 'gdtools needs 'systemfonts' to be installed. Install the packages
# in this order
#install.packages("systemfonts")
#install.packages("gdtools")
#install.packages("flextable")


library(systemfonts)
library(gdtools)
library(flextable)
library(dplyr)
library(tibble)
library(tidyverse)
library(xlsx)

# set working directory
#wd00 <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6/Results_from_ABI7500_for_MONIS6_species_specific_assays"
wd00 <- getwd()
setwd(wd00)
#define directory with output flies
wd14 <- "output14_sisterspecies_table"
# define full path for input directory
outdir10 <- paste(wd00,wd14, sep="/")
# remove previous versions of the 'outdir10'
unlink(outdir10, force = T, recursive = T)
# create the 'outdir10' in a new version
dir.create(outdir10)
# define directory with input files to read in
wddata <- "data"
#paste together a path for the directory that holds the inputs files
wd00_wddata <- paste(wd00,wddata,sep="/")

#list all files in wd - all the xls-files for which you want to 
# prepare plots from 
ls.fl01 <- list.files(wd00_wddata)
#make a variable with the element you want to search for
id1 <- "xls"
#grep for this variable in the list -  see this example: 
# https://stackoverflow.com/questions/35880242/r-selecting-element-from-list
ls.fl01.xls <- ls.fl01[grep(paste0(id1), ls.fl01)]
ls.fl01.xls <- ls.fl01.xls[grep("TS_list", ls.fl01.xls)]


infl01 <- paste0(wd00_wddata,"/",ls.fl01.xls)
# read in xls spreadsheet tab with raw flourescense data from ABI7500
rdt <- readxl::read_xlsx(infl01,
                         sheet =  1,skip = 1)
# first row contains the table legend
tbl.lgnd <- as.character(names(unlist(rdt[1,]))[1])
# change the column nmames to the first row
colnames(rdt) <- rdt[1,]
# exclude the first row
rdt <- rdt[-1,]
# substitute in column names
colnames(rdt) <- gsub(" ","_",colnames(rdt))


# Now also read in the excel file with the list of target species
#list all files in wd - all the xls-files for which you want to 
# prepare plots from 
ls.fl01 <- list.files(wd00_wddata)
#make a variable with the element you want to search for
id1 <- "xls"
#grep for this variable in the list -  see this example: 
# https://stackoverflow.com/questions/35880242/r-selecting-element-from-list
ls.fl01.xls <- ls.fl01[grep(paste0(id1), ls.fl01)]
ls.fl01.xls <- ls.fl01.xls[grep("Artsprio", ls.fl01.xls)]
ls.fl01.xls <- ls.fl01.xls[grep("v05", ls.fl01.xls)]

infl01 <- paste0(wd00_wddata,"/",ls.fl01.xls)
# read in xls spreadsheet tab with raw flourescense data from ABI7500
df_ts <- readxl::read_xlsx(infl01,
                         sheet =  2)
# replace the column names
colnames(df_ts) <- df_ts[2,]
# exclude rows
df_ts <- df_ts[-c(1,2),]
# substitute in column names
colnames(df_ts) <- gsub(" ","_",colnames(df_ts))
# exclude row if 'Genus_og_species_og_author' is NA
df_ts <- df_ts[!(is.na(df_ts$Genus_og_species_og_author)),]

# split string to get genus and species
Gnsp.splt <- strsplit(df_ts$Genus_og_species_og_author, " ")
# get the genus and species
genus.Nm <- sapply(Gnsp.splt, "[[", 1)
species.Nm <- sapply(Gnsp.splt, "[[", 2)
# paste together genus and species to get the full name
gsNm <- paste0(genus.Nm," ", species.Nm)
# modify the column with DK_NIS to have 'TS' for the target species
# so that it matches that genus_species names in the table
# with the prioritized TS species
rdt$DK_NIS[(rdt$Latinsk_artsnavn  %in%  gsNm)] <- "TS"
# re order the rows by the columns
rdt <- rdt %>% dplyr::arrange(Phylum,
                              Klasse, 
                              Orden,
                              Familie,
                              Latinsk_artsnavn)
# exclude row if phylum is NA
rdt <- rdt[!(is.na(rdt$Phylum)),]
# get the column names
clNMs <- colnames(rdt)
# split string to get genus and species
Gnsp.splt <- strsplit(rdt$Latinsk_artsnavn, " ")
# get the genus and species
genus.Nm <- sapply(Gnsp.splt, "[[", 1)
species.Nm <- sapply(Gnsp.splt, "[[", 2)
# get the genus and species and authors
gnsp.Nm <- paste(genus.Nm,species.Nm,sep=" ")
# get unique abbreviations for TS NTS and for geographical regions
catDKNIS <- unique(rdt$DK_NIS)
catLGeoRg <- unique(rdt$Lokal_geografi)
catOGeoRg <- unique(rdt$overordnet_geografi)
# order the categories alphabetically
catDKNIS <- catDKNIS[order(catDKNIS)]
catLGeoRg <- catLGeoRg[order(catLGeoRg)]
catOGeoRg <- catOGeoRg[order(catOGeoRg)]

# 
abbrev_txtx_for_lgnd <- c("Forkortelserne for de forskellige havområder er:  
                            EurW: European waters (ERMS scope) ,
                            Hol: Holland, Natl: NAtl Ocean , 
                            NAtl: North,  Atlantic , 
                            NEAtl: North East Atlantic, 
                            MW Atl: North west Atlantic Ocean, 
                            NS: North Sea, 
                            NTS: Nontarget species, 
                            TS: Target species, 
                            
                            Ger: Germany, 
                            BS: Baltic Sea, 
                            UK: United Kingdom;  
                            Fra: France, 
                            Bel: Belgium, 
                            DK: Denmark, 
                            
                            DK NIS: 
                            Angiver om det eftersøgt art blandt de ikke hjemmehørende arter (dvs en TS) 
                            eller om det ikke er at betragte som en ikke-hjemmehørende art (dvs. NTS)")

#_______
#https://stackoverflow.com/questions/8753531/repeat-rows-of-a-data-frame-n-times#8753732
clNMstr <- c("Phylum",
             "Klasse", 
             "Orden",
             "Familie")

ttcls <- ncol(rdt)
rdt.tx <- rdt[clNMstr]
rdt.tx <- rdt.tx %>% dplyr::distinct(Phylum,Klasse, Orden, Familie)
idxN <- ncol(rdt.tx)
Fcls <- rdt.tx[,idxN]
nFrw <- nrow(Fcls)
clAAfill<- rep("AA",nFrw)
Fcls2 <- unlist(as.vector(Fcls))
Fcls2 <- paste0(clAAfill,"_",Fcls2)
ncF <-(ttcls-idxN)
Fcln <- do.call("cbind", replicate(ncF, Fcls2, simplify = FALSE))
Fcln <- as.data.frame(Fcln)
rdt.txfF <- rdt.tx[,1:(idxN)]
Frws <- cbind(rdt.txfF,Fcln)
ncol(Frws)
colnames(Frws) <- colnames(rdt)
family.rows <- Frws

ttcls <- ncol(rdt)
rdt.tx <- rdt[clNMstr]
rdt.tx <- rdt.tx %>% dplyr::distinct(Phylum,Klasse, Orden)
idxN <- ncol(rdt.tx)
Fcls <- rdt.tx[,idxN]
nFrw <- nrow(Fcls)
clAAfill<- rep("AA",nFrw)
Fcls2 <- unlist(as.vector(Fcls))
Fcls2 <- paste0(clAAfill,"_",Fcls2)
ncF <-(ttcls-idxN)
Fcln <- do.call("cbind", replicate(ncF, Fcls2, simplify = FALSE))
Fcln <- as.data.frame(Fcln)
rdt.txfF <- rdt.tx[,1:(idxN)]
Frws <- cbind(rdt.txfF,Fcln)
ncol(Frws)
colnames(Frws) <- colnames(rdt)
order.rows <- Frws


ttcls <- ncol(rdt)
rdt.tx <- rdt[clNMstr]
rdt.tx <- rdt.tx %>% dplyr::distinct(Phylum,Klasse)
idxN <- ncol(rdt.tx)
Fcls <- rdt.tx[,idxN]
nFrw <- nrow(Fcls)
clAAfill<- rep("AA",nFrw)
Fcls2 <- unlist(as.vector(Fcls))
Fcls2 <- paste0(clAAfill,"_",Fcls2)
ncF <-(ttcls-idxN)
Fcln <- do.call("cbind", replicate(ncF, Fcls2, simplify = FALSE))
Fcln <- as.data.frame(Fcln)
rdt.txfF <- rdt.tx[,1:(idxN)]
Frws <- cbind(rdt.txfF,Fcln)
ncol(Frws)
colnames(Frws) <- colnames(rdt)
class.rows <- Frws


ttcls <- ncol(rdt)
rdt.tx <- rdt[clNMstr]
rdt.tx <- rdt.tx %>% dplyr::distinct(Phylum)
idxN <- ncol(rdt.tx)
Fcls <- rdt.tx[,idxN]
nFrw <- nrow(Fcls)
clAAfill<- rep("AA",nFrw)
Fcls2 <- unlist(as.vector(Fcls))
Fcls2 <- paste0(clAAfill,"_",Fcls2)
ncF <-(ttcls-idxN)
Fcln <- do.call("cbind", replicate(ncF, Fcls2, simplify = FALSE))
Fcln <- as.data.frame(Fcln)
rdt.txfF <- rdt.tx[,1:(idxN)]
Frws <- cbind(rdt.txfF,Fcln)
ncol(Frws)
colnames(Frws) <- colnames(rdt)
phyla.rows <- Frws


# combine the list of data frames into a data frame
df_r01 <- dplyr::bind_rows(phyla.rows,
                           class.rows,
                           order.rows,
                           family.rows,
                           rdt)
df_r02 <- df_r01 %>% dplyr::arrange(Phylum, Klasse, Orden, Familie)
#View(df_r02)
# replace across all columns:
#  https://stackoverflow.com/questions/29271549/replace-all-occurrences-of-a-string-in-a-data-frame
df_r04 <- df_r02 %>%
  mutate( across(
    .cols = everything(),
    ~str_replace( ., "AA_", "" )
  ) )
# identify the rows in 'Lokal_geografi' where AA_ appears, and get the index number
# for the rows
idxn_tax <- which(grepl("AA_",df_r02$Lokal_geografi))
#View(df_r02)
#_______________________________________________________________________________
# begin function 'unfill_vec'
#_______________________________________________________________________________
# using an opposite version of the 'fill' function from 'tidyr'
#https://github.com/tidyverse/tidyr/issues/250
unfill_vec <- function(x) {
  same <- x == dplyr::lag(x)
  ifelse(!is.na(same) & same, "AA", x)
}
#_______________________________________________________________________________
# end function 'unfill_vec'
#_______________________________________________________________________________
# use the unfill function by applying it to all columns
# # https://stackoverflow.com/questions/7303322/apply-function-to-each-column-in-a-data-frame-observing-each-columns-existing-da
# this returns each column as a list, so these needs to be binded back
# together in to a data frame using 'dplyr::bind_rows'
# but this must only be performed on the taxonomical category columns
# as the other columns with numbers will be filled with 'AA's 
df_r04.2 <- df_r04 %>% select(Phylum, Klasse, Orden, Familie) %>% 
  lapply(unfill_vec) %>% 
  dplyr::bind_rows()
# select colunms that are not taxonomical categories, and which have not been modified
df_r04.3 <- df_r04 %>% select(-Phylum, -Klasse, -Orden, -Familie)
# combine data frames again
df_r05 <- cbind(df_r04.2, df_r04.3)

#df_r05 <- as.data.frame(t(df_r05))
oldNms.f.clmns <- (colnames(df_r05))
# define a set of column names that are shorter
replc2.f.clmns <- c("Phyl", "Klas", "Orde", "Famil",
                    "LtN", 
                    "lG", "oG", "DKN")
# replace column names so that column names are even shorter
colnames(df_r05) <- replc2.f.clmns
# paste together to get abbreviations explanations
abbrxexpl  <- paste(replc2.f.clmns,": ",oldNms.f.clmns,sep="")
abbrxexpl <- paste(abbrxexpl, collapse = ", ")
#View(df_r05)

# replace across all columns:
#  https://stackoverflow.com/questions/29271549/replace-all-occurrences-of-a-string-in-a-data-frame
df_r05 <- df_r05 %>%
  mutate( across(
    .cols = everything(),
    ~str_replace( ., "AA", "" )
  ) )

# make it a flextable to be able to merge
ft_r02 <- flextable(df_r05)
# only merge for the rows in the 'idxn_tax' identified above
# as only these rows have the 'AA_' in the 'Obs_Dk' column
ft_r02 <- merge_h(x = ft_r02, i=idxn_tax)
ft_r02 <- merge_v(ft_r02, j = c("Phyl", "Klas", "Orde", "Famil"),
                  combine = T)

# define a filename to store the flextable in
tf2 <- paste0(wd00,"/",wd14,"/Table14_v01_TS_and_NTS_indented.html")
# store the flextable as a html file
save_as_html(
  'subheader.title2' = ft_r02,
  path = tf2,
  title = paste0("Table_11: ", tbl.lgnd)
)
# see the flextable
ft_r02

# make a filename
tf2 <- paste0(wd00,"/",wd14,"/Table14_v02_NIS_priority_indented.xlsx")
# or write it out as an excel file, as it will be copied into a word document
xlsx::write.xlsx(df_r05, tf2)
# or use another package ....

#install.packages("taxlist")
#library(taxlist)
# #
# 
# ## Show taxonomy of papyrus
# ?indented_list(Easplist, "papyrus")
 
# ## Include synonyms and taxon views
# indented_list(Easplist, "papyrus", level = TRUE, synonyms = TRUE,
#               secundum = "secundum")
# #

#