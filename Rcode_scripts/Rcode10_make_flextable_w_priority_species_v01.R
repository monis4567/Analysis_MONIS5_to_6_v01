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
wd10 <- "output10_species_priority_table"
# define full path for input directory
outdir10 <- paste(wd00,wd10, sep="/")
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
ls.fl01.xls <- ls.fl01.xls[grep("Artsprio", ls.fl01.xls)]
ls.fl01.xls <- ls.fl01.xls[grep("v05", ls.fl01.xls)]

infl01 <- paste0(wd00_wddata,"/",ls.fl01.xls)
# read in xls spreadsheet tab with raw flourescense data from ABI7500
rdt <- readxl::read_xlsx(infl01,
                        sheet =  2)
# replace the column names
colnames(rdt) <- rdt[2,]
# exlude rows
rdt <- rdt[-c(1,2),]
# exlude columns
rdt <- rdt[,-c(22,23)]
# define columns to keep
ctkeep <- c("Phylum", "Klasse", "Orden", "Familie", "Genus og species og author", 
  "Obs Dk", "Nseq NCBI for IkkHjmM Art", 
  "Mngl Nseq NCBI for tbsl art IkkHjmM Art", "Mngl Nseq NCBI for tbsl art i fam IkkHjmM Art", 
  "Mngl Nseq NCBI for tbsl art i ord IkkHjmM Art", "Mngl Nseq NCBI for tbsl art i gen IkkHjmM Art", 
  "PrNr", "Reference")
# limit to only columns specifeid in vector
rdt <-rdt[ctkeep]
# substitute in column names
colnames(rdt) <- gsub(" ","_",colnames(rdt))
# re order the rows by the columns
rdt <- rdt %>% dplyr::arrange(Phylum,
                Klasse, 
                Orden, 
                Familie, 
                Genus_og_species_og_author)
# exclude row if phylum is NA
rdt <- rdt[!(is.na(rdt$Phylum)),]


clNMs <- colnames(rdt)
# identify the rows where there is a reference
df_rdf <- rdt[!(is.na(rdt$Reference)),]
# split string to get genus and species
Gnsp.splt <- strsplit(df_rdf$Genus_og_species_og_author, " ")
# get the genus and species
genus.Nm <- sapply(Gnsp.splt, "[[", 1)
species.Nm <- sapply(Gnsp.splt, "[[", 2)
# get the authors and year for the publication with the detection assay
Ref.splt <- strsplit(df_rdf$Reference, ", ")
# get the authors
auth.ref <- sapply(Ref.splt, "[[", 1)
year.ref <- sapply(Ref.splt, "[[", 2)
# get the genus and species and authors
gnsp.Nm <- paste(genus.Nm,species.Nm,sep=" ")
# paste together genus species names and author and reference
dtct.ass.spc <- paste0("Af ",auth.ref," (",year.ref,") mod '",genus.Nm," ", species.Nm,"'")
# make the the refenreces for the detection assays one single string
dtct.ass.spc <- paste(dtct.ass.spc, collapse = ". ")
# make a text string  that can be used for the table legend
txt.dtct.ass.spc <- paste0("Tidligere detektionssystemer er udviklet: ", dtct.ass.spc)
# identify the column index number for the column named 'Reference'
idxclnmb.ref<- which(grepl("Reference",colnames(rdt)))
# exclude the the column named 'Reference'
rdt <- rdt[,-idxclnmb.ref]

#View(rdt)
#_______
#https://stackoverflow.com/questions/8753531/repeat-rows-of-a-data-frame-n-times#8753732

clNMstr <- c("Phylum", "Klasse", "Orden", "Familie")

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
# replace across all columns:
#  https://stackoverflow.com/questions/29271549/replace-all-occurrences-of-a-string-in-a-data-frame
df_r04 <- df_r02 %>%
  mutate( across(
    .cols = everything(),
    ~str_replace( ., "AA_", "" )
  ) )
# identify the rows in 'Obs_DK' where AA_ appears, and get the index number
# for the rows
idxn_tax <- which(grepl("AA_",df_r02$Obs_Dk))
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
df_r04.2 <- df_r04 %>% select(Phylum, Klasse, Orden, Familie, Genus_og_species_og_author) %>% 
  lapply(unfill_vec) %>% 
  dplyr::bind_rows()
# select that are not taxonomical categories, and which have not been modified
df_r04.3 <- df_r04 %>% select(Obs_Dk,Nseq_NCBI_for_IkkHjmM_Art, Mngl_Nseq_NCBI_for_tbsl_art_IkkHjmM_Art, 
                  Mngl_Nseq_NCBI_for_tbsl_art_i_fam_IkkHjmM_Art, Mngl_Nseq_NCBI_for_tbsl_art_i_ord_IkkHjmM_Art, 
                  Mngl_Nseq_NCBI_for_tbsl_art_i_gen_IkkHjmM_Art)
# combine data frames again
df_r05 <- cbind(df_r04.2, df_r04.3)

#df_r05 <- as.data.frame(t(df_r05))
oldNms.f.clmns <- (colnames(df_r05))
# define a set of column names that are shorter
replc2.f.clmns <- c("Phyl", "Klas", "Orde", "Famil", "GeSpAu", 
  "ODK", "nS", "mS", 
  "mSs", "mSo", 
  "mSg")
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
# define headers and subheaders for the flextable
subheader.title <- paste0("Liste over ikke-hjemmehørende arter for danske farvande,",
"som det ville være fornuftigt at prioritere i en overvågning af miljø-DNA. ",
"Alle arter er listet indenfor ‘phylum’, ‘klasse’, ‘orden’, ‘familie’ og ‘genus’ og",
"art, med author navn for hvem der har beskrevet arten og årstallet for beskrivelsen. ",
txt.dtct.ass.spc,". Forkortelser for kolonner: ",abbrxexpl) 

# define a filename to store the flextable in
tf2 <- paste0(wd00,"/",wd10,"/Table11_v01_NIS_priority_indented.html")
# store the flextable as a html file
save_as_html(
  'subheader.title2' = ft_r02,
  path = tf2,
  title = paste0("Table_11: ", subheader.title)
)
# see the flextable
ft_r02

# make a filename
tf2 <- paste0(wd00,"/",wd10,"/Table11_v02_NIS_priority_indented.xlsx")
# or write it out as an excel file, as it will be copied into a word document
xlsx::write.xlsx(df_r05, tf2)
# or use another package ....

#install.packages("taxlist")
library(taxlist)
# #
# 
# ## Show taxonomy of papyrus
# ?indented_list(Easplist, "papyrus")
# 
# ## Include synonyms and taxon views
# indented_list(Easplist, "papyrus", level = TRUE, synonyms = TRUE,
#               secundum = "secundum")
# #

#