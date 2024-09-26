
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
# set working directory
#wd00 <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6/Results_from_ABI7500_for_MONIS6_species_specific_assays"
wd00 <- getwd()
setwd (wd00)
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

clNmtorep <- "Klasse"
#clNmtorep <- "Phylum"
indxf.txcl <- which(clNmtorep==colnames(rdt))
nclrdt <- ncol(rdt)
txg.in.rdt <- as.vector(unique(as.data.frame(rdt)[,indxf.txcl]))
n.txg <- length(txg.in.rdt)
sq.txg <- seq(1,n.txg,1)
rwtadd.txg <- list()
for (i in sq.txg)
{
  Nm <- txg.in.rdt[i]
  print(Nm)
  rtadd <- rep(Nm,nclrdt)
  idxf <- indxf.txcl
  ctbblnk <- seq(1,idxf,1)
  df_txrep <- as.data.frame(t(rtadd))
  # if (length(ctbblnk)>1){
  # df_txrep[,c(ctbblnk)] <- ""
  #   }
  rwtadd.txg[[i]] <- df_txrep
}

# combine the list of data frames into a data frame
df_r01 <- dplyr::bind_rows(rwtadd.txg)
df_r01

colnames(df_r01) <- colnames(rdt)
df_r02 <- rbind(df_r01,rdt)

ft_r02 <- flextable(df_r02)
ft_r02 <- merge_h(x = ft_r02)
ft_r02

# or use another package ....

#install.packages("taxlist")
library(taxlist)
#

## Show taxonomy of papyrus
?indented_list(Easplist, "papyrus")

## Include synonyms and taxon views
indented_list(Easplist, "papyrus", level = TRUE, synonyms = TRUE,
              secundum = "secundum")
#

#