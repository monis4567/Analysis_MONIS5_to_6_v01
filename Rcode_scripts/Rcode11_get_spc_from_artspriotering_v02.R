#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-


#On the remote server. Ensure you have a directory where your 
#R packages are stored
# make one here:
# /home/sknu003/R/R_packages_for_Rv4_3
# load the R module
# $ module load R/4.3.2-foss-2023a

# # First make sure no unneeded modules have been loaded
# $ module purge
# You can use 'module spider R' to see which version of R is available on the remote server
# $ module spider R
# # I will try installing for Rv4_0_2
# # Start out by making a directory where all the packages can be placed inside
# $ mkdir R_packages_for_Rv4_0_2
# $ cd R_packages_for_Rv4_0_2/

# # start up R by typing R

# $ R

# ## In R 
# ## Run the lines below without the # sign. The other lines with 2 # signs are helpful comments.
# ## You also need to run each of these lines one by one individually
# # In R you now first need to specify a path to the directory where you want your packages
# # to be available for your R-code
# # Run these lines - changing the path your own directory for where the packages need to be 
# # replace my library path to your own library path
# # You will have to run each line one at a time

# # Run this line in R to specify the path to where you want the packages to placed:

# lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_3"

# # Continue by running these lines, one by one in R
# Sys.setenv(R_LIBS_USER="lib_path01")
# .libPaths("lib_path01")
# # change the path to where the packages should be installed from # see this website: https://stackoverflow.com/questions/15170399/change-r-default-library-path-using-libpaths-in-rprofile-site-fails-to-work
# .libPaths( c( lib_path01 , .libPaths() ) )
# .libPaths()
# .libPaths( c( lib_path01) )

## Or try pasting one long line with all commands, and installation of the 'taxizedb' library #
# lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_3"; Sys.setenv(R_LIBS_USER="lib_path01"); .libPaths("lib_path01"); .libPaths( c( lib_path01 , .libPaths() ) ); .libPaths(); .libPaths( c( lib_path01) ); install.packages(c("taxizedb", "taxize", "tidyverse", "readxl", "worms", "stringr", "dplyr"))

## Use this line here below to install a lot of packages needed for doing population genetic analysis
# lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_3"; Sys.setenv(R_LIBS_USER="lib_path01"); .libPaths("lib_path01"); .libPaths( c( lib_path01 , .libPaths() ) ); .libPaths(); .libPaths( c( lib_path01) ); if(!require("gaston")){install.packages("gaston", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require("hierfstat")){install.packages("hierfstat", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require("pegas")){install.packages("pegas")};if(!require("ips")){install.packages("ips", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require("tidyverse")){install.packages("tidyverse", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require("pals")){install.packages("pals", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require(adegenet)){install.packages("adegenet", repos='http://cran.us.r-project.org') };if(!require(apex)){install.packages("apex", repos='http://cran.us.r-project.org') };if(!require(mmod)){install.packages("mmod", repos='http://cran.us.r-project.org') };if(!require(tidyverse)){install.packages("tidyverse", repos='http://cran.us.r-project.org') };if(!require(pals)){install.packages("pals", repos='http://cran.us.r-project.org') };if(!require(ape)){install.packages("ape", repos='http://cran.us.r-project.org') };if(!require(RColorBrewer)) {install.packages("RColorBrewer", repos='http://cran.us.r-project.org') };if(!require(stringi)){install.packages("stringi", repos='http://cran.us.r-project.org') };if(!require(poppr)){install.packages("poppr", repos='http://cran.us.r-project.org') };if(!require(vegan)){install.packages("vegan", repos='http://cran.us.r-project.org')};if(!require(adegenet)){install.packages("adegenet", repos='http://cran.us.r-project.org')};if(!require(biogeo)){install.packages("biogeo", repos='http://cran.us.r-project.org') }
## Use this line here below to install a lot of packages needed for doing mapping
# lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_3"; Sys.setenv(R_LIBS_USER="lib_path01"); .libPaths("lib_path01"); .libPaths( c( lib_path01 , .libPaths() ) ); .libPaths(); .libPaths( c( lib_path01) ); if(!require(scales)){  install.packages("scales")};if(!require(fields)){  install.packages("fields")};if(!require(marmap)){  install.packages("marmap")};if(!require(TeachingDemos)){  install.packages("TeachingDemos")};if(!require(rworldmap)){  install.packages("rworldmap")};if(!require(rworldxtra)){  install.packages("rworldxtra")};require(rworldxtra);if(!require(readxl)){  install.packages("readxl")};if(!require(plyr)){  install.packages("plyr")};if(!require(mapdata)){  install.packages("mapdata")};if(!require(maps)){  install.packages("maps")};if(!require(mapplots)){  install.packages("mapplots")}; if(!require(purrr)){  install.packages("purrr")}

# ## In R 
# ## Run the lines below without the # sign. The other lines with 2 # signs are helpful comments.
# ## You also need to run each of these lines one by one individually


# install.packages(c("latticeExtra")) # Not available for R v3.4.1
# install.packages(c("robustbase")) # Not available for R v3.4.1
# install.packages(c("lessR")) # Not available for R v3.4.1
# install.packages(c("readr"))
# install.packages(c("dada2")) # Not available for R v3.4.1 and not for R v4.0.2
# install.packages(c("dada2", "readr", "dplyr", "taxize", "tidyr"))
# install.packages(c("taxizedb"))
# install.packages(c("ShortRead")) #Not available for R v3.4.1
# install.packages(c("taxize")) # Dependent on 'wikitaxa' which is Not available for R v3.4.1
# install.packages(c("tidyverse"))

# ## Note that this might exit with an error message . Read on for the next instructions below:

# ## As I had difficulties getting the 'dada2' package installed on the remote path
# ## I tried to look for solutions 
# ## I then looked up the 'dada2' package for R on the internet here:
# ## https://www.bioconductor.org/packages/release/bioc/html/dada2.html
# ## and this webpage recommended that I ran these two commands (again running one line at a time)
# ## You will "ShortRead" installed before you can install "dada2". Use "BiocManager" to install both
# ## First start with: "BiocManager"
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("ShortRead")
# BiocManager::install("dada2")

# ## This should get your 'dada2' package installed
# ## Check that R finishes without errors on packages 
# ## It might be that it instead finishes with a warning message. This is okay. As long as you are not getting an error message.
# ## It might be that some packages fail, but as long as you get "dada2" and "readr" and "BiocManager" installed without
# ## error I think it is just fine
# ## Now quit R 

# ## Now that you have installed all packages in your remote path : "/home/sknu003/R/R_packages_for_Rv4_3"
# # You are almost ready to run parts of the code together with their matching slurm sbatch submission scripts

# specify a path to where you have all your packages on your remote node
# NOTICE ! In order to have your packages available on the remote path, you will need to logon to your node, 
# and make a directory 
# called -e.g. : R_packages_for_Rv4_0_2
# Then load the R module
# module load R/v4.0.2
# Then start R byt typing:
# R
# Once R is started run the lines here below in section 01
# #_______________start section 01__________________________________________ 
# replace my library path to your own library path
#lib_path01 <- "/groups/hologenomics/phq599/data/R_packages_for_Rv3_6"
lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_3"
Sys.setenv(R_LIBS_USER="lib_path01")
.libPaths("lib_path01")
# change the path to where the packages should be installed from # see this website: https://stackoverflow.com/questions/15170399/change-r-default-library-path-using-libpaths-in-rprofile-site-fails-to-work
.libPaths( c( lib_path01 , .libPaths() ) )
.libPaths()
# #_______________end section 01__________________________________________
## You  will need to specify this path again later on when you are to run this R-script
## Before you can start this R-script on the remote server you will need to install the pacakges here

# library(c("dada2")) # Not available for R v3.4.1 and not for R v4.0.2
# library(c("dada2", "readr", "dplyr", "taxize", "tidyr"))
# library(c("latticeExtra")) # Not available for R v3.4.1
# library(c("lessR")) # Not available for R v3.4.1
# library(c("readr"))
# library(c("robustbase")) # Not available for R v3.4.1
# library(c("ShortRead")) #Not available for R v3.4.1
# library(c("taxize")) # Dependent on 'wikitaxa' which is Not available for R v3.4.1
# library("BiocManager")
# library("adegenet")
# library("ape")
# library("fields")
# library("mapdata")
# library("mapplots")
# library("maps")
# library("marmap")
# library("pillar")
# library("plyr")
# library("poppr")
 library("purrr")
# library("RColorBrewer")
 library("readxl")
# library("rworldmap")
# library("rworldxtra")
# library("scales")
# library("stringi")
 library("stringr")
# library("TeachingDemos")
 library("tibble")
 library("tidyverse")
# library("vegan")
 library("worrms")


wd00 <- getwd()
# wd00 <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6/Analysis_MONIS5_to_6_v01"
# setwd(wd00)

#define directory with output flies
wdout <- "output11_get_species_from_priority_table"
# define full path for input directory
outdir10 <- paste(wd00,wdout, sep="/")
# remove previous versions of the 'outdir10'
#unlink(outdir10, force = T, recursive = T)
# create the 'outdir10' in a new version
#dir.create(outdir10)
# define directory with input files to read in
wddata <- "data"
#paste together a path for the directory that holds the inputs files
wd00_wddata <- paste(wd00,wddata,sep="/")

## in section 02 here below 
## installing these packages is commented out, as they are not needed when you are to run this R-script
## But I have left them here in section 02 to be used when you install them the first time in your local 
## path
# if(!require(pillar)){
#   install.packages("pillar")
# }
library(pillar)
# if(!require(tibble)){
#   install.packages("tibble")
# }
library(tibble)
# if(!require(readxl)){
#   install.packages("readxl")
# }

# if(!require(worrms)){
#   install.packages("worrms")
# }
# if(!require(stringr)){
#   install.packages("stringi")
#   install.packages("stringr")
# }

library(worrms)
library(readxl)
library(stringr)

#define paths to working directories
#set the working dir
#setwd(wd00)
#wd00 <- getwd()
# paste path and directory together


wddata <- "data"
wd00_wddata <- paste(wd00,wddata,sep="/")
# get a lst of filesin the data directory
lfdd <- list.files(wd00_wddata,full.names = T) 
#
pthinf01 <-lfdd[grepl("Artsprioritering",lfdd)]
pthinf01 <- pthinf01[grepl("_v05",pthinf01)]
# read in excel file as tibble, skip 2 rows
tibl_inx01 <- readxl::read_xlsx(pthinf01,sheet = 2)
# get the tibble but not the first row
tibl_inx01 <- tibl_inx01[-1,]
colnames(tibl_inx01) <- tibl_inx01[1,]
tibl_inx01 <- tibl_inx01[-1,]
clNms <- colnames(tibl_inx01)
idxNmb.spc <- grep("pecies",clNms)
# use hte index number
tibl_inx01$Species <- tibl_inx01[,idxNmb.spc]
# subset to exclude if NA
tibl_inx01 <- tibl_inx01[!(is.na(tibl_inx01$Species)),]
spcNms <- as.character(unlist(tibl_inx01$Species))
spcNms <- gsub(" ", "_", spcNms)
# split the string by delimiter
spcNmssplt  <- strsplit(as.character(spcNms), "_")
# and use the splitted string to get each element and add back in a column
gnnms  <- sapply(spcNmssplt, "[[", 1)
spnms  <- sapply(spcNmssplt, "[[", 2)
# add back to tibble
tibl_inx01$genusNm <- gnnms
tibl_inx01$speciesNm <- spnms
tibl_inx01$Species <- paste0(tibl_inx01$genusNm," ",tibl_inx01$speciesNm)
# remove all non alphanumeric characters: https://stackoverflow.com/questions/10294284/remove-all-special-characters-from-a-string-in-r
tspE <- stringr::str_replace_all(tibl_inx01$Species, "[^[:alnum:]]", " ")
# split string and get lists nested in a list
gl <- strsplit(as.character(tspE), " ")
#get first and second element of nested list
gnl <- sapply(gl, "[[", 1)
spl <- sapply(gl, "[[", 2)
# paste genus and species name together
gspl <- paste0(gnl," ",spl)
# mkae list of unique genera names
uspl <- unique(gnl)
uspl <- unique(gspl)
# order the genus and species names alphabetically
uspl <- uspl[order(uspl)]

# if(!require(tidyverse)){
#   install.packages("tidyverse")
# }  
library(tidyverse)

#install packages
#get worms package
#get worms package
# if(!require(worrms)){
#   install.packages("worrms")
# }  
library(worrms)
# get the extra codes with functions
source(paste0(wd00,"/Rcode_scripts/worms_safe.R"))



# #https://stackoverflow.com/questions/39285570/error-in-curlcurl-fetch-memoryurl-handle-handle-couldnt-connect-to-ser
# library(httr)    
# set_config(use_proxy(url="10.3.100.207",port=8080))

# subset the list of species names
#uspll <- uspl[c(1,2,4:7)]
uspll <- uspl
#uspll <- uspl[19]
# example species
#uspll <-c("Petricolaria pholadiformis","Streblospio benedicti" )
#uspll <-c("Petricolaria pholadiformis","Streblospio benedicti" )
#make an empty list to add family to
lstap.f <- list()
#make an empty list to add genera to
lstap.g <- list()
#uspll <- uspll[42]
AphiaID <- 128900
#iterate over elements in list of species'
for (c in uspll){
  # get the matching number in the vector of elements
  i <- match(c,uspll)
  # get Aphia info for the species
  aphiaInfo<-GetSpeciesID(c)
  # get the AphiaID for the species
  aphiaID <- aphiaInfo$AphiaID
  # get the family name
  familyNm <- GetAphiaParameter(aphiaID,"family")
  # get  alist with info on tha family
  aphiaInfoFamily<-GetSpeciesID(familyNm)
  # get the AphiaID for the family
  af2 <- aphiaInfoFamily$AphiaID
  # get children of the AhpiaID for the family
  af3 <- GetAphiaChildren(af2)
  
  # remove results where genus is NA
  af3 <- af3 %>% dplyr::filter(!is.na(genus))
  # remove results where family is NA
  af3 <- af3 %>% dplyr::filter(!is.na(family))
  # remove results where order is NA
  af3 <- af3 %>% dplyr::filter(!is.na(order))
  # get a vector of genera in the family
  gn4 <- af3$genus
  
  #iterate over genera in family
  for (g2 in gn4)
  {
    #get number for genus in vector of genera
    j <- match(g2,gn4)
  # get the species ID info
  apInf <- GetSpeciesID(g2)
  # from this info get the Aphia ID number
  aphID <- apInf$AphiaID
  # use the Aphia ID number to get sub AphiaIDs
  atmp <- GetAphiaChildren(aphID)
  # append this data frame to a list of data frames
  lstap.g[[j]] <- as.data.frame(atmp)
  # end iteration over genera in the family
  }
  # make the list of data frames obtained for the genera per family a data frame
  df_ac01 <- as.data.frame(do.call(rbind,lstap.g))
  #limit to only include hits within the family searched for
  # if the genus name also appears in other families
  df_ac01[df_ac01$family==familyNm,]
    #make the list of data frames a single data frame
    lstap.f[[i]] <- as.data.frame(df_ac01)
  # end iteration over species and their families
  }
  #make the list of data fr ames a single data frame
df_ac02 <- as.data.frame(do.call(rbind,lstap.f))
df_ac02 <- df_ac02[!is.na(df_ac02$AphiaID),]
# paste the working directory and the output directory together
wd00_wd10 <- paste(wd00,wdout, sep="/")
# save the data frame to a csv file
flNm<-"priority_spc.csv"
folder_out <- wd00_wd10
write.table(df_ac02,file=paste0(folder_out,"/",flNm),row.names=F,col.names=T,sep=";",quote=F)
#read in the table
df_ac02.2 <- read.csv2(file=paste0(folder_out,"/",flNm),
                        #row.names=F,col.names=T,
                        sep=";",header=T)
# unique(df_ac02$class)
#unique(df_ac02$valid_name)
#View(df_ac02)
nAID <- length(df_ac02$AphiaID)
nAID <- seq(1,nAID,1)
#nAID <- seq(1,33,1)
lst_AID <- list()
#
for (AID in nAID)
  {
  A <- df_ac02$AphiaID[AID]
  #A <- 1507114
  AIDD <- GetSpeciesDistributions(A)
  lst_AID[[AID]] <- AIDD 
}
#bind the rows in each list in to one data frame
df_l01 <- data.table::rbindlist(lst_AID, fill=T)
df_l01 <- as.data.frame(df_l01)
# limit to only aline species
df_l01A <- df_l01[(df_l01$establishmentMeans=="Alien"),]
 unique(df_l01A$ScientificName)

locatNms <- df_l01 %>% dplyr::distinct( locality) 
locatNms <- locatNms[order(locatNms$locality),]

NEAlocs <- c("Belgian Exclusive Economic Zone",
"Belgium" ,
"Atlantic Europe",                        
"European waters (ERMS scope)" ,
"Germany",
"Netherlands",
"North Sea",
"United Kingdom" ,
"West Coast of Scotland" )
NEAlocs <- paste(NEAlocs, collapse = "|")
df_l02 <- df_l01[grepl(NEAlocs,df_l01$locality),]
NEAspc <- unique(df_l02$ScientificName)

# z04 <- df_ac02[df_ac02$class=="Actinopteri",]
# View(z04)
df_ac03 <- df_ac02[grep(" ",df_ac02$valid_name),]
#unique(df_ac03$class)
#
flNm<-"priority_spc.csv"
folder_out <- wd00_wd10
write.table(df_ac02,file=paste0(folder_out,"/",flNm),row.names=F,col.names=T,sep=";",quote=F)

#