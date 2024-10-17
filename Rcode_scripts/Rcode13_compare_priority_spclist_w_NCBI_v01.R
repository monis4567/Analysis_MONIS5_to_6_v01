#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-


# #On the remote server. Ensure you have a directory where your 
# #R packages are stored
# # make one here:
# # /home/sknu003/R/R_packages_for_Rv4_1
# # load the R module
# # $ module load R/4.1.0-gimkl-2020a
# 
# # # First make sure no unneeded modules have been loaded
# # $ module purge
# # You can use 'module spider R' to see which version of R is available on the remote server
# # $ module spider R
# # # I will try installing for Rv4_3
# # # Start out by making a directory where all the packages can be placed inside
# # $ mkdir R_packages_for_Rv4_3
# # $ cd R_packages_for_Rv4_3/
# 
# # # start up R by typing R
# 
# # $ R
# 
# # ## In R 
# # ## Run the lines below without the # sign. The other lines with 2 # signs are helpful comments.
# # ## You also need to run each of these lines one by one individually
# # # In R you now first need to specify a path to the directory where you want your packages
# # # to be available for your R-code
# # # Run these lines - changing the path your own directory for where the packages need to be 
# # # replace my library path to your own library path
# # # You will have to run each line one at a time
# 
# # # Run this line in R to specify the path to where you want the packages to placed:
# 
# # lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_1"
# 
# # # Continue by running these lines, one by one in R
# # Sys.setenv(R_LIBS_USER="lib_path01")
# # .libPaths("lib_path01")
# # # change the path to where the packages should be installed from # see this website: https://stackoverflow.com/questions/15170399/change-r-default-library-path-using-libpaths-in-rprofile-site-fails-to-work
# # .libPaths( c( lib_path01 , .libPaths() ) )
# # .libPaths()
# # .libPaths( c( lib_path01) )
# 
# ## Or try pasting one long line with all commands, and installation of the 'taxizedb' library #
# # lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_1"; Sys.setenv(R_LIBS_USER="lib_path01"); .libPaths("lib_path01"); .libPaths( c( lib_path01 , .libPaths() ) ); .libPaths(); .libPaths( c( lib_path01) ); install.packages(c("taxizedb", "taxize", "tidyverse", "readxl", "worms", "stringr", "dplyr"))
# 
# ## Use this line here below to install a lot of packages needed for doing population genetic analysis
# # lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_1"; Sys.setenv(R_LIBS_USER="lib_path01"); .libPaths("lib_path01"); .libPaths( c( lib_path01 , .libPaths() ) ); .libPaths(); .libPaths( c( lib_path01) ); if(!require("gaston")){install.packages("gaston", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require("hierfstat")){install.packages("hierfstat", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require("pegas")){install.packages("pegas")};if(!require("ips")){install.packages("ips", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require("tidyverse")){install.packages("tidyverse", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require("pals")){install.packages("pals", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require(adegenet)){install.packages("adegenet", repos='http://cran.us.r-project.org') };if(!require(apex)){install.packages("apex", repos='http://cran.us.r-project.org') };if(!require(mmod)){install.packages("mmod", repos='http://cran.us.r-project.org') };if(!require(tidyverse)){install.packages("tidyverse", repos='http://cran.us.r-project.org') };if(!require(pals)){install.packages("pals", repos='http://cran.us.r-project.org') };if(!require(ape)){install.packages("ape", repos='http://cran.us.r-project.org') };if(!require(RColorBrewer)) {install.packages("RColorBrewer", repos='http://cran.us.r-project.org') };if(!require(stringi)){install.packages("stringi", repos='http://cran.us.r-project.org') };if(!require(poppr)){install.packages("poppr", repos='http://cran.us.r-project.org') };if(!require(vegan)){install.packages("vegan", repos='http://cran.us.r-project.org')};if(!require(adegenet)){install.packages("adegenet", repos='http://cran.us.r-project.org')};if(!require(biogeo)){install.packages("biogeo", repos='http://cran.us.r-project.org') }
# ## Use this line here below to install a lot of packages needed for doing mapping
# # lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_1"; Sys.setenv(R_LIBS_USER="lib_path01"); .libPaths("lib_path01"); .libPaths( c( lib_path01 , .libPaths() ) ); .libPaths(); .libPaths( c( lib_path01) ); if(!require(scales)){  install.packages("scales")};if(!require(fields)){  install.packages("fields")};if(!require(marmap)){  install.packages("marmap")};if(!require(TeachingDemos)){  install.packages("TeachingDemos")};if(!require(rworldmap)){  install.packages("rworldmap")};if(!require(rworldxtra)){  install.packages("rworldxtra")};require(rworldxtra)if(!require(readxl)){  install.packages("readxl")};if(!require(plyr)){  install.packages("plyr")};if(!require(mapdata)){  install.packages("mapdata")};if(!require(maps)){  install.packages("maps")};if(!require(mapplots)){  install.packages("mapplots")}; if(!require(purrr)){  install.packages("purrr")}
# 
# # ## In R 
# # ## Run the lines below without the # sign. The other lines with 2 # signs are helpful comments.
# # ## You also need to run each of these lines one by one individually
# 
# 
# # install.packages(c("latticeExtra")) # Not available for R v3.4.1
# # install.packages(c("robustbase")) # Not available for R v3.4.1
# # install.packages(c("lessR")) # Not available for R v3.4.1
# # install.packages(c("readr"))
# # install.packages(c("dada2")) # Not available for R v3.4.1 and not for R v4.0.2
# # install.packages(c("dada2", "readr", "dplyr", "taxize", "tidyr"))
# # install.packages(c("taxizedb"))
# # install.packages(c("ShortRead")) #Not available for R v3.4.1
# # install.packages(c("taxize")) # Dependent on 'wikitaxa' which is Not available for R v3.4.1
# # install.packages(c("tidyverse"))
# 
# # ## Note that this might exit with an error message . Read on for the next instructions below:
# 
# # ## As I had difficulties getting the 'dada2' package installed on the remote path
# # ## I tried to look for solutions 
# # ## I then looked up the 'dada2' package for R on the internet here:
# # ## https://www.bioconductor.org/packages/release/bioc/html/dada2.html
# # ## and this webpage recommended that I ran these two commands (again running one line at a time)
# # ## You will "ShortRead" installed before you can install "dada2". Use "BiocManager" to install both
# # ## First start with: "BiocManager"
# # if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# # BiocManager::install("ShortRead")
# # BiocManager::install("dada2")
# 
# # ## This should get your 'dada2' package installed
# # ## Check that R finishes without errors on packages 
# # ## It might be that it instead finishes with a warning message. This is okay. As long as you are not getting an error message.
# # ## It might be that some packages fail, but as long as you get "dada2" and "readr" and "BiocManager" installed without
# # ## error I think it is just fine
# # ## Now quit R 
# 
# # ## Now that you have installed all packages in your remote path : "/home/sknu003/R/R_packages_for_Rv4_1"
# # # You are almost ready to run parts of the code together with their matching slurm sbatch submission scripts
# 
# # specify a path to where you have all your packages on your remote node
# # NOTICE ! In order to have your packages available on the remote path, you will need to logon to your node, 
# # and make a directory 
# # called -e.g. : R_packages_for_Rv4_3
# # Then load the R module
# # module load R/v4.0.2
# # Then start R byt typing:
# # R
# # Once R is started run the lines here below in section 01
# #_______________start section 01__________________________________________ 
# # replace my library path to your own library path
# #lib_path01 <- "/groups/hologenomics/phq599/data/R_packages_for_Rv3_6"
# lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_1"
# Sys.setenv(R_LIBS_USER="lib_path01")
# .libPaths("lib_path01")
# # change the path to where the packages should be installed from # see this website: https://stackoverflow.com/questions/15170399/change-r-default-library-path-using-libpaths-in-rprofile-site-fails-to-work
# .libPaths( c( lib_path01 , .libPaths() ) )
# .libPaths()
# #_______________end section 01__________________________________________
# ## You  will need to specify this path again later on when you are to run this R-script
# ## Before you can start this R-script on the remote server you will need to install the pacakges here
# 
# ## in section 02 here below 
# ## installing these packages is commented out, as they are not needed when you are to run this R-script
# ## But I have left them here in section 02 to be used when you install them the first time in your local 
# ## path

#load required libraries
library(purrr)
library(tidyverse)
# required function GetSpeciesDistributions(AphiaID)
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
# library("purrr")
# library("RColorBrewer")
# library("readxl")
# library("readxl")
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
library(jsonlite)



wd00 <- getwd()
# wd00 <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6/Analysis_MONIS5_to_6_v01"
# setwd(wd00)

#define directory with output flies
wdout <- "output12_compare_species_priority_table"
# define full path for input directory
outdir12 <- paste(wd00,wdout, sep="/")
# remove previous versions of the 'outdir12'
unlink(outdir12, force = T, recursive = T)
# create the 'outdir12' in a new version
dir.create(outdir12)
# define directory with input files to read in
wddata <- "data"
#paste together a path for the directory that holds the inputs files
wd00_wddata <- paste(wd00,wddata,sep="/")

library(worrms)
# get the extra codes with functions
source(paste0(wd00,"/Rcode_scripts/worms_safe.R"))


#define working directory
#wd00 <- "/home/hal9000/Documents/shrfldubuntu18/MONIS6_getNCBIseq"
#setwd(wd00)
# wd00 <- getwd()
# setwd(wd00)
wd00_wd01 <- outdir12

#Define input file name
inpfnm1 <- "priority_spc.csv"
inpfnm2 <- "MRG_matches.csv"
# paste path and file name together
wd00_inpfnm1 <- paste(wd00,"/data/",inpfnm1,sep="")
#wd00_inpfnm1 <- paste(wd00,"/",inpfnm1,sep="")
wd00_inpfnm2 <- paste(wd00,"/data/",inpfnm2,sep="")
# read in the csv file with species from the worms data base
df_spp01 <- read.csv(wd00_inpfnm1, sep=";", header=T,fill=T)
# read in the csv file with Marien regions close to Scandinavia
df_MRG01 <- read.csv(wd00_inpfnm2, sep=";", header=T,fill=T)
# grep for a space in the valid name -  i.e. the binomial species names
# since the family and genus names are irrellevant
df_spp02 <- df_spp01[grep(" ",df_spp01$valid_name),]
# only retain rows for unique valid names  - the following steps is to search 
# worms database once again to get the area of distribution for each of the
# species. This needs only to be done for the valid species
df_spp03 <- df_spp02[!duplicated(df_spp02[ , c("valid_name")]),]
#count the number of rows left
nrow(df_spp03)
unique(df_spp01$class)
#match(144714,df_spp03$AphiaID)
# limit to the first 10 rows  -use for trying out the code below on 
# a small scale
#df_spp03 <- df_spp03[1:10,]
#View(df_spp03)

library(purrr)
library(tidyverse)
# required function GetSpeciesDistributions(AphiaID)
#source("worms.R")
# Get list of aphia IDs e.g. from the column of an existing dataframe
AphiaList <- df_spp03$AphiaID
#ensure the AphiaID are numeric to enable you to search the worms database
AphiaList <- as.numeric(AphiaList)
# exclude any NA's from the AphiaList
AphiaList <- AphiaList[!is.na(AphiaList)]
# get the number of elements in the AphiaList -  i.e get the length of the vector 
noofAphIDs <-  length(AphiaList)
testNmbs<- c(338185,
             152270,
             314064,
             307607)
# iterate over test numbers
for (i in testNmbs){
tstNm <- AphiaList[grepl(i,AphiaList)]
print(tstNm)}

# try with shortened AphiaList
#AphiaList <- c(353452, 364147, 359400, 359401, 1473395, 576829, 326383, 1460751, 338526, 338528, 1460758, 558835, 871991, 1448034, 331106, 867786, 559058, 1421293, 338527, 829025, 559059, 131151, 1460759, 209866, 871993, 338529, 766742, 131152, 867787, 1379353, 1460756, 131153, 338530, 558837, 867788, 331131, 331109, 1460757, 829101, 157542, 1321265, 558834, 131154, 338531, 559063, 131155, 131156, 558843, 558844, 428899, 131157, 131158, 332986, 869237, 338447, 326470, 335302, 155278, 326471, 337494, 157179, 326472, 867785, 131107, 869236, 326473, 326112, 329535, 329536, 130353, 329537, 869882, 329538, 337183, 337184, 334135, 329539, 329540, 157499, 334136, 334137, 209805, 329541, 130357, 334138, 329542, 329543, 329544, 209804, 157500, 337194, 334139, 337198, 173701, 334140, 329545, 329546, 757997, 334141, 130359, 334142, 130360, 334143, 334144, 130362, 156263, 594747, 130363, 594751, 130364, 209802, 334146, 329548, 154982, 1039009, 1252770, 1039010, 1039008, 328633, 174800, 598899, 598793, 131126, 328634, 131127, 328635, 328636, 328637, 328638, 131128, 872662, 328639, 872661, 1518243, 328640, 1518242, 333774, 389407, 574377, 328641, 1518244, 328642, 328643, 872663, 328644, 598794, 574379, 598795, 1454317, 1454343, 1454344, 328646, 131129, 328647, 328648, 574380, 325653, 152336, 334742, 340637, 334743, 882753, 331747, 334744, 331748, 131171, 331749, 331750, 331751, 863020, 155543, 338556, 331753, 1258309, 868182, 131173, 155542, 862728, 155128, 863012, 331754, 867870, 334746, 131174, 334747, 331755, 562502, 209865, 131175, 882754, 331756, 331758, 131176, 863019, 331760, 334749, 334750, 331761, 331762, 872730, 131109, 478336, 871997, 131110, 598908, 333060, 871996, 327029, 327030, 414035, 703807, 327240, 327241, 327242, 327243, 327244, 327245, 333117, 327246, 414036, 414037, 155466, 327247, 703809, 333118, 327248, 333119, 131111, 327249, 327250, 327251, 598851, 131113, 414038, 333120, 872571, 703808, 327252, 327254, 327354, 327355, 333173, 327356, 327357, 327358, 327359, 327360, 327361, 333174, 1456099, 327363, 326703, 513934, 329168, 333952, 598916, 156232, 333953, 737814, 329169, 390265, 882799, 333955, 131134, 606830, 326983, 606811, 598850, 606834, 131115, 606828, 332753, 332754, 606831, 606810, 332755, 1455940, 131116, 131117, 332756, 606835, 598852, 332757, 131118, 131119, 606812, 606833, 326984, 131120, 606472, 326985, 606813, 606832, 869234, 606814, 606809, 606829, 606815, 131121, 131123, 131124, 606816, 606474, 606475, 606473, 606817, 1047230, 888478, 888479, 326987, 1346270, 326988, 1347212, 882752, 888480, 888481, 332759, 326989, 326990, 888482, 326991, 326992, 131125, 332050, 334825, 131180, 596180, 338557, 332051, 332052, 338558, 332053, 332054, 332055, 335298, 338190, 339371, 332056, 340211, 558894, 339372, 131181, 152314, 131183, 131184, 332057, 338561, 707360, 332058, 597008, 334826, 131185, 152439, 332059, 338564, 131186, 337929, 174805, 1288528, 334827, 869121, 701702, 131191, 131192, 1421113, 1418669, 332348, 1346502, 1346514, 1292026, 1339469, 328671, 459112, 326143, 338501, 338502, 157535, 330902, 338504, 330903, 330904, 224519, 330906, 1541186, 334571, 334572, 1423926, 335284, 390183, 330912, 330913, 330915, 330916, 872669, 328844, 388034, 328845, 389493, 329222, 181523, 131136, 1547885, 329018, 1547883, 1379645, 329019, 329020, 329021, 1547884, 329022, 334018, 558896, 888618, 598912, 329023, 329024, 157490, 329025, 329026, 329027, 329028, 332082, 1447298, 608083, 608080, 332083, 608087, 334829, 332085, 332086, 1415952, 332087, 332088, 1446449, 608082, 332089, 608084, 210026, 332090, 608086, 332091, 332092, 332093, 608081, 332095, 567182, 332096, 332097, 332098, 478335, 608085, 174808, 338185, 152270, 332100, 332101, 331410, 1468693, 1544542, 334675, 331411, 331412, 1456120, 331413, 1449818, 598854, 331414, 331415, 1456130, 131167, 1544543, 1544544, 1544545, 1521915, 334677, 331416, 331417, 576831, 131169, 331419, 1450107, 331420, 1449846, 1449844, 1449847, 331421, 174819, 329872, 330375, 1469962, 1470571, 1310099, 1470553, 334391, 1469959, 1469961, 131140, 330376, 334392, 330377, 330425)

# find distributions from the list
#df_distr01 <- purrr::map_dfr(AphiaList,GetSpeciesDistributions)
# lave en vector af FALSE værdier med samme længde som AphiaList
SearchSynonyms <- rep(FALSE,length(AphiaList))
# find distributions from the list – uden at søge efter synonymer
df_distr02 <- purrr::map2_dfr(AphiaList,SearchSynonyms,GetSpeciesDistributions)

folder_out <- wd00_wd01
flNm01 <- "df_distr01.csv"
flNm02 <- "df_distr02.csv"
flout01<- paste0(folder_out,"/",flNm01)
flout02<- paste0(folder_out,"/",flNm02)
#write.table(df_distr01,file=flout01,col.names=T,row.names=F,sep=";",na="",fileEncoding="UTF-8")
write.table(df_distr02,file=flout02,col.names=T,row.names=F,sep=";",na="",fileEncoding="UTF-8")

# subset the list to only comprise the last 100 elements
#AphiaList <-  AphiaList[(noofAphIDs-500):noofAphIDs]
# find distributions from the list
#df_distr01 <- purrr::map_dfr(AphiaList,GetSpeciesDistributions)

# reorder the data frame by the column with 'Name'
df_MRG01 <- df_MRG01[order(df_MRG01$Name),]
#limit to the 'higherGeography' names that appear in the list regions
# close to Scandinavia
df_distr03 <- df_distr02[df_distr02$higherGeography %in%  df_MRG01$Name,]
#write out the data frame afterwards
flNm03<-"lst_NE_Europe_priority_spc.csv"
folder_out <- wd00_wd01
flout03<- paste0(folder_out,"/",flNm03)
write.table(df_distr03,file=flout03,row.names=F,col.names=T,sep=";",quote=F)
