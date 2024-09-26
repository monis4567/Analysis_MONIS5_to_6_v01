#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

########################################################################################
# R-code for making sigmoid curve plots from ABI7500 qPCR probe plots

# Make sigmoid curve qPCR plots from excel files with raw data 
# exported from ABI 7500
########################################################################################
########################################################################################
#remove everything in the working environment, without a warning!!
#rm(list=ls())
########################################################################################
#install packages
# #get readxl package
# if(!require(readxl)){
#   install.packages("readxl")
# }  
library(readxl)
# #get ggplot package
# if(!require(ggplot2)){
#   install.packages("ggplot2")
# }  
library(ggplot2)
# #get pdp package
# if(!require(pdp)){
#   install.packages("pdp")
# }  
library(pdp)
library(dplyr)
##########################################################################################
# begin - install packages to be able to do the ggplot below
##########################################################################################
# #get tidyverse package
# if(!require(tidyverse)){
#   install.packages("tidyverse")
# }  
library(tidyverse)
#get broom package
# if(!require(broom)){
#   install.packages("broom")
# }  
library(broom)
# #get mgcv package
# if(!require(mgcv)){
#   install.packages("mgcv")
# }  
library(mgcv)
# #get tibble package
# if(!require(tibble)){
#   install.packages("tibble")
# }  
library(tibble)
library(tidyverse)
# set working directory
#wd00 <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6/Results_from_ABI7500_for_MONIS6_species_specific_assays"
wd00 <- getwd()
setwd (wd00)
#define directory with output flies
wd01 <- "output01_qPCR_plots_for_qpcrruns_for_test_of_specificity"
# define full path for input directory
outdir01 <- paste(wd00,wd01, sep="/")
# remove previous versions of the 'outdir01'
unlink(outdir01, force = T, recursive = T)
# create the 'outdir01' in a new version
dir.create(outdir01)
# define directory with input files to read in
wddata <- "data"
#paste together a path for the directory that holds the inputs files
wd00_wddata <- paste(wd00,wddata,sep="/")
#library(broom)
#library(mgcv)  #For the gam model
##########################################################################################
# end - install packages to be able to do the ggplot below
##########################################################################################

# make a range of colours for the geom_points in the ggplots
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
clBpalt1 <- c("darkolivegreen3",
              "#8B3E2F",
              "darkorchid2",
              "darkslategray1",
              "#EE1289",
              "#458B00",
              "cornsilk4",
              "gold",
              "dodgerblue4",
              "#FFA54F",
              "yellow",
              "slateblue4")

##########################################################################################
# Note about the input files for this code
##########################################################################################
# The excel files are prepared in the ABI7500 software as 
# individual raw data needed to reproduce the amplification plots
# Each file can be exported from ABI7500 software as individual 
# excel spreadsheets , remember to export ALL data, both raw data, 
# setup, and sample data
# all spreadsheets can then be zipped together, and 
# transfered to your own computer
# unzip the zip file with all xls-files in a folder that also works as
# working directory
##########################################################################################

#list all files in wd - all the xls-files for which you want to 
# prepare plots from 
ls.fl01 <- list.files(wd00_wddata)
#make a variable with the element you want to search for
id1 <- "xls"
#grep for this variable in the list -  see this example: 
# https://stackoverflow.com/questions/35880242/r-selecting-element-from-list
ls.fl01.xls <- ls.fl01[grep(paste0(id1), ls.fl01)]

files <- ls.fl01.xls
# get the raw data files
qr.dtf <- files[grepl("^qpcr",files)]
qr.dtf <- qr.dtf[grepl("1020",qr.dtf)]
# get the filename withoutthe xls ending
fNm.qr <- gsub("\\.xls","",qr.dtf)
# get the setup  files
st.f <- files[grepl("^setup",files)]
st.f <- st.f[grepl("qpcr1020",st.f)]
# get the qpcr numbers
qpcrNos <- gsub("^qpcr([0-9]+).*","\\1",qr.dtf)
# get the qpcr rundate
qpcrRundt <- gsub("^qpcr([0-9]+).*_rundate([0-9]+)_.*","\\2",qr.dtf)
# get the setup no 
setupNo <- gsub("setup_qpcr([0-9]+).*","\\1",st.f)

# read in xls spreadsheet tab with raw flourescense data from ABI7500
rdt <- readxl::read_xls(paste0(wd00_wddata,"/",qr.dtf),
                        sheet="Raw Data",
                        skip=6)
# In the manual for ABI7500 in chapter 1, on page 2-3
# https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_050334.pdf
# the filters are specified as monitoring the following colors
# About the Filters The 7500/7500 Fast system uses the following filters:
## Filter 	|1						      |2		  	|3	    			|4				  |	5
## Dye	   	|	 FAM  dye		      | JOE  dye| TAMRA  dye	| ROX  dye	|	Cy5 dye
##          |	 SYBR  Green dye	| VIC dye	| NED dye 		| Texas Red	|
##          |						        |		    	| Cy3 dye 		|				    |

# The present setup makes use of the FAM dye probe monitored dye
# and ROX dye as background dye
# To get these colors in separate columns
rdt$FAM <- rdt$`1`
rdt$ROX <- rdt$`4`
# read in the xlsx file with the setup of the qpcr
stu.dt <- readxl::read_xlsx(paste0(wd00_wddata,"/",st.f))
# get the row number where the setup of well starts
rwn.wtbs <- which("Well"==stu.dt[,1])
# get a seq from where 'Well' appears and the next 96 rows
setuprws <- seq(rwn.wtbs,(rwn.wtbs+96))
# limit the setup file to only comprise these rows with the setup
qstup <- stu.dt[setuprws,]
# make the setup data frame a data frame instead of a tibble 
df_qs <- as.data.frame(qstup)
# use the 1st row as column names 
colnames(df_qs) <- df_qs[1,]
# get the data frame without the 1st row
df_qs <- df_qs[-1,]
# make a column that has 'WellNo' as column name, in the qpcr setup data frame
df_qs$WellNo <- df_qs$Well
# make a column that has 'WellNo' as column name, in the qpcr rundata data frame
rdt$WellNo  <- rdt$Well
# exclude columns that have NA as column name
ctk <- colnames(df_qs)[!is.na(colnames(df_qs))]
df_qs <- df_qs[ctk]
ctk <- colnames(rdt)[!is.na(colnames(rdt))]
rdt <- rdt[ctk]
# load the library that allows for combining data frames with left_join
library(dplyr)
# use left_join to combine the data frames
# https://cmdlinetips.com/2020/10/4-ways-to-select-columns-from-a-dataframe-with-dplyrs-select/
dfb01 <- dplyr::left_join(rdt,
                          # select all columns to include           
                          df_qs %>% dplyr::select(everything()),
                          by = "WellNo")
# get the difference between probe and the background dye
dfb01$dR <- dfb01$FAM-dfb01$ROX
# use the 'scale' function to normalize the data
# the scaled difference in flourescense represents the 
# difference in the dye monitored
dfb01$ddR <- scale(dfb01$dR)


# copy the column with primer and probe combination under a different
# name
dfb01$FRP.comb <- dfb01$`Primer and probe combination`
# exclude the rows where the 'FRP.comb' is NA#
# since these are empty wells without reagents added 
dfb01 <- dfb01[!is.na(dfb01$FRP.comb),]
# also exclude if the well name is NA, as these also are empty wells
dfb01 <- dfb01[!is.na(dfb01$`Well Name`),]

#get the unique assays - to use for facet wrap
unq.FRP.comb <- unique(dfb01$FRP.comb)
#make a table of the unique assays, and turn in to a data frame
tu_df <- as.data.frame(table(unq.FRP.comb))
#count the elements
cul <- length(unq.FRP.comb)
#make a sequence of numbers and append to the data frame
tu_df$cul <- 1:cul
#get the number of elements
seq.cul <- tu_df$cul
# copy columns to have the same column names as used for making the MxPro
# plots
dfb01$wllnm2 <- dfb01$`Well Name`
dfb01$Cycles <- dfb01$Cycle
dfb01$well <- dfb01$WellNo

# make a function that can make the text in the legend in italics
# https://stackoverflow.com/questions/59554096/ggplot2-italics-in-the-legend
toexpr <- function(x, plain = NULL) {
  getfun <- function(x) {
    ifelse(x == plain, "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}


dfb01$wllnm3 <- dfb01$wllnm2
Nms <- dfb01$wllnm3
# Find third occurrence of a special character and drop everything before that in R
# https://stackoverflow.com/questions/35088337/find-third-occurrence-of-a-special-character-and-drop-everything-before-that-in
Nms <- gsub('^(?:[^_]*_){1}','',Nms)
# also substitute underscore with space
Nms <- gsub('_',' ',Nms)
dfb01$wllnm3 <-  Nms
# re order the data frame by species names
dfb01 <- dfb01 %>% dplyr::arrange(wllnm3)
#https://stackoverflow.com/questions/31751022/a-function-to-create-multiple-plots-by-subsets-of-data-frame
library(ggplot2)
# Make plots of the amplification
plot01 <- ggplot(dfb01, aes(
  x = Cycles,
  y = ddR, 
  group= well, 
  color = wllnm3)) +
  geom_point() + 
  # use a different color scale -  check this webpage for examples : https://sjspielman.github.io/introverse/articles/color_fill_scales.html
  # scale_color_viridis_d(option = "inferno", direction = -1,
  #                       labels = toexpr(unique(dfb01$wllnm3), plain = 'Wt')) +
  scale_color_manual(values=clBpalt1,
                     labels = toexpr(unique(dfb01$wllnm3), plain = 'Wt')) +  
  
  facet_wrap(~FRP.comb, nrow = 3) + #'facet_wrap' subsets by column value in dataframe
  geom_line() + #add lines
  #labs(color='Extracted sample') + # change the label for the legend
  labs(color='Ekstraheret prøve fra') + # change the label for the legend
  labs(x = "qPCR cyklusser", y = "ddR") +
  ggtitle(fNm.qr) # add a title to the plot - here the filename is used
# modify the plot background : https://www.statology.org/ggplot-background-color/
#plt01 <- plot01 + theme_minimal()
plt01 <- plot01 + theme_bw()
plt01
# get the qpcr no to include in the output figure
qpcrNo <- qpcrNos
#set variable to define if figures are to be saved
bSaveFigures<-T
# save the figure if the above 'bSaveFigures' is TRUE
if(bSaveFigures==T){
  ggsave(plot = plt01, 
         # define the output filenmae by pasting together 
         # qpcrrunno and qpcrrundate
         filename = paste0(outdir01,"/Fig01_v01_qpcrrun",
                           qpcrNo,"rundate",qpcrRundt,".png"),
         width=210,height=297,
         #width=297,height=210,
         units="mm",dpi=300)
}

# Try again but make individual plots for each assay
# First identify the FRP combinations
uFRP.comb <- unique(dfb01$FRP.comb)
nuFRP.comb <- length(uFRP.comb)
# make a sequence of numbers for the FRP combinations
sqFRP.comb <- seq(1,nuFRP.comb,1)
# make an empty list that 
lst_plts_qPCR <- list()
# iterate over this sequence
for (e in sqFRP.comb)
{print(e)
  #}
  # subset the data frame to only comprise one FRP combinations 
  s.dfb01  <- dfb01[(dfb01$FRP.comb==uFRP.comb[e]),]
  # get the genus and species name for the samples that amplify
  dfb01abov2 <- s.dfb01[(s.dfb01$ddR>=2),]
  wellNm <- dfb01abov2$`Well Name`
  wellNm <- gsub(" ","_",wellNm)
  NmGenus <- stringr::str_split_i(wellNm,"_",2)
  NmSpeci <- stringr::str_split_i(wellNm,"_",3)
  GenusNm <- unique(NmGenus)
  SpeciNm <- unique(NmSpeci)
  GenusNm <- GenusNm[1]
  SpeciNm <- SpeciNm[1]
  # make a plot 
  plot01 <- ggplot(s.dfb01, aes(
    x = Cycles,
    y = ddR, 
    group= well, 
    color = wllnm3)) +
    geom_point() + 
    # use a different color scale -  check this webpage for examples : https://sjspielman.github.io/introverse/articles/color_fill_scales.html
    #scale_color_viridis_d(option = "plasma") +
    # scale_color_viridis_d(option = "plasma", direction = -1,
    #                       labels = toexpr(unique(dfb01$wllnm3), plain = 'Wt')) +
    # 
    scale_color_manual(values=clBpalt1,
                       labels = toexpr(unique(dfb01$wllnm3), plain = 'Wt')) +  
    #facet_wrap(~FRP.comb, nrow = 3) + #'facet_wrap' subsets by column value in dataframe
    geom_line() + #add lines
    #labs(color='Extracted sample') + # change the label for the legend
    labs(color='Ekstraheret prøve fra') + # change the label for the legend
    labs(x = "qPCR cyklusser", y = "ddR") +
    # modify the plot background : https://www.statology.org/ggplot-background-color/
    #plt01 <- plot01 + theme_minimal()
    theme_bw()
    #ggtitle(fNm.qr) # add a title to the plot - here the filename is used
    # or use the assay oligos for title
    #ggtitle(uFRP.comb[e])
  
  plt01 <- plot01
  # 3 store the plot in the list
  lst_plts_qPCR[[e]] <- plt01
  # get only the first part of the FRP combination name
  uFRP.comb.sh <- gsub("(.*),.*,.*","\\1",uFRP.comb[e])
  # evalute if the figure should be saved - just set the TRUE to FALSE
  # to avoid getting a diagram stored
  if(bSaveFigures==T){
    ggsave(plot = plt01, 
           # define the output filenmae by pasting together 
           # qpcrrunno and qpcrrundate
           filename = paste0(outdir01,"/Fig01_v0",e,"_",GenusNm,"_",SpeciNm,"_",uFRP.comb.sh,"_qpcrrun",
                             qpcrNo,"rundate",qpcrRundt,".png"),
           width=210,height=297*0.4,
           #width=297,height=210,
           units="mm",dpi=300)
  }
  # end of iteration over FRP combinations
}


#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________


# #install packages
library(readxl)
library(ggplot2)
library(pdp)
library(tidyverse)
library(broom)
library(mgcv)
library(tibble)


# make a range of colours for the geom_points in the ggplots
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
clBpalt1 <- c("darkolivegreen3",
              "#8B3E2F",
              "darkorchid2",
              "darkslategray1",
              "#EE1289",
              "#458B00",
              "cornsilk4",
              "gold",
              "dodgerblue4",
              "#40E0D0",
              "#FFA54F",
              "slateblue4",
              "#76EE00",
              "#00008B")

##########################################################################################
# Note about the input files for this code
##########################################################################################
# The excel files are prepared in the ABI7500 software as 
# individual raw data needed to reproduce the amplification plots
# Each file can be exported from ABI7500 software as individual 
# excel spreadsheets , remember to export ALL data, both raw data, 
# setup, and sample data
# all spreadsheets can then be zipped together, and 
# transfered to your own computer
# unzip the zip file with all xls-files in a folder that also works as
# working directory
##########################################################################################

#list all files in wd - all the xls-files for which you want to 
# prepare plots from 
ls.fl01 <- list.files(wd00_wddata)
#make a variable with the element you want to search for
id1 <- "xls"
#grep for this variable in the list -  see this example: 
# https://stackoverflow.com/questions/35880242/r-selecting-element-from-list
ls.fl01.xls <- ls.fl01[grep(paste0(id1), ls.fl01)]

files <- ls.fl01.xls
# get the raw data files
qr.dtf <- files[grepl("^qpcr",files)]
qr.dtf <- qr.dtf[grepl("1067",qr.dtf)]
# get the filename withoutthe xls ending
fNm.qr <- gsub("\\.xls","",qr.dtf)
# get the setup  files
st.f <- files[grepl("^setup",files)]
st.f <- st.f[grepl("1067",st.f)]
# get the qpcr numbers
qpcrNos <- gsub("^qpcr([0-9]+).*","\\1",qr.dtf)
# get the qpcr rundate
qpcrRundt <- gsub("^qpcr([0-9]+).*_rundate([0-9]+)_.*","\\2",qr.dtf)
qpcrRundt <- gsub("qpcr([0-9]+).*_rundate([0-9]+)_.*","\\2",st.f)
# get the setup no 
setupNo <- gsub("setup_qpcr([0-9]+).*","\\1",st.f)

# read in xls spreadsheet tab with raw flourescense data from ABI7500
rdt <- readxl::read_xls(paste0(wd00_wddata,"/",qr.dtf),
                        sheet="Raw Data",
                        skip=6)
# In the manual for ABI7500 in chapter 1, on page 2-3
# https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_050334.pdf
# the filters are specified as monitoring the following colors
# About the Filters The 7500/7500 Fast system uses the following filters:
## Filter 	|1						      |2		  	|3	    			|4				  |	5
## Dye	   	|	 FAM  dye		      | JOE  dye| TAMRA  dye	| ROX  dye	|	Cy5 dye
##          |	 SYBR  Green dye	| VIC dye	| NED dye 		| Texas Red	|
##          |						        |		    	| Cy3 dye 		|				    |

# The present setup makes use of the FAM dye probe monitored dye
# and ROX dye as background dye
# To get these colors in separate columns
rdt$FAM <- rdt$`1`
rdt$ROX <- rdt$`4`
# read in the xlsx file with the setup of the qpcr
stu.dt <- readxl::read_xlsx(paste0(wd00_wddata,"/",st.f))
# get the row number where the setup of well starts
rwn.wtbs <- which("Well"==stu.dt[,1])
# get a seq from where 'Well' appears and the next 96 rows
setuprws <- seq(rwn.wtbs,(rwn.wtbs+96))
# limit the setup file to only comprise these rows with the setup
qstup <- stu.dt[setuprws,]
# make the setup data frame a data frame instead of a tibble 
df_qs <- as.data.frame(qstup)

# use the 1st row as column names 
colnames(df_qs) <- df_qs[1,]
# get the data frame without the 1st row
df_qs <- df_qs[-1,]
# make a column that has 'WellNo' as column name, in the qpcr setup data frame
df_qs$WellNo <- df_qs$Well
# make a column that has 'WellNo' as column name, in the qpcr rundata data frame
rdt$WellNo  <- rdt$Well
# exclude columns that have NA as column name
ctk <- colnames(df_qs)[!is.na(colnames(df_qs))]
df_qs <- df_qs[ctk]
ctk <- colnames(rdt)[!is.na(colnames(rdt))]
rdt <- rdt[ctk]
# load the library that allows for combining data frames with left_join
library(dplyr)
# use left_join to combine the data frames
# https://cmdlinetips.com/2020/10/4-ways-to-select-columns-from-a-dataframe-with-dplyrs-select/
dfb01 <- dplyr::left_join(rdt,
                          # select all columns to include           
                          df_qs %>% dplyr::select(everything()),
                          by = "WellNo")
# get the difference between probe and the background dye
dfb01$dR <- dfb01$FAM-dfb01$ROX
# use the 'scale' function to normalize the data
# the scaled difference in flourescense represents the 
# difference in the dye monitored
dfb01$ddR <- scale(dfb01$dR)


# copy the column with primer and probe combination under a different
# name
dfb01$FRP.comb <- dfb01$`Primer and probe combination`
# exclude the rows where the 'FRP.comb' is NA#
# since these are empty wells without reagents added 
dfb01 <- dfb01[!is.na(dfb01$FRP.comb),]
# also exclude if the well name is NA, as these also are empty wells
dfb01 <- dfb01[!is.na(dfb01$`Well Name`),]

# find negative control wells
dfNKwlls <- dfb01[grepl("NK",dfb01$`Well Name`),]
# among the  negative control wells, find those that gave rise to 
# amplification above 2 in 'ddR'
dfNKwlls <- dfNKwlls[(dfNKwlls$ddR>2),]
# find unique well nos for the negative control wells,  that gave rise to 
# amplification above 2 in 'ddR'
NKwellNos.pos.ampl <- unique(dfNKwlls$Well.x)
#subset to exclude these wells
dfb01 <- dfb01[(dfb01$Well.x!=NKwellNos.pos.ampl),]
#get the unique assays - to use for facet wrap
unq.FRP.comb <- unique(dfb01$FRP.comb)
#make a table of the unique assays, and turn in to a data frame
tu_df <- as.data.frame(table(unq.FRP.comb))
#count the elements
cul <- length(unq.FRP.comb)
#make a sequence of numbers and append to the data frame
tu_df$cul <- 1:cul
#get the number of elements
seq.cul <- tu_df$cul
# copy columns to have the same column names as used for making the MxPro
# plots
dfb01$wllnm2 <- dfb01$`Well Name`
dfb01$Cycles <- dfb01$Cycle
dfb01$well <- dfb01$WellNo

# make a function that can make the text in the legend in italics
# https://stackoverflow.com/questions/59554096/ggplot2-italics-in-the-legend
toexpr <- function(x, plain = NULL) {
  getfun <- function(x) {
    ifelse(x == plain, "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}


dfb01$wllnm3 <- dfb01$wllnm2
Nms <- dfb01$wllnm3
# Find third occurrence of a special character and drop everything before that in R
# https://stackoverflow.com/questions/35088337/find-third-occurrence-of-a-special-character-and-drop-everything-before-that-in
Nms <- gsub('^(?:[^_]*_){1}','',Nms)
# also substitute underscore with space
Nms <- gsub('_',' ',Nms)
# also substitute Crefor with Crepidula fornicata
Nms <- gsub("Crefor[0-9]{3}(.*)","Crepidula fornicata\\1",Nms)
dfb01$wllnm3 <-  Nms

# re order the data frame by species names
dfb01 <- dfb01 %>% dplyr::arrange(wllnm3)
# get the genus and species name for the samples that amplify
dfb01abov2 <- dfb01[(dfb01$ddR>=2),]
wellNm <- dfb01abov2$`Well Name`
wellNm <- gsub(" ","_",wellNm)
NmGenus <- stringr::str_split_i(wellNm,"_",2)
NmSpeci <- stringr::str_split_i(wellNm,"_",3)
GenusNm <- unique(NmGenus)
SpeciNm <- unique(NmSpeci)
GenusNm <- GenusNm[1]
SpeciNm <- SpeciNm[1]
uFRP.comb <- unique(dfb01$FRP.comb)
#https://stackoverflow.com/questions/31751022/a-function-to-create-multiple-plots-by-subsets-of-data-frame
library(ggplot2)
# Make plots of the amplification
plot01 <- ggplot(dfb01, aes(
  x = Cycles,
  y = ddR, 
  group= well, 
  color = wllnm3)) +
  geom_point() + 
  # use a different color scale -  check this webpage for examples : https://sjspielman.github.io/introverse/articles/color_fill_scales.html
  #scale_color_viridis_d(option = "inferno") +
  scale_color_manual(values=clBpalt1,
                     labels = toexpr(unique(dfb01$wllnm3), plain = 'Wt')) +  
  
  #facet_wrap(~FRP.comb, nrow = 3) + #'facet_wrap' subsets by column value in dataframe
  geom_line() + #add lines
  #labs(color='Extracted sample') + # change the label for the legend
  
  #labs(color='Extracted sample') + # change the label for the legend
  labs(color='Ekstraheret prøve fra') + # change the label for the legend
  labs(x = "qPCR cyklusser", y = "ddR")
  #ggtitle(fNm.qr) # add a title to the plot - here the filename is used
  #ggtitle(uFRP.comb)
# modify the plot background : https://www.statology.org/ggplot-background-color/
plt01 <- plot01 + theme_minimal()
plt01 <- plot01 + theme_bw()
plt01
# get the qpcr no to include in the output figure
qpcrNo <- qpcrNos
#set variable to define if figures are to be saved
bSaveFigures<-T
# save the figure if the above 'bSaveFigures' is TRUE
if(bSaveFigures==T){
  ggsave(plot = plt01, 
         # define the output filenmae by pasting together 
         # qpcrrunno and qpcrrundate
         filename = paste0(outdir01,"/Fig02_v01_Crepidula_fornicata_","qpcrrun",
                           qpcrNo,"rundate",qpcrRundt,".png"),
         width=210,height=297*0.4,
         #width=297,height=210,
         units="mm",dpi=300)
}

#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________

#list all files in wd - all the xls-files for which you want to 
# prepare plots from 
ls.fl01 <- list.files(wd00_wddata)
#make a variable with the element you want to search for
id1 <- "xls"
#grep for this variable in the list -  see this example: 
# https://stackoverflow.com/questions/35880242/r-selecting-element-from-list
ls.fl01.xls <- ls.fl01[grep(paste0(id1), ls.fl01)]

files <- ls.fl01.xls
# get the raw data files
qr.dtf <- files[grepl("^qpcr",files)]
qr.dtf <- qr.dtf[grepl("1068",qr.dtf)]
# get the filename withoutthe xls ending
fNm.qr <- gsub("\\.xls","",qr.dtf)
# get the setup  files
st.f <- files[grepl("^setup",files)]
st.f <- st.f[grepl("1068",st.f)]
# get the qpcr numbers
qpcrNos <- gsub("^qpcr([0-9]+).*","\\1",qr.dtf)
# get the qpcr rundate
qpcrRundt <- gsub("^qpcr([0-9]+).*_rundate([0-9]+)_.*","\\2",qr.dtf)
qpcrRundt <- gsub("qpcr([0-9]+).*_rundate([0-9]+)_.*","\\2",st.f)
# get the setup no 
setupNo <- gsub("setup_qpcr([0-9]+).*","\\1",st.f)

# read in xls spreadsheet tab with raw flourescense data from ABI7500
rdt <- readxl::read_xls(paste0(wd00_wddata,"/",qr.dtf),
                        sheet="Raw Data",
                        skip=6)
# In the manual for ABI7500 in chapter 1, on page 2-3
# https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_050334.pdf
# the filters are specified as monitoring the following colors
# About the Filters The 7500/7500 Fast system uses the following filters:
## Filter 	|1						      |2		  	|3	    			|4				  |	5
## Dye	   	|	 FAM  dye		      | JOE  dye| TAMRA  dye	| ROX  dye	|	Cy5 dye
##          |	 SYBR  Green dye	| VIC dye	| NED dye 		| Texas Red	|
##          |						        |		    	| Cy3 dye 		|				    |

# The present setup makes use of the FAM dye probe monitored dye
# and ROX dye as background dye
# To get these colors in separate columns
rdt$FAM <- rdt$`1`
rdt$ROX <- rdt$`4`
# read in the xlsx file with the setup of the qpcr
stu.dt <- readxl::read_xlsx(paste0(wd00_wddata,"/",st.f))
# get the row number where the setup of well starts
rwn.wtbs <- which("Well"==stu.dt[,1])
# get a seq from where 'Well' appears and the next 96 rows
setuprws <- seq(rwn.wtbs,(rwn.wtbs+96))
# limit the setup file to only comprise these rows with the setup
qstup <- stu.dt[setuprws,]
# make the setup data frame a data frame instead of a tibble 
df_qs <- as.data.frame(qstup)

# use the 1st row as column names 
colnames(df_qs) <- df_qs[1,]
# get the data frame without the 1st row
df_qs <- df_qs[-1,]
# make a column that has 'WellNo' as column name, in the qpcr setup data frame
df_qs$WellNo <- df_qs$Well
# make a column that has 'WellNo' as column name, in the qpcr rundata data frame
rdt$WellNo  <- rdt$Well
# exclude columns that have NA as column name
ctk <- colnames(df_qs)[!is.na(colnames(df_qs))]
df_qs <- df_qs[ctk]
ctk <- colnames(rdt)[!is.na(colnames(rdt))]
rdt <- rdt[ctk]
# load the library that allows for combining data frames with left_join
library(dplyr)
# use left_join to combine the data frames
# https://cmdlinetips.com/2020/10/4-ways-to-select-columns-from-a-dataframe-with-dplyrs-select/
dfb01 <- dplyr::left_join(rdt,
                          # select all columns to include           
                          df_qs %>% dplyr::select(everything()),
                          by = "WellNo")
# get the difference between probe and the background dye
dfb01$dR <- dfb01$FAM-dfb01$ROX
# use the 'scale' function to normalize the data
# the scaled difference in flourescense represents the 
# difference in the dye monitored
dfb01$ddR <- scale(dfb01$dR)


# copy the column with primer and probe combination under a different
# name
dfb01$FRP.comb <- dfb01$`Primer and probe combination`
# exclude the rows where the 'FRP.comb' is NA#
# since these are empty wells without reagents added 
dfb01 <- dfb01[!is.na(dfb01$FRP.comb),]
# also exclude if the well name is NA, as these also are empty wells
dfb01 <- dfb01[!is.na(dfb01$`Well Name`),]

#get the unique assays - to use for facet wrap
unq.FRP.comb <- unique(dfb01$FRP.comb)
#make a table of the unique assays, and turn in to a data frame
tu_df <- as.data.frame(table(unq.FRP.comb))
#count the elements
cul <- length(unq.FRP.comb)
#make a sequence of numbers and append to the data frame
tu_df$cul <- 1:cul
#get the number of elements
seq.cul <- tu_df$cul
# copy columns to have the same column names as used for making the MxPro
# plots
dfb01$wllnm2 <- dfb01$`Well Name`
dfb01$Cycles <- dfb01$Cycle
dfb01$well <- dfb01$WellNo

# make a function that can make the text in the legend in italics
# https://stackoverflow.com/questions/59554096/ggplot2-italics-in-the-legend
toexpr <- function(x, plain = NULL) {
  getfun <- function(x) {
    ifelse(x == plain, "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}


dfb01$wllnm3 <- dfb01$wllnm2
Nms <- dfb01$wllnm3
# Find third occurrence of a special character and drop everything before that in R
# https://stackoverflow.com/questions/35088337/find-third-occurrence-of-a-special-character-and-drop-everything-before-that-in
Nms <- gsub('^(?:[^_]*_){1}','',Nms)
# also substitute underscore with space
Nms <- gsub('_',' ',Nms)
# also substitute Crefor with Crepidula fornicata
Nms <- gsub("Crefor[0-9]{3}(.*)","Crepidula fornicata\\1",Nms)
Nms <- gsub('(^(?:[^ ]* ){2}).*','\\1',Nms)
Nms <- gsub(' $','',Nms)
dfb01$wllnm3 <-  Nms


# re order the data frame by species names
dfb01 <- dfb01 %>% dplyr::arrange(wllnm3)

uFRP.comb <- unique(dfb01$FRP.comb)

FRP.combs <- unique(dfb01$FRP.comb)
noofFRP.combs <- length(FRP.combs)
sqf_FRP.comb <- seq(1,noofFRP.combs,1)
sbfigl <- letters[sqf_FRP.comb]
FRP.comb <- FRP.combs
df_sbfigl <- as.data.frame(cbind(sbfigl,FRP.comb))

dfb01 <- dfb01 %>% dplyr::left_join(df_sbfigl,
                                    by="FRP.comb")

#https://stackoverflow.com/questions/31751022/a-function-to-create-multiple-plots-by-subsets-of-data-frame
library(ggplot2)
# Make plots of the amplification
plot01 <- ggplot(dfb01, aes(
  x = Cycles,
  y = ddR, 
  group= well, 
  color = wllnm3)) +
  geom_point() + 
  theme_minimal() +
  theme_bw() +
  # use a different color scale -  check this webpage for examples : https://sjspielman.github.io/introverse/articles/color_fill_scales.html
  #scale_color_viridis_d(option = "inferno") +
  scale_color_manual(values=clBpalt1,
                     labels = toexpr(unique(dfb01$wllnm3), plain = 'Wt')) +  
  
  facet_wrap(~sbfigl+FRP.comb, nrow = 2 ,
             labeller = label_bquote(col = bold(.(paste0(sbfigl,") "))) ~ italic(.(FRP.comb)) ) 
  ) + #'facet_wrap' subsets by column value in dataframe
  # # see : https://r-charts.com/ggplot2/facets/
  theme(strip.text = element_text(#face = "bold",
    color = "black",
    hjust = 0,
    size = 10),
    strip.background = element_rect(fill = c("white"),
                                    #linetype = "solid",
                                    color = "white",
                                    linewidth = 1)) +
  geom_line() + #add lines
  labs(color='Extracted sample') + # change the label for the legend
  labs(color='Ekstraheret prøve fra') + # change the label for the legend
  labs(x = "qPCR cyklusser", y = "ddR") 
#ggtitle(fNm.qr) # add a title to the plot - here the filename is used
#ggtitle("Cteide_Co1 systemer")
# modify the plot background : https://www.statology.org/ggplot-background-color/

plt01 <- plot01

plt01
# get the qpcr no to include in the output figure
qpcrNo <- qpcrNos
#set variable to define if figures are to be saved
bSaveFigures<-T
# save the figure if the above 'bSaveFigures' is TRUE
if(bSaveFigures==T){
  ggsave(plot = plt01, 
         # define the output filenmae by pasting together 
         # qpcrrunno and qpcrrundate
         filename = paste0(outdir01,"/Fig03_Ctenopharyngodon_idella_qpcrrun",
                           qpcrNo,"rundate",qpcrRundt,".png"),
         width=210*1.6,height=297*0.6,
         #width=297,height=210,
         units="mm",dpi=300)
}

