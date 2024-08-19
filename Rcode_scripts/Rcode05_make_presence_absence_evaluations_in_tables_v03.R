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
#install.packages("fields")
# if(!require(fields)){
#   install.packages("fields")
#   library(fields)
# }
library(fields)
## install the package 'gplots', to be able to translate colors to hex - function: col2hex
# #install.packages("gplots")
# if(!require(gplots)){
#   install.packages("gplots")
#   library(gplots)
# }
library(gplots)
library(ggplot2)
## install the package 'glad', to be able to color using the function 'myPalette'
#install.packages("glad")
#library(glad)
require(graphics)
#get package to do count number of observations that have the same value at earlier records:
# see this website: https://stackoverflow.com/questions/11957205/how-can-i-derive-a-variable-in-r-showing-the-number-of-observations-that-have-th
# #install.packages("plyr")
# if(!require(plyr)){
#   install.packages("plyr")
#   library(plyr)
# }
library(plyr)

# if(!require(xlsx)){
#   install.packages("xlsx")
#   library(xlsx)
# }
library(xlsx)

#get package to make maps - see this website: http://www.molecularecologist.com/2012/09/making-maps-with-r/
#install.packages("mapdata")
#library(mapdata)
#get package to make maps - see this website: http://www.molecularecologist.com/2012/09/making-maps-with-r/
#install.packages("maps")
#library(maps)
# #get package for shapefiles see this website: http://www.molecularecologist.com/2012/09/making-maps-with-r/
# install.packages(maptools)
# library(maptools)  #for shapefiles
# #get package for adding pies on the map
#install.packages("mapplots")
#library(mapplots)
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
#   library(envDocument)
# }
library(envDocument)
# First (before making maps)  make sure all the required packages are loaded
# https://uchicagoconsulting.wordpress.com/tag/r-ggplot2-maps-visualization/
# #install packages needed
# if(!require(maps)){
#   install.packages("maps")
# }
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
# }
library(cowplot)

# if(!require(googleway)){
#   install.packages("googleway")
# }
library(googleway)

# if(!require(ggrepel)){
#   install.packages("ggrepel")
# }
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
# }
library(sf)

# if(!require(rnaturalearth)){
#   install.packages("rnaturalearth")
# }
library(rnaturalearth)

# if(!require(rnaturalearthdata)){
#   install.packages("rnaturalearthdata")
# }
library(rnaturalearthdata)

# if(!require(ggforce)){
#   install.packages("ggforce")
#   
# }
library(ggforce)
# #get 'rnaturalearthhires' installed
# if(!require(rnaturalearthhires)){
#   #install.packages("rnaturalearthhires")
#   install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
#   
# }
library(rnaturalearthhires)
# # 
library("ggplot2")

library("sf")
library(ggforce)


#install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
# # 
library("rnaturalearth")
library("rnaturalearthdata")
library("rnaturalearthhires")

# define working directory
wd00 <- getwd()
#wd00 ="/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6/Analysis_MONIS5_to_6"
#wd00 <- "/Users/steenknudsen/Documents/Documents/NIVA_Ansaettelse_2020/NOVANA_proever_2018_2019/"
# define input directory
wd05 <- "output05_color_tables_from_Rcode_for_MST2017_2022samples"
# define out directory
wd06 <- "output06_presence_absence_evaluation_for_MST2017_2022_samples"
#paste dirs together
wd00_wd06 <- paste0(wd00,"/",wd06)
#Delete any previous versions of the output directory
unlink(wd00_wd06, recursive=TRUE)
#Create a directory to put resulting output files in
dir.create(wd00_wd06)
# set working directory
#setwd(wd00)
#getwd()
#paste dirs together
wd00_wd05 <- paste0(wd00,"/",wd05)

#_______________________________________________________________________________
# section 01 - start - read in and merge data frames 
#_______________________________________________________________________________
# # read in table prepared from previous code
dfg01 <- read.table(paste0(wd00_wd05,"/table07_out_smpls08.csv"), sep=",",
                  header = T)

# # use dplyr to only keep distinct rows, also add the '.keep_all = TRUE'
# # to retain all columns
dfg01 <- dfg01 %>% dplyr::distinct(qpcr.rplno,Lat_Species,MST.n_Dinds, 
                                   .keep_all = TRUE)


# Ensrue all data frames are data frames
dfg01 <- as.data.frame(dfg01,stringsAsFactors = FALSE)
dfg01[] <- lapply(dfg01, as.character)
# make a common column to use for matching columns up

dfg01$MSTD.LS.QN <- paste(dfg01$MST.n_Dinds,dfg01$Lat_Species,dfg01$qpcrno,sep="_")
dfg01.MST2021028 <- dfg01[grepl("MST2021028",dfg01$MST.n_Dinds),]


dfg03 <- dfg01

# check which detections of eDNA can be evaluated as being above 
# LOD, and evaluate on whether the organism was found
dfg03$NISfnd <- (dfg03$AbLOD_BeLOQ>0 | dfg03$aboveLOQ>0)
# evaluate if the organism was found 
dfg03$orgF <- dfg03$NISfnd


# subset to exclude Standards and NTC , so that only unknowns are retained
dfg03 <- dfg03[(dfg03$welltp=="Unknown"),]
#
ckeep <- c("yer","ssn","declon","declat","orgF","NISfnd",
           "MST.n_Dinds","Lat_Species","qpcrno")
#_______________________________________________________________________________
# section 01 - end - read in and merge data frames 
#_______________________________________________________________________________
#_______________________________________________________________________________
# section 02 - start - read in world map 
#_______________________________________________________________________________

#_______________________________________________________________________________
# First (before making maps)  make sure all the required packages are loaded
# https://uchicagoconsulting.wordpress.com/tag/r-ggplot2-maps-visualization/
# #install packages needed
# if(!require(maps)){
#   install.packages("maps")
# }
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
# }
library(cowplot)
# 
# if(!require(googleway)){
#   install.packages("googleway")
# }
library(googleway)

# if(!require(ggrepel)){
#   install.packages("ggrepel")
# }
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
# }
library(sf)

# if(!require(rnaturalearth)){
#   install.packages("rnaturalearth")
# }
library(rnaturalearth)

# if(!require(rnaturalearthdata)){
#   install.packages("rnaturalearthdata")
# }
library(rnaturalearthdata)

# if(!require(ggforce)){
#   install.packages("ggforce")
# 
# }
library(ggforce)
# #get 'rnaturalearthhires' installed
# if(!require(rnaturalearthhires)){
#   #install.packages("rnaturalearthhires")
#   install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
#   
# }
library(rnaturalearthhires)
# # 
library("ggplot2")

library("sf")
library(ggforce)


#install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
# # 
library("rnaturalearth")
library("rnaturalearthdata")
library("rnaturalearthhires")


theme_set(theme_bw())
# # Get a map, use a high number for 'scale' for a coarse resolution
# use a low number for scale for a high resolution
# if the map 'world' does not exist, then download it
world <- ne_countries(scale = 10, returnclass = "sf")
library(ggplot2)
#_______________________________________________________________________________
# section 02 - end - read in world map 
#_______________________________________________________________________________


#_______________________________________________________________________________
# section 03 - start - plot on map the number of detected species per year 
#_______________________________________________________________________________

# the ggplot with facet wrap will be without
# a map for spring 2017. To correct for this, I will add an empty row
# with NA for "spring 2017"
# copy the data frame

# but first remove the NAs
dfg05 <-dfg03[!(is.na(dfg03$yea)),]
dfg05$yer <- dfg05$yea
dfg05_MST0267_2020.06.16 <-  dfg05[grepl("MST0267_2020-06-16",dfg05$MST.n_Dinds),]
dfg05_MST0267_2020.06.16 <- dfg05_MST0267_2020.06.16[order(dfg05_MST0267_2020.06.16$Lat_Species),]

# View(dfg05_MST0267_2020.06.16)
# No samples were collected in spring 2017
# Which means that ggplot with facet wrap will be without
# a map for spring 2017. To correct for this, I will add an empty row
# with NA for "spring 2017"
# First get the columns name
clnmdfg <- colnames(dfg05)
#head(dfg02.nOF,3)
# get the number of columns
ncldfg <- length(clnmdfg)
# identify the column number for the columns named 'yer' and 'ssn'
nycl <- which(clnmdfg=="yer") 
nscl <- which(clnmdfg=="ssn")
nsNISfcl <- which(clnmdfg=="NISfnd")
# make an empty vector with only NAs
nwrw <- rep(0,ncldfg)
# then replace in the column numbers that will be year and season
nwrw[nycl] <- 2017
nwrw[nscl] <- "spri"
nwrw[nsNISfcl] <- F
# make the tibble a data frame , to be able to append with 'rbind'
dfg05 <- as.data.frame(dfg05)
# then add an empty row , that has been filled for year and season
dfg05 <- rbind(dfg05,nwrw)

# ##______________________________________________________________________________
# # subset to try out on only the year 2018
# dfg05 <- dfg05[(dfg05$yer==2018),]
# ##______________________________________________________________________________

# identify column index numbers
nlon<- which(grepl("declon",clnmdfg))
nlat<- which(grepl("declat",clnmdfg))
nyer<- which(grepl("yer",clnmdfg))
nlva<- which(grepl("vanda",clnmdfg))
# make sure all columns that need to be numeric are numeric
cols_to_convert <- clnmdfg[c(nlon,nlat,nyer,nlva)]
dfg05[cols_to_convert]  <- lapply(dfg05[cols_to_convert], as.numeric)
# make the column with organism found a boolean
#https://stackoverflow.com/questions/33930188/convert-dataframe-column-to-1-or-0-for-true-false-values-and-assign-to-dataf
dfg05$orgFnd <- as.logical(dfg05$orgF)

dfg05$orgFnd2 <- (dfg05$orgFnd*1)
# Make an empty column
dfg05$ssnno <- NA
# evalutate on the season column to add a number version of the season
dfg05$ssnno[(dfg05$ssn=="spri")] <- "1st"
dfg05$ssnno[(dfg05$ssn=="fall")] <- "2nd"
# paste together year and the numbered part of the year
dfg05$yer_ssn <- paste0(dfg05$yer,", ",dfg05$ssnno," halvår")

#unique(dfg05$yer_ssn)
# subset to only have the TRUE detections
dfg07 <- dfg05[(dfg05$orgFnd2>0),]

dfg07.MST2021028 <- dfg07[grepl("MST2021028",dfg07$MST.n_Dinds),]
#nrow(dfg07.MST2021028)
#unique(dfg07$orgFnd2)
dfg07 <- dfg07[!(is.na(dfg07$orgFnd2)),]
# sum a variable by group using dplyr
library(dplyr)
# https://stackoverflow.com/questions/1660124/how-to-sum-a-variable-by-group
# summarise specific variables, not all
dfg07.01 <- dfg07 %>% 
  dplyr::group_by(MST.n_Dinds,
                  declon,
                  declat ,
                  yer,
                  #lokalitet_vanda,
                  yer_ssn
  ) %>% 
  dplyr::summarise(across(c(orgFnd2), list(moF2 = mean, soF2 = sum)))

dfg07.01l7 <- dfg07.01[(dfg07.01$orgFnd2_soF2>=7),]
#View(dfg07.01l7)
# the dplyr sum per category makes a tibble, and makes columns contain 
# characters, these need to be numeric to be able to plot
dfg07.01 <- as.data.frame(dfg07.01)
clnmdfg <- colnames(dfg07.01)
# identify column index numbers
nlon<- which(grepl("declon",clnmdfg))
nlat<- which(grepl("declat",clnmdfg))
nyer<- which("yer"==clnmdfg)
nlva<- which(grepl("vanda",clnmdfg))
# make sure all columns that need to be numeric are numeric
cols_to_convert <- clnmdfg[c(nlon,nlat,nyer,nlva)]
dfg07.01[cols_to_convert]  <- lapply(dfg07.01[cols_to_convert], as.numeric)

# exclude rows with NAs  using tidyr
# https://stackoverflow.com/questions/4862178/remove-rows-with-all-or-some-nas-missing-values-in-data-frame
library(tidyr)
dfg07.01 <- dfg07.01 %>% tidyr::drop_na()

# Find the maximum number of sum of organisms found 
mxorgF <- max(dfg07.01$orgFnd2_soF2)
orgFrng <- seq(0,mxorgF,1)
# make a color range
CLscl <- colorRampPalette(RColorBrewer::brewer.pal(
  mxorgF, "BuPu"))(length(orgFrng)-1)

# make a color range
CLscl <- colorRampPalette(RColorBrewer::brewer.pal(
  mxorgF, "YlGnBu"))(length(orgFrng)-1)
# make another color range
CLscl <- colorRampPalette(c("white","yellow","seagreen2","deepskyblue","dodgerblue1"))(mxorgF)

# make the color range a data frame to match number of organisms found
dfCScl <- as.data.frame(cbind(orgFrng,CLscl))
# match the number of organisms found to get the color  
dfg07.01$CLscl <- dfCScl$CLscl[match(dfg07.01$orgFnd2_soF2,dfCScl$orgFrng)]


dfg05$yer_ssn2 <- dfg05$yer_ssn
dfg05$yer_ssn2 <- gsub("1st halvår","jan - jun",dfg05$yer_ssn2)
dfg05$yer_ssn2 <- gsub("2nd halvår","jul - nov",dfg05$yer_ssn2)

dfg07.01$yer_ssn2 <- dfg07.01$yer_ssn
dfg07.01$yer_ssn2 <- gsub("1st halvår","jan - jun",dfg07.01$yer_ssn2)
dfg07.01$yer_ssn2 <- gsub("2nd halvår","jul - nov",dfg07.01$yer_ssn2)

library(ggplot2)
# clear up memory
gc()
#make plot
p05 <- ggplot(data = world) +
  geom_sf(color = "black", fill = "azure3", lwd=0.1, stroke=0.1) +
  # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  theme_void() +
  # Plot points for TRUE find of NIS (i.e. above LOD) 
  # using 'geom_pointdensity' to get a heat color for overlaid points
  
  # geom_pointdensity(data=dfg07, 
  #                   mapping= aes(
  #                     x =declon, 
  #                     y = declat) ,
  #                   show.legend = T,
  #                   size=5) +
  # 
  
# set the color scale for the  'geom_pointdensity' for 
# the heat color for overlaid points
#scale_color_viridis_c() +

# plot all points with crosses, to show were sampling was performed
# NOTE that this part that plots the crosses as placed first, so that 
# the next layers with numbered points are placed on top, and will cover
# the crosses
geom_point(data=dfg05, 
           aes( x =declon, 
                y = declat),
           shape=3) +
  # add points for sampled locations with the number of of organisms found
  # on the location, NOTE that the variable needs to be modified with
  # 'as.factor()' 
  # see this question : https://stackoverflow.com/questions/43359050/error-continuous-value-supplied-to-discrete-scale-in-default-data-set-example
  geom_point(data= dfg07.01,
             aes( x=declon,
                  y=declat,
                  fill=as.factor(orgFnd2_soF2)),
             shape= 21,
             stroke = 0.1,
             size=3.2) +
  # use the color range prepared above
  # scale_fill_manual(values = alpha(c(CLscl),0.8)) +
  scale_fill_manual(values = c(CLscl)) +
  # add text labels on the point 
  geom_text(data= dfg07.01, 
            aes( x=declon,
                 y=declat,
                 label=orgFnd2_soF2),
            size=1.8) +
  #Arrange in facets
  # ggplot2::facet_wrap( ~ yer_ssn,
  #                      drop=FALSE,
  #                      dir="h",
  #                      ncol = 2,
  #                      labeller = label_bquote(cols =
  #                                                italic(.(as.character(yer_ssn))))  ) +
  # 
  ggplot2::facet_wrap( ~ yer_ssn2,
                       drop=FALSE,
                       dir="h",
                       ncol = 2) +
  
  # # see : https://r-charts.com/ggplot2/facets/
  theme(strip.text = element_text(#face = "bold",
    color = "black",
    hjust = 0,
    size = 8),
    strip.background = element_rect(fill = c("white"),
                                    #linetype = "solid",
                                    color = "white",
                                    linewidth = 1)) +
  
  # place legend on bottom
  #theme(legend.position = "bottom") +
  # or remove the legend
  #theme(legend.position = "none") +
  #define limits of the plot 
  ggplot2::coord_sf(xlim = c(6.6, 17.2),
                    ylim = c(54.2, 58.4), 
                    expand = FALSE) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # or remove them all completely
  #http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()) 
#change axis labels
#p05t <- p05 + xlab("longitude") + ylab("latitude")
#p05t <- p05 + xlab("længdegrad") + ylab("breddegrad")
p05t <- p05 + xlab(" ") + ylab(" ")
#change the header for the legend on the side, 
#this must be done for both 'fill', 'color' and 'shape', to avoid 
#getting separate legends
p05t <- p05t + labs(color='')
p05t <- p05t + labs(fill='Antal NIS')
p05t <- p05t + labs(shape='')
p05t <- p05t + labs(size='')
# https://github.com/tidyverse/ggplot2/issues/3492
#adjust tick marks on axis
p05t <- p05t + scale_y_continuous(breaks=seq(54.0, 58.4,1))
p05t <- p05t + scale_x_continuous(breaks=seq(6.0, 17.2,2.0))
# see the plot
p05t
bSaveFigures<-T
if(bSaveFigures==T){
  ggsave(plot = p05t, 
         filename = paste0(wd00_wd06,"/Fig05_v01_map_of_NIS_detected_2017_to_2023.png"),
         width=210*0.64,height=297,
         #width=210*0.8,height=297,
         #width=297,height=210,
         #width=297,height=210,
         #width=1.4*297,height=210,
         units="mm",dpi=300)
}


library(ggplot2)
#make plot
p05 <- ggplot(data = world) +
  geom_sf(color = "black", fill = "azure3", lwd=0.1, stroke=0.1) +
  # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  theme_void() +
# plot all points with crosses, to show were sampling was performed
# NOTE that this part that plots the crosses as placed first, so that 
# the next layers with numbered points are placed on top, and will cover
# the crosses
geom_point(data=dfg05, 
           aes( x =declon, 
                y = declat),
           shape=3) +
  # add points for sampled locations with the number of of organisms found
  # on the location, NOTE that the variable needs to be modified with
  # 'as.factor()' 
  # see this question : https://stackoverflow.com/questions/43359050/error-continuous-value-supplied-to-discrete-scale-in-default-data-set-example
  geom_point(data= dfg07.01,
             aes( x=declon,
                  y=declat,
                  fill=as.factor(orgFnd2_soF2)),
             shape= 21,
             stroke = 0.1,
             size=3.2) +
  # use the color range prepared above
  # scale_fill_manual(values = alpha(c(CLscl),0.8)) +
  scale_fill_manual(values = c(CLscl)) +
  # add text labels on the point 
  geom_text(data= dfg07.01, 
            aes( x=declon,
                 y=declat,
                 label=orgFnd2_soF2),
            size=1.8) +
  #Arrange in facets
  # ggplot2::facet_wrap( ~ yer_ssn,
  #                      drop=FALSE,
  #                      dir="h",
  #                      ncol = 2,
  #                      labeller = label_bquote(cols =
  #                                                italic(.(as.character(yer_ssn))))  ) +
  # 
  ggplot2::facet_wrap( ~ yer_ssn,
                       drop=FALSE,
                       dir="h") +
  
  # # see : https://r-charts.com/ggplot2/facets/
  theme(strip.text = element_text(#face = "bold",
    color = "black",
    hjust = 0,
    size = 8),
    strip.background = element_rect(fill = c("white"),
                                    #linetype = "solid",
                                    color = "black",
                                    linewidth = 1)) +
  # place legend on bottom
  #theme(legend.position = "bottom") +
  # or remove the legend
  #theme(legend.position = "none") +
  #define limits of the plot 
  ggplot2::coord_sf(xlim = c(6.6, 17.2),
                    ylim = c(54.2, 58.4), 
                    expand = FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # or remove them all completely
  #http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()) 
#change axis labels
#p05t <- p05 + xlab("longitude") + ylab("latitude")
p05t <- p05 + xlab("længdegrad") + ylab("breddegrad")
p05t <- p05 + xlab(" ") + ylab(" ")
#change the header for the legend on the side, 
#this must be done for both 'fill', 'color' and 'shape', to avoid 
#getting separate legends
p05t <- p05t + labs(color='')
p05t <- p05t + labs(fill='Antal NIS')
p05t <- p05t + labs(shape='')
p05t <- p05t + labs(size='')
# https://github.com/tidyverse/ggplot2/issues/3492
#adjust tick marks on axis
p05t <- p05t + scale_y_continuous(breaks=seq(54.0, 58.4,1))
p05t <- p05t + scale_x_continuous(breaks=seq(6.0, 17.2,2.0))
# see the plot
p05t
bSaveFigures<-T
if(bSaveFigures==T){
  ggsave(plot = p05t, 
         filename = paste0(wd00_wd06,"/Fig05_v02_map_of_NIS_detected_2017_to_2023.png"),
         #width=210*0.64,height=297,
         #width=210*0.8,height=297,
         #width=210,height=297,
         width=210,height=0.70*297,
         #width=297,height=210,
         #width=1.4*297,height=210,
         units="mm",dpi=300)
}

dfg05_2017_2020 <- dfg05[grepl("2017|2018|2019|2020",dfg05$yer_ssn),]
dfg07.01_2017_2020 <- dfg07.01[grepl("2017|2018|2019|2020",dfg07.01$yer_ssn),]
library(ggplot2)
#make plot
p05 <- ggplot(data = world) +
  geom_sf(color = "black", fill = "azure3", lwd=0.1, stroke=0.1) +
  # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  theme_void() +
  # plot all points with crosses, to show were sampling was performed
  # NOTE that this part that plots the crosses as placed first, so that 
  # the next layers with numbered points are placed on top, and will cover
  # the crosses
  geom_point(data=dfg05_2017_2020, 
             aes( x =declon, 
                  y = declat),
             shape=3) +
  # add points for sampled locations with the number of of organisms found
  # on the location, NOTE that the variable needs to be modified with
  # 'as.factor()' 
  # see this question : https://stackoverflow.com/questions/43359050/error-continuous-value-supplied-to-discrete-scale-in-default-data-set-example
  geom_point(data= dfg07.01_2017_2020,
             aes( x=declon,
                  y=declat,
                  fill=as.factor(orgFnd2_soF2)),
             shape= 21,
             stroke = 0.1,
             size=6.2) +
  # use the color range prepared above
  # scale_fill_manual(values = alpha(c(CLscl),0.8)) +
  scale_fill_manual(values = c(CLscl)) +
  # add text labels on the point 
  geom_text(data= dfg07.01_2017_2020, 
            aes( x=declon,
                 y=declat,
                 label=orgFnd2_soF2),
            size=3.8) +
  #Arrange in facets
  # ggplot2::facet_wrap( ~ yer_ssn,
  #                      drop=FALSE,
  #                      dir="h",
  #                      ncol = 2,
  #                      labeller = label_bquote(cols =
  #                                                italic(.(as.character(yer_ssn))))  ) +
  # 
  ggplot2::facet_wrap( ~ yer_ssn,
                       drop=FALSE,
                       dir="h",
                       ncol = 2) +
  
  # # see : https://r-charts.com/ggplot2/facets/
  theme(strip.text = element_text(#face = "bold",
    color = "black",
    hjust = 0,
    size = 8),
    strip.background = element_rect(fill = c("white"),
                                    #linetype = "solid",
                                    color = "black",
                                    linewidth = 1)) +
  # place legend on bottom
  #theme(legend.position = "bottom") +
  # or remove the legend
  #theme(legend.position = "none") +
  #define limits of the plot 
  ggplot2::coord_sf(xlim = c(6.6, 17.2),
                    ylim = c(54.2, 58.4), 
                    expand = FALSE) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # or remove them all completely
  #http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()) 
#change axis labels
#p05t <- p05 + xlab("longitude") + ylab("latitude")
#p05t <- p05 + xlab("længdegrad") + ylab("breddegrad")
p05t <- p05 + xlab(" ") + ylab(" ")
#change the header for the legend on the side, 
#this must be done for both 'fill', 'color' and 'shape', to avoid 
#getting separate legends
p05t <- p05t + labs(color='')
p05t <- p05t + labs(fill='Antal NIS')
p05t <- p05t + labs(shape='')
p05t <- p05t + labs(size='')
# https://github.com/tidyverse/ggplot2/issues/3492
#adjust tick marks on axis
p05t <- p05t + scale_y_continuous(breaks=seq(54.0, 58.4,1))
p05t <- p05t + scale_x_continuous(breaks=seq(6.0, 17.2,2.0))
# see the plot
p05t
bSaveFigures<-T
if(bSaveFigures==T){
  ggsave(plot = p05t, 
         filename = paste0(wd00_wd06,"/Fig05_v03_map_of_NIS_detected_2017_to_2020.png"),
         #width=210*0.64,height=297,
         #width=210*0.8,height=297,
         #width=210,height=297,
         width=210,height=297,
         #width=297,height=210,
         #width=1.4*297,height=210,
         units="mm",dpi=300)
}
dfg05_sbs <- dfg05[grepl("2021|2022|2023",dfg05$yer_ssn),]
dfg07.01_sbs <- dfg07.01[grepl("2021|2022|2023",dfg07.01$yer_ssn),]
library(ggplot2)
#make plot
p05 <- ggplot(data = world) +
  geom_sf(color = "black", fill = "azure3", lwd=0.1, stroke=0.1) +
  # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  theme_void() +
  # plot all points with crosses, to show were sampling was performed
  # NOTE that this part that plots the crosses as placed first, so that 
  # the next layers with numbered points are placed on top, and will cover
  # the crosses
  geom_point(data=dfg05_sbs, 
             aes( x =declon, 
                  y = declat),
             shape=3) +
  # add points for sampled locations with the number of of organisms found
  # on the location, NOTE that the variable needs to be modified with
  # 'as.factor()' 
  # see this question : https://stackoverflow.com/questions/43359050/error-continuous-value-supplied-to-discrete-scale-in-default-data-set-example
  geom_point(data= dfg07.01_sbs,
             aes( x=declon,
                  y=declat,
                  fill=as.factor(orgFnd2_soF2)),
             shape= 21,
             stroke = 0.1,
             size=6.2) +
  # use the color range prepared above
  # scale_fill_manual(values = alpha(c(CLscl),0.8)) +
  scale_fill_manual(values = c(CLscl)) +
  # add text labels on the point 
  geom_text(data= dfg07.01_sbs, 
            aes( x=declon,
                 y=declat,
                 label=orgFnd2_soF2),
            size=3.8) +
  #Arrange in facets
  # ggplot2::facet_wrap( ~ yer_ssn,
  #                      drop=FALSE,
  #                      dir="h",
  #                      ncol = 2,
  #                      labeller = label_bquote(cols =
  #                                                italic(.(as.character(yer_ssn))))  ) +
  # 
  ggplot2::facet_wrap( ~ yer_ssn,
                       drop=FALSE,
                       dir="h",
                       ncol = 2) +
  
  # # see : https://r-charts.com/ggplot2/facets/
  theme(strip.text = element_text(#face = "bold",
    color = "black",
    hjust = 0,
    size = 8),
    strip.background = element_rect(fill = c("white"),
                                    #linetype = "solid",
                                    color = "black",
                                    linewidth = 1)) +
  # place legend on bottom
  #theme(legend.position = "bottom") +
  # or remove the legend
  #theme(legend.position = "none") +
  #define limits of the plot 
  ggplot2::coord_sf(xlim = c(6.6, 17.2),
                    ylim = c(54.2, 58.4), 
                    expand = FALSE) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # or remove them all completely
  #http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()) 
#change axis labels
#p05t <- p05 + xlab("longitude") + ylab("latitude")
#p05t <- p05 + xlab("længdegrad") + ylab("breddegrad")
p05t <- p05 + xlab(" ") + ylab(" ")
#change the header for the legend on the side, 
#this must be done for both 'fill', 'color' and 'shape', to avoid 
#getting separate legends
p05t <- p05t + labs(color='')
p05t <- p05t + labs(fill='Antal NIS')
p05t <- p05t + labs(shape='')
p05t <- p05t + labs(size='')
# https://github.com/tidyverse/ggplot2/issues/3492
#adjust tick marks on axis
p05t <- p05t + scale_y_continuous(breaks=seq(54.0, 58.4,1))
p05t <- p05t + scale_x_continuous(breaks=seq(6.0, 17.2,2.0))
# see the plot
p05t
bSaveFigures<-T
if(bSaveFigures==T){
  ggsave(plot = p05t, 
         filename = paste0(wd00_wd06,"/Fig05_v04_map_of_NIS_detected_2021_to_2023.png"),
         #width=210*0.64,height=297,
         #width=210*0.8,height=297,
         #width=210,height=297,
         width=210,height=297*0.85,
         #width=297,height=210,
         #width=1.4*297,height=210,
         units="mm",dpi=300)
}

#_______________________________________________________________________________
# section 03 - end - plot on map the number of detected species per year 
#_______________________________________________________________________________
# Only keep unique rows for 'MST.n_Dinds', and 'lokalitet_vanda' and 
# 'lok_pos_lat' and 'lok_pos_lon'
dfg05$lok_pos_lat <- dfg05$declat
dfg05$lok_pos_lon <- dfg05$declon
dfg05$MST.nummer <- dfg05$MST.n_Dinds
dfg05$Dato_inds <- gsub("(^.*)_(.*$)","\\2",dfg05$MST.n_Dinds)
# This only retains the rows that represents unique sample locations per year 
dfg08 <- dfg05 %>% dplyr::distinct( MST.n_Dinds,
                                    lokalitet_vanda,
                                    declat,
                                    declon,
                                    volfilt_mL,
                                    .keep_all = TRUE)
# only keep the columns needed for presenting the sampled locations
dfg08 <- dfg08 %>%
  dplyr::select(MST.nummer,
                Dato_inds,
                volfilt_mL,
                lokalitet_vanda,
                lok_pos_lat,
                lok_pos_lon)
# select only the columns that are to be made numeric, and make them
# numeric
dfg08.ncl <- 
  dfg08 %>% dplyr::select(lokalitet_vanda,
                          volfilt_mL,
                          lok_pos_lat,
                          lok_pos_lon) %>% 
  dplyr::mutate_if(is.character, as.numeric)
# select the other columns
dfg08.ccl <- 
  dfg08 %>% dplyr::select(MST.nummer,
                          Dato_inds)
# bind together in columns the numeric and the non numeric columns
dfg08 <- cbind(dfg08.ncl,dfg08.ccl)
# Make a function that can round numbers and apply on a vector comprising
# the column names
rndNmb <- function(nTF){as.numeric(format(round(nTF, 4), nsmall = 4))}
clMtM <-c("lok_pos_lat", "lok_pos_lon")
# apply the function prepared above
dfg08[, clMtM] <- lapply(clMtM, function(x) rndNmb(dfg08[[x]]))
# get the column names
clNm <- colnames(dfg08)
# identify column index numbers
ixcMST<- which(grepl("MST",clNm))
ixcDI <- which(grepl("Dato_inds",clNm))
# identify column numbers not yet identified
otcl <- seq(1:length(clNm)) [!(seq(1:length(clNm)) %in% c(ixcMST,ixcDI))]
nOrd <- c(ixcMST,ixcDI,otcl)
# reorder the data frame
dfg08 <- dfg08[,nOrd]
# reorder the data frame
dfg08 <- dfg08[order(dfg08$Dato_inds),]
# do NOT include if 'Dato_inds' is 0
dfg08 <- dfg08[!(dfg08$Dato_inds==0),]
# for der hvor vandvolumne er ikke er kendt, indsæt da istedet
# 'ukendt'
dfg08$volfilt_mL[is.na(dfg08$volfilt_mL)] <- "ukendt"
dfg08$volfilt_mL[dfg08$volfilt_mL==0] <- "ukendt"
#
# change the colunm names
colnames(dfg08)[grepl("pos_lon",colnames(dfg08))] <- "længdegrad" 
colnames(dfg08)[grepl("pos_lat",colnames(dfg08))] <- "breddegrad" 
colnames(dfg08)[grepl("Dato_inds",colnames(dfg08))] <- "dato indsamlet"
colnames(dfg08)[grepl("volfilt",colnames(dfg08))] <- "volumen vand filtreret (mL)"
# write a copy of the table
write.table(dfg08,paste0(wd00_wd06,"/table09_v01_all_locations_sampled.csv"),
           row.names = F, col.names = T, sep = ";")
# subset data frame to only comprise rows where there was a detection of eDNA
# from NonIndeginous Species
tblOF <- dfg05[(dfg05$orgFnd==T),]
# subset data frame to only comprise rows where there was NO detection of eDNA
tblNF <- dfg05[(dfg05$orgFnd==F),]

# only keep the columns needed for presenting the sampled locations
dfg09 <- tblOF %>%
  dplyr::select(MST.nummer,
                Dato_inds,
                lokalitet_vanda,
                Lat_Species,
                ssn,
                yea)
# indsæt måneder for sæsoner 
dfg09$ssn[(dfg09$ssn=="spri")] <- "jan-jun"
dfg09$ssn[(dfg09$ssn=="fall")] <- "jul-nov"
# reorder the data frame to be sorted by first species
# then by year, and season, and then by MST number
dfg09 <- dfg09 %>% dplyr::arrange(Lat_Species,
                                  yea,
                                  ssn,
                                  MST.nummer)
# change the colunm names
colnames(dfg09)[grepl("Lat_Species",colnames(dfg09))] <- "latinsk artsnavn" 
colnames(dfg09)[grepl("ssn",colnames(dfg09))] <- "sæson" 
colnames(dfg09)[grepl("Dato_inds",colnames(dfg09))] <- "dato indsamlet"
colnames(dfg09)[grepl("yea",colnames(dfg09))] <- "år"
# write a copy of the table
write.table(dfg09,paste0(wd00_wd06,"/table10_v01_all_pos_detections_abLOD.csv"),
          row.names = F, col.names = T, sep = ";")


#_______________________________________________________________________________
#_______________________________________________________________________________
# only retain specific columns
dfg08 <- dfg05 %>% dplyr::select(Lat_Species,MST.n_Dinds,orgFnd2)
# ensure it is a data frame
dfg08 <- as.data.frame(dfg08)
# exclude rows where the 'Lat_Species' is 0 or is NA
dfg08 <- dfg08[(dfg08$Lat_Species!="0"),]
dfg08 <- dfg08[(!is.na(dfg08$Lat_Species)),]
# check which rows are not unique
dfg08.dpl  <- dfg08 %>% 
  dplyr::group_by(Lat_Species,
                  MST.n_Dinds) %>% 
  dplyr::filter(n()>1)
# Only keep unique rows for 'MSTno_di', and 'Lat speciesr'
dfg08 <- dfg08 %>% dplyr::distinct( Lat_Species,
                                    MST.n_Dinds,
                                    .keep_all = TRUE)
# ensure it is a data frame
dfg08 <- as.data.frame(dfg08)
# rearrange the data frame to a wider version
dfg08 <- dfg08 %>% tidyr::pivot_wider(names_from = Lat_Species,
                                      values_from = orgFnd2)
# ensure it is a data frame
dfg08 <- as.data.frame(dfg08)
# get the column names
clNm <- colnames(dfg08)
# order the columns alphabetically
dfg08 <- dfg08[,order(clNm)]
# get the column names again in the new order
clNm <- colnames(dfg08)
# identify column index numbers
ixcMST <- which(grepl("MST",clNm))
# identify column numbers not yet identified
otcl <- seq(1:length(clNm)) [!(seq(1:length(clNm)) %in% c(ixcMST))]
nOrd <- c(ixcMST,otcl)
# reorder the data frame
dfg08 <- dfg08[,nOrd]
# make it a data frame
dfg08 <- as.data.frame(dfg08)
# # make all columns numeric
dfg08[] <- lapply(dfg08, as.character)
# # Replacing NULL values in a data.frame
# library(dplyr)
dfg08 <- dfg08 %>% replace(.=="NULL", "N")
# replace NAs with only N
dfg08[is.na(dfg08)] <- "N"
# write out the table
write.table(dfg08,paste0(wd00_wd06,"/table11_v01_all_pos_detections_abLOD.csv"),
            row.names = F, col.names = T, sep = ";")
#



library(dplyr)
# summarize by counting per year_season and by location
dfg05.01 <- dfg05 %>% dplyr::group_by(yer_ssn,MST.n_Dinds) %>% 
  dplyr::summarize(
    count = n(),
    mean = mean(orgFnd2 , na.rm = TRUE), 
    sd = sd(orgFnd2 , na.rm = TRUE)) 
# get the count of species screened for per year and season
dfg05.01 <- dfg05.01 %>% dplyr::distinct(yer_ssn,count)
# exclude the observations of locations that did not find eDNA from the organism 
dfg05.02 <- dfg05[(dfg05$orgFnd2>0),]
#dfg05.02 <- dfg05
dfg05.03<- dfg05.02[!is.na(dfg05.02$Lat_Species),]
# summarize by counting per year_season nad by Latin name species
# dfg05.03 <- dfg05.03 %>% dplyr::group_by(Lat_Species,yer_ssn) %>% 
#   summarise(Freq_OrgF = sum(orgFnd2))
# dfg05.03<- dfg05.03[!is.na(dfg05.03$Freq_OrgF),]


dfg05.03 <- dfg05.03 %>%
  dplyr::select(Lat_Species,yer_ssn,orgFnd2) %>%
  dplyr::group_by(Lat_Species,yer_ssn) %>%
  dplyr::summarise(tot_sum = sum(orgFnd2)) 
#
dfg05.03$yer_ssn2 <- dfg05.03$yer_ssn
dfg05.03$yer_ssn2 <- gsub("1st halvår","jan - jun",dfg05.03$yer_ssn2)
dfg05.03$yer_ssn2 <- gsub("2nd halvår","jul - nov",dfg05.03$yer_ssn2)
# make a range of colours for the geom_points in the ggplots
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                "#ffa4ff", "#ff209f", "#8c0020", "#06167F",
                "#B2B23E", "#A56608", "khaki2", "#97FF59",
                "#065d74")

library(ggplot2)
library(tidyverse)
# https://stackoverflow.com/questions/59554096/ggplot2-italics-in-the-legend
toexpr <- function(x, plain = NULL) {
  getfun <- function(x) {
    ifelse(x == plain, "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

# make stacked bar plot
p05 <- ggplot(dfg05.03, aes(fill=Lat_Species, 
               y=tot_sum, 
               x=yer_ssn2)) + 
  # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  theme_void() +
  # make the x-axis labels rotate 90 degrees
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  # alter the labels along the x- and y- axis
  labs(y = "positive miljø-DNA detektioner", x = "årstal og sæson") + 
  guides(fill=guide_legend(title="Latinsk artsnavn")) +
  # use a manual  fill color scale, and make use of the function above for
  # making the species names in italics
  scale_fill_manual(values=cbbPalette,
                     labels = toexpr(unique(dfg05.03$Lat_Species), plain = 'Wt')) +
  # geom_text(data=dfg05.01,aes(label = count), vjust = -0.2) +
  geom_bar(position='stack', stat='identity')


bSaveFigures<-T
if(bSaveFigures==T){
  ggsave(plot = p05, 
         filename = paste0(wd00_wd06,"/Fig06_v01_stckbar_plot_NIS_detected_on_each_site_2021_to_2023.png"),
         #width=210*0.64,height=297,
         #width=210*0.8,height=297,
         #width=210,height=297,
         width=210,height=297*0.7,
         #width=297,height=210,
         #width=1.4*297,height=210,
         units="mm",dpi=300)
}

#_______________________________________________________________________________
# section 03 - start - make individual plots w detections on map per species
#_______________________________________________________________________________

spcsNMs <- unique(dfg05$Lat_Species)
n_spcsNMs<- length(spcsNMs)
sq_spcsNMs <- seq(1,n_spcsNMs,1)
dfg05.04 <- dfg05
dfg07.04 <- dfg07.01
# find unique year and seasons
uyer_ssn<- unique(dfg05$yer_ssn)
library(ggplot2)
# iterate over species
for (ns in sq_spcsNMs)
{
  #
  sbs_spcNm <- spcsNMs[ns]
  print(paste0("making plot ",ns," for ", sbs_spcNm))
  dfg05.04 <- dfg05[(dfg05$Lat_Species==sbs_spcNm),]
  sbs_spcNm_wu <- gsub(" ","_",sbs_spcNm)
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  ins <-stringr::str_pad(ns, 2, pad = "0")
  # iterate over year and season
  for (ys in uyer_ssn)
  {
    # check if there is no entry row for the year and season
    # if this is missing, then add it 
    if (!any(grepl(ys,dfg05.04$yer_ssn)))
    { 
      ncld05.04 <- ncol(dfg05.04)
      #dfg05.04 <- rbind(dfg05.04,c(rep(NA,ncld05.04)))
      dfg05.04 <- rbind(dfg05.04,NA)
      #View(dfg05.04)
      rwnNA <- which(is.na(dfg05.04$yer_ssn))
      dfg05.04$yer_ssn[rwnNA] <- ys
    }
    # end iteration over year and season
  }
  
  dfg05.04$yer_ssn2 <- dfg05.04$yer_ssn
  dfg05.04$yer_ssn2 <- gsub("1st halvår","jan - jun",dfg05.04$yer_ssn2)
  dfg05.04$yer_ssn2 <- gsub("2nd halvår","jun - nov",dfg05.04$yer_ssn2)
  
  dfg05.04$or
  # see this example for getting the plot title in italics
  #https://stackoverflow.com/questions/32555531/how-to-italicize-part-one-or-two-words-of-an-axis-title
  Titl_W_italSpcNm <- bquote(' '~italic(.(sbs_spcNm)))
  #make plot
  p05 <- ggplot(data = world) +
    geom_sf(color = "black", fill = "azure3", lwd=0.1, stroke=0.1) +
    # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
    theme_void() +
    # Plot points for TRUE find of NIS (i.e. above LOD) 
    # using 'geom_pointdensity' to get a heat color for overlaid points
    
    # geom_pointdensity(data=dfg07, 
    #                   mapping= aes(
    #                     x =declon, 
    #                     y = declat) ,
    #                   show.legend = T,
    #                   size=5) +
    # 
    
    # set the color scale for the  'geom_pointdensity' for 
    # the heat color for overlaid points
    #scale_color_viridis_c() +
    #View(dfg05.04)
    # plot all points with crosses, to show were sampling was performed
    # NOTE that this part that plots the crosses as placed first, so that 
    # the next layers with numbered points are placed on top, and will cover
    # the crosses
    geom_point(data=dfg05.04, 
               aes( x =declon, 
                    y = declat,
                    size= orgFnd,
                    shape=orgFnd,
                    fill= orgFnd)) +
    scale_shape_manual(values = c(3,21)) +
    scale_fill_manual(values = c("black",alpha("red",0.7))) +
    scale_size_manual(values = c(2.6,3.0)) +
    #Arrange in facets
    ggplot2::facet_wrap( ~ yer_ssn2,
                         drop=FALSE,
                         dir="h",
                         ncol = 2) +
    
    # ggplot2::facet_wrap( ~ yer_ssn,
    #                      drop=FALSE,
    #                      dir="h",
    #                      ncol = 2,
    #                      labeller = label_bquote(cols =
    #                                                italic(.(as.character(yer_ssn))))  ) +
    # 
    # # see : https://r-charts.com/ggplot2/facets/
    theme(strip.text = element_text(#face = "bold",
      color = "black",
      hjust = 0,
      size = 10),
      strip.background = element_rect(fill = c("white"),
                                      #linetype = "solid",
                                      color = "white",
                                      linewidth = 1)) +
    # place legend on bottom
    #theme(legend.position = "bottom") +
    # or remove the legend
    # theme(legend.position = "none") +
    #define limits of the plot 
    ggplot2::coord_sf(xlim = c(6.6, 17.2),
                      ylim = c(54.2, 58.4), 
                      expand = FALSE) +
    # or remove them all completely
    #http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()) 
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #change axis labels
  #p05t <- p05 + xlab("longitude") + ylab("latitude")
  #p05t <- p05 + xlab("længdegrad") + ylab("breddegrad")
  p05t <- p05 + xlab(" ") + ylab(" ")
  #change the header for the legend on the side, 
  #this must be done for both 'fill', 'color' and 'shape', to avoid 
  #getting separate legends
  #p05t <- p05t + labs(color='miljø-DNA detekteret')
  p05t <- p05t + labs(fill='miljø-DNA detekteret')
  p05t <- p05t + labs(title=Titl_W_italSpcNm) # add a title to the plot that is the species name
  p05t <- p05t + labs(shape='miljø-DNA detekteret')
  p05t <- p05t + labs(size='miljø-DNA detekteret')
  # https://github.com/tidyverse/ggplot2/issues/3492
  #adjust tick marks on axis
  p05t <- p05t + scale_y_continuous(breaks=seq(54.0, 58.4,1))
  p05t <- p05t + scale_x_continuous(breaks=seq(6.0, 17.2,2.0))
  # see the plot
  #p05t
  bSaveFigures<-T
  if(bSaveFigures==T){
    ggsave(plot = p05t, 
           filename = paste0(wd00_wd06,"/Fig07_v",ins,"_map_of_",sbs_spcNm_wu,"_detected_2017_to_2023.png"),
           width=210*0.60,height=297,
           #width=210*0.8,height=297,
           #width=297,height=210,
           #width=297,height=210,
           #width=1.4*297,height=210,
           units="mm",dpi=300)
  }
  # end iteration over latin species names
}
#_______________________________________________________________________________
# section 03 - end - make individual plots w detections on map per species
#_______________________________________________________________________________

#
#_______________________________________________________________________________
# section 04 - start - make maps with sampled locations
#_______________________________________________________________________________
# subset data frame to only comprise
# 
dfg08 <- dfg05 %>% dplyr::select(
                        declat,
                        declon,
                        ssn,
                        yer_ssn,
                        mnt,
                        yer,
                        yea) %>% 
              dplyr::distinct()
dfg08$ssn2 <- dfg08$ssn
dfg08$ssn2[grepl("spr",dfg08$ssn2)] <- "A. jan - jun"
dfg08$ssn2[grepl("fall",dfg08$ssn2)] <- "B. jul - nov"
# exclude if not sampled in 20- something
dfg08 <- dfg08[grepl("20",dfg08$yea),]
# reorder by year
dfg08 <- dfg08 %>% dplyr::arrange(yer) 
#make plot
p05 <- ggplot(data = world) +
  geom_sf(color = "black", fill = "azure3", lwd=0.1, stroke=0.1) +
  # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  theme_void() +
  # plot all points with colors to indicate sample year
  geom_point(data=dfg08, 
             aes( x =declon, 
                  y = declat,
                  fill= yea),
             shape=21,
             size=4,
             position = position_jitter(w = 0.04, h = 0.04)) +
  # color the fills by sample year
  scale_fill_viridis_d(option = 'viridis', direction=1) +   
  #Arrange in facets
  ggplot2::facet_wrap( ~ ssn2 ,
                       drop=FALSE,
                       dir="h",
                       ncol = 2) +
  # # see : https://r-charts.com/ggplot2/facets/
  theme(strip.text = element_text(#face = "bold",
    color = "black",
    hjust = 0,
    size = 12),
    strip.background = element_rect(fill = c("white"),
                                    #linetype = "solid",
                                    color = "white",
                                    linewidth = 1)) +
  # move the legend
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  #change the header for the legend on the side, 
  #this must be done for both 'fill', 'color' and 'shape', to avoid 
  #getting separate legends
  labs(color='') +
  labs(fill='indsamlingsår', nrow=1) +
  #labs(title=Titl_W_italSpcNm) + # add a title to the plot that is the species name
  #labs(shape='') +
  #labs(size='') +
  #define limits of the plot 
  ggplot2::coord_sf(xlim = c(6.6, 17.2),
                    ylim = c(54.2, 58.4), 
                    expand = FALSE) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # or remove them all completely
  #http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()) 
#change axis labels
#p05t <- p05 + xlab("longitude") + ylab("latitude")
#p05t <- p05 + xlab("længdegrad") + ylab("breddegrad")
p05t <- p05 + xlab(" ") + ylab(" ")
# https://github.com/tidyverse/ggplot2/issues/3492
#adjust tick marks on axis
p05t <- p05t + scale_y_continuous(breaks=seq(54.0, 58.4,1))
p05t <- p05t + scale_x_continuous(breaks=seq(6.0, 17.2,2.0))
# see the plot
p05t
bSaveFigures<-T
if(bSaveFigures==T){
  ggsave(plot = p05t, 
         filename = paste0(wd00_wd06,"/Fig08_v01_sampled_locations_2017_to_2023.png"),
         width=210,height=297*0.40,
         #width=210*0.8,height=297,
         #width=297,height=210,
         #width=297,height=210,
         #width=1.4*297,height=210,
         units="mm",dpi=300)
}
# save the data frame
write.csv(dfg05,file=paste0(wd00_wd06,"/table_12_v01_results_for_plotting_on_map.csv"))

#_______________________________________________________________________________
# section 04 - end - make maps with sampled locations
#_______________________________________________________________________________


#_______________________________________________________________________________
# section 05 - start - make linear models
#_______________________________________________________________________________

# see: https://stackoverflow.com/questions/25752909/multiple-ggplot-linear-regression-lines
require(ggplot2)
require(reshape2)
# identify the number of years and seasons
ord.smplye_ssn <- unique(dfg05.03$yer_ssn)[order(unique(dfg05.03$yer_ssn))]
nyessn <- length(ord.smplye_ssn)
sqyessn <- seq(1,nyessn,1)
# combine year_and_season and the number for year_and_season
df_yessn <- as.data.frame(cbind(ord.smplye_ssn,sqyessn))
# change the column names
colnames(df_yessn) <- c("yer_ssn","no.fyssn")
# join the data frames by yer and season
dfg05.04 <- dfg05.03 %>% dplyr::left_join( 
  df_yessn,
  by="yer_ssn")
# ensure the variable is numeris
dfg05.04$no.fyssn <- as.numeric(dfg05.04$no.fyssn)
# make the plot
p05 <- ggplot(dfg05.04) +
  geom_point(aes(no.fyssn,tot_sum, 
                  colour=Lat_Species),) + 
  geom_smooth(aes(no.fyssn,tot_sum, 
                  colour=Lat_Species), method=lm, 
              se=T, 
              fill="grey67",
              level=0.90) +
  # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  theme_void() +
  #Arrange in facets
  ggplot2::facet_wrap( ~ Lat_Species,
                       drop=FALSE,
                       dir="h",
                       #ncol = 2,
                       labeller = label_bquote(cols =
                                                 italic(.(as.character(Lat_Species))))  ) +
  
  # # see : https://r-charts.com/ggplot2/facets/
  theme(strip.text = element_text(#face = "bold",
    color = "black",
    hjust = 0,
    size = 8),
    strip.background = element_rect(fill = c("white"),
                                    #linetype = "solid",
                                    color = "white",
                                    linewidth = 1)) +
  labs(x = "år og sæson", y = "antal miljø DNA detektioner") +
  #guides(fill=guide_legend(title="Latinsk artsnavn")) +
  scale_color_manual(values= cbbPalette) +
  #scale_fill_manual(values= cbbPalette) +
  # or remove the legend
  theme(legend.position = "none") +
  
  scale_x_continuous(breaks=seq(1,nyessn,1),labels=ord.smplye_ssn) + 
  # make the x-axis labels rotate 90 degrees
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))
# or remove them all completely
# #http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
# theme(
#   axis.text.x = element_blank(),
#   axis.text.y = element_blank(),
#   axis.ticks = element_blank()) +

#p05

bSaveFigures<-T
if(bSaveFigures==T){
  ggsave(plot = p05, 
         filename = paste0(wd00_wd06,"/Fig09_v01_linear_models_over_years.png"),
         width=210,height=297*0.6,
         #width=210*0.8,height=297,
         #width=297,height=210,
         #width=297,height=210,
         #width=1.4*297,height=210,
         units="mm",dpi=300)
}

#_______________________________________________________________________________
# section 05 - end - make linear models
#_______________________________________________________________________________


#_______________________________________________________________________________
# section 06 - start - make linear models split by seasons
#_______________________________________________________________________________

# see: https://stackoverflow.com/questions/25752909/multiple-ggplot-linear-regression-lines
require(ggplot2)
require(reshape2)
# get the season
dfg05.03$ssn <- gsub(".*, (.*)","\\1",dfg05.03$yer_ssn)
dfg05.03$ssn2 <- as.numeric(substr(dfg05.03$ssn,1,1))
dfg05.03$yer <- gsub("(.*), (.*)","\\1",dfg05.03$yer_ssn)
dfg05.03$yer <- as.numeric(dfg05.03$yer)
# identify the number of years and seasons
ord.smplye_ssn <- unique(dfg05.03$yer_ssn)[order(unique(dfg05.03$yer_ssn))]
nyessn <- length(ord.smplye_ssn)
sqyessn <- seq(1,nyessn,1)
# combine year_and_season and the number for year_and_season
df_yessn <- as.data.frame(cbind(ord.smplye_ssn,sqyessn))
# change the column names
colnames(df_yessn) <- c("yer_ssn","no.fyssn")
# join the data frames by yer and season
dfg05.04 <- dfg05.03 %>% dplyr::left_join( 
  df_yessn,
  by="yer_ssn")
# ensure the variable is numeris
dfg05.04$no.fyssn <- as.numeric(dfg05.04$no.fyssn)
# make the plot
p05 <- ggplot(dfg05.04) +
  geom_point(aes(no.fyssn,tot_sum, 
                  colour=Lat_Species),) + 
  geom_smooth(aes(no.fyssn,tot_sum, 
                  colour=Lat_Species), method=lm, 
              se=T, 
              fill="grey67",
              level=0.90) +
  # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  theme_void() +
  #Arrange in facets
  ggplot2::facet_wrap( ~ Lat_Species ~ ssn2,
                       drop=FALSE,
                       dir="h",
                       ncol = 2,
                       labeller = label_bquote(cols =
                                                 italic(.(as.character(Lat_Species))))  ) +
  
  # # see : https://r-charts.com/ggplot2/facets/
  theme(strip.text = element_text(#face = "bold",
    color = "black",
    hjust = 0,
    size = 8),
    strip.background = element_rect(fill = c("white"),
                                    #linetype = "solid",
                                    color = "white",
                                    linewidth = 1)) +
  labs(x = "år og sæson", y = "antal miljø DNA detektioner") +
  #guides(fill=guide_legend(title="Latinsk artsnavn")) +
  scale_color_manual(values= cbbPalette) +
  #scale_fill_manual(values= cbbPalette) +
  # or remove the legend
  theme(legend.position = "none") +
  
  scale_x_continuous(breaks=seq(1,nyessn,1),labels=ord.smplye_ssn) + 
  # make the x-axis labels rotate 90 degrees
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

#p05

bSaveFigures<-T
if(bSaveFigures==T){
  ggsave(plot = p05, 
         filename = paste0(wd00_wd06,"/Fig09_v02_linear_models_over_years.png"),
         width=210,height=297*2.2,
         #width=210*0.8,height=297,
         #width=297,height=210,
         #width=297,height=210,
         #width=1.4*297,height=210,
         units="mm",dpi=300)
}

#_______________________________________________________________________________
# section 06 - end - make linear models split by seasons
#_______________________________________________________________________________



#_______________________________________________________________________________
# section 07 - start - make linear models split by seasons
#_______________________________________________________________________________

# see: https://stackoverflow.com/questions/25752909/multiple-ggplot-linear-regression-lines
require(ggplot2)
require(reshape2)
uLsp <- unique(dfg05.03$Lat_Species)
nSPc <-length(uLsp)

fTs <- ceiling(nSPc/4)
sfTs <- seq(1,fTs,1)

# get the season
dfg05.03$ssn <- gsub(".*, (.*)","\\1",dfg05.03$yer_ssn)
dfg05.03$ssn2 <- as.numeric(substr(dfg05.03$ssn,1,1))
dfg05.03$yer <- gsub("(.*), (.*)","\\1",dfg05.03$yer_ssn)
dfg05.03$yer <- as.numeric(dfg05.03$yer)
dfg05.03$ssn3 <- dfg05.03$ssn2
dfg05.03$ssn3[(grepl(1,dfg05.03$ssn2))] <- "jan - jun"
dfg05.03$ssn3[(grepl(2,dfg05.03$ssn2))] <- "jul - nov"
# identify the number of years and seasons
ord.smplye_ssn <- unique(dfg05.03$yer_ssn)[order(unique(dfg05.03$yer_ssn))]
nyessn <- length(ord.smplye_ssn)
sqyessn <- seq(1,nyessn,1)
# combine year_and_season and the number for year_and_season
df_yessn <- as.data.frame(cbind(ord.smplye_ssn,sqyessn))
# change the column names
colnames(df_yessn) <- c("yer_ssn","no.fyssn")
# join the data frames by yer and season
dfg05.04 <- dfg05.03 %>% dplyr::left_join( 
  df_yessn,
  by="yer_ssn")
# ensure the variable is numeris
dfg05.04$no.fyssn <- as.numeric(dfg05.04$no.fyssn)

# itersate over portions of 4's -
for (e in sfTs)
{
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  Pe <- stringr::str_pad(e, 2, pad = "0")
  # start the iteration from 1 below the first
  e <- e-1
  # then get every 4th element from the vector
  sta.e <- e*4+1
  stp.e <- e*4+4
  print(paste0(sta.e," to ",stp.e))
  lssSpc <- uLsp[sta.e:stp.e]
  #}
  sptGf<- paste(lssSpc,collapse = "|")
  # subset to only comprise the species listed
  dfg05.05 <- dfg05.04[grepl(sptGf,dfg05.04$Lat_Species),]
  
  # make the plot
  p05 <- ggplot(dfg05.05) +
    geom_point(aes(no.fyssn,tot_sum, 
                   colour=Lat_Species),) +
    # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
    theme_void() +
    geom_smooth(aes(no.fyssn,tot_sum, 
                    colour=Lat_Species), method=lm, 
                se=T, 
                fill="grey67",
                level=0.90) +
    #Arrange in facets
    ggplot2::facet_wrap( ~ Lat_Species ~ ssn3,
                         drop=FALSE,
                         dir="h",
                         ncol = 2,
                         labeller = label_bquote(cols =
                                                   italic(.(as.character(Lat_Species))) ~ "," ~ .(as.character(ssn3))
                         ) ) +
    
    # # see : https://r-charts.com/ggplot2/facets/
    theme(strip.text = element_text(#face = "bold",
      color = "black",
      hjust = 0,
      size = 8),
      strip.background = element_rect(fill = c("white"),
                                      #linetype = "solid",
                                      color = "white",
                                      linewidth = 1)) +
    labs(x = "år og sæson", y = "antal miljø DNA detektioner") +
    #guides(fill=guide_legend(title="Latinsk artsnavn")) +
    scale_color_manual(values= cbbPalette[sta.e:stp.e]) +
    #scale_fill_manual(values= cbbPalette) +
    # or remove the legend
    theme(legend.position = "none") +
    
    scale_x_continuous(breaks=seq(1,nyessn,1),labels=ord.smplye_ssn) + 
    # make the x-axis labels rotate 90 degrees
    theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))
  
  #p05
  
  bSaveFigures<-T
  if(bSaveFigures==T){
    ggsave(plot = p05, 
           filename = paste0(wd00_wd06,"/Fig10_v",Pe,"_linear_models_over_years.png"),
           width=210,height=297*0.8,
           #width=210*0.8,height=297,
           #width=297,height=210,
           #width=297,height=210,
           #width=1.4*297,height=210,
           units="mm",dpi=300)
  }
  
}
#_______________________________________________________________________________
# section 07 - end - make linear models split by seasons
#_______________________________________________________________________________
