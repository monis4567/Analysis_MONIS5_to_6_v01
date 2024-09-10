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
# #see this
#website
#on how to only install required packages
# #https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
# if (!require("pacman")) install.packages("pacman")
# pacman::p_load(
#   scales, 
#   fields, 
#   gplots,
#   plyr)#,
# #ReporteRs)
library(plyr)
library(scales)
library(gplots)
library(fields)
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
# # devtools::install_github("davidgohel/ReporteRs")
# # devtools::install_github("davidgohel/officer")
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
# # install package if required
# if(!require(envDocument)){
#   install.packages("envDocument")
#   library(envDocument)
# }
library(envDocument)


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
# }
library(ggforce)
#get 'rnaturalearthhires' installed
# if(!require(rnaturalearthhires)){
#   #install.packages("rnaturalearthhires")
#   install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
# }
library(rnaturalearthhires)
# # 
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

#https://www.eleanor-jackson.com/post/searching-for-spring/
options(stringsAsFactors = F)
#get spocc package
# if(!require(spocc)){
#   install.packages("spocc")
# 
# }  
library(spocc)
# #get rinat package
# if(!require(rinat)){
#   remotes::install_github("ropensci/rinat")
#   install.packages("rinat")
# }  
library(rinat)
library("tidyverse")
library("httr")
library("jsonlite")
library("dplyr") 
#load libraries
library(readr)
library(rgbif)
library(dismo)
library(tidyr)

# define working directory

#wd00 <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6/Analysis_MONIS5_to_6"
#wd00 <- "/Users/steenknudsen/Documents/Documents/NIVA_Ansaettelse_2020/NOVANA_proever_2018_2019/"
wd00 <- getwd()
# define input directory
wd06 <- "output06_presence_absence_evaluation_for_MST2017_2022_samples"
# define out directory
wd08 <- "output08_stckbarplot_w_GBIF_and_iNat_for_MST2017_2023_records"
#paste dirs together
wd00_wd08 <- paste0(wd00,"/",wd08)
#Delete any previous versions of the output directory
unlink(wd00_wd08, recursive=TRUE)
#Create a directory to put resulting output files in
dir.create(wd00_wd08)
# set working directory
#setwd(wd00)
#getwd()
#paste dirs together
wd00_wd06 <- paste0(wd00,"/",wd06)


#
df_A08 <- df_A07
#
df_A08$orgFnd2 <- 1
df_A08$orgFnd2[grepl("ingen",df_A08$reccat)] <- 0

df_A08 <- df_A08 %>%
  dplyr::select(Lat_Species,yer_ssn.per,orgFnd2, source) %>%
  dplyr::group_by(Lat_Species,yer_ssn.per,source) %>%
  dplyr::summarise(tot_sum = sum(orgFnd2)) 

# make a range of colours for the geom_points in the ggplots
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                "#ffa4ff", "#ff209f", "#8c0020", "#06167F",
                "#B2B23E", "#A56608","grey21", "khaki2", "#97FF59",
                "#065d74","grey56","orange","seagreen",
                "white","red","orchid3","blue","yellow",
                "limegreen","cyan","magenta","tomato3")

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
p05 <- ggplot(df_A08, aes(fill=Lat_Species, 
                            y=tot_sum, 
                            x=yer_ssn.per)) + 
  # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  theme_void() +
  # make the x-axis labels rotate 90 degrees
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  # alter the labels along the x- and y- axis
  labs(y = "positive miljø-DNA detektioner", x = "årstal og sæson") + 
  guides(fill=guide_legend(title="Latinsk artsnavn")) +
  #Arrange in facets
  ggplot2::facet_wrap( ~ source,
                       drop=FALSE,
                       dir="h",
                       ncol = 1,
                       labeller = label_bquote(cols =
                                                 .(as.character(source))
                       ) ) +
  # use a manual  fill color scale, and make use of the function above for
  # making the species names in italics
  scale_fill_manual(values=cbbPalette,
                    labels = toexpr(unique(df_A08$Lat_Species), plain = 'Wt')) +
  # geom_text(data=dfg05.01,aes(label = count), vjust = -0.2) +
  geom_bar(position='stack', stat='identity')

#p05

bSaveFigures<-T
if(bSaveFigures==T){
  ggsave(plot = p05, 
         filename = paste0(wd00_wd08,"/Fig13_v01_stckbar_plot_NISrecords_from_arterdk_iNat_MONIS6_2021_to_2023.png"),
         width=210,height=297,
         units="mm",dpi=300)
}
