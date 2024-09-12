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

#install.packages("tableHTML")
# #https://cran.r-project.org/web/packages/tableHTML/vignettes/tableHTML.html
# if(!require(tableHTML)){
#   install.packages("tableHTML")
#   library(tableHTML)
# }
library(tableHTML)
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

# #install.packages("worrms")
# if(!require(worrms)){
#   install.packages("worrms")
#   library(worrms)
# }
# #install.packages("taxize")
# if(!require(taxize)){
#   install.packages("taxize")
#   library(taxize)
# }
library(worrms)
library(taxize)

# define working directory

#wd00 <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6/Analysis_MONIS5_to_6"
#wd00 <- "/Users/steenknudsen/Documents/Documents/NIVA_Ansaettelse_2020/NOVANA_proever_2018_2019/"
wd00 <- getwd()
# define input directory
wd06 <- "output06_presence_absence_evaluation_for_MST2017_2022_samples"
wd07 <- "output07_map_eDNA_detections_with_GBIF_and_iNat_for_MST2017_2023_samples"
# define out directory
wd08 <- "output08_stckbarplot_w_GBIF_and_iNat_for_MST2017_2023_records"
#paste dirs together
wd00_wd07 <- paste0(wd00,"/",wd07)
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

# define file name for input file
inpflNm <- paste0(wd00_wd07,"/Table06_iNat_arterdk_and_MONIS6_records_2017_2023.csv")
# read in data frame
df_A07 <- read.table(inpflNm,sep = ";")
#
df_A08 <- df_A07

# make a column that evaluates the season sampled
df_A07$ssn.no <- NA
# evaluate on the season column to add a number version of the season
df_A07$ssn.no[(df_A07$mnt<=6)] <- "1st"
df_A07$ssn.no[(df_A07$mnt>6)]  <- "2nd"

#
df_A08$orgFnd2 <- 1
df_A08$orgFnd2[grepl("ingen",df_A08$reccat)] <- 0
#
df_A08 <- df_A08 %>%
  dplyr::select(Lat_Species,yer_ssn.per,ssn.per,yer,orgFnd2, source) %>%
  dplyr::group_by(Lat_Species,yer_ssn.per,yer,ssn.per,source) %>%
  dplyr::summarise(tot_sum = sum(orgFnd2)) 
# make a range of colours for the geom_points in the ggplots
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                "#ffa4ff", "#ff209f", "#8c0020", "#06167F",
                "#B2B23E", "#A56608","grey21", "khaki3", "#97FF59",
                "#065d74","grey56","orange","seagreen",
                "bisque","red","orchid3","blue","yellow",
                "limegreen","cyan","magenta","tomato3")
# check this website for help on how to retrieve taxonomical information
# https://www.r-bloggers.com/2017/01/extracting-and-enriching-ocean-biogeographic-information-system-obis-data-with-r/
# Check taxonomy of a single taxon
# Get the WoRMS ID for a single species – here, Atlantic cod, Gadus morhua:
  m_sp_aphia <- worrms::wm_name2id(name = "Gadus morhua")
# Then get the full WoRMS record, as a list by Aphia ID:
  m_sp_taxo <- worrms::wm_record(id = m_sp_aphia)
# Or as a tibble by name - here specifying exact match only (fuzzy = FALSE)
# and restricting to marine species (marine_only = TRUE)
  m_sp_taxo <- worrms::wm_records_names(name = "Gadus morhua", fuzzy = FALSE, marine_only = TRUE)
  # substitute on the incorrect species names
  df_A08$Lat_Species <- gsub("Oncorhyncus gorbuscha","Oncorhynchus gorbuscha",
       df_A08$Lat_Species)
  df_A08$Lat_Species <- gsub("Prorocentrum minimum","Prorocentrum cordatum",
                             df_A08$Lat_Species)
  # get all the species in the data frame
  LTSpc <- unique(df_A08$Lat_Species)
  # Get taxonomy for multiple species
# Start with a data frame of species names:
  m_sp <- data_frame(sciname = c(LTSpc))
# Then get the WoRMS records for each:
  m_sp_taxo <- worrms::wm_records_names(name = m_sp$sciname, fuzzy = FALSE, marine_only = TRUE)
# For 'n' species this returns a list of 'n' tables. Convert these into a 
  # single table with 'n' rows:
  m_sp_taxo <- bind_rows(m_sp_taxo)
# exclude taxonomically unaccepted names
  m_sp_taxo <- m_sp_taxo[!grepl("unaccept",m_sp_taxo$status),]
# define columns to keep
ctkeep <- c( "AphiaID", "scientificname", "authority", "status", 
             "valid_AphiaID", "valid_name", 
             "valid_authority", "parentNameUsageID", "kingdom", "phylum", 
             "class", "order", "family", "genus")
# only keep specified columns
m1_sp_taxo <- m_sp_taxo[ctkeep]
# copy and rename column
m1_sp_taxo$Lat_Species <- m1_sp_taxo$scientificname
# join the data frames
df_A08 <- df_A08 %>% dplyr::left_join(m1_sp_taxo,
                            by="Lat_Species")
# prepare taxonomical group_catagories for the species
# to be able to assign gradient colors to the stacked bar plots
#https://stackoverflow.com/questions/49818271/stacked-barplot-with-colour-gradients-for-each-bar
library(ggplot2)
# make a column with class and order combined
df_A08$class_order <- paste0(df_A08$class,"_" ,df_A08$order)
# make a column with phylum, class and order combined
df_A08$phylum_class_order <- paste0(df_A08$phylum,"_",df_A08$class,"_" ,df_A08$order)
df_A08$phylum_class_order_latspc <- paste0(df_A08$phylum,"_",df_A08$class,"_" ,df_A08$order,"_",df_A08$Lat_Species)
df_A08$class_order_latspc <- paste0(df_A08$class,"_" ,df_A08$order,"_",df_A08$Lat_Species)
# make a column with kingdom, phylum, and class
df_A08$kingdom_phylum_class <- paste0(df_A08$kingdom,"_",df_A08$phylum,"_",df_A08$class)
# combine categories
df_A08$group <- paste0(df_A08$phylum, "-", 
                       df_A08$class_order,
                        sep = "")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# start - function 'ColourPalleteMulti'
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ColourPalleteMulti <- function(df, group, subgroup){
  # Find how many colour categories to create and the number of colours in each
  categories <- aggregate(as.formula(paste(subgroup,
                                           group, sep="~" )), 
                          df, function(x) length(unique(x)))
  category.start <- (scales::hue_pal(l = 100)(nrow(categories))) # Set the top of the colour pallete
  category.end  <- (scales::hue_pal(l = 30)(nrow(categories))) # set the bottom
  
  # Build Colour pallette
  colours <- unlist(lapply(1:nrow(categories),
                           function(i){
                             colorRampPalette(colors = c(category.start[i], 
                                                         category.end[i]))(categories[i,2])}))
  return(colours)
}
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# end - function 'ColourPalleteMulti'
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Build the colour pallete using the function
col.phycla <-ColourPalleteMulti(df_A08, "phylum", "class_order_latspc")

# combinde columns to a data frame
class_order_latspc <- (unique(df_A08$class_order_latspc))
df_cfphy <- as.data.frame(cbind(class_order_latspc,col.phycla))

library(ggplot2)
library(tidyverse)
# https://stackoverflow.com/questions/59554096/ggplot2-italics-in-the-legend
toexpr <- function(x, plain = NULL) {
  getfun <- function(x) {
    ifelse(x == plain, "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}
# get the years
yer <- unique(df_A08$yer)
nyer <- length(yer)
emptyer <- rep(" ",nyer)
sqfyer <- seq(1,nyer,1)
#yer <- c(rep(yer,2))
yer <- yer[(order(yer))]
yer_and_empt <- c(emptyer,yer)
seqfyeandempt <- c(rep(sqfyer,2))
df_sqfyer <- as.data.frame(cbind(yer_and_empt,seqfyeandempt))
df_sqfyer <- df_sqfyer[order(df_sqfyer$seqfyeandempt),]
yere <- df_sqfyer$yer_and_empt
# substitute in the 'source'
df_A08$source <- gsub("iNat" ,"iNaturalist" ,df_A08$source)
df_A08$source <- gsub("arter\\.dk" ,"www.arter.dk" ,df_A08$source)
# reorder the data frame by multiple columns
# to ensure the colors match the 
df_A08 <- df_A08[with(df_A08, order(df_A08$phylum_class_order_latspc, df_A08$Lat_Species)), ]
# make stacked bar plot
p05 <- ggplot(df_A08, aes(#fill=Lat_Species, 
                          fill = phylum_class_order_latspc,
                            y=tot_sum, 
                            x=yer_ssn.per)) + 
  # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  #theme_void() +
  theme_bw() +
  # change the header on the facet wrap
  # https://stackoverflow.com/questions/41631806/change-facet-label-text-and-background-colour
  theme(strip.background =element_rect(colour = 'white',fill="white"))+
  theme(strip.text = element_text(colour = 'black', face="bold", hjust=0.1) ) +
  # Remove grid, background color, and top and right borders from ggplot2
  # https://stackoverflow.com/questions/10861773/remove-grid-background-color-and-top-and-right-borders-from-ggplot2
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # make the x-axis labels rotate 90 degrees
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  scale_x_discrete(label = yere) +
  # alter the labels along the x- and y- axis
  labs(y = "positive miljø-DNA detektioner", x = "årstal og sæson") + 
  guides(fill=guide_legend(title="Latinsk artsnavn",
                           ncol=1)) +
  
  #Arrange in facets
  ggplot2::facet_wrap( ~ source+ssn.per,
                       drop=FALSE,
                       dir="h",
                       ncol = 2,
                       labeller = label_bquote(cols =
                          .(as.character(paste0(source,", ",ssn.per)))
                       ) ) +
  # change the x axis
  #scale_x_discrete(name="year", breaks=c(as.character(yer))) +
  # use a manual  fill color scale, and make use of the function above for
  # making the species names in italics
  # set the colors for the fill manually
  scale_fill_manual(values=c(col.phycla),
                    labels = toexpr(unique(df_A08$Lat_Species), 
                                    plain = 'Wt')) +
  # scale_fill_manual(values=cbbPalette,
  #                   labels = toexpr(unique(df_A08$Lat_Species), plain = 'Wt')) +
  # geom_text(data=dfg05.01,aes(label = count), vjust = -0.2) +
  geom_bar(position='stack', stat='identity',  width = 1.8)

#p05

bSaveFigures<-T
if(bSaveFigures==T){
  ggsave(plot = p05, 
         filename = paste0(wd00_wd08,"/Fig13_v01_stckbar_plot_NISrecords_from_arterdk_iNat_MONIS6_2021_to_2023.png"),
         width=210,height=297*0.7,
         units="mm",dpi=300)
}

# join the data frames
df_A07 <- df_A07 %>% dplyr::left_join(m1_sp_taxo,
                                      by="Lat_Species")
# make a column with class and order combined
df_A07$class_order <- paste0(df_A07$class,"_" ,df_A07$order)
# make a column with phylum, class and order combined
df_A07$phylum_class_order <- paste0(df_A07$phylum,"_",df_A07$class,"_" ,df_A07$order)
df_A07$phylum_class_order_latspc <- paste0(df_A07$phylum,"_",df_A07$class,"_" ,df_A07$order,"_",df_A07$Lat_Species)
df_A07$class_order_latspc <- paste0(df_A07$class,"_" ,df_A07$order,"_",df_A07$Lat_Species)
# make a column with kingdom, phylum, and class
df_A07$kingdom_phylum_class <- paste0(df_A07$kingdom,"_",df_A07$phylum,"_",df_A07$class)
# combine categories
df_A07$group <- paste0(df_A07$phylum, "-", 
                       df_A07$class_order,
                       sep = "")
# combine data frames using lefy join to get colors added to the data frame
df_A09 <- df_A07 %>% dplyr::left_join(df_cfphy,
                            by="class_order_latspc")

# combine data frames using lefy join to get colors added to the data frame
df_A08 <- df_A08 %>% dplyr::left_join(df_cfphy,
                                      by="class_order_latspc")
# find the minimum and maximum for the years sampled
mn.yer <- min(df_A08$yer)
mx.yer <- max(df_A08$yer)
# exclude rows if zero
# No need to plot points that have zero detections
# and no need to evaluate on the zero detections
df_A08 <- df_A08[!(df_A08$tot_sum==0),]
# identify which rows have detections, that occurs in how many rows
cntsA08 <- plyr::ddply(df_A08, .(df_A08$ssn.per, 
                           df_A08$source,
                           df_A08$Lat_Species), nrow)
names(cntsA08) <- c("ssn.per","source", "Lat_Species", "Freq_pYer_p_src")
# join the data frame with all observations with the frequency count 
# per season per species per source
df_A08 <- dplyr::left_join(df_A08, 
                 cntsA08, by = c("ssn.per","source", "Lat_Species"))
# exclude rows if the 'Freq_pYer_p_src' count is 2 or less, as it is not 
# possible to infer linear regression on less than 3 points, so there is no
# point in including such points
df_A08 <- df_A08[(df_A08$Freq_pYer_p_src>=3),]
# ensure the package required for plotting is loaded
library(ggplot2)
# subset to only comprise one group
phylum.rec <- unique(df_A08$phylum)
# make a sequence for this range
nfPhy<- seq(1,length(phylum.rec),1)
# iterate over numbers for phylum
for (phrc in nfPhy)
{
  print(phrc)
  #}
  # get the phylum name with an underscore
  phylum_wu <- gsub(" ","_",phylum.rec[phrc])
  # subset the datafram
  df_A10 <- df_A08[(df_A08$phylum==phylum_wu),]
  # exclude rows if zero
  # No need to plot points that have zero detections
  df_A10 <- df_A10[!(df_A10$tot_sum==0),]
  # work out an upper level for the plot
  mx_totsum <- max(df_A10$tot_sum)
  # use the plyr function to get the upper 10th level see this question:
  # https://stackoverflow.com/questions/6461209/how-to-round-up-to-the-nearest-10-or-100-or-x
  max_no_of_detct <- plyr::round_any(mx_totsum, 10, f = ceiling)  
  
    # check if the data frame is empty
  if (!empty(df_A10)==T)
    # if it is not empty , then make a plot
  {
    #create plot to visualize fitted linear regression model
    p05 <- ggplot(df_A10,aes(x=yer,
                             y=tot_sum,
                             color=phylum_class_order_latspc,
                             fill=phylum_class_order_latspc),
                  shape=21,
                  size=2)+
      geom_point() +
      # add a linear regression model with standard error that has a
      # level of 0.90
      geom_smooth(method = "lm", se=T, level = 0.90) +
      #theme(legend.position="none") +
      #Arrange in facets
      ggplot2::facet_wrap( ~ source+ssn.per,
                           drop=FALSE,
                           dir="h",
                           ncol = 2,
                           labeller = label_bquote(cols =
                                                     .(as.character(paste0(source,", ",ssn.per)))
                           ) ) +
      # use previous defined scale of colors
      scale_fill_manual(values=c(cbbPalette),
                        labels = toexpr(unique(df_A10$Lat_Species), 
                                        plain = 'Wt') ) +
      scale_color_manual(values=c(cbbPalette),
                         labels = toexpr(unique(df_A10$Lat_Species), 
                                         plain = 'Wt')) +
      # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
      theme_bw() +
      # change the header on the facet wrap
      # https://stackoverflow.com/questions/41631806/change-facet-label-text-and-background-colour
      theme(strip.background =element_rect(colour = 'white',fill="white"))+
      theme(strip.text = element_text(colour = 'black', face="bold", hjust=0.1) ) +
      #theme_void() +
      # Remove grid, background color, and top and right borders from ggplot2
      # https://stackoverflow.com/questions/10861773/remove-grid-background-color-and-top-and-right-borders-from-ggplot2
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) +
      # set the lower and upper limit for the axis
      coord_cartesian(ylim = c(0, max_no_of_detct))  +
      
      # make the x-axis labels rotate 90 degrees
      theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
      # make the x axis have breaks that are represented every 1 increment      
      scale_x_continuous(breaks=seq(mn.yer, mx.yer, 1)) +
      # both 'fill' and 'color' needs to be changed
      guides(fill=guide_legend(title="Latinsk artsnavn",
                               ncol=1),
             color=guide_legend(title="Latinsk artsnavn",
                                ncol=1)) +
      
      # alter the labels along the x- and y- axis
      labs( x = "årstal og sæson") +
      # break the label across more lines, this website recommends
      # more complicated approaches, 
      # https://www.gangofcoders.net/solution/ggplot2-two-line-label-with-expression/
      #Use atop to fake a line break
      ylab(expression(atop("fund eller detektioner", 
                             paste("per sæson"))))
    #p05
    #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
    ins <- stringr::str_pad(phrc, 2, pad = "0")
    # find out how many 'source' that contribute to the plot with facet wrap plots
    # and use this factor to calculate a factor that can be used for 
    # setting the height of the plot
    nrw_hfplt <- length(unique(df_A10$source))
    fct.hght <- 1/(3/nrw_hfplt)
    
    # evaluate whether to store the plot
    bSaveFigures<-T
    if(bSaveFigures==T){
      ggsave(plot = p05, 
             filename = paste0(wd00_wd08,"/Fig14_v",ins,"_linear_regr_",phylum_wu,".png"),
             width=210,height=297*0.6*fct.hght,
             units="mm",dpi=300)
# end check for whether to store plot
    }
# end 'if check' if data frame is empty
  }
# end iteration over phylum
}
#_______________________________________________________________________________


# paste together to have a variable that reflects the species, the source and
# the season sampled
df_A08$LtSp.ssn.src <- paste(df_A08$Lat_Species,
                            df_A08$ssn.per,
                            df_A08$source,
                            sep="_")
# find the unique combinations for species, seasdon and source
LtSp.ssn.src <- unique(df_A08$LtSp.ssn.src)
# and count them, to make a sequence of numbers
nLtSp.ssn.src <- length(LtSp.ssn.src)
sqLtSp.ssn.src <- seq(1,nLtSp.ssn.src,1)
#

lst_lm_spc.ssn.src <- list()
for (i in sqLtSp.ssn.src)
  {
  #i <- 59
  #}
  oLtSp_sn_src <- LtSp.ssn.src[i]
  print(oLtSp_sn_src)
  # split the string
  lngNmspl <- strsplit(as.character(oLtSp_sn_src), "_")
  # and use the splitted string to get each element
  LtsNm <- sapply(lngNmspl, "[[", 1)
  ssnNm <- sapply(lngNmspl, "[[", 2)
  srcNm <- sapply(lngNmspl, "[[", 3)
  
  
  df_A11 <- df_A08[(df_A08$LtSp.ssn.src==oLtSp_sn_src),]
  
  yer.tsm.lm <- lm(tot_sum~ yer, data = df_A11)
  #cor(x, y)  # calculate correlation between x and y
  # calculate correlation between tot_sum and yer 
  cor_yer.tsm.lm <- cor(df_A11$tot_sum, df_A11$yer)  
  summ.lm_yer.tsm  <-summary(yer.tsm.lm)
  # get the R2 and adjusted R2 values 
  adj.r.squared <- summ.lm_yer.tsm$adj.r.squared
  adj.r.squared <- round(adj.r.squared, digits = 2)
  r.squared <- summ.lm_yer.tsm$r.squared
  r.squared <- round(r.squared, digits = 2)
  # get the intercept and the increment and p-values
  intcpt <- summ.lm_yer.tsm$coefficients[1,1]
  intcpt <- round(intcpt, digits = 2)
  incrm <- summ.lm_yer.tsm$coefficients[2,1]
  incrm <- round(incrm, digits = 2)
  pval_lm_yer.tsm <- summ.lm_yer.tsm$coefficients[2,4]
  pval_lm_yer.tsm <- round(pval_lm_yer.tsm, digits = 2)
  
  # combine values in to a data frame
  df_lm_yer.tsm <- as.data.frame(cbind(intcpt,
                                        incrm,
                                        r.squared,
                                        adj.r.squared,
                                        pval_lm_yer.tsm,
                                       oLtSp_sn_src,
                                       LtsNm,
                                       ssnNm,
                                       srcNm
                                       )) 
  # collect the data frame in a list
  lst_lm_spc.ssn.src[[i]] <- df_lm_yer.tsm
  # end iteration over seq numbers for latinspc_sssn_source
  }

# combine the list of data frames into a data frame
df_A11_lm_sum <- dplyr::bind_rows(lst_lm_spc.ssn.src, .id = "column_label")
# exclude if  'r.squared' is NaN
df_A11_lm_sum <- df_A11_lm_sum[!(grepl("NaN",df_A11_lm_sum$r.squared)),]

# re-order the data frame by species, by season sampled, and then by
# source
df_A11_lm_sum <- df_A11_lm_sum %>% dplyr::arrange(LtsNm, ssnNm, srcNm)
# Check if the 
df_A11_lm_sum$sgnfpval <- df_A11_lm_sum$pval_lm_yer.tsm<=0.05
df_A11_lm_sum$sgnf_eval <- "ikke sign"
df_A11_lm_sum$sgnf_eval[(df_A11_lm_sum$sgnfpval==T)] <- "sign"
df_lm_sum <- df_A11_lm_sum[(df_A11_lm_sum$sgnfpval==T),]
# define columns to keep
ctkeep <- c("LtsNm", "ssnNm", "srcNm", 
            "intcpt", "incrm", "r.squared", "adj.r.squared", 
"pval_lm_yer.tsm","sgnf_eval")
# only keep specified columns
df_lm_sum02 <- df_A11_lm_sum[ctkeep]
# change column headers
colnames(df_lm_sum02) <- c( "Latinsk artsnavn",
                            "sæson",
                            "kilde for fund",
                            "skæring",
                            "hældning",
                            "R2",
                            "jR2",
                            "p-værdi",
                            "sigifikans vurdering")

library(htmlTable)
# make a table caption
capt_tbl02 <-        paste0(
  "Tabel 1. Værdier for lineær regressionsmodeller for antallet af fund af ikke-hjemmehørende",
  " arter per år per sæson for MONIS6 projektet og to hjemmesider (arter.dk og iNaturalist.org). ",
  "Skæring og hældning er for den estimerede lineære model. Vurderingen om",
  " den lineære sammenhæng mellem antallet af fund og året er gjort en p-værdi",
  " på 0,05 eller mindre. For de arter hvor der er en en sammenhæng er angiver kolonnen ",
  "'singifikans vurdering' om den lineære model er signifikant eller ej.",
  " Ikke alle arter, sæsoner og kilder er inkluderet i tabellen, da der for nogle arter og nogle",
  " kilder ikke var tre eller flere datapunkter at udføre lineær regression på."
)
# show the table
t.HTML06 <- df_lm_sum02 %>%
  addHtmlTableStyle(align = "r") %>%
  htmlTable(caption = capt_tbl02, rnames = FALSE)
t.HTML06

# Make a filename to store the html table with  values
filNm.for_html <- paste0(
  wd00_wd08,
  "/Table01_v01",
  "_html_table_linear_regress_model",
  ".html"
)
# save the html table
htmltools::save_html(t.HTML06, file = filNm.for_html)

  #
  
  #
  
  