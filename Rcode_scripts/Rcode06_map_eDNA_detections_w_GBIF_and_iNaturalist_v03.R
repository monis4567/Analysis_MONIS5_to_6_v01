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
wd07 <- "output07_map_eDNA_detections_with_GBIF_and_iNat_for_MST2017_2023_samples"
#paste dirs together
wd00_wd07 <- paste0(wd00,"/",wd07)
#Delete any previous versions of the output directory
unlink(wd00_wd07, recursive=TRUE)
#Create a directory to put resulting output files in
dir.create(wd00_wd07)
# set working directory
#setwd(wd00)
#getwd()
#paste dirs together
wd00_wd06 <- paste0(wd00,"/",wd06)

#_________________________
# 01 start - make trycatch function, in case the genus name is not on GBIF 
#_________________________
# https://www.statology.org/r-trycatch/
try.c.iNat <- function(tx, boundslim, endval){
  tryCatch(
    {
      #get iNaturalist records
      g <- rinat::get_inat_obs(
        taxon_name = tx,
        quality = "research",
        geo=T, #only include geo referenced results
        bounds = boundslim,
        maxresults = endval)
      return(g)
    },
    error=function(e) {
      message('An Error Occurred')
      print(e)
    },
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
      return(NA)
    }
  )
}
#_________________________
# 01 end - make trycatch function, in case the genus name is not on GBIF 
#_________________________


# https://github.com/cran/rinat
#_______________________________________________________________________________
# nelng cuts on the eastern boundary
# nelat cut on the northern border
#                             |
#                       nelng |
#             N         nelat___    |
#             |                     | y-axis is lat
#             |                     |
#             |                     |
#   W____________________E          |
#             |                     |
#             |                     |
#             |                     |
# ___  swlat  S                     |
#     |swlng
#     |
#_____________________________
#             x-axis is lon
# swlng cuts on the western boundary
# swlat cut on the southern border
#try defining your own bounding box
set_nelat= 58
set_nelng= 15.4
set_swlat= 54.4
set_swlng= 8
#try defining your own bounding box
set_nelat= 59.5
set_nelng= 17
set_swlat= 53
set_swlng= 7
# needs to be in the format 
#  'min_lat','min_lon','max_lat', 'max_lon'
# whcih equals
#  'min_y','min_x','max_y', 'max_x'
boundslim <- c(set_swlat, set_swlng, set_nelat, set_nelng)
## Search using the boundslimits defined above
df_iNat01 <- get_inat_obs(taxon_name = "Mya arenaria",
                          quality = "research",
                          bounds = boundslim,
                          maxresults = 500)
# make the ggplot
plt_amp01 <- ggplot(data = df_iNat01, aes(x = longitude,
                                          y = latitude,
                                          colour = scientific_name)) +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "grey95",
               color = "gray40",
               linewidth = 0.1) +
  geom_point(size = 0.7, alpha = 0.5) +
  coord_fixed(xlim = range(df_iNat01$longitude, na.rm = TRUE),
              ylim = range(df_iNat01$latitude, na.rm = TRUE)) +
  theme_bw()
# see the plot
plt_amp01



#_______________________________________________________________________________

#make a list of organisms
lst_orgnsm <- c(
  "Mya arenaria",
  "Callinectes sapidus")
#_______________________________________________________________________________

# read in the saved  data frame
dfg05<- read.csv(file=paste0(wd00_wd06,"/table_12_v01_results_for_plotting_on_map.csv"))
# vmake a alist of the species names
latSpcNms  <- unique(dfg05$Lat_Species)
# limit to species names that have a 'sapce' included
latSpcNms <- latSpcNms[(grepl(" ",latSpcNms))]
#_______________________________________________________________________________
# copy the list into a new object
lst_orgnsm <- latSpcNms

lst_spcs <- lst_orgnsm
#_______________________________________________________________________________

# 01 start - make trycatch function, in case the genus name is not on GBIF 
#_______________________________________________________________________________
# https://www.statology.org/r-trycatch/
try.c.iNat <- function(tx, boundslim, endval){
  tryCatch(
    {
      #get iNaturalist records
      g <- rinat::get_inat_obs(
        taxon_name = tx,
        quality = "research",
        geo=T, #only include geo referenced results
        bounds = boundslim,
        maxresults = endval)
      return(g)
    },
    error=function(e) {
      message('An Error Occurred')
      print(e)
    },
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
      return(NA)
    }
  )
}
#_________________________
# 01 end - make trycatch function, in case the genus name is not on GBIF 
#_________________________
kee <- c("scientific_name",
         "datetime",
         "place_guess",
         "latitude",
         "longitude",
         "tag_list",
         "common_name",
         "url",
         "image_url",
         "species_guess",
         "iconic_taxon_name",
         "taxon_id",
         "num_identification_agreements",
         "num_identification_disagreements",
         "observed_on_string",
         "observed_on",
         "time_observed_at",
         "time_zone")

#make an empty list to use for collecting data frame
lst_tx_gobs <- list()
#start a growing number
i <- 1
#iterate over taxon names in list 
for (tx in lst_spcs)
{
  #  print(tx)}
  print(tx)
  # substitute the underscore
  tx <- gsub("_"," ",tx)
  ## Search for the  species using the boundslimits defined above
  g <- try.c.iNat(tx, boundslim, 5000)
  # limit to only specific columns, otherwise it ends up being too much data
  g <- g[kee]
  # check if there are no data, and in that case, add NAs for the 'kee' columns
  if(!is.null(colnames(g)))
  {g <- g} else {
    df_tmp <- as.data.frame(t(as.matrix(kee)))
    df_tmp <- rbind(df_tmp,rep(NA,length(kee)))
    colnames(df_tmp) <- df_tmp[1,]
    df_tmp <- df_tmp[-1,]
    g <- df_tmp
  }
  # add the taxon name that was used for making the search
  g$txNmsrch <- tx
  # make the entire data frame characters
  g[] <- lapply(g, as.character)
  # append the data frame to the list of data frames
  # store it as the i'th element 
  lst_tx_gobs[[i]] <- g
  # increase the count of i by one
  i <- i+1
  # end iteration over amphibian species in the list
}
#bind the rows in each list in to one data frame
df_g03 <- data.table::rbindlist(lst_tx_gobs, fill=T)
df_g03 <- as.data.frame(df_g03)
# if there is no latitude, then omit the row
df_g03 <- df_g03[!is.na(df_g03$lat),]
# copy the column with the scientific name
df_g03$scientific_name2 <- df_g03$scientific_name



#_______________________________________________________________________________

# https://uchicagoconsulting.wordpress.com/tag/r-ggplot2-maps-visualization/
# #install packages needed
# if(!require(maps)){
#   install.packages("maps")
# }
# 
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
#   library(rnaturalearth)
# }
library(rnaturalearth)
# if(!require(rnaturalearthdata)){
#   install.packages("rnaturalearthdata")
# }
library(rnaturalearthdata)
# #install rgeos
# if(!require(rgeos)){
#   install.packages("rgeos")
#   library(rgeos)
# }
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
theme_set(theme_bw())
library("sf")
library(ggforce)


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
# copy the data frame
df_A05 <- df_g03
# copy columns
df_A05$Latitude <- df_A05$latitude
df_A05$Longitude <- df_A05$longitude
# make sure they are numeric
df_A05$Longitude <- as.numeric(df_A05$Longitude)
df_A05$Latitude <- as.numeric(df_A05$Latitude)
# add an evaluation category based on the finding
df_A05$eval <- "iNat"
dfg05$eval <- "no eDNA"
# evaluate based on whether there was any eDNA recorded
dfg05$eval[(dfg05$orgFnd2>0)] <- "eDNA found"
# define columns to keep
keep <- c(
  "scientific_name", 
  "datetime", 
  #"place_guess", 
  "latitude", 
  "longitude", 
  #"tag_list", "common_name", "url", "image_url", "species_guess", 
  #"iconic_taxon_name", "taxon_id", "num_identification_agreements", 
  #"num_identification_disagreements", "observed_on_string", "observed_on", 
  #"time_observed_at", "time_zone", 
  "txNmsrch", 
  "scientific_name2", 
  "eval")
# keep only selected columns
df_A05 <- df_A05[keep]
# rename columns to have common names
df_A05$Dato_inds <- df_A05$datetime
df_A05$declat <- df_A05$latitude
df_A05$declon <- df_A05$longitude
df_A05$Lat_Species <- df_A05$txNmsrch
# get the year from the date
# see: https://stackoverflow.com/questions/36568070/extract-year-from-date
df_A05$yer<- format(as.Date(df_A05$Dato_inds),"%Y")
# get the month  from the date
df_A05$mnt<- format(as.Date(df_A05$Dato_inds),"%m")
# make years numeric 
df_A05$mnt <- as.numeric(df_A05$mnt)
df_A05$yer <- as.numeric(df_A05$yer)
df_A05$ssnno <- NA
# evalutate on the season column to add a number version of the season
df_A05$ssnno[(df_A05$mnt<=6)] <- "1st"
df_A05$ssnno[(df_A05$mnt>6)] <- "2nd"
# paste together year and the numbered part of the year
#df_A05$yer_ssn <- paste0(df_A05$yer,", ",df_A05$ssnno," halvår")
# get the years and seasons,
# 
uyssn <- unique(dfg05$yer_ssn2)
noss<- length(uyssn)
# Repeat rows making each repeated rows following the original rows and assign new variables for each row [duplicate]
# https://stackoverflow.com/questions/69160320/repeat-rows-making-each-repeated-rows-following-the-original-rows-and-assign-new?noredirect=1&lq=1
library(dplyr)
library(tidyr)
# use tidyr to replicate
df_A06  <- df_A05 %>%
  tidyr::uncount(noss) %>%
  dplyr::mutate(yer_ssn2 = rep(c(uyssn), length.out = n()))



# define columns to keep
keepcol <- c( "Dato_inds",
              "declat",
              "declon",
              "yer" , 
              "mnt",
              "yer_ssn2",
              "eval",
              "Lat_Species")
# only keep specific columns
df_A06 <- df_A06[keepcol]
dfg06 <- dfg05[keepcol]
#bind rows together in columns
df_A06 <- rbind(df_A06,dfg06)

# exclude rows if it does not contain a space, and if it contains a zero
df_A06 <- df_A06[grepl(" ",df_A06$Lat_Species),]
df_A06 <- df_A06[!grepl("0",df_A06$Lat_Species),]


# make color range
cl06 <- c("firebrick3","steelblue1","white")
nspo2 <- length(unique(df_A06$eval))
# # make a datra frame to sort monitoring categories
# df_mcat <- as.data.frame(cbind(unique(df_A06$eval), c(2,3,1)))
# # change column names
# colnames(df_mcat) <- c("eval","ordcat")
# # match back to main data frame
# df_A06$ordcat <- df_mcat$ordcat[match(df_A06$eval,df_mcat$eval)]
# ordcat1 <- unique(df_A06$ordcat)
# reorder the data frame
# df_A06 <- df_A06[order(df_A06$family, df_A06$ordcat), ]
# replace in the evaluation category
df_A06$eval <- gsub("eDNA found","eDNA detected",df_A06$eval)
df_A06$eval <- gsub("no_eDNA","no eDNA",df_A06$eval)
df_A06$eval <- gsub("gbif","GBIF",df_A06$eval)
# opy the columns with latitude and longitude
df_A06$Latitude  <- as.numeric(df_A06$declat)
df_A06$Longitude <- as.numeric(df_A06$declon)
# get the number of latin species
LSpc <- unique(df_A06$Lat_Species)
unique(df_A06$yer_ssn2)
# make a sequence for this range
nfLSpc<- seq(1,length(LSpc),1)
# iterate over species
for (i in nfLSpc)
  {print(paste0(i, " making plot for ",LSpc[i]))
  #}
   # i <- 7
    #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
    ins <- stringr::str_pad(i, 2, pad = "0")
    # get the species name with an underscore
    sbs_spcNm_wu <- gsub(" ","_",LSpc[i])
    #}
  # make a function to use for making facet wrap headers
  # https://ggplot2.tidyverse.org/reference/as_labeller.html
  # https://stackoverflow.com/questions/63857833/changing-the-facet-wrap-labels-using-labeller-in-ggplot2
  #appender <- function(nm1) paste0(df_tx01$class[match(nm1,df_tx01$family )],": ", nm1)
  
    df_A06.1 <- df_A06[grepl(LSpc[i],df_A06$Lat_Species),]
  # add a level for jittering points
  jitlvl <- 0.03
  #make plot
  
  # identify the missing seasons for sampled period
  usss <- unique(df_A06.1$yer_ssn2) 
  msssn <- uyssn[!(uyssn %in% usss) ]
  nmsssn <- length(msssn)
  nclA06.1 <- ncol(df_A06.1)
  eRw <- rep(0,nclA06.1)
  clNmA06.1 <- colnames(df_A06.1)
  dfeR <- as.data.frame(t((eRw)))
  colnames(dfeR) <- clNmA06.1
  dfeR <- dfeR[rep(seq_len(nrow(dfeR)), each = nmsssn), ]
  ixcN <- which(grepl("yer_ssn2",colnames(df_A06.1)))
  dfeR[,ixcN] <- msssn
  df_A06.2 <- rbind(df_A06.1,dfeR)
  # make 'iNat' categories for the zeroes for the added years
  df_A06.2$eval[(df_A06.2$eval==0)] <- "iNat"
  
  #df_A06.2 <- df_A06.1
  # subset data frame
  df_A07 <- df_A06.2[(df_A06.2$eval=="iNat"),]
  
  # make a data frame 
  df_A07.ne <- df_A06.1[(df_A06.2$eval=="no eDNA"),]
  df_A07.ed <- df_A06.1[(df_A06.2$eval=="eDNA detected"),]
  
  #make plot
  p05 <- ggplot(data = world) +
    geom_sf(color = "black", fill = "azure3", lwd=0.1) +
    # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
    theme_void() +
    #https://ggplot2.tidyverse.org/reference/position_jitter.html
    #https://stackoverflow.com/questions/15706281/controlling-the-order-of-points-in-ggplot2
    # use 'geom_jitter' instead of 'geom_point' 
    geom_jitter(data = df_A07 ,
                aes(x = Longitude, y = Latitude),
                color=alpha(c("black"),c(0.6)),
                size=1.6,
                fill=alpha(c("yellow"),c(0.6)),
                shape=22,
                width = jitlvl, #0.07, jitter width 
                height = jitlvl) + #, #0.07, # jitter height
    #Arrange in facets
    ggplot2::facet_wrap( ~ yer_ssn2,
                         drop=FALSE,
                         dir="h",
                         ncol = 2,
                         labeller = label_bquote(cols =
                               .(as.character(yer_ssn2))
                         ) ) +
    # alter the them strip above
    theme(strip.text = element_text(#face = "bold",
      color = "black",
      hjust = 0,
      size = 8),
      strip.background = element_rect(fill = c("white"),
                                      #linetype = "solid",
                                      color = "white",
                                      linewidth = 1)) +
    # ggplot2::facet_wrap( ~ Lat_Species,
    #                      drop=FALSE,
    #                      ncol = 4,
    #                      labeller = as_labeller(appender))+
    #define limits of the plot 
    ggplot2::coord_sf(xlim = c(6.6, 17.2),
                      ylim = c(54.2, 58.4), 
                      expand = FALSE) +
    # # Add points for sampling locations -  this has to be from a separate data frame
    # The 'shape=3' makes crosses for sampled locations without eDNA
    geom_point(data=df_A07.ne,aes(x=Longitude,y=Latitude),
               shape=3,colour=alpha("#000000",1),size=1.6) +
    # The 'shape=21' makes circular points for sampled locations with eDNA present
    geom_point(data=df_A07.ed,aes(x=Longitude,y=Latitude),
               shape=21,colour=alpha("#000000",1),fill=alpha("firebrick3",0.6),
               size=2.4) +
    labs(title = LSpc[i]) +
    # make the  labels turn 90 degrees on the x-axis
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    # or remove them all completely
    #http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()) +
    # legend on bottom
    theme(legend.position="bottom")  + 
    # make the title in italic
    theme(plot.title = element_text(face = "italic", size =10))  
  
  #change axis labels
  #p05t <- p05 + xlab("longitude") + ylab("latitude")
  #p05t <- p05 + xlab("længdegrad") + ylab("breddegrad")
  p05t <- p05 + xlab(" ") + ylab(" ")
  #change the header for the legend on the side, 
  #this must be done for both 'fill', 'color' and 'shape', to avoid 
  #getting separate legends
  # p05t <- p05t + labs(color='')
  # p05t <- p05t + labs(fill='')
  # p05t <- p05t + labs(shape='')
  # p05t <- p05t + labs(size='')
  #get the number of species
  ncat <- length(unique(df_A07$eval))
  # https://github.com/tidyverse/ggplot2/issues/3492
  #repeat 'black' a specified number of times
  filltxc = rep("black", ncat)
  #adjust tick marks on axis
  p05t <- p05t + scale_y_continuous(breaks=seq(54.2, 58.4,2))
  p05t <- p05t + scale_x_continuous(breaks=seq(6.6, 17.2,4))
  
  # see the plot
  p05t
  
  #
  bSaveFigures<-T
  if(bSaveFigures==T){
    ggsave(plot = p05t, 
           filename = paste0(wd00_wd07,"/Fig11_v",ins,"_map_of_",sbs_spcNm_wu,"_detected_2017_to_2023.png"),
           width=210*0.60,height=297,
           #width=210*0.8,height=297,
           #width=297,height=210,
           #width=297,height=210,
           #width=1.4*297,height=210,
           units="mm",dpi=300)
  }
  # end iteration over Latin species names
}

#



uyssn <- unique(dfg05$yer_ssn2)
noss<- length(uyssn)
# Repeat rows making each repeated rows following the original rows and assign new variables for each row [duplicate]
# https://stackoverflow.com/questions/69160320/repeat-rows-making-each-repeated-rows-following-the-original-rows-and-assign-new?noredirect=1&lq=1
library(dplyr)
library(tidyr)
#ensure the months are as numeric values
dfg05$mnt <- as.numeric(gsub("(^.*)-(.*)-(.*$)","\\2",dfg05$Dato_inds))
df_A05$ssnno2 <- df_A05$ssnno

df_A05$ssnno2 <- gsub("1st","jan - jun",df_A05$ssnno)
df_A05$ssnno2 <- gsub("2nd","jul - nov",df_A05$ssnno2)

df_A05$yer_ssn2 <- paste0(df_A05$yer,", ",df_A05$ssnno2)
df_A06 <- df_A05
df_A06 <- df_A06[(df_A06$yer>=2017),]
df_A06 <- df_A06[(df_A06$yer<=2023),]
# define columns to keep
keepcol <- c( "Dato_inds",
              "declat",
              "declon",
              "yer" , 
              "mnt",
              "yer_ssn2",
              "eval",
              "Lat_Species")
# only keep specific columns
df_A06 <- df_A06[keepcol]
dfg06 <- dfg05[keepcol]

#bind rows together in columns
df_A06 <- rbind(df_A06,dfg06)

unique(df_A06$mnt)
# exclude rows if it does not contain a space, and if it contains a zero
df_A06 <- df_A06[grepl(" ",df_A06$Lat_Species),]
df_A06 <- df_A06[!grepl("0",df_A06$Lat_Species),]
# add a column with source cagtegories
df_A06$source <- NA
df_A06$source[grepl("iNat",df_A06$eval)] <- "iNat"
df_A06$source[grepl("eDNA",df_A06$eval)] <- "MONIS6"
# make a data frame that can be used for matching short evaluation and record
# categories to obtain long category names
evc01 <- c("iNat","no eDNA","eDNA found")
evc02 <- c("iNaturalist","ingen DNA","miljø DNA fundet")
df_evc <- as.data.frame(cbind(evc01,evc02))

# define directory path for data files
# read file with records from www.arter.dk
wd00_wddata <- paste0(wd00,"/data")
lst.flwddat <- list.files(wd00_wddata, full.names = T)
inf01 <- lst.flwddat[grepl("fund_arter_dk",lst.flwddat)]
#read the xls file with 
df_fa01  <- readxl::read_excel(inf01, sheet = 2)
# copy and modify column names to make it match the other data frame 
# with MONIS6 records
# make a column that has the record category information
df_fa01$reccat <- "www.arter.dk"
df_fa01$source <- "arter.dk"
# get the longitude and latitude for the sampled locations
df_fa01$declat <- df_fa01$Lat
df_fa01$declon <- df_fa01$Long
# get sampling year and month
df_fa01$yer <- as.numeric(gsub(".*-(.*)","\\1",df_fa01$Observationsdato))
df_fa01$mnt <- as.numeric(gsub(".*\\/(.*)-(.*)","\\1",df_fa01$Observationsdato))
# get the Latin species names
df_fa01$Lat_Species <- df_fa01$`Taxon latinsk navn`
# ensure the latitude and longitude ar numeric
df_A06$declat <- as.numeric(df_A06$declat)
df_A06$declon <- as.numeric(df_A06$declon)
unique(df_A06$eval)
# match to get the long record category for each sample
df_A06$reccat <- df_evc$evc02[match(df_A06$eval,df_evc$evc01)]
# define a vector with columns name to keep
ctkeep <- c("Lat_Species",
            "reccat",
            "yer",
            "mnt",
            "declat",
            "declon",
            "source")
# subset data frame by limiting to only the columns required for the mapping
df_A07 <- df_A06[ctkeep]
df_fa02 <- df_fa01[ctkeep]
# combine tthe data frames by adding rows, since they now have the
# same column names
df_A07 <- rbind(df_A07,df_fa02)
# make sure the months are as numeric values to be able to evaluate on the
# values
df_A07$mnt <- as.numeric(df_A07$mnt)
df_A07$yer <- as.numeric(df_A07$yer)

# substitute in the species names
df_A07$Lat_Species <-  gsub("Crassostrea gigas" ,"Magallana gigas",df_A07$Lat_Species)
df_A07$Lat_Species <- gsub("Colpomenia peregrine" , "Colpomenia peregrina" ,df_A07$Lat_Species)
# make a column that evaluates the season sampled
df_A07$ssn.per <- NA
# evalutate on the season column to add a number version of the season
df_A07$ssn.per[(df_A07$mnt<=6)] <- "jan - jun"
df_A07$ssn.per[(df_A07$mnt>6)]  <- "jul - nov"
# make columns that holds both year and sampled period
df_A07$yer_ssn.per <- paste0(df_A07$yer,", ",df_A07$ssn.per)
df_A07[is.na(df_A07$ssn.per),]

# limit the years to only cover the MONIS6 years sampled
df_A07 <- df_A07[(df_A07$yer<=2023),]
df_A07 <- df_A07[(df_A07$yer>=2017),]
# make color range
fill.f.pnt <- c("firebrick3","steelblue1","khaki3","white")
symb.f.pnt <- c(21,22,25,3)
size.f.pnt <- c(2.2,1.6,1.6,2.8)
orde.f.pnt <- c(1,3,4,2)
cate.f.pnt <- c("miljø DNA fundet","iNaturalist","www.arter.dk","ingen DNA")
colo.f.pnt <- c(rep("black",4))
alph.f.pnt <- c(0.6,1,0.6,0.6)
# combine vectors to a data frame
df_pntcats <- as.data.frame(cbind(fill.f.pnt,
                                  symb.f.pnt,
                                  size.f.pnt,
                                  orde.f.pnt,
                                  colo.f.pnt,
                                  alph.f.pnt,
                                  cate.f.pnt))
# get the number of latin species
LSpc <- unique(df_A07$Lat_Species)
LSpc <- LSpc[order(LSpc)]

# make a sequence for this range
nfLSpc<- seq(1,length(LSpc),1)
# iterate over species
for (i in nfLSpc)
{print(paste0(i, " making plot for ",LSpc[i]))
  # }
  #  i <- 3
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  ins <- stringr::str_pad(i, 2, pad = "0")
  # get the species name with an underscore
  sbs_spcNm_wu <- gsub(" ","_",LSpc[i])
  sbs_spcNm_wu
  #}
  # make a function to use for making facet wrap headers
  # https://ggplot2.tidyverse.org/reference/as_labeller.html
  # https://stackoverflow.com/questions/63857833/changing-the-facet-wrap-labels-using-labeller-in-ggplot2
  #appender <- function(nm1) paste0(df_tx01$class[match(nm1,df_tx01$family )],": ", nm1)
  
  df_A07.1 <- df_A07[grepl(LSpc[i],df_A07$Lat_Species),]
  # add a level for jittering points
  jitlvl <- 0.03
  #make plot
  uyssn <- unique(df_A07$yer_ssn.per)
  noss<- length(uyssn)
  # identify the missing seasons for sampled period
  usss <- unique(df_A07.1$yer_ssn.per) 
  msssn <- uyssn[!(uyssn %in% usss) ]
  nmsssn <- length(msssn)
  nclA07.1 <- ncol(df_A07.1)
  eRw <- rep(0,nclA07.1)
  clNmA07.1 <- colnames(df_A07.1)
  dfeR <- as.data.frame(t((eRw)))
  colnames(dfeR) <- clNmA07.1
  dfeR <- dfeR[rep(seq_len(nrow(dfeR)), each = nmsssn), ]
  ixcN <- which(grepl("yer_ssn.per",colnames(df_A07.1)))
  dfeR[,ixcN] <- msssn
  df_A07.2 <- rbind(df_A07.1,dfeR)
  # make 'iNat' categories for the zeroes for the added years
  df_A07.2$reccat[(df_A07.2$reccat=="0")] <- "iNaturalist"
  # df_A07.2[is.na(df_A07.2$reccat),]
  # # unique(df_A07.2$reccat)
  # # subset data frame
  # df_A08 <- df_A07.2[(df_A07.2$eval=="iNat"),]
  # # make a data frame 
  # df_A08.ne <- df_A07.2[(df_A07.2$reccat=="no eDNA"),]
  # df_A08.ed <- df_A07.2[(df_A07.1$reccat=="eDNA detected"),]
  # 
  # # change the labels to Danish
  # df_A06.2$eval[(df_A06.2$eval=="eDNA detected")] <- "miljø DNA fundet"
  # df_A06.2$eval[(df_A06.2$eval=="no eDNA")] <- "ingen DNA"
  # df_A06.2$eval[(df_A06.2$eval=="iNat")] <- "iNaturalist"
  # # find  the categories for the points to plot 
  # and use for the manual scales
  catf_plt <- unique(df_A07.2$reccat)
  # match to get the fill, the pch symbol, the size and the order
  sz_f_mnscl <- as.numeric(df_pntcats$size.f.pnt[match(catf_plt,df_pntcats$cate.f.pnt)])
  sy_f_mnscl <- as.numeric(df_pntcats$symb.f.pnt[match(catf_plt,df_pntcats$cate.f.pnt)])
  fi_f_mnscl <- df_pntcats$fill.f.pnt[match(catf_plt,df_pntcats$cate.f.pnt)]
  co_f_mnscl <- df_pntcats$colo.f.pnt[match(catf_plt,df_pntcats$cate.f.pnt)]
  al_f_mnscl <- as.numeric(df_pntcats$alph.f.pnt[match(catf_plt,df_pntcats$cate.f.pnt)])
  # get the order for plotting the point
  df_A07.2$or_f_pnt <- df_pntcats$orde.f.pnt[match(df_A07.2$reccat,df_pntcats$cate.f.pnt)]
  df_A07.2$or_f_pnt <- as.numeric(df_A07.2$or_f_pnt)
  # reorder the data frame to be able to plot the  points in
  # an order that allows the iNaturalist and www.arter.dk records to be plottted
  # below the MONIS6 eDNA results
  df_A07.2 <- df_A07.2[order(df_A07.2$or_f_pnt, decreasing = T),]
  #dev.off()
  #make plot
  p05 <- ggplot(data = world) +
    geom_sf(color = "black", fill = "azure3", lwd=0.1) +
    # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
    theme_void() +
    #https://ggplot2.tidyverse.org/reference/position_jitter.html
    #https://stackoverflow.com/questions/15706281/controlling-the-order-of-points-in-ggplot2
    # use 'geom_jitter' instead of 'geom_point' 
    geom_point(data = df_A07.2 ,
                aes(x = declon, 
                    y = declat,
                    #color=reccat,
                    size=reccat,
                    fill=reccat,
                    shape=reccat) ) + #,
                
                # width = jitlvl, #0.07, jitter width 
                # height = jitlvl) + #, #0.07, # jitter height
    #Arrange in facets
    ggplot2::facet_wrap( ~ yer_ssn.per,
                         drop=FALSE,
                         dir="h",
                         ncol = 2,
                         labeller = label_bquote(cols =
                                                   .(as.character(yer_ssn.per))
                         ) ) +
    # alter the them strip above
    theme(strip.text = element_text(#face = "bold",
      color = "black",
      hjust = 0,
      size = 8),
      strip.background = element_rect(fill = c("white"),
                                      #linetype = "solid",
                                      color = "white",
                                      linewidth = 1)) +
    
    # ggplot2::facet_wrap( ~ Lat_Species,
    #                      drop=FALSE,
    #                      ncol = 4,
    #                      labeller = as_labeller(appender))+
    #define limits of the plot 
    ggplot2::coord_sf(xlim = c(6.6, 17.2),
                      ylim = c(54.2, 58.4), 
                      expand = FALSE) +
    # # # Add points for sampling locations -  this has to be from a separate data frame
    # # The 'shape=3' makes crosses for sampled locations without eDNA
    # geom_point(data=df_A07.ne,aes(x=Longitude,y=Latitude),
    #            shape=3,colour=alpha("#000000",1),size=1.6) +
    # # The 'shape=21' makes circular points for sampled locations with eDNA present
    # geom_point(data=df_A07.ed,aes(x=Longitude,y=Latitude),
    #            shape=21,colour=alpha("#000000",1),fill=alpha("firebrick3",0.6),
    #            size=2.4) +
    # Set the colors of the outline of the points manually
    scale_color_manual(values=c(alpha( co_f_mnscl,al_f_mnscl) ), breaks = catf_plt ) +
    # # Set the colors of the fill of the points manually
    #scale_fill_manual(values=c(fi_f_mnscl), breaks = catf_plt) +
    scale_fill_manual(values=c(alpha( fi_f_mnscl,al_f_mnscl) ), breaks = catf_plt ) +
    # Set the shape of the points manually
    scale_shape_manual(values=c(sy_f_mnscl),breaks = catf_plt) +
    # Set the size of the points manually
    scale_size_manual(values=c(sz_f_mnscl),breaks = catf_plt) +
    #
    labs(title = LSpc[i]) +
    # make the  labels turn 90 degrees on the x-axis
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    
    # or remove them all completely
    #http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()) +
    # legend on bottom
    #theme(legend.position="bottom")  + 
    
    # make the title in italic
    theme(plot.title = element_text(face = "italic", size =10))  
  
  #change axis labels
  #p05t <- p05 + xlab("longitude") + ylab("latitude")
  #p05t <- p05 + xlab("længdegrad") + ylab("breddegrad")
  p05t <- p05 + xlab(" ") + ylab(" ")
  #change the header for the legend on the side, 
  #this must be done for both 'fill', 'color' and 'shape', to avoid 
  #getting separate legends
  p05t <- p05t + labs(color='')
  p05t <- p05t + labs(fill='')
  p05t <- p05t + labs(shape='')
  p05t <- p05t + labs(size='')
  #get the number of species
  ncat <- length(unique(df_A07$reccat))
  # https://github.com/tidyverse/ggplot2/issues/3492
  #repeat 'black' a specified number of times
  filltxc = rep("black", ncat)
  #adjust tick marks on axis
  p05t <- p05t + scale_y_continuous(breaks=seq(54.2, 58.4,2))
  p05t <- p05t + scale_x_continuous(breaks=seq(6.6, 17.2,4))
  
  # see the plot
  p05t
  
  #
  bSaveFigures<-T
  if(bSaveFigures==T){
    ggsave(plot = p05t, 
           filename = paste0(wd00_wd07,"/Fig12_v",ins,"_map_of_",
                             sbs_spcNm_wu,"_detected_2017_to_2023.png"),
           width=210*0.70,height=297,
           #width=210*0.8,height=297,
           #width=297,height=210,
           #width=297,height=210,
           #width=1.4*297,height=210,
           units="mm",dpi=300)
  }
  # end iteration over latin species names
}


#
#
