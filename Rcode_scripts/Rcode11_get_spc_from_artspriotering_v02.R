#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-


#On the remote server. Ensure you have a directory where your 
#R packages are stored
# make one here:
# /home/sknu003/R/R_packages_for_Rv4_1
# load the R module
# $ module load R/4.1.0-gimkl-2020a

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

# lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_1"

# # Continue by running these lines, one by one in R
# Sys.setenv(R_LIBS_USER="lib_path01")
# .libPaths("lib_path01")
# # change the path to where the packages should be installed from # see this website: https://stackoverflow.com/questions/15170399/change-r-default-library-path-using-libpaths-in-rprofile-site-fails-to-work
# .libPaths( c( lib_path01 , .libPaths() ) )
# .libPaths()
# .libPaths( c( lib_path01) )

## Or try pasting one long line with all commands, and installation of the 'taxizedb' library #
# lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_1"; Sys.setenv(R_LIBS_USER="lib_path01"); .libPaths("lib_path01"); .libPaths( c( lib_path01 , .libPaths() ) ); .libPaths(); .libPaths( c( lib_path01) ); install.packages(c("taxizedb", "taxize", "tidyverse", "readxl", "worms", "stringr", "dplyr"))

## Use this line here below to install a lot of packages needed for doing population genetic analysis
# lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_1"; Sys.setenv(R_LIBS_USER="lib_path01"); .libPaths("lib_path01"); .libPaths( c( lib_path01 , .libPaths() ) ); .libPaths(); .libPaths( c( lib_path01) ); if(!require("gaston")){install.packages("gaston", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require("hierfstat")){install.packages("hierfstat", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require("pegas")){install.packages("pegas")};if(!require("ips")){install.packages("ips", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require("tidyverse")){install.packages("tidyverse", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require("pals")){install.packages("pals", dependencies = TRUE, INSTALL_opts = '--no-lock')};if(!require(adegenet)){install.packages("adegenet", repos='http://cran.us.r-project.org') };if(!require(apex)){install.packages("apex", repos='http://cran.us.r-project.org') };if(!require(mmod)){install.packages("mmod", repos='http://cran.us.r-project.org') };if(!require(tidyverse)){install.packages("tidyverse", repos='http://cran.us.r-project.org') };if(!require(pals)){install.packages("pals", repos='http://cran.us.r-project.org') };if(!require(ape)){install.packages("ape", repos='http://cran.us.r-project.org') };if(!require(RColorBrewer)) {install.packages("RColorBrewer", repos='http://cran.us.r-project.org') };if(!require(stringi)){install.packages("stringi", repos='http://cran.us.r-project.org') };if(!require(poppr)){install.packages("poppr", repos='http://cran.us.r-project.org') };if(!require(vegan)){install.packages("vegan", repos='http://cran.us.r-project.org')};if(!require(adegenet)){install.packages("adegenet", repos='http://cran.us.r-project.org')};if(!require(biogeo)){install.packages("biogeo", repos='http://cran.us.r-project.org') }
## Use this line here below to install a lot of packages needed for doing mapping
# lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_1"; Sys.setenv(R_LIBS_USER="lib_path01"); .libPaths("lib_path01"); .libPaths( c( lib_path01 , .libPaths() ) ); .libPaths(); .libPaths( c( lib_path01) ); if(!require(scales)){  install.packages("scales")};if(!require(fields)){  install.packages("fields")};if(!require(marmap)){  install.packages("marmap")};if(!require(TeachingDemos)){  install.packages("TeachingDemos")};if(!require(rworldmap)){  install.packages("rworldmap")};if(!require(rworldxtra)){  install.packages("rworldxtra")};require(rworldxtra)if(!require(readxl)){  install.packages("readxl")};if(!require(plyr)){  install.packages("plyr")};if(!require(mapdata)){  install.packages("mapdata")};if(!require(maps)){  install.packages("maps")};if(!require(mapplots)){  install.packages("mapplots")}; if(!require(purrr)){  install.packages("purrr")}

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

# ## Now that you have installed all packages in your remote path : "/home/sknu003/R/R_packages_for_Rv4_1"
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
# # replace my library path to your own library path
# #lib_path01 <- "/groups/hologenomics/phq599/data/R_packages_for_Rv3_6"
# lib_path01 <- "/home/sknu003/R/R_packages_for_Rv4_1"
# Sys.setenv(R_LIBS_USER="lib_path01")
# .libPaths("lib_path01")
# # change the path to where the packages should be installed from # see this website: https://stackoverflow.com/questions/15170399/change-r-default-library-path-using-libpaths-in-rprofile-site-fails-to-work
# .libPaths( c( lib_path01 , .libPaths() ) )
# .libPaths()
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
library("BiocManager")
library("adegenet")
library("ape")
library("fields")
library("mapdata")
library("mapplots")
library("maps")
library("marmap")
library("pillar")
library("plyr")
library("poppr")
library("purrr")
library("RColorBrewer")
library("readxl")
library("readxl")
library("rworldmap")
library("rworldxtra")
library("scales")
library("stringi")
library("stringr")
library("TeachingDemos")
library("tibble")
library("tidyverse")
library("vegan")
library("worrms")



#===============================================================================
# worms functions - start
#===============================================================================
#Use the libraries
library(stats)
library(tidyverse)
library(jsonlite) #https://cran.r-project.org/web/packages/jsonlite/
library(httr)

wd00 <- getwd()
# wd00 <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6/Analysis_MONIS5_to_6_v01"
# setwd(wd00)

#define directory with output flies
wdout <- "output11_get_species_from_priority_table"
# define full path for input directory
outdir10 <- paste(wd00,wdout, sep="/")
# remove previous versions of the 'outdir10'
unlink(outdir10, force = T, recursive = T)
# create the 'outdir10' in a new version
dir.create(outdir10)
# define directory with input files to read in
wddata <- "data"
#paste together a path for the directory that holds the inputs files
wd00_wddata <- paste(wd00,wddata,sep="/")

#ErrorList <- read.table("badnames.txt",sep=";",header=F,stringsAsFactors=F,encoding="Windows-1251")
#ErrorList <- ErrorList[,1]

BadName<-function(name){
  #browser()
  
  ErrorList <- c("",
                 "Andre grønalger",
                 "Blomsterplanter, uspec.",
                 "Blågrøn",
                 "Blågrønalger og bakterier",
                 "Blegrxnalge",
                 "Blegrxnalgelignende",
                 "Blågrønalge",
                 "Brown bush",
                 "Brown crust",
                 "Brunskorper",
                 "Enårige trådformede grøn-",
                 "brunalger",
                 "Div. blegrxnalger",
                 "Div. blågrønalger",
                 "Div. centriske kiselalger",
                 "Div. furealger",
                 "Div. kiselalger",
                 "Diverse",
                 "Ebria-lign. flagellat",
                 "fasthæftede oprette alger",
                 "Fedtmøg - Pilayella, Ectocarpus",
                 "Fingrenet grøn busk",
                 "Green endophyte",
                 "Green endozoic",
                 "Grønalger kalkborende",
                 "Gulgrøn",
                 "Hjuldyr fg",
                 "Hjuldyr æg",
                 "Kolonidannende chlorococcale grxnalger",
                 "Kolonidannende chlorococcale grønalger",
                 "minibakterier",
                 "Nxgen furealge",
                 "Nxgne furealger",
                 "Nøgen furealge",
                 "Nøgen furealger",
                 "Ovale blegrxnalgeceller",
                 "Ovale blågrønalgeceller",
                 "Ovale blågrønalgeceller",
                 "Ovale chlorococcale grxnalger,",
                 "Ovale chlorococcale grønalger,",
                 "Ingen arter registreret",
                 "Kiselalger",
                 "løstliggende alger",
                 "muslinge ånderør",
                 "Red bush",
                 "Red calcified crust",
                 "Red crust",
                 "Rødalger på lavt vand",
                 "Rødkødskorp, tyk",
                 "Rødkødskorp, tynd",
                 "Små makrofytter",
                 "svovlbakterier",
                 "Terebellidae",
                 "Trådalger-rø,br,gr - generelt",
                 "Trådform.enår.rødalg",
                 "Trådformede enårige brunalger",
                 "Trådformede enårige grønalger",
                 "Trådformede enårige rødalger",
                 "Bakterier",
                 "Ciliater (fra fytoplanktonundersøgelser)",
                 "Diverse (fra fytoplanktonundersøgelser)",
                 "Flagellate",
                 "Flagellate x",
                 "Monader",
                 "Nanoflagellates, unidentified",
                 "Nanoplankton, unidentified",
                 "Picoplankton, unidentified",
                 "Ultrananoplankton, unidentified",
                 "Unidentified",
                 "Unidentified flagellates",
                 "Unidentified flagellates",
                 "Unidentifred ciliates",
                 "Unidentifred oligotrich ciliates",
                 "Unidentifred tintinnids",
                 "Nauplie",
                 "Parasitisk copepod",
                 "Unidentified",
                 "Unidentifred ciliates",
                 "Virus")
  
  #ErrorList <- iconv(ErrorList, from="Windows-1252", to="UTF-8")
  
  if(tolower(name) %in% tolower(ErrorList)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

FixNames<-function(name){
  # fixes to name used before searching in WoRMS
  name<-gsub("  "," ",name) # replace two spaces with one
  name<-gsub(" spp\\.","",name)
  name<-gsub(" sp\\.","",name)
  name<-gsub(" gr\\.","",name)
  name<-gsub(" indet\\.","",name)
  name<-gsub(" indet","",name)
  name<-gsub(" s\\.l\\.","",name)
  name<-gsub("Acrchaetium","Acrochaetium",name)
  name<-gsub(" tetrasporophyte","",name)
  name<-gsub("tetrasporophyte","",name)
  name<-gsub(" \\(æg\\)","",name)
  name<-gsub(" \\(egg\\)","",name)
  name<-gsub(" æg","",name)
  name<-gsub(" egg","",name)
  name<-gsub(" eggs","",name)
  name<-gsub(" nauplier","",name)
  name<-gsub(", cysts","",name)
  name<-gsub(" cysts","",name)
  name<-gsub(" cysts","",name)
  name<-gsub("pungitus","pungitius",name)
  name<-gsub("spioniidae","spionidae",name)
  name<-gsub("f. scrubsolei","",name)
  name<-gsub("f. semispina","",name)
  name<-gsub("f. truncata","",name)
  name<-gsub("v. variabile","",name)
  name<-gsub("-group","",name)
  name<-gsub(" group","",name)
  name<-gsub(" <celle>","",name)
  name<-gsub("enkeltceller","",name)
  name<-gsub("filament","",name)
  name<-gsub(", Børsteorme","",name)
  name<-gsub("Harmothoë","Harmothoe",name)
  name<-gsub("Tetraëdron","Tetraedron",name)
  name<-gsub("HYDRACARINA I, Vandmider","Hydracarina",name)
  name<-gsub("MYSIDACEA, Mysider, Kårer","Mysidacea",name)
  name<-gsub("OSTRACODA, Muslingekrebs","Ostracoda",name)
  name<-gsub("ANTHOZOA, Koraldyr","Anthozoa",name)
  name<-gsub("ECHINODERMATA II, Ophiuroidea Slangestje","Echinodermata",name)
  name<-gsub("GASTROPODA, PROSOBRANCHIA, NUDIBRANCHIA","Nudibranchia",name)
  name<-gsub("HYDROZOA, athecata","Athecata",name)
  name<-gsub("anguste-tabulatum","angustitabulatum",name)
  name<-gsub("flos-aquae","flosaquae",name)
  name<-gsub("BIVALVIA, muslinger","Bivalvia",name)
  name<-gsub("NEMATODA, Rundorme","Nematoda",name)
  name<-gsub("POLYCHAETA, Havbørsteorme","Polychaeta",name)
  name<-gsub("NEMERTINI Slimbændler","Nemertini",name)
  name<-gsub("HYDROZOA, Polypdyr","Hydrozoa",name)
  name<-gsub("Chaetoceros, phaeoceros-group","Chaetoceros (Phaeoceros)",name)
  name<-gsub(", delicatissima"," delicatissima",name)
  name<-gsub(", seriata"," seriata",name)
  name<-gsub("på veg/sten","",name)
  name<-gsub("på alger","",name)
  name<-gsub("gul på sten","",name)
  name<-gsub("på sten","",name)
  name<-gsub("på sediment","",name)
  name<-gsub("på veg","",name)
  name<-gsub("spp., båndformer","",name)
  name<-gsub("5-10 µm","",name)
  name<-gsub("400-500 µm","",name)
  name<-gsub("0,5-0,75 mm","",name)
  name<-gsub("Centriske kiselalge.*","",name)
  name<-gsub("Ovale.*","",name)
  name<-gsub("Blågrøn.*","",name)
  name<-gsub("Amphithoe rubricata","Ampithoe rubricata",name)
  name<-gsub("Abranchus microstomus","Abranchus microstoma",name)
  name<-gsub("Acrchaetium savianum","Acrochaetium savianum",name)
  #name<-gsub("","",name)
  
  #name<-gsub("","",name)
  return(name)
}

# <- 369346
# Anisocladella: AphiaID=369346 (Valid ID)
#searchtext <-  "Anisocladella"
GetSpeciesIDFuzzy<-function(searchtext){
  # ---------- get the AphiaID from the search text using fuzzy match  -----------------------------------------------
  #Build the URL to get the data from
  searchtext2 <- gsub(" ","%20",searchtext)
  url<-sprintf("http://marinespecies.org/rest/AphiaRecordsByMatchNames?scientificnames[]=%s&marine_only=false",searchtext2)
  x<-http_status(GET(url))
  ID<-""
  if(x$reason!="OK"){
    cat(paste0(searchtext,": ",x$reason,"\n"))
    return(ID)
  }
  AphiaRecord <- fromJSON(url) 
  if(length(AphiaRecord)>0){
    AphiaRecord <- as.data.frame(AphiaRecord[[1]])
    if(nrow(filter(AphiaRecord,status=="Accepted"))>0){
      AphiaRecord <- AphiaRecord %>% 
        filter(status=="Accepted")
    }
    if(!is.na(AphiaRecord$scientificname[1])){
      ID<-AphiaRecord$AphiaID[1]
      cat(paste0(searchtext,": AphiaID=",ID,"\n"))
    }
  }
  return(ID)
  
}

GetSpeciesID<-function(searchtext){
  # ---------- get the AphiaID from the search text -----------------------------------------------
  #Build the URL to get the data from
  
  df <- list(Species=searchtext,
             AphiaID=NA,
             ScientificName=NA,
             Rank=NA,
             ParentID=NA)
  
  if(BadName(searchtext)==TRUE){
    # Species is in a list of names that won't return anything useful
    # e.g. "Blågrønalger og bakterier","Brown bush",
    cat(paste0("'",searchtext,"' is on the list of bad names. No search made\n"))
    return(df)
  }
  
  if(grepl("\\/", searchtext)==TRUE){
    # invalid character
    return(df)
  }
  
  searchtext2 <- gsub(" ","%20",searchtext)
  url<-sprintf("http://marinespecies.org/rest/AphiaIDByName/%s?marine_only=false",searchtext2)
  
  x<-http_status(GET(url))
  if(x$reason!="OK"){
    cat(paste0(searchtext,": ",x$reason," (trying fuzzy search instead)\n"))
    
    AphiaID <- GetSpeciesIDFuzzy(searchtext)
    if(AphiaID==""){
      return(df)
    }
  }else{
    #Get the AphiaID
    AphiaID <- fromJSON(url)
    cat(paste0(searchtext,": AphiaID=",AphiaID))
  }
  
  # ---------- get the Aphia record from the AphiaID -----------------------------------------------
  
  
  url<-sprintf("http://marinespecies.org/rest/AphiaRecordByAphiaID/%d",AphiaID)
  url
  AphiaRecord <- fromJSON(url)
  
  
  validID<-AphiaRecord$valid_AphiaID
  if(is.null(validID)){
    cat(paste0(" (No validID in record)\n"))
  }else{
    if(validID != AphiaID){
      cat(paste0(" (Using AphiaID=",validID,")\n"))
      
      AphiaIDorig <- AphiaID
      AphiaID <- validID
      AphiaRecordOrig<-AphiaRecord
      
      # get the correct record
      url<-sprintf("http://marinespecies.org/rest/AphiaRecordByAphiaID/%d",AphiaID)
      AphiaRecord <- fromJSON(url)
    }else{
      cat(paste0(" (Valid ID)\n"))
    }
  }
  
  df$AphiaID<-AphiaID
  df$ScientificName<-AphiaRecord$scientificname
  df$Rank<-AphiaRecord$rank
  df$ParentID<-AphiaRecord$parentNameUsageID
  return(df)
}

#AphiaID=206343
#AphiaID=1496728
#GetSpeciesDistributions(1496728)
GetSpeciesDistributions<-function(AphiaID){
  
  url<-sprintf("http://marinespecies.org/rest/AphiaRecordByAphiaID/%d",AphiaID)
  AphiaRecord <- fromJSON(url)
  cat(paste0(AphiaRecord$scientificname," [",AphiaID,"]\n"))
  
  # ---------- get the Aphia synonyms -----------------------------------------------
  
  url<-sprintf("http://marinespecies.org/rest/AphiaSynonymsByAphiaID/%d?offset=1",AphiaID)
  
  x<-http_status(GET(url))
  if(x$reason=="OK"){
    dfSynonyms <- fromJSON(url)
    synonyms<-paste0(dfSynonyms$scientificname,collapse=", ")
    cat(paste0("  Synonyms: ",synonyms,"\n"))
    bSynonyms = TRUE
  }else{
    cat(paste0("  No synonyms found\n"))
    bSynonyms = FALSE
  }
  
  # ---------- Get all distributions for a given AphiaID -----------------------------------------------
  
  # loop through AphiaID for synonyms
  cat("  Get distributions:\n")
  nc <- nchar(AphiaRecord$scientificname)
  if(bSynonyms==TRUE){
    ncs <- max(nchar(dfSynonyms$scientificname),na.rm=T)
    if(ncs>nc){
      nc<-ncs
    }
  }
  
  url <- sprintf("http://marinespecies.org/rest/AphiaDistributionsByAphiaID/%d", AphiaID);
  x<-http_status(GET(url))
  if(x$reason=="OK"){
    Synonym<-AphiaRecord$scientificname
    spaces<-paste(replicate(nc-nchar(Synonym), " "), collapse = "")
    
    bFoundDistribution=TRUE
    distribution <- fromJSON(url)
    distribution$AphiaID <- AphiaID
    distribution$SynonymID <- AphiaID
    distribution$Synonym <- Synonym
    
    cat(paste0("   ",Synonym," [",AphiaID,"]: ",spaces, nrow(distribution)," records \n"))
  }else{
    bFoundDistribution=FALSE
  }
  
  if(bSynonyms==TRUE){
    
    for(id in dfSynonyms$AphiaID){
      Synonym<-dfSynonyms$scientificname[dfSynonyms$AphiaID==id]
      spaces<-paste(replicate(nc-nchar(Synonym), " "), collapse = "")
      cat(paste0("   ",Synonym," [",id,"]: ",spaces))
      url <- sprintf("http://marinespecies.org/rest/AphiaDistributionsByAphiaID/%d", id);
      x<-http_status(GET(url))
      if(x$reason=="OK"){
        distributionSynonym <- fromJSON(url)
        
        cat(paste0(nrow(distributionSynonym)," records \n"))
        
        distributionSynonym$AphiaID <- AphiaID
        distributionSynonym$SynonymID <- id
        distributionSynonym$Synonym <- Synonym
        if(bFoundDistribution==TRUE){
          distribution <- bind_rows(distribution,distributionSynonym)
        }else{
          distribution <- distributionSynonym
          bFoundDistribution=TRUE
        }
        
      }else{
        cat(paste0("no records \n"))
      }
    }
  }
  
  if(!exists("distribution")){
    
    # check for null values in the record
    if(is.null(AphiaRecord$kingdom)){
      kingdom<-NA
    }else{
      kingdom<-AphiaRecord$kingdom
    }
    if(is.null(AphiaRecord$phylum)){
      phylum<-NA
    }else{
      phylum<-AphiaRecord$phylum
    }
    if(is.null(AphiaRecord$class)){
      class<-NA
    }else{
      class<-AphiaRecord$class
    }
    if(is.null(AphiaRecord$order)){
      order<-NA
    }else{
      order<-AphiaRecord$order
    }
    if(is.null(AphiaRecord$family)){
      family<-NA
    }else{
      family<-AphiaRecord$family
    }
    if(is.null(AphiaRecord$genus)){
      genus<-NA
    }else{
      genus<-AphiaRecord$genus
    }
    
    distribution <- data.frame(
      AphiaID=AphiaID,
      SynonymID=AphiaID,
      ScientificName=AphiaRecord$scientificname,
      Synonym=AphiaRecord$scientificname,
      Kingdom=kingdom,
      Phylum=phylum,
      Class=class,
      Order=order,
      Family=family,
      Genus=genus,
      stringsAsFactors=FALSE)
    
  }else{
    distribution$ScientificName <- AphiaRecord$scientificname
    distribution$Kingdom <- AphiaRecord$kingdom
    distribution$Phylum <- AphiaRecord$phylum
    distribution$Class <- AphiaRecord$class
    distribution$Order <- AphiaRecord$order
    distribution$Family <- AphiaRecord$family
    distribution$Genus <- AphiaRecord$genus
  }
  
  
  nameslist <- c("ScientificName","AphiaID","Synonym","SynonymID","Kingdom","Phylum","Class","Order","Family","Genus")
  
  nameslist2 <- names(distribution)[!names(distribution) %in% nameslist ]
  nameslist <- c(nameslist,nameslist2)
  nameslist<-nameslist[nameslist %in% names(distribution)]
  distribution <- distribution[,nameslist]
  
  return(distribution)
  
  
  #distribution <- distribution %>%
  #  filter(!is.na(establishmentMeans))
  
  # http://marinespecies.org/aphia.php?p=manual#topic22
  #
  # Images, specimens, vernaculars, notes, distributions, taxa, … carry quality indicator icons.
  # 
  # Checked: verified by a taxonomic editor
  # Trusted: edited by a thematic editor
  # Unreviewed: has not been verified by a taxonomic editor
  # 
  
}



GetAphiaRecord<-function(AphiaID){
  url<-sprintf("http://marinespecies.org/rest/AphiaRecordByAphiaID/%d",AphiaID)
  x<-http_status(GET(url))
  if(x$reason=="OK"){
    AphiaRecord <- fromJSON(url)
    cat(paste0(AphiaRecord$scientificname," [",AphiaID,"]\n"))
    return(unlist(AphiaRecord))
  }else{
    cat(paste0("no records \n"))
    return(c(rank="",parentNameUsageID=NA))
  }
}

# return the parent AphiaID for a given AphiaID
GetParentAphiaID<-function(AphiaID){
  record<-GetAphiaRecord(AphiaID)
  if(is.null(record)){
    return(NA)
  }else{
    return(as.numeric(record["parentNameUsageID"]))
  }
}

# return the parameter from an Aphia record for a given AphiaID
GetAphiaParameter<-function(AphiaID,parameter){
  record<-GetAphiaRecord(AphiaID)
  if(is.null(record)){
    return(NA)
  }else{
    return(record[[parameter]])
  }
}

# return list of child records for an AphiaID
GetAphiaChildren <- function(AphiaID){
  n<-1
  urlbase<-sprintf("http://marinespecies.org/rest/AphiaChildrenByAphiaID/%d?offset=",AphiaID)
  # max 50 rows can be retreived at one time
  repeat{
    url <- paste0(urlbase,n)
    x<-http_status(GET(url))
    if(x$reason!="OK"){
      break
    }else{
      tempAphiaRecord <- fromJSON(url)
      if(n==1){
        AphiaRecord <- tempAphiaRecord
      }else{
        AphiaRecord <- bind_rows(AphiaRecord,tempAphiaRecord)
      }
    }
    n<-n+500
  }
  if(n>1){
    return(AphiaRecord)
  }else{
    cat(paste0("no records \n"))
    return(NULL)
  }
}

GetAphiaSynonymIDs <- function(AphiaID){
  n<-1
  urlbase<-sprintf("http://marinespecies.org/rest/AphiaSynonymsByAphiaID/%d?offset=",AphiaID)
  # max 50 rows can be retreived at one time
  repeat{
    url <- paste0(urlbase,n)
    x<-http_status(GET(url))
    if(x$reason!="OK"){
      break
    }else{
      tempAphiaRecord <- fromJSON(url)
      if(n==1){
        AphiaRecord <- tempAphiaRecord
      }else{
        AphiaRecord <- bind_rows(AphiaRecord,tempAphiaRecord)
      }
    }
    n<-n+50
  }
  if(n>1){
    return(AphiaRecord)
  }else{
    cat(paste0("no records \n"))
    return(NULL)
  }
}

GetVernacularNames<-function(AphiaID,language_select=c("dan","eng"),include_lang_name=F){
  df <- GetAphiaSynonymIDs(AphiaID)
  ids <-  c(AphiaID,df$AphiaID)
  res <- data.frame(AphiaID=integer(),SynonymID=integer(),vernacular=character(),language_code=character(),language=character())
  for(id in ids){
    url<-sprintf("http://marinespecies.org/rest/AphiaVernacularsByAphiaID/%d",id)
    x<-http_status(GET(url))
    if(x$reason=="OK"){
      restemp <- fromJSON(url)
      restemp$AphiaID<-AphiaID
      restemp$SynonymID<-id
      res <- bind_rows(res,restemp)
    }
  }
  res <- res %>%
    filter(language_code %in% language_select) 
  res <- res %>%
    select(AphiaID,SynonymID,language_code,language,vernacular)
  if(!include_lang_name){
    res <- res %>%
      select(-language)
  }
  return(res)
}

GetVernacularNamesList <- function(AphiaList){
  res <- data.frame(AphiaID=integer(),SynonymID=integer(),language_code=character(),vernacular=character())
  for(i in AphiaList){
    temp <- GetVernacularNames(i)
    cat(paste0(i,": ",nrow(temp)," rows found\n"))
    res <- bind_rows(res,temp)
  }
  return(res)
}

#===============================================================================
# worms functions - end
#===============================================================================
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
source(paste0(wd00,"/Rcode_scripts/worms.R"))



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
#make the list of data frames a single data frame
df_ac02 <- as.data.frame(do.call(rbind,lstap.f))
df_ac02 <- df_ac02[!is.na(df_ac02$AphiaID),]
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