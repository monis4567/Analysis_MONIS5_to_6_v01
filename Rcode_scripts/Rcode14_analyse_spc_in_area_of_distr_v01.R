#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# set path to object
#wd00 <- "/home/hal9000/Documents/shrfldubuntu18/MONIS6_getNCBIseq"
#
#setwd(wd00)
#getwd()
# if(!require(worms)){
#   install.packages('worms')
# }
# library("worms")

#load required libraries
library(purrr)
library(tidyverse)
library(ape)

#define working directory
#wd00 <- "/home/hal9000/Documents/shrfldubuntu18/MONIS6_getNCBIseq"
#setwd(wd00)
wd00 <- getwd()
setwd(wd00)
# required function GetSpeciesDistributions(AphiaID)

library(worrms)
# get the extra codes with functions
# required function GetSpeciesDistributions(AphiaID)
source(paste0(wd00,"/Rcode_scripts/worms_safe.R"))
#
wd13 <- "output13_spc_in_area_of_distribution"
wd00_wd13 <- paste0(wd00,"/",wd13)

#define directory with output flies
wdout <- wd13
# define full path for input directory
outdir <- paste(wd00,wdout, sep="/")
# remove previous versions of the 'outdir10'
unlink(outdir, force = T, recursive = T)
# create the 'outdir10' in a new version
dir.create(outdir)
#
wd_data <- "data"
# paste together a path for the directory with input data in
wd00_wddata <- paste0(wd00,"/",wd_data)
#Define input file name
inpfnm1 <- "lst_NE_Europe_priority_spc.csv"
inpfnm2 <- "df_distr02.csv"
# paste path and file name together
wd00_inpfnm1 <- paste(wd00_wddata,"/",inpfnm1,sep="")
#wd00_inpfnm1 <- paste(wd00,"/",inpfnm1,sep="")
wd00_inpfnm2 <- paste(wd00_wddata,"/",inpfnm2,sep="")
# read in the csv file with species from the worms data base
df_NEEps01 <- read.csv(wd00_inpfnm1, sep=";", header=T,fill=T)
# read in the csv file with Marien regions close to Scandinavia
df_ds02 <- read.csv(wd00_inpfnm2, sep=";", header=T,fill=T)
#nrow(df_ds02)
# also read in the list of species that are prioritized for monitoring
# define input filename
pthinf01 <- "Artsprioritering_v01_2022may.xlsx"
pthinf01 <- paste0(wd00_wddata,"/",pthinf01)
# read in excel file as tibble, skip 2 rows
tibl_inx01 <- readxl::read_xlsx(pthinf01)
# get the tibble but not the first row
tibl_inx01 <- tibl_inx01[-1,]
# remove all non alphanumeric characters: https://stackoverflow.com/questions/10294284/remove-all-special-characters-from-a-string-in-r
tspE <- stringr::str_replace_all(tibl_inx01$Species, "[^[:alnum:]]", " ")
# split string and get lists nested in a list
gl <- strsplit(as.character(tspE), " ")
#get first and second element of nested list
gnl <- sapply(gl, "[[", 1)
spl <- sapply(gl, "[[", 2)
# paste genus and species name together
gspl <- paste0(gnl," ",spl)
# add back to tibble
tibl_inx01$NIS_SpeciesNm <- gspl
tibl_inx01$NIS_GenusNm <- gnl
# make list of unique genera_species names
uspl <- unique(gspl)
# also get related species for some extra NIS 
extNIS <- c("Chattonella subsalsa",
            "Chattonella",
            "Ostreopsidoideae",
            "Alexandrium",
            "Alexandrium ostenfeldii",
            "Ctenopharyngodon idella",
            "Prymneisum parvum",
            "Prymnesiaceae",
            "Prymnesiales",
            "Prymnesiophycidae",
            "Prymnesium")
# combine them into a vector
uspl <- c(extNIS,uspl)

spc <- uspl[1]
# make an empty list to add to
lst_uAphIDs <- list()
#iterate over list of species names to get Aphia ID numbers
for (spc in uspl){
# get Aphia info for the species
  aphiaInfo<-GetSpeciesID(spc)
  # get the AphiaID for the species
  aphiaID <- aphiaInfo$AphiaID
  # 
  lst_uAphIDs[[spc]] <- as.numeric(aphiaID)
  # 
}

AphiaList <- lst_uAphIDs
#make the list a data frame
df_Aphlst <- do.call(rbind.data.frame, AphiaList)
colnames(df_Aphlst) <- "AphiaIDno"
# get the genus and species names
nms.splt <- strsplit(names(AphiaList), " ")
# get the first name
refNm01 <- sapply(nms.splt, "[[", 1)
# get the index number for the genus-species names that are without a space
# these are probably families and orders, and the 'GetSpeciesDistributions'
# function requires genus-species names, and the idea is to get species 
# distributions for the species, not for the families and orders
idxn.nnm <- which(!(grepl(" ",names(AphiaList))))
# subset the 'AphiaList' to exclude the names that are without a space
AphiaList <- AphiaList[-idxn.nnm]
# get the index number for Chaetoceros
idxn.nnm_Chaetoceros <- which(grepl("Chaetoceros peruvianus",names(AphiaList)))
AphiaList <- AphiaList[idxn.nnm_Chaetoceros]
# make a vector of FALSE værdier med samme længde som AphiaList
SearchSynonyms <- rep(FALSE,length(AphiaList))
#df_us_distr01 <- purrr::map2_dfr(AphiaList,SearchSynonyms,GetSpeciesDistributions)
distributions <- purrr::map_dfc(AphiaList, GetSpeciesDistributions) %>% 
  filter(AphiaID == SynonymID)
# copy the vector into a new vector
df_us_distr01 <- distributions
# limit to only include valid records
df_us_distr01 <- df_us_distr01[df_us_distr01$recordStatus=="valid",]
df_NEEps01 <- df_NEEps01[df_NEEps01$recordStatus=="valid",]
df_ds02 <- df_ds02[df_ds02$recordStatus=="valid",]
# identify the unique family names and order names
fam_nis<- unique(df_us_distr01$Family)
ord_nis<- unique(df_us_distr01$Order)
sciNm_nis <- unique(df_us_distr01$ScientificName)
# limit the distribution of the NIS to only comprise the species listed by
# MST as prioritized
df_us_distr02 <- df_us_distr01[df_us_distr01$ScientificName %in% uspl,]
# keep only rows that are unique for the species name
df_us_distr01 <- df_us_distr01[!duplicated(df_us_distr01[ , c("ScientificName")]), ]

#df_NEEps01[grepl("Pry",df_NEEps01$Class),]

# only retian the part of the data frame where the family name 
# is the same as for the nonindegenous species  
df_NEEps02 <- df_NEEps01[df_NEEps01$Family %in% fam_nis,]
df_NEEps02 <- df_NEEps01[df_NEEps01$Order %in% ord_nis,]
# check the column names
clNMNE02 <- colnames(df_NEEps02)
clNMud01 <- colnames(df_us_distr01)
clNMud01  %in% clNMNE02
# subset the data frame to the same columns that are in the other
# data frame
df_us_distr02 <- df_us_distr01[clNMNE02]
# bind the data frames together by row. To include the nonindegenous species 
df_NEEps03 <- rbind(df_NEEps02,df_us_distr02)

Fam_nis <- df_NEEps03$Family[grepl("Prymnes",df_NEEps03$ScientificName)]
#Ord_nis <- df_NEEps03$Order[grepl("Prymnes",df_NEEps03$ScientificName)]
df_NEEps03$ScientificName[which(Fam_nis==df_NEEps03$Family)]
#df_NEEps03$ScientificName[which(Ord_nis==df_NEEps03$Family)]

# keep only rows that are unique for the species name
df_NEEps03 <- df_NEEps03[!duplicated(df_NEEps03[ , c("ScientificName")]), ]

#unique(df_NEEps03[c('Family', 'ScientificName')])
df_uNEEps03 <- df_NEEps03[!duplicated(df_NEEps03[c('Family', 'ScientificName')]),]
# define a list of North Eastern Atlantic location names
NEAloc<- c("European waters (ERMS scope)",
  "North West Atlantic",
  "Belgian Eclusive Economic Zone",
  "North Sea",
  "Baie de la Seine",
  "Zeebrugge",
  "English Channel",
  "Estonian Eclusive Economic Zone",
  "Devon",
  "Swedish Eclusive Economic Zone",
  "Danish part of the North Sea",
  "Roscoff",
  "Dutch Eclusive  Economic Zone",
  "Nieuwpoort",
  "Oostende",
  "Gulf of Maine",
  "Baltic sea",
  "East Coast of England",
  "St-Jean-de-Luz",
  "Neeltje Jans, Buitenhaven",
  "France",
  "East Gulf of Finland",
  "North Atlantic",
  "Dutch part of the North Sea",
  "Germany",
  "Oostduinkerke",
  "Polish Eclusive Economic Zone",
  "North East Atlantic",
  "German Bight",
  "Southeastern North Atlantic",
  "East North Atlantic",
  "French part of the Bay of Biscay",
  "Denmark",
  "Belgian part of the North Sea",
  "Atlantic Ocean")
# make a column for non target species
df_uNEEps03$DK_NIS <- "nontargetspecies"
# replace the non target species  that are "targetspecies"
df_uNEEps03$DK_NIS[(df_uNEEps03$AphiaID  %in% df_Aphlst$AphiaIDno)] <- "targetspecies"
#check if the location name is on the list
# and limit to only comprise these species
df_uNEEps04 <- df_uNEEps03[(df_uNEEps03$locality %in% NEAloc),]

# define columns to keep
clke <- c("ScientificName","AphiaID",
          "Phylum","Class","Order","Family",
          "Genus","locality",
          "higherGeography" ,"DK_NIS")
# only keep defined columns
df_uNEEps05 <- df_uNEEps04[clke]

#View(df_uNEEps05)
# get the family names and the corresponding target species names
tspc.fam <- df_uNEEps05$Family[(df_uNEEps05$DK_NIS=="targetspecies")]
tspc.scNm <- df_uNEEps05$ScientificName[(df_uNEEps05$DK_NIS=="targetspecies")]
tspc.gnNm <- df_uNEEps05$Genus[(df_uNEEps05$DK_NIS=="targetspecies")]
# bind  as columns into a data frame
df_ts_FsciNm <- as.data.frame(cbind(tspc.fam,tspc.gnNm,tspc.scNm))

# rearrange the data frame by taxonomy
df_uNEEps05.1 <- df_uNEEps05 %>% dplyr::arrange(Phylum, 
                                                Class, 
                                                Order, 
                                                Family, 
                                                Genus, 
                                                ScientificName)
outflNm01 <- paste0(wd00_wd13,"/list_of_target_spc_and_nontarget_spc.csv")
write.table(df_uNEEps05.1,
            file=outflNm01,
            quote=F,
            sep=","
            )
#
#df_ts_FsciNm$tspc.fam %in% df_uNEEps05.1$Family
#https://cran.r-project.org/web/packages/tableHTML/vignettes/tableHTML.html
# if(!require(tableHTML)){
#   install.packages("tableHTML")
#   library(tableHTML)
# }
library(tableHTML)
# # see : https://ardata-fr.github.io/flextable-book/index.html
# if(!require("flextable")){
#   install.packages("flextable", dependencies = TRUE, INSTALL_opts = '--no-lock')
# }
library("flextable")

# # IN ubuntu in a terminal, first install
# sudo add-apt-repository -y ppa:cran/imagemagick
# sudo apt-get update
# sudo apt-get install -y libmagick++-dev
# see : https://davidgohel.github.io/flextable/reference/as_image.html
# if(!require("magick")){
#   install.packages("magick", dependencies = TRUE, INSTALL_opts = '--no-lock')
# }
library("magick")

# # Insalling "xlsx" requires that you first run in a terminal:
# # sudo apt install libbz2-dev
# if(!require("xlsx")){
#   install.packages("xlsx", dependencies = TRUE, INSTALL_opts = '--no-lock')
# }
library("xlsx")
# # install "remotes"
# if(!require("remotes")){
#   install.packages("remotes", dependencies = TRUE, INSTALL_opts = '--no-lock')
# }
library("remotes")
library("patchwork")

# make a table caption, using the list of abbreviations for locations
table_capt01 <- paste("Tabel 1. Artsliste over invasive arter og tætte slægtninge.")
# show the table
t.HTML03 <- df_uNEEps05 %>%
  htmlTable::addHtmlTableStyle(align = c(rep("l",ncol(df_uNEEps05)))) %>%
  htmlTable::htmlTable(caption = table_capt01, rnames = FALSE)
t.HTML03

#df_uNEEps05
# Jeg og Kim har lige haft en snak om prioriteringen og ønsker følgende arter
#prioriteret i den anførte rækkefølge.
# 
# 1. Beroe ovata (Højeste harmonia score af de mulige arter, kan ikke fanges i NOVANA overvågning derfor høj prioritet, dog højt PriNmb)
# 2. Anguillicola crassus (Middelhøj harmonia score, kan ikke fanges i NOVANA overvågning, samfundsmæssig særlig interessant pga. den negative effekt på ål)
# 3. Mytilopsis leucophaeata (Høj harmonia score, interessant da den ikke har været obs. i NOVANA, lavt PriNmb, dog brakvandsart og derfor begrænset hvor meget den kan påvirke det marine miljø)
# 4. Crepidula fornicata (Middelhøj harmonia score, lavt PriNmb, dog observeres den allerede i NOVANA 261 obs inde på Miljøportalen)
# 5. Crassostrea virginica (hvis der er tid/midler til) (Lavt harmonia score, ikke observeret i NOVANA overvågningenm lavt PriNmb)

imp_NIS<- c("Beroe ovata", 
"Anguillicola crassus", 
"Mytilopsis leucophaeata", 
"Crepidula fornicata", 
"Crassostrea virginica")
# get 5 main target NIS genera in a list
library(stringr)
gnl1 <- sapply(stringr::str_split(imp_NIS," "), "[[", 1)
# get only unique genus names
df_uNEEps05.unq.gen <- df_uNEEps05[(!duplicated(df_uNEEps05$Genus)),]
# define columns to keep
clke <- c("Family",
          "Genus")
# only keep defined columns
df_uNEEps05.unq.gen <- df_uNEEps05.unq.gen[clke]

fml1 <- df_uNEEps05.unq.gen$Family[(df_uNEEps05.unq.gen$Genus %in% gnl1)]
df_uNEEps05.1 <- df_uNEEps05[(df_uNEEps05$Family %in% fml1),]


# make a table caption, using the list of abbreviations for locations
table_capt02 <- paste("Tabel 2. Artsliste over 5 udvalgte invasive arter og tætte slægtninge.")
# show the table
t.HTML04 <- df_uNEEps05.1 %>%
  htmlTable::addHtmlTableStyle(align = c(rep("l",ncol(df_uNEEps05.1)))) %>%
  htmlTable::htmlTable(caption = table_capt02, rnames = FALSE)
t.HTML04

# count species per family
spcpfam <- df_NEEps03 %>% dplyr::count(Family)
# count species per Order
spcpord <- df_NEEps03 %>% dplyr::count(Order)
spcpgen <- df_NEEps03 %>% dplyr::count(Genus)
# add back the count of species per family and order
df_us_distr02$spcperfam <- spcpfam$n[match(df_us_distr02$Family,spcpfam$Family)]
df_us_distr02$spcperord <- spcpord$n[match(df_us_distr02$Order,spcpord$Order)]
df_us_distr02$spcpergen <- spcpgen$n[match(df_us_distr02$Genus,spcpgen$Genus)]


# View(df_us_distr02)
# nrow(df_NEEps03)
# head(df_us_distr01,4)
# head(df_ds02,4)
# head(df_NEEps03,4)
# nrow(df_ds02)
# nrow(df_NEEps01)
# # limit records to the valid records
# nrow(df_ds02)
# nrow(df_NEEps03)
# get a list of related nonindeginous species
rnis <- df_NEEps03$ScientificName
rnis_table <- rnis[1:3]
rnis_table <- rnis[1:300]
#exclude 'Acipenser ruthenus' , as this species makes th code crash
rnis_table <- rnis[!grepl("Acipenser ruth",rnis)]
#rnis_table <- rnis[grepl("Acipenser ruth",rnis)]
tableID <- rnis_table
rnis_table[1]

library(rentrez)
library(ape)
#https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
#https://docs.ropensci.org/rentrez/

##  Provide ENTREZ API key
## options(ENTREZ_KEY="7ba07daba5fe26e44ea50deda6506d977a09")
set_entrez_key("7ba07daba5fe26e44ea50deda6506d977a09")
Sys.getenv("ENTREZ_KEY")

# /////////////////////////////////////////////////////////////////////////////////////////
# function get_no_of_accno_per_spc
# /////////////////////////////////////////////////////////////////////////////////////////
get_no_of_accno_per_spc <- function(tableID,
                                    Start_from=1,
                                    sleep_interval=60,
                                    timeout_limit=60, rmxval=200){
  # The "rmxval" value is used for the rentrez::entrez_search function, that searches NCBI
  # it is the maximum of sequences that will be fetched. This can be inreased of course.
  # However , since the main purpose of this search is to check if there is any sequences at all 
  # available on NCBI it should be sufficient to check first of 200 sequences are available

  require(taxizedb)
  require(ape)
  require(rentrez)
  

  
  all_spcs <- names(table(tableID))
  all_classifications <- list() # prepare list for taxize output
  lst_nac <- list()
  lst_acc <- list()
  o=length(all_spcs) # number of taxids
  
  # define max bytes
  max_bytes <- (2^30)-1
  #Start_from = 1 # change if loop needs to be restarted due to time-out
  # sleep_interval = 60 # time to sleep (sec) before new attempt after failed 
  wrong_taxid_matches <- c()
  remove_entries <- c()
  print(paste0("Step 1 of 3: processing: ", o , " species")) 
  #Get ncbi classification of each entry
  for (cl in Start_from:o) {
    #print(cl)}
    attempts<- 0
    max_attempts <- 10
    curr_spc <- all_spcs[cl]
    cat(paste0("Working on species[",cl,"/",o,"] ", curr_spc,"\n"))
    
    repeat{
      attempts <- attempts + 1
      cat(paste0("     attempt ",attempts))
      
      start_time <- Sys.time()
      fasta_rec_match <- NULL
      tryCatchResult <- tryCatch(                     # Using tryCatch() function
        
        expr = {                    # Setting up the expression
          #fetch the classification information
          seterm <- paste0(curr_spc,"[Organism]")
          se_res <- rentrez::entrez_search(db="nuccore", term=seterm, retmax=rmxval, use_history = TRUE)
          sesid <- se_res$ids
          fasta_summary <- rentrez::entrez_summary(db="nuccore",id=sesid)
          seq_lengths <- extract_from_esummary(fasta_summary,"slen")
          data_mb <- sum(seq_lengths) / (1000^2)
          cat(" found ",length(sesid)," sequences, data =", data_mb, "MB\n" )
          # get cumulative sum of lengths
          bytes_cum <- cumsum(seq_lengths)
          # get list of as many ids as possible 
          # without exceeding max_bytes
          sesid <- sesid[bytes_cum<max_bytes]
          
          fasta_rec_match <- entrez_fetch(db="nuccore", id=sesid, rettype="fasta")
          all_recs <- fasta_rec_match
          wd00 <- getwd()
          #fasta_rec_match <- taxize::classification(curr_spc, db = "ncbi")   # AGR
          temp <- tempfile(pattern= "tmpfasta",tmpdir = wd00, fileext = ".fas")
          write(all_recs, temp)
          p_recs <- ape::read.dna(temp,format="fasta")
          lprec <- labels(p_recs)
          # delete the temporary file, as there is no need to store them
          #Check its existence
          if (file.exists(temp)) {
            #Delete file if it exists
            file.remove(temp)
          }
          # exclude sequences deposited with "satellite"
          # as these are most likely microsatellite, and cannot be used for 
          # species specific detection
          lprec <- lprec[!grepl("satellite",lprec)]
          # also remove the unverified theoretical sequences, as these are fairly useless
          lprec <- lprec[!grepl("UNVERIFIED",lprec)]
          # collapse the list of sequence headers in to one single string.
          # Notice that delimiter is double semi colons
          allseqNms <- paste(unlist(lprec),collapse=";;")
          # split string and get lists nested in a list
          lbls01 <- strsplit(as.character(lprec), " ")
          # get only NCBI accession numbers
          lbls02 <- sapply(lbls01, "[[", 1)
          ""
        },
        
        error = function(e){        # Error message that should be returned
          cat(paste0(" ...error: ",e$message,"\n"))
          return(e$message)
        }) # end of tryCatch
      
      # # just to see the TryCatch result in case there is an error
      # if (!is.null(fasta_rec_match)){
      #       print("Here is the result from TrCatch")
      #       print(tryCatchResult)
      #       }
      
      #The if test to check whether the repeat function should be stopped
      # if a tax match was found, then change success to TRUE
      
      doSleep <- FALSE
      if (!is.null(fasta_rec_match)) {
        # the function was successful
        cat(" ...successful!\n")
        break
      }else if (attempts >= max_attempts) {
        # too many attempts - stop and add this taxid to the "wrong" list
        cat(paste0("tried finding fasta_rec_match ", attempts, " times, but does not exist\n"))
        break
      }else{
        # unsuccessful attempt
        time_elapsed <- Sys.time() - start_time
        if(time_elapsed<timeout_limit){
          # the error occurred within timeout limit
          if(grepl("HTTP",tryCatchResult)==TRUE){
            # HTTP error e.g. 429
            cat("HTTP error\n")
            doSleep <- TRUE
          }else{
            # it was likely a "1770542 error" so no need to retry
            cat("Skipping this species\n")
            break
          }
        }else{
          # the error was likely a timeout
          cat("Looked like a timeout\n")
          doSleep <- TRUE
        }
      }
      if(doSleep==TRUE){
        #- sleep before trying again
        cat(paste0("Trying again in ", sleep_interval, " sec.\n"))
        Sys.sleep(sleep_interval)
      }
    } # end of repeat
    #check if fasta_rec_match is wrong
    if (is.null(fasta_rec_match)) {
      wrong_taxid_matches <- c(wrong_taxid_matches,curr_spc)
    } else {
      all_classifications[length(all_classifications)+1] <- fasta_rec_match
      #lst_nac[length(lst_nac)+1] <- length(lbls02)
      lst_nac[curr_spc] <- length(lbls02)
      lst_acc[curr_spc] <- allseqNms
    }        
  } # end of loop (cl in Start_from:o)
  # make the list a data frame
  df_ac01 <- as.data.frame(do.call(rbind,lst_nac))
  df_acNm <- as.data.frame(do.call(rbind,lst_acc))
  # make the row names a column
  df_ac01$genus_species <- row.names(df_ac01)
  df_acNm$genus_species <- row.names(df_acNm)
  #change the column names
  colnames(df_ac01) <- c("no_of_accn","genus_species")
  colnames(df_acNm) <- c("seq_name_for_accn","genus_species")
  #match between data frame to one data frame that has sequence counts and the sequence headers
  df_acNm$no_of_accn <- df_ac01$no_of_accn[match(df_ac01$genus_species,df_acNm$genus_species)]
  return(df_acNm)  
  #delete tmp files afterwards
  lst_of_fls <- list.files(wd00, full.names = TRUE)
  lst_of_fasfls <- lst_of_fls[grepl(".fas",lst_of_fls)]
  lst_of_tmpfasfls <- lst_of_fasfls[grepl("tmp",lst_of_fasfls)]
  do.call(file.remove, list(lst_of_tmpfasfls))
}
# /////////////////////////////////////////////////////////////////////////////////////////
# use the function to check how many useful sequences can be found on NCBI for each
# species -  this might take a while to run for all species
df_ac <- get_no_of_accno_per_spc(rnis_table,Start_from=1,sleep_interval=10,timeout_limit=60, rmxval=300)
# add back the number of sequences available on NCBI back to the data frame that lists all the
# closely related species that  potentially could give rise to false positive detection,
# due to similarity in sequences
df_NEEps03$No_of_seq_onNBCI <- df_ac$no_of_accn[match(df_NEEps03$ScientificName,df_ac$genus_species)]
# # match back the NIS target species by comparing Orders
# df_NEEps03$NIS_target_spc <- df_us_distr02$ScientificName[match(df_NEEps03$Order,df_us_distr02$Order)]
# # add back target species for the non target species that were missing an order level
# df_NEEps03$NIS_target_spc[is.na(df_NEEps03$NIS_target_spc)] <- df_us_distr02$ScientificName[match(df_NEEps03$Class[is.na(df_NEEps03$NIS_target_spc)],df_us_distr02$Class)]
# # also get the NIS target species family 
# df_NEEps03$NIS_target_spc_Family <- df_us_distr02$Family[match(df_NEEps03$NIS_target_spc,df_us_distr02$ScientificName)]
# df_NEEps03$NIS_target_spc_Order <- df_us_distr02$Order[match(df_NEEps03$NIS_target_spc,df_us_distr02$ScientificName)]
# # and match back the number of non target species per target species
df_NEEps03$NIS_target_spc_perfam <- df_us_distr02$spcperfam[match(df_NEEps03$NIS_target_spc,df_us_distr02$ScientificName)]
df_NEEps03$NIS_target_spc_perord <- df_us_distr02$spcperord[match(df_NEEps03$NIS_target_spc,df_us_distr02$ScientificName)]
df_NEEps03$NIS_target_spc_pergen <- df_us_distr02$spcpergen[match(df_NEEps03$NIS_target_spc,df_us_distr02$ScientificName)]
# 
# # subset data frame 
# df_NEEps04 <- df_NEEps03[df_NEEps03$NIS_target_spc_perfam<10,]
# # count NAs within group. see: https://stackoverflow.com/questions/24477748/r-count-na-by-group
# df_NEEps05 <- aggregate(No_of_seq_onNBCI ~ NIS_target_spc, data=df_NEEps03, function(x) {sum(is.na(x))}, na.action = NULL)
# # change column names
# colnames(df_NEEps05) <- c("NIS_target_spc","No_of_NTspc_miss_seq_onNBCI")
# # add back number of non target species per NIS within NIS TS family
# df_NEEps05$NTspc_perfam <- df_NEEps03$NIS_target_spc_perfam[match(df_NEEps05$NIS_target_spc,df_NEEps03$NIS_target_spc)]
# df_NEEps05$NTspc_perord <- df_NEEps03$NIS_target_spc_perord[match(df_NEEps05$NIS_target_spc,df_NEEps03$NIS_target_spc)]
# df_NEEps05$NTspc_pergen <- df_NEEps03$NIS_target_spc_pergen[match(df_NEEps05$NIS_target_spc,df_NEEps03$NIS_target_spc)]
# 
# 
#
# add back number of available NCBI sequences for the NIS target species
tibl_inx01$No_of_seq_onNBCI_for_NIST <- df_NEEps03$No_of_seq_onNBCI[match(tibl_inx01$NIS_SpeciesNm,df_NEEps03$ScientificName)]
# add back taxonomic categories
tibl_inx01$family_for_NIST <- df_NEEps03$Family[match(tibl_inx01$NIS_SpeciesNm,df_NEEps03$ScientificName)]
tibl_inx01$order_for_NIST <- df_NEEps03$Order[match(tibl_inx01$NIS_SpeciesNm,df_NEEps03$ScientificName)]
tibl_inx01$class_for_NIST <- df_NEEps03$Class[match(tibl_inx01$NIS_SpeciesNm,df_NEEps03$ScientificName)]
tibl_inx01$phylum_for_NIST <- df_NEEps03$Phylum[match(tibl_inx01$NIS_SpeciesNm,df_NEEps03$ScientificName)]
# match back the NIS target species by comparing Orders
df_NEEps03$NIS_target_spc <- tibl_inx01$NIS_SpeciesNm[match(df_NEEps03$Order,tibl_inx01$order_for_NIST)]
# count NAs within group. see: https://stackoverflow.com/questions/24477748/r-count-na-by-group
df_NEEps05 <- aggregate(No_of_seq_onNBCI ~ NIS_target_spc, data=df_NEEps03, function(x) {sum(is.na(x))}, na.action = NULL)
# change column names
colnames(df_NEEps05) <- c("NIS_target_spc","No_of_NTspc_miss_seq_onNBCI")
# split string and get lists nested in a list
gl05 <- strsplit(as.character(df_NEEps05$NIS_target_spc), " ")
#get first and second element of nested list
gnl05 <- sapply(gl05, "[[", 1)
# add back genus name to data frame
df_NEEps05$NIS_target_gen <- gnl05 
# match back family name for target NIS
df_NEEps05$NIS_target_fam<- df_NEEps03$Family[match(df_NEEps05$NIS_target_gen,df_NEEps03$Genus)]
df_NEEps05$NIS_target_ord<- df_NEEps03$Order[match(df_NEEps05$NIS_target_gen,df_NEEps03$Genus)]
df_NEEps05$NIS_target_cla<- df_NEEps03$Class[match(df_NEEps05$NIS_target_gen,df_NEEps03$Genus)]
# Match back number of non target species species missing sequences from NCBI
tibl_inx01$No_of_NTspc_miss_seq_onNBCI <- df_NEEps05$No_of_NTspc_miss_seq_onNBCI[match(tibl_inx01$NIS_SpeciesNm,df_NEEps05$NIS_target_spc)]
tibl_inx01$No_of_NTspc_miss_seq_onNBCI[is.na(tibl_inx01$No_of_NTspc_miss_seq_onNBCI)] <- df_NEEps05$No_of_NTspc_miss_seq_onNBCI[match(tibl_inx01$family_for_NIST[is.na(tibl_inx01$No_of_NTspc_miss_seq_onNBCI)],df_NEEps05$NIS_target_fam)]
tibl_inx01$No_of_NTspc_miss_seq_onNBCI[is.na(tibl_inx01$No_of_NTspc_miss_seq_onNBCI)] <- df_NEEps05$No_of_NTspc_miss_seq_onNBCI[match(tibl_inx01$order_for_NIST[is.na(tibl_inx01$No_of_NTspc_miss_seq_onNBCI)],df_NEEps05$NIS_target_ord)]
tibl_inx01$No_of_NTspc_miss_seq_onNBCI[is.na(tibl_inx01$No_of_NTspc_miss_seq_onNBCI)] <- df_NEEps05$No_of_NTspc_miss_seq_onNBCI[match(tibl_inx01$class_for_NIST[is.na(tibl_inx01$No_of_NTspc_miss_seq_onNBCI)],df_NEEps05$NIS_target_cla)]

# # add back number of non target species per NIS within NIS TS family
df_NEEps05$NTspc_perfam <- df_NEEps03$NIS_target_spc_perfam[match(df_NEEps05$NIS_target_spc,df_NEEps03$NIS_target_spc)]
df_NEEps05$NTspc_perord <- df_NEEps03$NIS_target_spc_perord[match(df_NEEps05$NIS_target_spc,df_NEEps03$NIS_target_spc)]
df_NEEps05$NTspc_pergen <- df_NEEps03$NIS_target_spc_pergen[match(df_NEEps05$NIS_target_spc,df_NEEps03$NIS_target_spc)]
# # add back number of non target species per NIS within NIS TS family
tibl_inx01$NTspc_perfam <- df_NEEps05$NTspc_perfam[match(tibl_inx01$NIS_SpeciesNm,df_NEEps05$NIS_target_spc)]
tibl_inx01$NTspc_perfam[is.na(tibl_inx01$NTspc_perfam)] <- df_NEEps05$NTspc_perfam[match(tibl_inx01$NIS_GenusNm[is.na(tibl_inx01$NTspc_perfam)],df_NEEps05$NIS_target_gen)]
tibl_inx01$NTspc_perfam[is.na(tibl_inx01$NTspc_perfam)]  <- df_NEEps05$NTspc_perfam[match(tibl_inx01$family_for_NIST[is.na(tibl_inx01$NTspc_perfam)],df_NEEps05$NIS_target_fam)]
tibl_inx01$NTspc_perfam[is.na(tibl_inx01$NTspc_perfam)]  <- df_NEEps05$NTspc_perfam[match(tibl_inx01$order_for_NIST[is.na(tibl_inx01$NTspc_perfam)],df_NEEps05$NIS_target_ord)]
# # add back number of non target species per NIS within NIS TS order
tibl_inx01$NTspc_perord <- df_NEEps05$NTspc_perord[match(tibl_inx01$NIS_SpeciesNm,df_NEEps05$NIS_target_spc)]
tibl_inx01$NTspc_perord[is.na(tibl_inx01$NTspc_perord)] <- df_NEEps05$NTspc_perord[match(tibl_inx01$NIS_GenusNm[is.na(tibl_inx01$NTspc_perord)],df_NEEps05$NIS_target_gen)]
tibl_inx01$NTspc_perord[is.na(tibl_inx01$NTspc_perord)]  <- df_NEEps05$NTspc_perord[match(tibl_inx01$family_for_NIST[is.na(tibl_inx01$NTspc_perord)],df_NEEps05$NIS_target_fam)]
tibl_inx01$NTspc_perord[is.na(tibl_inx01$NTspc_perord)]  <- df_NEEps05$NTspc_perord[match(tibl_inx01$order_for_NIST[is.na(tibl_inx01$NTspc_perord)],df_NEEps05$NIS_target_ord)]
# # add back number of non target species per NIS within NIS TS genus
tibl_inx01$NTspc_pergen <- df_NEEps05$NTspc_pergen[match(tibl_inx01$NIS_SpeciesNm,df_NEEps05$NIS_target_spc)]
tibl_inx01$NTspc_pergen[is.na(tibl_inx01$NTspc_pergen)] <- df_NEEps05$NTspc_pergen[match(tibl_inx01$NIS_GenusNm[is.na(tibl_inx01$NTspc_pergen)],df_NEEps05$NIS_target_gen)]
tibl_inx01$NTspc_pergen[is.na(tibl_inx01$NTspc_pergen)]  <- df_NEEps05$NTspc_pergen[match(tibl_inx01$family_for_NIST[is.na(tibl_inx01$NTspc_pergen)],df_NEEps05$NIS_target_fam)]
tibl_inx01$NTspc_pergen[is.na(tibl_inx01$NTspc_pergen)]  <- df_NEEps05$NTspc_pergen[match(tibl_inx01$order_for_NIST[is.na(tibl_inx01$NTspc_pergen)],df_NEEps05$NIS_target_ord)]
# sum up columns to get a priority number
tibl_inx01$prioritNmb <- rowSums(tibl_inx01[ , c("No_of_NTspc_miss_seq_onNBCI","NTspc_perfam","NTspc_perord","NTspc_pergen")], na.rm=TRUE)

# re order the data frame by th summed up column - a low value
# is a NIS target species to prioritize
tibl_inx01 <- tibl_inx01[order(tibl_inx01$prioritNmb),]

capt_tbl02 <- " Tabel 1. Prioriteret list over NIS arter som der kan udvikles specifikke eDNA systemer imod."
# show the table
t.HTML05 <- tibl_inx01 %>%
  htmlTable::addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  #htmlTable::addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable::htmlTable(caption = capt_tbl02)
t.HTML05
# write the resulting table as an excel file
require("xlsx")
write.xlsx(tibl_inx01, "Artsprioritering_v02.xlsx")
#View(df_NEEps05)

