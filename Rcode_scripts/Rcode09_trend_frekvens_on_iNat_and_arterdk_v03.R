#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
# R-code provided for the project:
# “MONIS6”
library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
library(purrr)


# place the working directory in a variable
wd00 <- getwd()
# define directory with inputfile in
wd07 <- "output07_map_eDNA_detections_with_GBIF_and_iNat_for_MST2017_2023_samples"
file <- "Table06_iNat_arterdk_and_MONIS6_records_2017_2023.csv"
wd00_wd07_file <- paste(wd00,wd07,file,sep="/")
# read in the csv-file
df_obs <- read.table(wd00_wd07_file, sep=";", header=T)
# filter to only comprise the 'reccat' that are "iNaturalist"
df_obs <- df_obs %>% dplyr::filter((reccat=="iNaturalist")) 

# exclude if there is no Latin species name 
df_obs <- df_obs  %>%
  dplyr::rename(Latspecies=Lat_Species)   %>%
  dplyr::filter(!is.na(Latspecies)) %>%
  dplyr::filter(!(Latspecies==0)) 

# get the sampling year and sampling month, and make the values numeric
df_obs <- df_obs %>%
  dplyr::rename(year=yer) %>%
  dplyr::rename(month=mnt) %>%
  dplyr::mutate(year=as.numeric(year)) %>%
  dplyr::mutate(month=as.numeric(month))
# split the string by delimiter
ssnmnth  <- strsplit(as.character(df_obs$yer_ssn.per), ", ")
# and use the splitted string to get each element and add back in a column
df_obs$ssnmnth <- sapply(ssnmnth, "[[", 2)
# make it a factor
df_obs <- df_obs %>%
  dplyr::mutate(season=as.factor(ssnmnth))

# substitute to remove the space in the seasons 
df_obs$season <- gsub(" ","",df_obs$season)
# make a column with a variable that can represent the sampling location
df_obs$loc.pos <- paste(df_obs$declon,df_obs$declat,sep="_")
df_obs$stn <-  df_obs$loc.pos
# count up the number of attempts to search for the Latspecies, per year, 
# per season
df_srchfor <- df_obs %>%
  dplyr::group_by(Latspecies, year, season) %>%
  dplyr::summarise(n_srchfor=n(), .groups="drop")
  
#  dplyr::rename(stn=loc.pos)
# count up the number of actual detections for the Latspecies, per year, 
# per season, per sampling location
df_finds <- df_obs %>%
  #dplyr::filter(!(orgFnd2==0)) %>%
  dplyr::group_by(Latspecies, year, season) %>%
  dplyr::summarise(n_fnd=n(), .groups="drop")  %>%
  dplyr::mutate(locat.srch=1)

#
df <- dplyr::distinct(df_obs, Latspecies,stn) %>%
  merge(dplyr::distinct(df_obs, stn), all=T)  %>%
  merge(dplyr::distinct(df_obs, year), all=T)  %>%
  merge(dplyr::distinct(df_obs, season), all=T) 

#______
df_finds$locat.srch==df_finds$n_fnd

#use the data frame with counts of attempts to search for the organism
# and join it with the data frame with the counts of the detections of the 
# organism, use 'dplyr::mutate' to replace any eventual NAs with '0'
df_srchfnd <- df_srchfor %>% 
  left_join(df_finds, by=c("Latspecies", "year", "season")) %>%  
  dplyr::mutate(n_fnd=ifelse(is.na(n_fnd),0,n_fnd))
#
dff <- df_srchfnd %>%
  dplyr::group_by(Latspecies, year,season) %>%
  #dplyr::summarise(n=n(), f=sum(detected)/n(), .groups="drop") %>%
  dplyr::mutate(Latspecies_season=paste0(Latspecies,"_",season))
# make a frequency element 'f', that is based on the number of detections
# per Latspecies per year per season, where the 'f' is the number of detections
# per attempts made to find the organism
dff$f <- dff$n_fnd/dff$n_srchfor
dff$f <- dff$n_fnd
#unique(dff$f)

# if there are less than 3 observations of the species then, it is not possible
# to do a linear model , so start by counting the number of observations per #
# species per season
df_cntprseas <- dff %>% 
  dplyr::group_by(Latspecies_season) %>% 
  dplyr::summarise(n_perseason = n())
# join the main data frame with the data frame that counts up
# detections of species per season, and use filter to exclude
# when there is less than 3 observations
dff <- dff %>% 
  left_join(df_cntprseas, by=c("Latspecies_season")) %>%
  dplyr::filter((n_perseason>=3)) 

# since the jul-nov season is half a season later, the year can be added 0.5 
dff <- dff %>% dplyr::mutate(year=ifelse(season=="jan-jun",year,year+0.5))

dff_split <- dff %>%
  split(.$Latspecies_season)

lm_list <- purrr::map(dff_split, lm, formula= f ~ year)

# lm_list[1] %>%
#   purrr::map(summary.lm) %>%
#   purrr::map(c("coefficients")) %>%
#   purrr::map_dbl(8)  %>%

pvalues <- lm_list %>% 
  purrr::map(summary.lm) %>% 
  purrr::map(c("coefficients")) %>% 
  purrr::map_dbl(8)  %>% # 8th element is the p-value 
  broom::tidy() %>% 
  dplyr::arrange(desc(x)) %>% 
  dplyr::rename(Latspecies_season=names, p_val = x)

slopes <- lm_list %>% 
  purrr::map(summary.lm) %>% 
  purrr::map(c("coefficients")) %>% 
  map_dbl(2)  %>% # 8th element is the p-value 
  broom::tidy() %>% 
  dplyr::arrange(desc(x)) %>% 
  dplyr::rename(Latspecies_season=names, slope = x)


get_fitted_vals <- function(lm, conf_int=0.9){
  df <- lm %>% predict.lm(
    interval="confidence", 
    level=conf_int) %>%
    as.data.frame() %>%
    rename(pred=fit)
  yr <- lm$model %>%
    select(year)
  df <- yr %>%
    bind_cols(df) %>%
    ungroup()
  return(df)
}

df_fit <- lm_list %>% purrr::map(get_fitted_vals) %>%
  bind_rows(.id="Latspecies_season")

df_fit <- df_fit %>%
  dplyr::mutate(Latspecies = stringr::str_split_i(Latspecies_season,"_",1),
                season = stringr::str_split_i(Latspecies_season,"_",2))

df_fit <- df_fit %>%
  left_join(pvalues, by="Latspecies_season") 

df_fit <- df_fit %>%
  dplyr::mutate(pred=ifelse(is.na(p_val),NA,pred)) %>%
  dplyr::mutate(style=ifelse(p_val<=0.1, "1", "2"))

df_fit <- df_fit %>%
  left_join(slopes, by="Latspecies_season") %>% 
  dplyr::mutate(year=ifelse(season=="jan-jun",year,year+0.5))


p_sig <- 0.1

df_note <-  df_fit %>%
  dplyr::mutate(slope=ifelse(is.na(p_val),NA, slope)) %>%
  dplyr::distinct(Latspecies, season, slope, p_val) %>%
  dplyr::mutate(p_text=ifelse(is.na(slope),NA_character_, ifelse(p_val<0.001,
                                                                 "p<0.001",
                                                                 ifelse(p_val> p_sig, "n.s.",
                                                                        paste0("p=",round(p_val,3)))))) %>%
  dplyr::mutate(note=ifelse(is.na(slope),NA_character_,
                            paste0(ifelse(slope>0,"+",""),
                                   round(100*slope,1),"%/år (", p_text, ")"))) %>%
  dplyr::group_by(Latspecies) %>%
  dplyr::mutate(n=sum(!is.na(slope))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(offset=ifelse(n>1,15,0)) %>%
  dplyr::mutate(y=max(dff$f,na.rm=T), year=min(dff$year,na.rm=T)) %>%
  dplyr::mutate(y=ifelse(season=="jan-jun",y,y-offset))

#View(df_note)

#dff <- dff %>% dplyr::mutate(year=ifelse(season=="jan-jun",year,year+0.5))
mn.yer <- round(min(dff$year),0)
mx.yer <- ceiling(max(dff$year))

pal_season <- c("#0000FF", "#FF0000")
pal_season <- c("blueviolet", "darkorange4")
pal_season <- c("slateblue4", "goldenrod4")

#View(dff)

p <- ggplot(dff) +
  geom_ribbon(data=df_fit,
              aes(x=year, ymin=lwr, ymax=upr, 
                  group=season, fill=season),
              colour=NA, alpha=0.1) +
  geom_line(data=df_fit, aes(x=year, y=pred, group=season, colour=season, 
                             linetype = style, alpha=style)) +
  geom_point(aes(x=year, y=f, group=season, colour=season)) +
  geom_text(data=df_note, aes(x=year, y=y, group=season, 
                              colour=season, label=note),
            show.legend = F,
            hjust=0, size=3) +
  facet_wrap(.~Latspecies, ncol=4,
             labeller = label_bquote(col = italic(.(Latspecies))) ) +
  # make the x axis have breaks that are represented every 1 increment      
  scale_x_continuous(breaks=seq(mn.yer, mx.yer, 1)) +
  
  scale_colour_manual(values=pal_season, name="sæson") +
  scale_alpha_manual(values=c(0.8, 0.6), guide="none") +
  scale_fill_manual(values=pal_season, name="sæson", guide="none") +
  scale_linetype_manual(values=c("solid", "longdash"), guide="none") +
  xlab("år") +
  ylab("antal af fund") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        legend.position="bottom",
        strip.text = element_text(hjust = 0))
# see the plot
p
#ggsave(p, file="cjm/trend_frekvens.png", height=20, width=20, units="cm", dpi=300, bg="white")
# define output directory
wd08 <- "output08_stckbarplot_w_GBIF_and_iNat_for_MST2017_2023_records"
outflNm <- "Fig20_trend_frekvens_iNat.png"
wd00wd08_outflNm <- paste(wd00,wd08,outflNm, sep="/")
ggsave(p, 
       file=wd00wd08_outflNm, 
       height=20, width=20, units="cm", dpi=300, bg="white")

#


#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
#____________________________________________________________________________#
# place the working directory in a variable
wd00 <- getwd()
# define directory with inputfile in
wd07 <- "output07_map_eDNA_detections_with_GBIF_and_iNat_for_MST2017_2023_samples"
file <- "Table06_iNat_arterdk_and_MONIS6_records_2017_2023.csv"
wd00_wd07_file <- paste(wd00,wd07,file,sep="/")
# read in the csv-file
df_obs <- read.table(wd00_wd07_file, sep=";", header=T)
#unique(df_obs$reccat)
# filter to only comprise the 'reccat' that are "iNaturalist"
#df_obs <- df_obs %>% dplyr::filter((reccat=="iNaturalist")) 
df_obs <- df_obs %>% dplyr::filter((reccat=="www.arter.dk")) 
# exclude if there is no Latin species name 
df_obs <- df_obs  %>%
  dplyr::rename(Latspecies=Lat_Species)   %>%
  dplyr::filter(!is.na(Latspecies)) %>%
  dplyr::filter(!(Latspecies==0)) 

# get the sampling year and sampling month, and make the values numeric
df_obs <- df_obs %>%
  dplyr::rename(year=yer) %>%
  dplyr::rename(month=mnt) %>%
  dplyr::mutate(year=as.numeric(year)) %>%
  dplyr::mutate(month=as.numeric(month))
# split the string by delimiter
ssnmnth  <- strsplit(as.character(df_obs$yer_ssn.per), ", ")
# and use the splitted string to get each element and add back in a column
df_obs$ssnmnth <- sapply(ssnmnth, "[[", 2)
# make it a factor
df_obs <- df_obs %>%
  dplyr::mutate(season=as.factor(ssnmnth))

# substitute to remove the space in the seasons 
df_obs$season <- gsub(" ","",df_obs$season)
# make a column with a variable that can represent the sampling location
df_obs$loc.pos <- paste(df_obs$declon,df_obs$declat,sep="_")
df_obs$stn <-  df_obs$loc.pos
# count up the number of attempts to search for the Latspecies, per year, 
# per season
df_srchfor <- df_obs %>%
  dplyr::group_by(Latspecies, year, season) %>%
  dplyr::summarise(n_srchfor=n(), .groups="drop")

#  dplyr::rename(stn=loc.pos)
# count up the number of actual detections for the Latspecies, per year, 
# per season, per sampling location
df_finds <- df_obs %>%
  #dplyr::filter(!(orgFnd2==0)) %>%
  dplyr::group_by(Latspecies, year, season) %>%
  dplyr::summarise(n_fnd=n(), .groups="drop")  %>%
  dplyr::mutate(locat.srch=1)

#
df <- dplyr::distinct(df_obs, Latspecies,stn) %>%
  merge(dplyr::distinct(df_obs, stn), all=T)  %>%
  merge(dplyr::distinct(df_obs, year), all=T)  %>%
  merge(dplyr::distinct(df_obs, season), all=T) 

#______
df_finds$locat.srch==df_finds$n_fnd

#use the data frame with counts of attempts to search for the organism
# and join it with the data frame with the counts of the detections of the 
# organism, use 'dplyr::mutate' to replace any eventual NAs with '0'
df_srchfnd <- df_srchfor %>% 
  left_join(df_finds, by=c("Latspecies", "year", "season")) %>%  
  dplyr::mutate(n_fnd=ifelse(is.na(n_fnd),0,n_fnd))
#
dff <- df_srchfnd %>%
  dplyr::group_by(Latspecies, year,season) %>%
  #dplyr::summarise(n=n(), f=sum(detected)/n(), .groups="drop") %>%
  dplyr::mutate(Latspecies_season=paste0(Latspecies,"_",season))
# make a frequency element 'f', that is based on the number of detections
# per Latspecies per year per season, where the 'f' is the number of detections
# per attempts made to find the organism
dff$f <- dff$n_fnd/dff$n_srchfor
dff$f <- dff$n_fnd
#unique(dff$f)

# if there are less than 3 observations of the species then, it is not possible
# to do a linear model , so start by counting the number of observations per #
# species per season
df_cntprseas <- dff %>% 
  dplyr::group_by(Latspecies_season) %>% 
  dplyr::summarise(n_perseason = n())
# join the main data frame with the data frame that counts up
# detections of species per season, and use filter to exclude
# when there is less than 3 observations
dff <- dff %>% 
  left_join(df_cntprseas, by=c("Latspecies_season")) %>%
  dplyr::filter((n_perseason>=3)) 

# since the jul-nov season is half a season later, the year can be added 0.5 
dff <- dff %>% dplyr::mutate(year=ifelse(season=="jan-jun",year,year+0.5))

dff_split <- dff %>%
  split(.$Latspecies_season)

lm_list <- purrr::map(dff_split, lm, formula= f ~ year)

# lm_list[1] %>%
#   purrr::map(summary.lm) %>%
#   purrr::map(c("coefficients")) %>%
#   purrr::map_dbl(8)  %>%

pvalues <- lm_list %>% 
  purrr::map(summary.lm) %>% 
  purrr::map(c("coefficients")) %>% 
  purrr::map_dbl(8)  %>% # 8th element is the p-value 
  broom::tidy() %>% 
  dplyr::arrange(desc(x)) %>% 
  dplyr::rename(Latspecies_season=names, p_val = x)

slopes <- lm_list %>% 
  purrr::map(summary.lm) %>% 
  purrr::map(c("coefficients")) %>% 
  map_dbl(2)  %>% # 8th element is the p-value 
  broom::tidy() %>% 
  dplyr::arrange(desc(x)) %>% 
  dplyr::rename(Latspecies_season=names, slope = x)


get_fitted_vals <- function(lm, conf_int=0.9){
  df <- lm %>% predict.lm(
    interval="confidence", 
    level=conf_int) %>%
    as.data.frame() %>%
    dplyr::rename(pred=fit)
  yr <- lm$model %>%
    dplyr::select(year)
  df <- yr %>%
    bind_cols(df) %>%
    dplyr::ungroup()
  return(df)
}

df_fit <- lm_list %>% purrr::map(get_fitted_vals) %>%
  bind_rows(.id="Latspecies_season")

df_fit <- df_fit %>%
  dplyr::mutate(Latspecies = stringr::str_split_i(Latspecies_season,"_",1),
                season = stringr::str_split_i(Latspecies_season,"_",2))

df_fit <- df_fit %>%
  left_join(pvalues, by="Latspecies_season") 

df_fit <- df_fit %>%
  dplyr::mutate(pred=ifelse(is.na(p_val),NA,pred)) %>%
  dplyr::mutate(style=ifelse(p_val<=0.1, "1", "2"))

df_fit <- df_fit %>%
  left_join(slopes, by="Latspecies_season") %>% 
  dplyr::mutate(year=ifelse(season=="jan-jun",year,year+0.5))


p_sig <- 0.1

df_note <-  df_fit %>%
  dplyr::mutate(slope=ifelse(is.na(p_val),NA, slope)) %>%
  dplyr::distinct(Latspecies, season, slope, p_val) %>%
  dplyr::mutate(p_text=ifelse(is.na(slope),NA_character_, ifelse(p_val<0.001,
                                                                 "p<0.001",
                                                                 ifelse(p_val> p_sig, "n.s.",
                                                                        paste0("p=",round(p_val,3)))))) %>%
  dplyr::mutate(note=ifelse(is.na(slope),NA_character_,
                            paste0(ifelse(slope>0,"+",""),
                                   round(100*slope,1),"%/år (", p_text, ")"))) %>%
  dplyr::group_by(Latspecies) %>%
  dplyr::mutate(n=sum(!is.na(slope))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(offset=ifelse(n>1,5,0)) %>%
  dplyr::mutate(y=max(dff$f,na.rm=T), year=min(dff$year,na.rm=T)) %>%
  dplyr::mutate(y=ifelse(season=="jan-jun",y,y-offset))

#View(df_note)

#dff <- dff %>% dplyr::mutate(year=ifelse(season=="jan-jun",year,year+0.5))
mn.yer <- round(min(dff$year),0)
mx.yer <- ceiling(max(dff$year))

pal_season <- c("#0000FF", "#FF0000")
pal_season <- c("blueviolet", "darkorange4")
#View(dff)

p <- ggplot(dff) +
  geom_ribbon(data=df_fit,
              aes(x=year, ymin=lwr, ymax=upr, 
                  group=season, fill=season),
              colour=NA, alpha=0.1) +
  geom_line(data=df_fit, aes(x=year, y=pred, group=season, colour=season, 
                             linetype = style, alpha=style)) +
  geom_point(aes(x=year, y=f, group=season, colour=season)) +
  geom_text(data=df_note, aes(x=year, y=y, group=season, 
                              colour=season, label=note),
            show.legend = F,
            hjust=0, size=3) +
  facet_wrap(.~Latspecies, ncol=4,
             labeller = label_bquote(col = italic(.(Latspecies))) ) +
  # make the x axis have breaks that are represented every 1 increment      
  scale_x_continuous(breaks=seq(mn.yer, mx.yer, 1)) +
  
  scale_colour_manual(values=pal_season, name="sæson") +
  scale_alpha_manual(values=c(0.8, 0.6), guide="none") +
  scale_fill_manual(values=pal_season, name="sæson", guide="none") +
  scale_linetype_manual(values=c("solid", "longdash"), guide="none") +
  xlab("år") +
  ylab("antal af fund") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        legend.position="bottom",
        strip.text = element_text(hjust = 0))
# see the plot
p
#ggsave(p, file="cjm/trend_frekvens.png", height=20, width=20, units="cm", dpi=300, bg="white")
# define output directory
wd08 <- "output08_stckbarplot_w_GBIF_and_iNat_for_MST2017_2023_records"
outflNm <- "Fig21_trend_frekvens_arterdk.png"
wd00wd08_outflNm <- paste(wd00,wd08,outflNm, sep="/")
ggsave(p, 
       file=wd00wd08_outflNm, 
       height=20, width=20, units="cm", dpi=300, bg="white")

#
#____________________________________________________________________________#
#____________________________________________________________________________#

# read in the csv-file
df_obs <- read.table(wd00_wd07_file, sep=";", header=T)
# ma a vector with categories to filter in dplyr
cat.trg <- c("www.arter.dk","iNaturalist")
# Use dplyr::filter to subset the data frame , to only
# comprise 'reccat' that matches the elements in the 'cat.trg'
df_obs <- df_obs %>% dplyr::filter((reccat %in% cat.trg)) 
#unique(df_obs$reccat)

# exclude if there is no Latin species name 
df_obs <- df_obs  %>%
  dplyr::rename(Latspecies=Lat_Species)   %>%
  dplyr::filter(!is.na(Latspecies)) %>%
  dplyr::filter(!(Latspecies==0)) 

# get the sampling year and sampling month, and make the values numeric
df_obs <- df_obs %>%
  dplyr::rename(year=yer) %>%
  dplyr::rename(month=mnt) %>%
  dplyr::mutate(year=as.numeric(year)) %>%
  dplyr::mutate(month=as.numeric(month))
# split the string by delimiter
ssnmnth  <- strsplit(as.character(df_obs$yer_ssn.per), ", ")
# and use the splitted string to get each element and add back in a column
df_obs$ssnmnth <- sapply(ssnmnth, "[[", 2)
# make it a factor
df_obs <- df_obs %>%
  dplyr::mutate(season=as.factor(ssnmnth))

# substitute to remove the space in the seasons 
df_obs$season <- gsub(" ","",df_obs$season)
# make a column with a variable that can represent the sampling location
df_obs$loc.pos <- paste(df_obs$declon,df_obs$declat,sep="_")
df_obs$stn <-  df_obs$loc.pos
# count up the number of attempts to search for the Latspecies, per year, 
# per season, per recording category
df_srchfor <- df_obs %>%
  dplyr::group_by(Latspecies, year, season, reccat) %>%
  dplyr::summarise(n_srchfor=n(), .groups="drop")

# count up the number of actual detections for the Latspecies, per year, 
# per season, per recording category
df_finds <- df_obs %>%
  #dplyr::filter(!(orgFnd2==0)) %>%
  dplyr::group_by(Latspecies, year, season, reccat) %>%
  dplyr::summarise(n_fnd=n(), .groups="drop")  %>%
  dplyr::mutate(locat.srch=1)

#
df <- dplyr::distinct(df_obs, Latspecies,stn,reccat) %>%
  merge(dplyr::distinct(df_obs, stn,reccat), all=T)  %>%
  merge(dplyr::distinct(df_obs, year,reccat), all=T)  %>%
  merge(dplyr::distinct(df_obs, season,reccat), all=T) 

#______
df_finds$locat.srch==df_finds$n_fnd

#use the data frame with counts of attempts to search for the organism
# and join it with the data frame with the counts of the detections of the 
# organism, use 'dplyr::mutate' to replace any eventual NAs with '0'
df_srchfnd <- df_srchfor %>% 
  left_join(df_finds, by=c("Latspecies", "year", "season", "reccat")) %>%  
  dplyr::mutate(n_fnd=ifelse(is.na(n_fnd),0,n_fnd))
#
dff <- df_srchfnd %>%
  dplyr::group_by(Latspecies, year,season, reccat) %>%
  #dplyr::summarise(n=n(), f=sum(detected)/n(), .groups="drop") %>%
  #dplyr::mutate(Latspecies_season=paste0(Latspecies,"_",season))
  dplyr::mutate(Latspecies_season_reccat=paste0(Latspecies,"_",season,"_",reccat))
# make a frequency element 'f', that is based on the number of detections
# per Latspecies per year per season, where the 'f' is the number of detections
# per attempts made to find the organism
dff$f <- dff$n_fnd/dff$n_srchfor
dff$f <- dff$n_fnd
#unique(dff$f)

# if there are less than 3 observations of the species then, it is not possible
# to do a linear model , so start by counting the number of observations per #
# species per season
df_cntprseas <- dff %>% 
  dplyr::group_by(Latspecies_season_reccat) %>% 
  dplyr::summarise(n_perseason = n())
# join the main data frame with the data frame that counts up
# detections of species per season, and use filter to exclude
# when there is less than 3 observations
dff <- dff %>% 
  left_join(df_cntprseas, by=c("Latspecies_season_reccat")) %>%
  dplyr::filter((n_perseason>=3)) 

# since the jul-nov season is half a season later, the year can be added 0.5 
dff <- dff %>% dplyr::mutate(year=ifelse(season=="jan-jun",year,year+0.5))
# paste together season and record category
dff$season_reccat <- paste0(dff$season,"_",dff$reccat)

dff_split <- dff %>%
  split(.$Latspecies_season_reccat,3)

lm_list <- purrr::map(dff_split, lm, formula= f ~ year)

# lm_list[1] %>%
#   purrr::map(summary.lm) %>%
#   purrr::map(c("coefficients")) %>%
#   purrr::map_dbl(8)  %>%

pvalues <- lm_list %>% 
  purrr::map(summary.lm) %>% 
  purrr::map(c("coefficients")) %>% 
  purrr::map_dbl(8)  %>% # 8th element is the p-value 
  broom::tidy() %>% 
  dplyr::arrange(desc(x)) %>% 
  dplyr::rename(Latspecies_season_reccat=names, p_val = x)

slopes <- lm_list %>% 
  purrr::map(summary.lm) %>% 
  purrr::map(c("coefficients")) %>% 
  map_dbl(2)  %>% # 8th element is the p-value 
  broom::tidy() %>% 
  dplyr::arrange(desc(x)) %>% 
  dplyr::rename(Latspecies_season_reccat=names, slope = x)


get_fitted_vals <- function(lm, conf_int=0.9){
  df <- lm %>% predict.lm(
    interval="confidence", 
    level=conf_int) %>%
    as.data.frame() %>%
    dplyr::rename(pred=fit)
  yr <- lm$model %>%
    dplyr::select(year)
  df <- yr %>%
    bind_cols(df) %>%
    dplyr::ungroup()
  return(df)
}

df_fit <- lm_list %>% purrr::map(get_fitted_vals) %>%
  bind_rows(.id="Latspecies_season_reccat")

df_fit <- df_fit %>%
  dplyr::mutate(Latspecies = stringr::str_split_i(Latspecies_season_reccat,"_",1),
                season_reccat = sapply(stringr::str_split(Latspecies_season_reccat,"_",n=2), "[[", 2) )

df_fit <- df_fit %>%
  left_join(pvalues, by="Latspecies_season_reccat") 

df_fit <- df_fit %>%
  dplyr::mutate(pred=ifelse(is.na(p_val),NA,pred)) %>%
  dplyr::mutate(style=ifelse(p_val<=0.1, "1", "2"))
# split string by delimeter
seas_rect <- strsplit(as.character(df_fit$season_reccat), "_")
# and use the splitted string to get each element
seass <- sapply(seas_rect, "[[", 1)
# add back as column on data frame
df_fit$season <- seass
#
df_fit <- df_fit %>%
  left_join(slopes, by="Latspecies_season_reccat") %>% 
  dplyr::mutate(year=ifelse(season=="jan-jun",year,year+0.5))


p_sig <- 0.1

df_note <-  df_fit %>%
  dplyr::mutate(slope=ifelse(is.na(p_val),NA, slope)) %>%
  dplyr::distinct(Latspecies, season_reccat, slope, p_val) %>%
  dplyr::mutate(p_text=ifelse(is.na(slope),NA_character_, ifelse(p_val<0.001,
                                                                 "p<0.001",
                                                                 ifelse(p_val> p_sig, "n.s.",
                                                                        paste0("p=",round(p_val,3)))))) %>%
  dplyr::mutate(note=ifelse(is.na(slope),NA_character_,
                            paste0(ifelse(slope>0,"+",""),
                                   round(100*slope,1),"%/år (", p_text, ")"))) %>%
  dplyr::group_by(Latspecies) %>%
  dplyr::mutate(n=sum(!is.na(slope))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(offset=ifelse(n>1,12,0)) %>%
  dplyr::mutate(y=max(dff$f,na.rm=T), year=min(dff$year,na.rm=T)) 

df_note <- df_note %>%
  dplyr::mutate(season= stringr::str_split_i(season_reccat,"_",1),
                reccat= stringr::str_split_i(season_reccat,"_",2)) %>%
  dplyr::mutate(y=ifelse(season=="jan-jun",y,y-offset)) %>%
  dplyr::mutate(y=ifelse(reccat=="iNaturalist",y,y-24))

#View(df_note)

#dff <- dff %>% dplyr::mutate(year=ifelse(season=="jan-jun",year,year+0.5))
mn.yer <- round(min(dff$year),0)
mx.yer <- ceiling(max(dff$year))

pal_season <- c("#0000FF", "#FF0000")
pal_season <- c( "slateblue4","firebrick",
                 "orchid3","darkorange3")
# get the individual unique seasons and record categories
reccat <- unique(dff$season_reccat)
cf.seas <- rbind(reccat,pal_season)
# make the color categories for the names
colnames(cf.seas) <- cf.seas[1,]
cf.seas <- cf.seas[-1,]
# get an upper limit for the y-axis
mx.fnd <- max(df_fit$upr)
# substitute in the column with season and record category, to hvae something
# more readable
# dff$season_reccat <- gsub("_",",",dff$season_reccat)
# make a plot
p <- ggplot(dff) +
  geom_ribbon(data=df_fit,
              aes(x=year, ymin=lwr, ymax=upr, 
                  group=season_reccat, fill=season_reccat),
              colour=NA, alpha=0.1) +
  geom_line(data=df_fit, aes(x=year, y=pred, group=season_reccat, 
                             colour=season_reccat, 
                             linetype = style, alpha=style)) +
  geom_point(aes(x=year, y=f, group=season_reccat, colour=season_reccat)) +
  geom_text(data=df_note, aes(x=year, y=y, group=season_reccat, 
                              colour=season_reccat, label=note),
            show.legend = F,
            hjust=0, size=3) +
  facet_wrap(.~Latspecies, ncol=4,
             labeller = label_bquote(col = italic(.(Latspecies))) ) +
  # make the x axis have breaks that are represented every 1 increment      
  scale_x_continuous(breaks=seq(mn.yer, mx.yer, 1)) +
  scale_y_continuous(limits = c(0, mx.fnd)) +
  
  scale_colour_manual(values=c(cf.seas), name="sæson og hjemmeside") +
  scale_alpha_manual(values=c(0.8, 0.6), guide="none") +
  scale_fill_manual(values=c(cf.seas), name="sæson og hjemmeside", guide="none") +
  scale_linetype_manual(values=c("solid", "longdash"), guide="none") +
  xlab("år") +
  ylab("antal af fund") +
  theme_minimal() + 
  guides(colour = guide_legend(nrow = 2),
         fill = guide_legend(nrow = 2)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        legend.position="bottom",
        strip.text = element_text(hjust = 0))
# see the plot
p
# define output directory
wd08 <- "output08_stckbarplot_w_GBIF_and_iNat_for_MST2017_2023_records"
outflNm <- "Fig22_trend_frekvens_arterdk_og_iNat.png"
wd00wd08_outflNm <- paste(wd00,wd08,outflNm, sep="/")
ggsave(p, 
       file=wd00wd08_outflNm, 
       height=20, width=20, units="cm", dpi=300, bg="white")
