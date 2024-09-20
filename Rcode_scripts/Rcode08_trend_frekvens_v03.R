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
# file <-"output06_presence_absence_evaluation_for_MST2017_2022_samples/table11_v01_all_pos_detections_abLOD.csv"
file <-"output06_presence_absence_evaluation_for_MST2017_2022_samples/table10_v01_all_pos_detections_abLOD.csv"
file <-"output06_presence_absence_evaluation_for_MST2017_2022_samples/table_12_v01_results_for_plotting_on_map.csv"
# read in the csv-file
df_obs <- read.table(file, sep=",", header=T)

# exclude if there is no Latin species name 
df_obs <- df_obs  %>%
  dplyr::rename(Latspecies=Lat_Species)   %>%
  dplyr::filter(!is.na(Latspecies)) %>%
  dplyr::filter(!(Latspecies==0)) %>%
  dplyr::filter(!(lokalitet_vanda==0)) %>%
  dplyr::rename(stn=lokalitet_vanda)
# get the sampling year and sampling month, and make the values numeric
df_obs <- df_obs %>%
  dplyr::mutate(year=as.numeric(substr(Dato_inds,1,4)),
         month=as.numeric(substr(Dato_inds,6,7)))
# split the string by delimiter
ssnmnth  <- strsplit(as.character(df_obs$yer_ssn2), ", ")
# and use the splitted string to get each element and add back in a column
df_obs$ssnmnth <- sapply(ssnmnth, "[[", 2)
# make it a factor
df_obs <- df_obs %>%
  dplyr::mutate(season=as.factor(ssnmnth))

# substitute to remove the space in the seasons 
df_obs$season <- gsub(" ","",df_obs$season)

# count up the number of attempts to search for the Latspecies, per year, per season
df_srchfor <- df_obs %>%
  dplyr::group_by(Latspecies, year, season) %>%
  dplyr::summarise(n_srchfor=n(), .groups="drop")
# count up the number of actual detections for the Latspecies, per year, per season
# but exlude those rows where organism was not found - i.e. when it
# was 'orgFnd2==0'
df_finds <- df_obs %>%
  dplyr::filter(!(orgFnd2==0)) %>%
  dplyr::group_by(Latspecies, year, season) %>%
  dplyr::summarise(n_fnd=n(), .groups="drop") 
#
df <- dplyr::distinct(df_obs, Latspecies) %>%
  merge(dplyr::distinct(df_obs, stn), all=T)  %>%
  merge(dplyr::distinct(df_obs, year), all=T)  %>%
  merge(dplyr::distinct(df_obs, season), all=T) 
# use the data frame with counts of attempts to search for the organism
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

dff_split <- dff %>%
  split(.$Latspecies_season)
# since the jul-nov season is half a season later, the year can be added 0.5 
dff <- dff %>% dplyr::mutate(year=ifelse(season=="jan-jun",year,year+0.5))

lm_list <- purrr::map(dff_split, lm, formula= f ~ year)

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
  dplyr::mutate(offset=ifelse(n>1,0.2,0)) %>%
  dplyr::mutate(y=max(dff$f,na.rm=T), year=min(dff$year,na.rm=T)) %>%
  dplyr::mutate(y=ifelse(season=="jan-jun",y,y-offset))

#View(df_note)

#dff <- dff %>% dplyr::mutate(year=ifelse(season=="jan-jun",year,year+0.5))
mn.yer <- round(min(dff$year),0)
mx.yer <- ceiling(max(dff$year))

pal_season <- c("#0000FF", "#FF0000")

p <- ggplot(dff) +
  geom_ribbon(data=df_fit,
              aes(x=year, ymin=lwr, ymax=upr, 
                  group=season, fill=season),
              colour=NA, alpha=0.1) +
  geom_line(data=df_fit,
            aes(x=year, y=pred, group=season, colour=season, 
                linetype = style, alpha=style)) +
  geom_point(aes(x=year, y=f, group=season, colour=season)) +
  geom_text(data=df_note, aes(x=year, y=y, group=season, 
                              colour=season, label=note),
            show.legend = F,
            hjust=0, size=3) +
  facet_wrap(.~Latspecies, ncol=4,
             labeller = label_bquote(col = italic(.(Latspecies))) ) +
  scale_colour_manual(values=pal_season, name="Sæson") +
  scale_alpha_manual(values=c(0.8, 0.6), guide="none") +
  scale_linetype_manual(values=c("solid", "longdash"), guide="none") +
  scale_fill_manual(values=pal_season, name="Sæson", guide="none") +
  # make the x axis have breaks that are represented every 1 increment      
  scale_x_continuous(breaks=seq(mn.yer, mx.yer, 1)) +
  xlab("år") +
  ylab("frekvens af fund") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        legend.position="bottom",
        strip.text = element_text(hjust = 0))
# define output directory
wd08 <- "output08_stckbarplot_w_GBIF_and_iNat_for_MST2017_2023_records"
outflNm <- "Fig18_trend_frekvens.png"
wd00wd08_outflNm <- paste(wd00,wd08,outflNm, sep="/")
ggsave(p, 
       file=wd00wd08_outflNm, 
       height=20, width=20, units="cm", dpi=300, bg="white")
#_______________________________________________________________________________

#geom_ribbon(aes(ymin=data$lower, ymax=data$upper), linetype=2, alpha=0.1)
p

#_______________________________________________________________________________

