library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
library(purrr)

# file <-"output06_presence_absence_evaluation_for_MST2017_2022_samples/table11_v01_all_pos_detections_abLOD.csv"
file <-"output06_presence_absence_evaluation_for_MST2017_2022_samples/table10_v01_all_pos_detections_abLOD.csv"

df_obs <- read.table(file, sep=";", header=T) %>%
  rename(species=latinsk.artsnavn,
         stn=lokalitet_vanda)  %>%
  filter(!is.na(species))

df_obs <- df_obs %>%
  mutate(year=as.numeric(substr(dato.indsamlet,1,4)),
         month=as.numeric(substr(dato.indsamlet,6,7)))

 df_obs <- df_obs %>%
   mutate(season=as.factor(sæson))

df_finds <- df_obs %>%
  group_by(stn, species, year, season) %>%
  summarise(n=n(), .groups="drop") %>%
  mutate(detected=1)

df <- distinct(df_obs, species) %>%
  merge(distinct(df_obs, stn), all=T)  %>%
  merge(distinct(df_obs, year), all=T)  %>%
  merge(distinct(df_obs, season), all=T) 


df <- df %>%
  left_join(df_finds, by=c("stn", "species", "year", "season")) %>%
  mutate(detected=ifelse(is.na(detected),0,detected))

df <- df %>%
  mutate(stn=as.factor(stn))


dff <- df %>%
  group_by(species, year,season) %>%
  summarise(n=n(), f=sum(detected)/n(), .groups="drop") %>%
  mutate(species_season=paste0(species,"_",season))


dff_split <- dff %>%
  split(.$species_season)

lm_list <- purrr::map(dff_split, lm, formula= f ~ year)

pvalues <- lm_list %>% 
  map(summary.lm) %>% 
  map(c("coefficients")) %>% 
  map_dbl(8)  %>% # 8th element is the p-value 
  broom::tidy() %>% 
  dplyr::arrange(desc(x)) %>% 
  rename(species_season=names, p_val = x)

slopes <- lm_list %>% 
  map(summary.lm) %>% 
  map(c("coefficients")) %>% 
  map_dbl(2)  %>% # 8th element is the p-value 
  broom::tidy() %>% 
  dplyr::arrange(desc(x)) %>% 
  rename(species_season=names, slope = x)


get_fitted_vals <- function(lm){
  return(lm$fitted.values)
}

df_fit <- lm_list %>% purrr::map(get_fitted_vals) %>%
  bind_rows(.id="species_season")

df_fit <- df_fit %>%
  pivot_longer(cols=2:ncol(df_fit), names_to="year", values_to="pred")

df_fit <- df_fit %>%
  mutate(year=as.numeric(year))

df_fit <- df_fit %>%
  mutate(species = stringr::str_split_i(species_season,"_",1),
         season = stringr::str_split_i(species_season,"_",2))


df_fit$year <- sort(unique(df$year))[df_fit$year]

df_fit <- df_fit %>%
  left_join(pvalues, by="species_season") 

df_fit <- df_fit %>%
  mutate(pred=ifelse(is.na(p_val),NA,pred)) %>%
  mutate(style=ifelse(p_val<=0.1, "1", "2"))

df_fit <- df_fit %>%
  left_join(slopes, by="species_season") 

p_sig <- 0.1
df_note <-  df_fit %>%
  mutate(slope=ifelse(is.na(p_val),NA, slope)) %>%
  distinct(species, season, slope, p_val) %>%
  mutate(p_text=ifelse(is.na(slope),NA_character_, ifelse(p_val<0.001,
                                                          "p<0.001",
                                                          ifelse(p_val> p_sig, "n.s.",
                                                                 paste0("p=",round(p_val,3)))))) %>%
  mutate(note=ifelse(is.na(slope),NA_character_,
                     paste0(ifelse(slope>0,"+",""),round(100*slope,1),"%/år (", p_text, ")"))) %>%
  group_by(species) %>%
  mutate(n=sum(!is.na(slope))) %>%
  ungroup() %>%
  mutate(offset=ifelse(n>1,0.05,0)) %>%
  mutate(y=max(dff$f,na.rm=T), year=min(dff$year,na.rm=T)) %>%
  mutate(y=ifelse(season=="jan-jun",y,y-offset))

pal_season <- c("#0000FF", "#FF0000")

p <- ggplot(dff) +
  geom_line(data=df_fit, aes(x=year, y=pred, group=season, colour=season, 
                             linetype = style, alpha=style)) +
  geom_point(aes(x=year, y=f, group=season, colour=season)) +
  geom_text(data=df_note, aes(x=year, y=y, group=season, 
                              colour=season, label=note),
            show.legend = F,
            hjust=0, size=3) +
  facet_wrap(.~species, ncol=4) +
  scale_colour_manual(values=pal_season, name="Sæson") +
  scale_alpha_manual(values=c(0.8, 0.6), guide="none") +
  scale_linetype_manual(values=c("solid", "longdash"), guide="none") +
  xlab("år") +
  ylab("frekvens af fund") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        legend.position="bottom",
        strip.text = element_text(hjust = 0))

#ggsave(p, file="cjm/trend_frekvens.png", height=20, width=20, units="cm", dpi=300, bg="white")

p