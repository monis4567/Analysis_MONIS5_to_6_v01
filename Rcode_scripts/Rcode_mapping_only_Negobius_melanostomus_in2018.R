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
#_______________________________________________________________________________
# read in the csv-file
df_obs <- read.table(file, sep=",", header=T)
# exclude if there is no Latin species name 
df_obs <- df_obs  %>%
  dplyr::rename(Latspecies=Lat_Species)   %>%
  dplyr::filter(!is.na(Latspecies)) %>%
  dplyr::filter(!(Latspecies==0)) %>%
  dplyr::filter(!(lokalitet_vanda==0)) %>%
  dplyr::rename(stn=lokalitet_vanda)
# get the observations for "Neogobius melanostomus" in 2018
df_obsNeomel2018 <- df_obs %>% filter(Latspecies=="Neogobius melanostomus") %>%
  dplyr::filter(yea==2018)

library("rnaturalearth")
library("rnaturalearthdata")
library("rnaturalearthhires")
library("ggrepel")

# # Get a map, use a high number for 'scale' for a coarse resolution
# use a low number for scale for a high resolution
# if the map 'world' does not exist, then download it
world <- ne_countries(scale = 10, returnclass = "sf")
library(ggplot2)

#make plot
p05 <- ggplot(data = world) +
  geom_sf(color = "black", fill = "azure3", lwd=0.1) +
  # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  theme_void() +
  #https://ggplot2.tidyverse.org/reference/position_jitter.html
  #https://stackoverflow.com/questions/15706281/controlling-the-order-of-points-in-ggplot2
  # use 'geom_jitter' instead of 'geom_point' 
  geom_point(data = df_obsNeomel2018 ,
             aes(x = declon, 
                 y = declat,
                 color=orgFnd) ) + #,
  
  ggrepel::geom_text_repel(data = df_obsNeomel2018 ,
                           aes(x = declon, 
                               y = declat,
                               label=MST.nummer) ,
                           size=1.0,
                           alpha=0.8) +
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
                                    color = "white",
                                    linewidth = 1)) +
  #define limits of the plot 
  ggplot2::coord_sf(xlim = c(6.6, 17.2),
                    ylim = c(54.2, 58.4), 
                    expand = FALSE) +
  # Set the colors of the outline of the points manually
  #scale_color_manual(values=c(alpha( co_f_mnscl,al_f_mnscl) ), breaks = catf_plt ) +
  
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
p05t <- p05 + xlab(" ") + ylab(" ")
#change the header for the legend on the side, 
#this must be done for both 'fill', 'color' and 'shape', to avoid 
#getting separate legends
# p05t <- p05t + labs(color='')
# p05t <- p05t + labs(fill='')
# p05t <- p05t + labs(shape='')
# p05t <- p05t + labs(size='')
#adjust tick marks on axis
p05t <- p05t + scale_y_continuous(breaks=seq(54.2, 58.4,2))
p05t <- p05t + scale_x_continuous(breaks=seq(6.6, 17.2,4))

# see the plot
#p05t
# define output directory
wd08 <- "output08_stckbarplot_w_GBIF_and_iNat_for_MST2017_2023_records"
outflNm <- "Fig19_map_Neogobius_melanostomus_2018.png"
wd00wd08_outflNm <- paste(wd00,wd08,outflNm, sep="/")
ggsave(p05t, 
       file=wd00wd08_outflNm, 
       height=297*0.3, 
       width=210, 
       units="mm", dpi=300, bg="white")

#

df_obsNeomel2018.pos <- df_obsNeomel2018[(df_obsNeomel2018$orgF==TRUE),]
