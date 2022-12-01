library(ggplot2)
library(Rmisc)
library(dplyr)
library(gridExtra)
library(gtable)
library(ggpubr)
library(grid)
library("cowplot")
library(reshape2)
library(oceanmap)
library(sf)
library(ggmap)
library(data.table)
library(plyr)

setwd("/Users/Christopher Johns/Dropbox/Lith Adsorption/Spreadsheets")
adsorption <- read.csv("Lith_Adsorption_Mastersheet_updated_lith_protection.csv")
virus_wash <- read.csv("virus counts.csv")

lith_protection <- subset(adsorption, experiment.name == "Lith Protection 1" | experiment.name == "Lith Protection 2", select = c(1:31))

lith_protection$cell.concentration <- as.numeric(lith_protection$cell.concentration)
lith_protection$lith.concentration <- as.numeric(as.character(lith_protection$lith.concentration))

lith_protection_cells <- summarySE(lith_protection, measurevar = "cell.concentration",  
                      groupvars = c("experiment.name", "hour","treatment","liths"))

lith_protection_liths <- summarySE(lith_protection, measurevar = "lith.concentration",  
                      groupvars = c("experiment.name", "hour","treatment","liths"))



#plotting the growth curve 
  ggplot(lith_protection_cells, aes(x=hour, y=log10(cell.concentration)))+  
  geom_point(size=5)+theme_bw()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  theme_classic()+
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour")+  
  scale_y_continuous(name = bquote('Cell concentration; log10 ('*cells ~ mL^-1*')'))+
  # scale_colour_manual("plotID",
  #                     labels = c("OM control","IM control","OM + coccoliths", "IM + coccoliths"),
  #                     values = c("control im-control" = "#FFCC00",
  #                                "viral im-control" = "#FF3300",
  #                                "control im-cells+liths" = "#006699",  
  #                                "viral im-cells+liths" = "#3399FF"))  
  facet_wrap(~experiment.name)

ggplot(lith_protection, aes(x=as.factor(hour),y=percent.calcified, colour=treatment))+  
    geom_boxplot()+
    geom_point(aes(colour = treatment),
               size = 4,
               position = position_jitterdodge())+
    theme_bw()+
    theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme_classic()+
    # scale_colour_manual("plotID",
    #                     labels = c("OM control","IM control","OM + coccoliths", "IM + coccoliths"),
    #                     values = c("control im-control" = "#FFCC00",
    #                                "viral im-control" = "#FF3300",
    #                                "control im-cells+liths" = "#006699",  
    #                                "viral im-cells+liths" = "#3399FF"))+
    scale_y_continuous(name = "Percent calcified cells", limits = c(0,100), breaks = c(0,25,50,75,100))+
    scale_x_discrete(name = "Hour")+  
  facet_wrap(~experiment.name)

virus_wash %>% 
  filter(date == "11/11/2020") %>%
  arrange(percent.decrease) %>%
  mutate(sample = factor(sample, levels=c("filtrate",  
                                          "3",  
                                          "5"))) %>%  
ggplot(aes(x=sample,y=percent.decrease, colour=sample))+  
  geom_boxplot()+
  geom_point(aes(colour = sample),
             size = 4,
             position = position_jitterdodge())+
  theme_bw()+
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_classic()+
  # scale_colour_manual("plotID",
  #                     labels = c("OM control","IM control","OM + coccoliths", "IM + coccoliths"),
  #                     values = c("control im-control" = "#FFCC00",
  #                                "viral im-control" = "#FF3300",
  #                                "control im-cells+liths" = "#006699",  
  #                                "viral im-cells+liths" = "#3399FF"))+
  scale_y_continuous(name = "virus concentration")+
  scale_x_discrete(name = "Virolith preperation step")

virus_wash %>%  
  filter(date == "11/11/2020") %>%  
  filter(sample != "filtrate") %>%  
  filter(sample != "filtrate wash") %>%  
  filter(sample != "2") %>%
  arrange(total.viruses) %>%
  mutate(sample = factor(sample, levels=c("lysate", "3", "4", "5", "6"))) %>%  
  ggplot(aes(x=sample,y=log10(total.viruses), colour=sample))+  
  geom_boxplot()+
  geom_point(aes(colour = sample),
             size = 4,
             position = position_jitterdodge())+
  theme_bw()+
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_classic()+
  # scale_colour_manual("plotID",
  #                     labels = c("OM control","IM control","OM + coccoliths", "IM + coccoliths"),
  #                     values = c("control im-control" = "#FFCC00",
  #                                "viral im-control" = "#FF3300",
  #                                "control im-cells+liths" = "#006699",  
  #                                "viral im-cells+liths" = "#3399FF"))+
  scale_y_continuous(name = "virus concentration")+
  scale_x_discrete(name = "Virolith preperation step")+facet_wrap(~date)

virus_wash %>%  
  filter(date == "11/11/2020") %>%  
  filter(sample != "filtrate wash") %>%  
  filter(sample != "2") %>%  
  filter(sample != "4") %>%  
  filter(sample != "6") %>%
  arrange(ehv.to.liths) %>%
  mutate(sample = factor(sample, levels=c("filtrate", "3", "5"))) %>%  
  ggplot(aes(x=sample,y=ehv.to.liths, colour=sample))+  
  geom_boxplot()+
  geom_point(aes(colour = sample),
             size = 4,
             position = position_jitterdodge())+
  theme_bw()+
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_classic()+
  # scale_colour_manual("plotID",
  #                     labels = c("OM control","IM control","OM + coccoliths", "IM + coccoliths"),
  #                     values = c("control im-control" = "#FFCC00",
  #                                "viral im-control" = "#FF3300",
  #                                "control im-cells+liths" = "#006699",  
  #                                "viral im-cells+liths" = "#3399FF"))+
  scale_y_continuous(name = "Number of viruses")+
  scale_x_discrete(name = "Virolith preperation step")+facet_wrap(~date)
