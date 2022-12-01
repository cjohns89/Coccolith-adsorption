library(ggplot2)
library(Rmisc)
library(dplyr)
library(gridExtra)
library(gtable)
library(ggpubr)
library(grid)
library("cowplot")
library(reshape2)
library(data.table)
library(plyr)
library(expss)
library(FSA)

setwd("/Users/cjohns/Dropbox/Lith Adsorption/Spreadsheets")

mpn <- read.csv("./MPNs/virolith_mpn.csv")
traditional_mpn <- read.csv("./MPNs/20221003_liths_or_no_liths_traditional_mpn.csv")

#mpn time serie####################################################################################################
mpn$corrected_id <- as.factor(paste(mpn$row_letter, mpn$cell_concentration, mpn$virus_host_ratio, mpn$liths, sep="-"))
mpn$plotID <- as.factor(paste(mpn$virus_host_ratio, mpn$liths, sep="-"))

#conversion into cell concentration
blank_naked <- 0.040484375
blank_liths <- 0.062934375

mpn$od_blank <- rep (c("0.062934375", "0.040484375"), each=48)
mpn$od_blank <- as.numeric(mpn$od_blank)


mpn$od_blank_sub <- mpn$od - mpn$od_blank
mpn$est_cell_concentration <- mpn$od_blank_sub *5e7


mpn <- data.table::data.table(mpn, key= c("corrected_id"))
mpn[, change_od:= (od-(od[match("0", hour)])), by= c("corrected_id")]

mpn_cell_od_100000 <- mpn %>% filter(cell_concentration == "100000") %>%
  summarySE(measurevar = "change_od",  
            groupvars = c("hour","liths","virus_host_ratio"), na.rm = TRUE)

mpn_cells <- mpn %>% filter(cell_concentration == "100000") %>%
  summarySE(measurevar = "est_cell_concentration",  
            groupvars = c("hour","liths","virus_host_ratio"), na.rm = TRUE)

mpn_cell_od_10000 <- mpn %>% filter(cell_concentration == "10000") %>%
  summarySE(measurevar = "change_od",  
            groupvars = c("hour","liths","virus_host_ratio"), na.rm = TRUE)

#plotting as od values
no_liths <- mpn_cell_od_100000 %>% filter(liths == "control") %>% filter(hour < 241) %>% filter(virus_host_ratio != "control") %>% filter(virus_host_ratio != "0.001_1") %>%
  arrange(change_od) %>%
  mutate(virus_host_ratio = factor(virus_host_ratio, levels=c("control",  
                                                              "10_1",
                                                              "1_1",  
                                                              "0.1_1",  
                                                              "0.01_1",  
                                                              "0.001_1"))) %>%  
  ggplot(aes(x=hour, y=change_od, colour=virus_host_ratio, shape=virus_host_ratio))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 15, angle = 45, hjust = 0.75),
        axis.text.y = element_text(size = 15),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=change_od-se, ymax=change_od+se), width= 5, size=0.75)+  
  geom_line(aes(colour=virus_host_ratio,group=virus_host_ratio))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192,216,240,264,288,312))+  
  scale_y_continuous(name = "OD (750nm)", limits = c(0,0.4), breaks = c(0,0.1,0.2,0.3,0.4))+
  scale_colour_manual("V:H Ratio",
                      # labels = c("OM control","OM + coccoliths","IM control","IM + coccoliths"),
                      values = c("10_1" = "#009E73",
                                 "1_1" = "#009E73",
                                 "0.1_1" = "#009E73",
                                 "0.01_1" = "#009E73"), guide = "none")+  
  scale_shape_manual(values = c(15,16,17,18))

liths <- mpn_cell_od_100000 %>% filter(liths == "liths") %>% filter(hour < 241) %>% filter(virus_host_ratio != "control") %>% filter(virus_host_ratio != "0.001_1") %>%
  arrange(change_od) %>%
  mutate(virus_host_ratio = factor(virus_host_ratio, levels=c("control",  
                                                              "10_1",
                                                              "1_1",  
                                                              "0.1_1",  
                                                              "0.01_1",  
                                                              "0.001_1"))) %>%  
  ggplot(aes(x=hour, y=change_od, colour=virus_host_ratio, shape=virus_host_ratio))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 15, angle = 45, hjust = 0.75),
        axis.text.y = element_text(size = 15),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=change_od-se, ymax=change_od+se), width= 5, size=0.75)+  
  geom_line(aes(colour=virus_host_ratio,group=virus_host_ratio))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192,216,240,264,288,312))+  
  scale_y_continuous(name = "OD (750nm)", limits = c(0,0.4), breaks = c(0,0.1,0.2,0.3,0.4))+
  scale_colour_manual("V:H Ratio",
                      # labels = c("OM control","OM + coccoliths","IM control","IM + coccoliths"),
                      values = c("10_1" = "#56B4E9",
                                 "1_1" = "#56B4E9",
                                 "0.1_1" = "#56B4E9",
                                 "0.01_1" = "#56B4E9"), guide = "none")+  
  scale_shape_manual(values = c(15,16,17,18))

#plotting as the converted od values to cell concentration
mpn_cells %>% filter(hour < 241) %>%    
  filter(virus_host_ratio != "control") %>%
  arrange(est_cell_concentration) %>%
  mutate(virus_host_ratio = factor(virus_host_ratio, levels=c("control",  
                                                              "10_1",
                                                              "1_1",  
                                                              "0.1_1",  
                                                              "0.01_1",  
                                                              "0.001_1"))) %>%  
  ggplot(aes(x=hour, y=est_cell_concentration, colour=virus_host_ratio))+  
  geom_point(size=5)+theme_bw()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+
  geom_errorbar(aes(ymin=est_cell_concentration-se, ymax=est_cell_concentration+se), width= 5, size=0.75)+  
  geom_line(aes(colour=virus_host_ratio,group=virus_host_ratio))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168,192,216,240,264,288,312))+  
  scale_y_continuous()+
  # scale_colour_manual("plotID",
  #                     # labels = c("OM control","OM + coccoliths","IM control","IM + coccoliths"),
  #                     values = c("control im-control" = "#FFCC00",
  #                                "viral im-control" = "#FF3300",
  #                                "control im-cells+liths" = "#006699",  
  #                                "viral im-cells+liths" = "#3399FF"))  
  facet_wrap(~liths)

mpn_144 <- filter(mpn, hour == "144")
mpn_144 <- mpn_144 %>% filter(cell_concentration == "100000") %>% filter(virus_host_ratio == "10_1")

mpn_240 <- mpn %>% filter(hour == "240") %>% filter(virus_host_ratio != "control") %>% filter(virus_host_ratio != "0.001_1")
#subsetting all of the different groups
mpn_240_10_1 <- mpn_240 %>% filter(virus_host_ratio == "10_1")
mpn_240_1_1 <- mpn_240 %>% filter(virus_host_ratio == "1_1")
mpn_240_0.1_1 <- mpn_240 %>% filter(virus_host_ratio == "0.1_1")
mpn_240_0.01_1 <- mpn_240 %>% filter(virus_host_ratio == "0.01_1")


shapiro.test(mpn_240_10_1$change_od) #pvalue is less than 0.05 (0.0343), implying the data is not normal and I'll use a Mann-Whitney
wilcox.test(change_od ~ liths, data = mpn_240_10_1, paired = FALSE) #pvalue 0.06496
shapiro.test(mpn_240_1_1$change_od)  #pvalue is less than 0.05 (6e-4), implying the data is not normal and I'll use a Mann-Whitney
wilcox.test(change_od ~ liths, data = mpn_240_1_1, paired = FALSE) #pvalue 0.1949
shapiro.test(mpn_240_0.1_1$change_od) #pvalue is less than 0.05 (0.0036), implying the data is not normal and I'll use a Mann-Whitney
wilcox.test(change_od ~ liths, data = mpn_240_0.1_1, paired = FALSE) #pvalue 0.1949
shapiro.test(mpn_240_0.01_1$change_od) #pvalue is greater than 0.05 (0.7825), implying the data is normal and I'll use a t-test
t.test(change_od ~ liths, data = mpn_240_0.01_1, var.equal = TRUE) #pvalue 0.2286



mpn_box_plots <- rbind(mpn_144, mpn_240)
  

  
od_boxplot <- mpn %>% subset(hour == "240") %>%
  filter(cell_concentration == "100000") %>%
  filter(virus_host_ratio != "control") %>%  
  filter(virus_host_ratio != "0.001_1") %>%
  arrange(change_od) %>%
  mutate(plotID = factor(plotID, levels=c("10_1-control",
                                          "10_1-liths", 
                                          "1_1-control", 
                                          "1_1-liths",  
                                          "0.1_1-control", 
                                          "0.1_1-liths",
                                          "0.01_1-control", 
                                          "0.01_1-liths"))) %>%  
  ggplot(aes(x=plotID, y=log10(change_od), colour=plotID))+  
  geom_boxplot()+
  geom_point(aes(colour = plotID),
             size = 4,
             position = position_jitterdodge())+
  theme_bw()+
  theme_classic()+
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_colour_manual("plotID",
                      # labels = c("OM control","OM + coccoliths","IM control","IM + coccoliths"),
                      values = c("10_1-control" = "#009E73",
                                 "10_1-liths" = "#56B4E9",
                                 "1_1-control" = "#009E73",
                                 "1_1-liths" = "#56B4E9",
                                 "0.1_1-control" = "#009E73",
                                 "0.1_1-liths" = "#56B4E9",
                                 "0.01_1-control" = "#009E73",
                                 "0.01_1-liths" = "#56B4E9"))+
  scale_y_continuous()+
  scale_x_discrete(name = "Virus Host Ratio")

ggarrange(no_liths, liths, od_boxplot, virolith_cells, align = 'hv', ncol = 2, nrow = 2)


mpn_cell_od_10000 %>% filter(hour < 241) %>% 
  arrange(change_od) %>%
  mutate(virus_host_ratio = factor(virus_host_ratio, levels=c("control",  
                                                              "10_1",
                                                              "1_1",  
                                                              "0.1_1",  
                                                              "0.01_1",  
                                                              "0.001_1"))) %>%  
  ggplot(aes(x=hour, y=change_od, colour=virus_host_ratio))+  
  geom_point(size=5)+theme_bw()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+
  geom_errorbar(aes(ymin=change_od-se, ymax=change_od+se), width= 5, size=0.75)+  
  geom_line(aes(colour=virus_host_ratio,group=virus_host_ratio))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168,192,216,240,264,288,312))+  
  scale_y_continuous()+
  # scale_colour_manual("plotID",
  #                     # labels = c("OM control","OM + coccoliths","IM control","IM + coccoliths"),
  #                     values = c("control im-control" = "#FFCC00",
  #                                "viral im-control" = "#FF3300",
  #                                "control im-cells+liths" = "#006699",  
  #                                "viral im-cells+liths" = "#3399FF"))  
  facet_wrap(~liths)

mpn_cell_od <- mpn %>%
  summarySE(measurevar = "change_od",  
            groupvars = c("cell_concentration", "hour","liths","virus_host_ratio"), na.rm = TRUE)


mpn_cell_od %>% 
  arrange(change_od) %>%
  mutate(virus_host_ratio = factor(virus_host_ratio, levels=c("control",  
                                                              "10_1",
                                                              "1_1",  
                                                              "0.1_1",  
                                                              "0.01_1",  
                                                              "0.001_1"))) %>%  
  ggplot(aes(x=hour, y=change_od, colour=virus_host_ratio))+  
  geom_point(size=5)+theme_bw()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+
  geom_errorbar(aes(ymin=change_od-se, ymax=change_od+se), width= 5, size=0.75)+  
  geom_line(aes(colour=virus_host_ratio,group=virus_host_ratio))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168,192,216))+  
  scale_y_continuous()+
  # scale_colour_manual("plotID",
  #                     # labels = c("OM control","OM + coccoliths","IM control","IM + coccoliths"),
  #                     values = c("control im-control" = "#FFCC00",
  #                                "viral im-control" = "#FF3300",
  #                                "control im-cells+liths" = "#006699",  
  #                                "viral im-cells+liths" = "#3399FF"))  
  facet_grid(cell_concentration~liths)

#traditional mpn to measure infectiousnesss########################################################################
#removing the blanks
traditional_mpn <- filter(traditional_mpn, dilution_factor != "blank")

traditional_mpn$corrected_id <- as.factor(paste(traditional_mpn$sample_id, traditional_mpn$dilution_factor, traditional_mpn$liths, sep="-"))
traditional_mpn$plotID <- as.factor(paste(traditional_mpn$dilution_factor, traditional_mpn$liths, sep="-"))


traditional_mpn <- data.table(traditional_mpn, key= c("corrected_id"))
traditional_mpn[, change_od:= (od-(od[match("0", hour)])), by= c("corrected_id")]



traditional_mpn %>%  
  filter(hour > 2) %>%
  arrange(change_od) %>% 
  mutate(dilution_factor = factor(dilution_factor, levels=c("0",
                                                            "5", 
                                                            "25",  
                                                            "125", 
                                                            "625", 
                                                            "3125", 
                                                            "15625",
                                                            "78125", 
                                                            "390625", 
                                                            "1953125"))) %>%  
  ggplot(aes(x=as.factor(dilution_factor),y=change_od, colour=liths))+  
  geom_boxplot()+
  geom_point(aes(colour = liths),
             size = 4,
             position = position_jitterdodge())+
  theme_bw()+
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_classic()+
  # scale_colour_manual("plotID",
  #                     # labels = c("OM control","OM + coccoliths","IM control","IM + coccoliths"),
  #                     values = c("control im-control" = "#FFCC00",
  #                                "viral im-control" = "#FF3300",
  #                                "control im-cells+liths" = "#006699",  
  #                                "viral im-cells+liths" = "#3399FF"))+
  scale_y_continuous(name = "change in od")+
  scale_x_discrete(name = "dilution factor")+  
  facet_wrap(~as.factor(hour))

traditional_mpn_cell_od <- mpn %>%
  summarySE(measurevar = "change_od",  
            groupvars = c("cell_concentration", "hour","liths","virus_host_ratio"), na.rm = TRUE)