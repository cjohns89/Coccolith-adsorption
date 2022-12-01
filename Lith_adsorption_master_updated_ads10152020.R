setwd("/Users/cjohns/Dropbox/Lith Adsorption/Spreadsheets/adsorptive_exchange")
adsorption <- read.csv("Lith_Adsorption_Mastersheet_20191011.csv")
lithspikestress <- read.csv("Lith Spike Infection II 20181308.csv")
measurements <- read.csv("percent_calcified_measurements.csv")
lith_infection <- read.csv("Lith_infection_20182903.csv")
library(ggplot2)
library(dplyr)
library(Rmisc)
library(tidyverse)
library(gridExtra)
library(gtable)
library(ggpubr)
library(grid)
library("cowplot")
library(reshape2)
library(data.table)
library(plotly)
library(forcats)
library(rcompanion)
library("agricolae") #need to install questionr package prior, otherwise it will fail
library(multcomp)
library(multcompView)
library(lawstat)
library(FSA)
library(broom)
library(lemon)
library(doBy)
#percent calcified via lith adsorption######################################################################################################################################################################################
#filtering everything that I've given 374 liths to
growth <- subset(adsorption, liths == "control" | liths == "cells+liths", select = c(1:21))

#filtering out experiments not used in this analysis
just374 <- filter(growth, host.strain != "1516")
just374 <- filter(just374, host.strain != "phaeocystis")
just374 <- filter(just374, experiment.name != "Type O")
just374 <- filter(just374, experiment.name != "Morphotype Screen")
just374 <- filter(just374, experiment.name != "ViroLiths")
just374 <- filter(just374, experiment.name != "Lith Adsorption Dilution Series")
just374 <- filter(just374, experiment.name != "G_Oceanica Lith Adsorption")
just374 <- filter(just374, experiment.name != "calcein_dilution_series")
just374 <- filter(just374, treatment != "bleached 374 liths") 
just374 <- filter(just374, treatment != "374 proteinase liths")
just374 <- filter(just374, treatment != "374 lipase liths") 
just374 <- filter(just374, treatment != "374 amylase liths ")
just374 <- filter(just374, treatment != "374 amylase liths")
just374 <- filter(just374, experiment.name != "viroliths 2")
just374 <- filter(just374, experiment.name != "iso_tiso")
just374 <- filter(just374, experiment.name != "viroliths 3")
just374 <- filter(just374, experiment.name != "viroliths 4")

#uses percent calcified as detected by side scatter (SSC)
ggplot(just374,aes(x=as.factor(hour), y=percent.calcified, colour= liths))+ geom_boxplot()+
  geom_point(aes(colour = liths),  
             size = 4, 
             position = position_jitterdodge())+  
  theme_bw()+  
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),  
                   panel.grid.major = element_blank(),  
                   panel.grid.minor = element_blank())+  
  theme_classic()+  
  scale_color_manual(name = "Treatment", labels = c("control" = "Naked", "cells+liths" = "Cells + Liths"), values = c("control" = "#E69F00", "cells+liths" = "#0072B2"))+  
  scale_y_continuous(name = "Percent calcified cells (Side scatter; SSC)", limits = c(0,100), breaks = c(0,25,50,75,100))+  
  scale_x_discrete(name = "Hour")
            
#uses percent calcified as detected by 520nm fluorescence          
ggplot(just374,aes(x=as.factor(hour), y=percent.calcified.calcein, colour= liths))+ geom_boxplot()+
  geom_point(aes(colour = liths),  
             size = 4, 
             position = position_jitterdodge())+
  theme_bw()+  
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+  
  scale_color_manual(name = "Treatment", labels = c("control" = "Naked", "cells+liths" = "Cells + Liths"), values = c("control" = "#E69F00", "cells+liths" = "#0072B2"))+  
  scale_y_continuous(name = "Percent calcified cells (520nm fluorescence)", limits = c(0,100), breaks = c(0,25,50,75,100))+  
  scale_x_discrete(name = "Hour")

#summarzing the two different datasets
just374_sum <- summarySE(just374, measurevar = "percent.calcified", groupvars = c("hour","liths"))
just374_calcein_sum <- summarySE(just374, measurevar = "percent.calcified.calcein", groupvars = c("hour","liths"))

#comparing the percent calcified across different measurements
measurements %>%
  arrange(percent.calcified) %>%
  mutate(metric = factor(metric, levels=c("micro", "520", "ssc"))) %>%

ggplot(aes(x=as.factor(hour), y=percent.calcified, colour= metric))+ geom_boxplot()+
  geom_point(aes(colour = metric),  
             size = 4, 
             position = position_jitterdodge())+
  theme_bw()+  
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+
  scale_color_manual(name = "Measurement Type", labels = c("micro" = "Microscopy", "520" = "Flow Cytometry (520nm)", "ssc" = "Flow Cytometry (SSC)"), values = c("micro" = "#44AA99", "520" = "#999933", "ssc" = "#332288"))+
  scale_y_continuous(name = "Percent calcified cells", limits = c(0,100), breaks = c(0,25,50,75,100))+
  scale_x_discrete(name = "Hour")

#growth rates######################################################################################################################################################################################
#calculating growth rates
#plotting the specific growth rates after summarizing each replicate and finding the slope
just374$statID <- as.factor(paste(just374$liths, just374$replicate, just374$experiment.name, just374$treatment,sep="-"))

#need to convert it into a numeric value
just374$cell.concentration <- as.numeric(as.character(just374$cell.concentration))
slopes <- just374 %>% group_by(experiment.name, statID, liths) %>%
  do(model = lm(log(cell.concentration) ~ hour, data = .)) %>%
  mutate(coef=coef(model)["hour"])

#converts it into per day
slopes$growth_day <- slopes$coef *24

rates <- slopes %>%
  arrange(growth_day) %>%
  mutate(liths = factor(liths, levels=c("control", "cells+liths"))) %>%

ggplot(aes(x= liths, y=growth_day, colour = liths))+  
  geom_boxplot()+  
  geom_point(aes(colour = liths), size = 5,
             position = position_jitterdodge())+  
  theme_classic()+
  theme(legend.position = "none",  
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  scale_x_discrete(name = "Treatment", labels = c("Control","Non-calcified cells\n+coccoliths"))+
  scale_y_continuous(name = bquote('Specific Growth Rate ('*mu ~ day^-1*')'))+
  scale_fill_manual(values = c("#E69F00","#0072B2"))+
  scale_colour_manual(labels = c("Control",
                                 "Non-calcified cells\n+coccoliths"),
                      values = c("control" = "#E69F00",
                                 "cells+liths" = "#0072B2"))+  
  annotate("text", x=2, y= 0.95, label = "p-value = 0.06046")

#looking at the distribution of the data
ggplot(slopes, aes(x=growth_day, fill=liths)) +
  geom_density(alpha=0.4)
ggplot(slopes, aes(x=liths, y = growth_day)) +
  geom_violin()

#stats
ggplot(slopes, aes(growth_day)) +
  geom_histogram(fill = "white", color = "grey30")+facet_wrap(~liths)
#histogram and qqplots
hist(slopes$growth_day)
qqnorm(slopes$growth_day, pch = 1, frame = FALSE)
qqline(slopes$growth_day, col = "steelblue", lwd = 2)

shapiro.test(slopes$growth_day) #this suggests the data is not normally distributed, the pvalue is less than 0.05 (0.01829)
levene.test(slopes$growth_day, slopes$liths)#based on the levene test the variances are homogenous 
#data isn't normally distributed so I'm going to do non-parametric test, mann whitney because samples are independent from each other
wilcox_output <- wilcox.test(growth_day ~ liths, data = slopes, paired = FALSE) #not paired data so this is a Mann-Whitney for independent samples
#also looked at t.test to confirm
t.test(growth_day ~ liths, data = slopes, var.equal = TRUE)

#not significantly different p-value = 0.06046

#looking at the summary stats
sum_stats_growth <- summarySE(slopes, measurevar = "growth_day",  
                                   groupvars = c("liths"))
#other hosts#######################################################################################################################################################################################
#this shows any host that wasn't CCMP374
ehux_1516 <- subset(growth, host.strain == "1516")
ehux_1516 <- filter(ehux_1516, hour < 100)
haptophytes <- subset(growth, host.strain == "phaeocystis" | host.strain == "isochyrsis" | host.strain == "tisochyrsis", select = c(1:21))
haptophytes <- filter(haptophytes, hour < 100)

all_other_hosts <-subset(growth, host.strain == "1516" | host.strain == "phaeocystis" | host.strain == "isochyrsis" | host.strain == "tisochyrsis")
all_other_hosts <- filter(all_other_hosts, hour < 100)

ehux <- ggplot(ehux_1516,aes(x=as.factor(hour), y=percent.calcified, colour= liths))+ geom_boxplot()+
  geom_point(aes(colour = liths),  
             size = 4, 
             position = position_jitterdodge())+  
  theme_bw()+  
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+  
  scale_color_manual(name = "Treatment", labels = c("control" = "Naked", "cells+liths" = "Cells + Liths"), values = c("control" = "#E69F00", "cells+liths" = "#0072B2"))+  
  scale_y_continuous(name = "Percent calcified cells", limits = c(0,75), breaks = c(0,25,50,75))+  
  scale_x_discrete(name = "Hour")

other_hosts <- all_other_hosts %>% filter(liths != "control") %>% ggplot(aes(x=as.factor(hour), y=percent.calcified, colour= host.strain))+ geom_boxplot()+
  geom_point(size = 4, 
             position = position_jitterdodge())+  
  theme_bw()+
  theme_classic()+  
  theme(legend.position = c(0.2,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  scale_color_manual(name = "Treatment", values = c("1516" = "#CC6677",
                                                    "isochyrsis" = "#117733",
                                                    "phaeocystis" = "#332288",
                                                    "tisochyrsis" = "#44AA99"))+
  scale_y_continuous(name = "Percent calcified cells", limits = c(0,40), breaks = c(0,10,20,30,40))+  
  scale_x_discrete(name = "Hour", breaks = c(0,24,48,72,96))

#summary stats for the two different plots
# sum_stats_other_hosts <- summarySE(otherhosts, measurevar = "percent.calcified",  
#                                             groupvars = c("experiment.name","hour","liths"))

g_ocean <- measurements %>% filter(cell_type == "g_ocean") %>%
  arrange(percent.calcified) %>%
  mutate(metric = factor(metric, levels=c("520", "micro"))) %>%
  ggplot(aes(x=as.factor(hour), y=percent.calcified, colour= metric))+ geom_boxplot()+
  geom_point(aes(colour = metric),  
             size = 4, 
             position = position_jitterdodge())+
  theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.2,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  scale_color_manual(name = "Measurement Type", labels = c("micro" = "Microscopy", "520" = "Flow Cytometry (520nm)"), values = c("micro" = "#44AA99", "520" = "#999933"))+
  scale_y_continuous(name = "Percentage of calcein positive cells", limits = c(0,10), breaks = c(0,2,4,6,8,10))+
  scale_x_discrete(name = "Hour")

ggarrange(other_hosts, g_ocean, align = 'h', ncol = 2, nrow = 2, common.legend = FALSE, legend = "bottom")

 
#bleached liths######################################################################################################################################################################################
#plotting the data for just the bleached experiments
bleachedliths <- subset(adsorption, experiment.name == "Organic Removal and Calcein" | experiment.name == "UNCW_beast", select = c(1:21))
bleachedliths <- subset(bleachedliths, treatment == "control" | treatment == "374 liths" | treatment == "calcein"  | treatment == "bleached 374 liths", select = c(1:21))

#filtering out all of the timepoints greater than 24 hrs
bleachedliths_24hrs <- filter(bleachedliths, hour < 40)

bleachedliths_24hrs %>%
  arrange(percent.calcified) %>%
  mutate(liths = factor(liths, levels=c("control", "cells+liths", "cells+liths_bleach"))) %>%
ggplot(aes(x=as.factor(hour), y=percent.calcified, colour= liths))+ geom_boxplot()+
  geom_point(aes(colour = liths, shape = experiment.name),  
             size = 4, 
             position = position_jitterdodge())+  
  theme_bw()+  
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+  
  scale_color_manual(name = "Treatment", labels = c("control" = "Naked", "cells+liths" = "Cells + Liths", "cells+liths_bleach" = "Cells + Oxidized Liths"),  
                     values = c("control" = "#E69F00", "cells+liths" = "#0072B2", "cells+liths_bleach" = "#D55E00"))+  
  scale_y_continuous(name = "Percent calcified cells", limits = c(0,100), breaks = c(0,25,50,75))+  
  scale_x_discrete(name = "Hour")+  
  scale_shape_discrete(name = "Experiment", labels = c("Organic Removal and Calcein" = "Oxy Liths 1", "UNCW_beast" = "Oxy Liths 2"))

#histogram of the data distribution
hist(bleachedliths$percent.calcified)
ggplot(bleachedliths, aes(percent.calcified)) +
  geom_histogram(fill = "white", color = "grey30")+facet_wrap(~liths)


#subsetting the data into two time points
t0_bleach <- subset(bleachedliths, hour == "0", select = c(1:21))
ggplot(t0_bleach, aes(percent.calcified)) +
  geom_histogram(fill = "white", color = "grey30", binwidth = 0.5)
t24_bleach <- subset(bleachedliths, hour == "24", select = c(1:21))
ggplot(t24_bleach, aes(percent.calcified)) +
  geom_histogram(fill = "white", color = "grey30", binwidth = 0.5)


#data is not normally distributed therefore, will do a kruskal-wallis
#just the 0hr timepoint
krusk_t0 <- kruskal.test(percent.calcified~liths, data = t0_bleach)
#just the 24hr timepoint
krusk_t24 <- kruskal.test(percent.calcified~liths, data = t24_bleach)

#post-hoc dunn's test with bonferroni pvalue adj for 0hr
DT_t0 = dunnTest(percent.calcified ~ liths,
              data=t0_bleach,
              method="bonferroni")
DT_t0 = DT_t0$res
DT_t0
cldList(comparison = DT_t0$Comparison,
        p.value    = DT_t0$P.adj,
        threshold  = 0.05)
DT_t0letters <- cldList(P.adj ~ Comparison, data = DT_t0, threshold = 0.05)

#post-hoc dunn's test with bonferroni pvalue adj for 24hr
DT_t24 = dunnTest(percent.calcified ~ liths,
                 data=t24_bleach,
                 method="bonferroni")
DT_t24 = DT_t24$res
DT_t24
DT_t24letters <- cldList(P.adj ~ Comparison, data = DT_t24, threshold = 0.05)
DT_t24letters

#looking at the summary statistics
bleached_stats <- summarySE(bleachedliths_24hrs, measurevar = "percent.calcified",  
          groupvars = c("hour","liths"))


#stress markers ######################################################################################################################################################################################
#plotting the growth data from the experiment with the stress markers
stress <- subset(lithspikestress, treatment == "control" | treatment == "cells and liths", select = 1:23)

growthcurve <- summarySE(stress, measurevar = "cell.concentration", groupvars = c("treatment","time"))
#plotting the cell densities
#rearranges the order or the treatments so the control appears first
cells <- growthcurve %>% #using a pipe here
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control", "cells and liths"))) %>%

ggplot(aes(x=time, y=cell.concentration, colour = treatment))+  #with the pipe you don't need df
  geom_point(size = 3)+
  theme_bw()+theme(legend.text = element_text(lineheight = 2),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),  
                   axis.text.x = element_blank(),  
                   axis.title.x = element_blank())+  
  theme_classic()+
  geom_errorbar(aes(ymin=cell.concentration-se, ymax=cell.concentration+se),
                width= 5, size=0.75)+
  geom_line(aes(colour = treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96))+
  scale_y_continuous(name = bquote('Cell concentration ('*cells ~ mL^-1 ~ x10^6*')'),  
                     breaks = c(0,500000,1000000,1500000,2000000),  
                     labels = c(0,0.5,1,1.5,2))+
  labs(color = "Treatment")+
  scale_colour_manual(labels = c("Control",
                                 "Non-calcified cells\n+coccoliths"),
                      values = c("control" = "#E69F00",
                                 "cells and liths" = "#0072B2"))

#plotting the stress markers
#ros
#what's basal stress like
ros_data <- read.csv("growth data.csv")
#removing all of the strains that aren't CCMP374
ros_data <- filter(ros_data, strain == "CCMP374")
 
#removing all of the NA ros values
ros_data <- ros_data[!is.na(ros_data$normalized.ros),]
ros_data <- filter(ros_data, time < 100)

#looking for normality
ggplot(ros_data, aes(normalized.ros)) +
  geom_histogram(fill = "white", color = "grey30")

CI_ros <- wilcox.test(ros_data$normalized.ros,
            alternative="two.sided",
            correct=TRUE,
            conf.int=TRUE,
            conf.level=0.95)
median(ros_data$normalized.ros)


#95% confidence interval
upperros <- log10(447.4)
medianros <- log10(348.5832) 
lowerros <- log10(288.5)


#looking where stressed infected cells might fall
stressed_cells <- read.csv("Lith Spike Infection II 20181308.csv")
stressed_cells <- subset(stressed_cells, treatment == "cells liths virus" | treatment == "cells and virus", select = c(1:23))
stressed_cells <- filter(stressed_cells, time > 50)
ggplot(stressed_cells, aes(normalized.ros)) +
  geom_histogram(fill = "white", color = "grey30",binwidth = 100)

high_ros <- log10(median(stressed_cells$normalized.ros))

ros <- stress %>%
  arrange(normalized.ros) %>%
  mutate(treatment = factor(treatment, levels=c("control", "cells and liths"))) %>%

ggplot(aes(x= as.factor(time), y=log10(normalized.ros), colour = treatment))+  
  geom_boxplot()+
  geom_point(aes(colour = treatment),
             size = 4,
             position = position_jitterdodge())+  
  geom_hline(yintercept = medianros)+  
  geom_hline(yintercept = upperros, colour = "black", linetype = "dashed")+  
  geom_hline(yintercept = lowerros, colour = "black", linetype = "dashed")+
  geom_hline(yintercept = high_ros, colour = "red", linetype = "dashed")+
  theme_bw()+theme(legend.position = "none",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  theme_classic()+
  scale_x_discrete(name = "Hour")+
  scale_y_continuous(name = "Reactive oxygen species\n(median fluorescence units; log10")+
  scale_colour_manual(labels = c("Control",
                                 "Non-calcified cells\n+coccoliths"),
                      values = c("control" = "#E69F00",
                                 "cells and liths" = "#0072B2"))+
  scale_fill_manual(values = c("#E69F00","#0072B2"))

#sytox
#what's basal sytox like
sytox_data <- read.csv("growth data.csv")
sytox_data <- filter(sytox_data, strain == "CCMP374")

#removing all of the NA ros values
sytox_data <- sytox_data[!is.na(sytox_data$normalized.sytox),]
sytox_data <- filter(sytox_data, time < 100)


ggplot(combined_sytox_data, aes(normalized.sytox)) +
  geom_histogram(fill = "white", color = "grey30") #need to use the median, not normal

#95% confidence interval
CI_sytox <- wilcox.test(sytox_data$normalized.sytox,
                  alternative="two.sided",
                  correct=TRUE,
                  conf.int=TRUE,
                  conf.level=0.95)
uppersytox <- 5.850005
mediansytox <- 4.799947
lowersytox <- 3.999997


#average high sytox values
high_sytox <- median(stressed_cells$normalized.sytox)

death <- stress %>%
  arrange(normalized.sytox) %>%
  mutate(treatment = factor(treatment, levels=c("control", "cells and liths"))) %>%
  
  ggplot(aes(x= as.factor(time), y=normalized.sytox, colour = treatment))+  
  geom_point(aes(colour = treatment),
             size = 4,
             position = position_jitterdodge())+ 
  geom_boxplot()+  
  geom_hline(yintercept = mediansytox)+  
  geom_hline(yintercept = uppersytox, colour = "black", linetype = "dashed")+  
  geom_hline(yintercept = lowersytox, colour = "black", linetype = "dashed")+  
  geom_hline(yintercept = high_sytox, colour = "red", linetype = "dashed")+
  theme_bw()+theme(legend.position = "none",  
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),  
                   axis.text.x = element_blank(),  
                   axis.title.x = element_blank())+  
  theme_classic()+
  scale_x_discrete(name = "Hour")+
  scale_y_continuous(name = "Percent dead", limits = c(0,67))+ 
  scale_colour_manual(labels = c("Control",
                                 "Non-calcified cells\n+coccoliths"),
                      values = c("control" = "#E69F00",
                                 "cells and liths" = "#0072B2"))+
  scale_fill_manual(values = c("#E69F00","#0072B2"))

#daf
#looking at how this compares to average DAF in other experiments
daf_data <- read.csv("growth data.csv")
daf_data <- filter(daf_data, strain == "CCMP374")

#removing all of the NA ros values
daf_data <- daf_data[!is.na(daf_data$normalized.ros),]
daf_data <- filter(daf_data, time < 100)


#median
mediandaf <- log10(median(daf_data$normalized.daf)) #685.2319 is the mean of the values
#95% confidence interval
CI_daf <- wilcox.test(daf_data$normalized.daf,
                      alternative="two.sided",
                      correct=TRUE,
                      conf.int=TRUE,
                      conf.level=0.95)
upperdaf <- log10(6974)
mediandaf <- log10(6114)
lowerdaf <- log10(5285)

high_daf <- log(median(stressed_cells$normalized.daf))

#plot
daf <- stress %>%
  arrange(normalized.daf) %>%
  mutate(treatment = factor(treatment, levels=c("control", "cells and liths"))) %>%
  
  ggplot(aes(x= as.factor(time), y=log10(normalized.daf), colour = treatment))+  
  geom_point(aes(colour = treatment),
             size = 4,
             position = position_jitterdodge())+ 
  geom_boxplot()+ 
  geom_hline(yintercept = mediandaf)+  
  geom_hline(yintercept = upperdaf, colour = "black", linetype = "dashed")+  
  geom_hline(yintercept = lowerdaf, colour = "black", linetype = "dashed")+  
  geom_hline(yintercept = high_daf, colour = "red", linetype = "dashed")+
  theme_bw()+theme(legend.position = "none",  
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  theme_classic()+
  scale_x_discrete(name = "Hour")+
  scale_y_continuous(name = "Intracellular nitric oxide\n(median fluorescence; log10)")+ 
  scale_colour_manual(labels = c("Control",
                                 "Non-calcified cells\n+coccoliths"),
                      values = c("control" = "#E69F00",
                                 "cells and liths" = "#0072B2"))+
  scale_fill_manual(values = c("#E69F00","#0072B2"))


#photophys######################################################################################################################################################################################
#for this to run you need to run all of the code for the PI datasets, see Lith_adsorption_photophysiology.R
#plotting just the photophysiology data


#combining plots
ggarrange(death, ros, daf, pmax, align = 'hv',  
          label.x = 0, 
          nrow = 2, ncol = 2, hjust = -0.35, vjust = 26, 
          common.legend = TRUE, legend = "bottom")
#morphotypes######################################################################################################################################################################################
#plotting the morphotypes
morphotypes <- subset(adsorption, experiment.name == "Morphotype Screen" | experiment.name == "ViroLiths" | experiment.name == "G_Oceanica Lith Adsorption", select = 1:29)
morphotypes <- subset(morphotypes, treatment == "control " | treatment == "type O" | treatment == "type Over A" | treatment == "type R" | treatment == "type A" | treatment == "g_oceanica  liths", select = 1:29)

morpho <- morphotypes %>% #using a pipe here
  arrange(percent.calcified) %>%
  mutate(treatment = factor(treatment, levels=c("control ",  
                                                "type A",
                                                "type O",  
                                                "type Over A",  
                                                "type R", 
                                                "g_oceanica  liths"))) %>%
  
  ggplot(aes(x=as.factor(hour), y=percent.calcified, colour=treatment))+  
  geom_boxplot()+
geom_point(aes(colour = treatment),
           size = 4,
           position = position_jitterdodge())+theme_bw()+
  theme(legend.position = c(0.2,1), legend.text = element_text(lineheight = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
    theme_classic()+
scale_x_discrete(name = "Hour")+
scale_y_continuous(name = "Percent calcified cells", limits = c(0,100), breaks = c(0,25,50,75))+
scale_colour_manual("Lith Morphotype", values = c("control " = "#E69F00",
                                                  "type A" = "#CC6677",
                                                  "type O" = "#117733",
                                                  "type Over A" = "#332288",
                                                  "type R" = "#44AA99",
                                                  "g_oceanica  liths" = "#999933" ),
                    labels = c("control " = "Naked","type A" = "Type A","type O" = "Type O","type Over A" = "Type Over A","type R" = "Type R","g_oceanica  liths" = "G. Oceanica"))

morpho_48hr <- filter(morphotypes, hour < 49)

morpho_48hr %>%
  arrange(percent.calcified) %>%
  mutate(treatment = factor(treatment, levels=c("control ",  
                                                "type A",
                                                "type O",  
                                                "type Over A",  
                                                "type R", 
                                                "g_oceanica  liths"))) %>%
  ggplot(aes(x=as.factor(hour), y=percent.calcified, colour=treatment))+  
  geom_boxplot()+
  geom_point(aes(colour = treatment),
             size = 4,
             position = position_jitterdodge())+theme_bw()+
  theme(legend.position = c(0.2,1), legend.text = element_text(lineheight = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_classic()+
  scale_x_discrete(name = "Hour")+
  scale_y_continuous(name = "Percent calcified cells", limits = c(0,100), breaks = c(0,25,50,75,100))+
  scale_colour_manual("Lith Morphotype", values = c("control " = "#E69F00",
                                                    "type A" = "#CC6677",
                                                    "type O" = "#117733",
                                                    "type Over A" = "#332288",
                                                    "type R" = "#44AA99",  
                                                    "g_oceanica  liths" = "#999933"),
                      labels = c("control " = "Naked","type A" = "Type A","type O" = "Type O","type Over A" = "Type Over A","type R" = "Type R", "g_oceanica  liths" = "G. oceanica"))


morpho_nocontrol <- filter(morphotypes, liths != "control")
qqnorm(morphotypes$percent.calcified, pch = 1, frame = FALSE)
qqline(morphotypes$percent.calcified, col = "steelblue", lwd = 2)

model <- lm(percent.calcified~treatment, data = morphotypes)

morpho_48hr_stats <- filter(morphotypes, hour == "48")

ggplot(morpho_48hr_stats, aes(percent.calcified)) +
  geom_histogram(fill = "white", color = "grey30")+facet_wrap(~treatment)

KW_morpho <- kruskal.test(morpho_48hr_stats$percent.calcified, morpho_48hr_stats$treatment, method = "bonferroni")

DT_morpho <- dunnTest(percent.calcified ~ treatment, data = morpho_48hr_stats, method = "bonferroni")

DT_morpho_letters = DT_morpho$res

DT_morpho_letters <- cldList(P.adj ~ Comparison, data = DT_morpho_letters, threshold = 0.05)
DT_morpho_letters

morpho_sum <- summarySE(morphotypes, measurevar = "percent.calcified", groupvars = c("hour","date","liths","treatment"))


#linear regressions for change in lith concentrations over time########################  
#we are goingn to do the slopes for the lith dilution series first
lith_ads <- adsorption[!is.na(adsorption$lith.concentration), ]
lith_ads <- filter(lith_ads, lith.concentration != 0)
lith_ads <- filter(lith_ads, experiment.name != "1000yr old liths")
lith_ads <- filter(lith_ads, hour < 97)
lith_ads <- filter(lith_ads, host.strain != "1516")
lith_ads <- filter(lith_ads, experiment.name != "Organic Removal and Calcein")
lith_ads <- filter(lith_ads, experiment.name != "374 growth and adsorption")
lith_ads <- filter(lith_ads, treatment != "calcified")
lith_ads <- filter(lith_ads, number.of.liths.per.cell != "50:1 G. Oceanica; 50:1 Ehux")

ads_data_dilution <- subset(lith_ads, experiment.name == "Lith Adsorption Dilution Series" | experiment.name == "calcein_dilution_series", select = 1:31)

#removing unused columns
ads_data_dilution <- ads_data_dilution[,-(24:31),drop=FALSE] 

#make new IDs
lithdilution <- ads_data_dilution
lithdilution$ID <- as.factor(paste(lithdilution$number.of.liths.per.cell, lithdilution$replicate, lithdilution$target.density, lithdilution$treatment, sep="-"))
lithdilution$rephour <- as.factor(paste(lithdilution$replicate, lithdilution$hour, sep="-"))

#calc ratio loss and transforming to lognormal and finding the slope
lithdilution<- data.table (lithdilution, key= c("ID") )
lithdilution [, ratiolossln:= log((lith.concentration/(lith.concentration[match("0", hour)]))), by= c("ID")]
#making a new ID so I can plot out the specific treatments
lithdilution$plotID <- as.factor(paste(lithdilution$number.of.liths.per.cell, lithdilution$target.density, sep="-"))
#run the linear regressions
ratioloss_lith_slopes <- lithdilution %>% group_by(experiment.name, plotID, ID, treatment, liths, replicate, number.of.liths.per.cell) %>%
  do(model = lm(ratiolossln ~ hour, data = .)) %>%
  mutate(coef=coef(model)["hour"])
#need to make a new ID in order to average them by the correct treatment
ratioloss_lith_slopes$avgID <- as.factor(paste(ratioloss_lith_slopes$experiment.name, ratioloss_lith_slopes$treatment, sep="-"))
#convert to a data table first
ratioloss_lith_slopes <- data.table(ratioloss_lith_slopes, key = c("avgID"))
#averages the treatments
ratioloss_lith_slopes [, blanksavg:=   mean(coef), by = "avgID"]
#new ID for subtracting out the averaged blank
ratioloss_lith_slopes$rephour <- as.factor(paste(ratioloss_lith_slopes$experiment.name, ratioloss_lith_slopes$replicate, ratioloss_lith_slopes$hour, sep="-"))
#removes the blanks  
ratioloss_lith_slopes [, blanksubs:= (coef-(blanksavg[match ("blank", treatment)])), by=c("rephour")]
#removes the blanks from the data frame for plotting  
ratioloss_lith_slopes <- as.data.table(filter(ratioloss_lith_slopes, treatment != "blank"))


#plots as boxplot 
dilution <- ggplot(ratioloss_lith_slopes, aes(x= plotID, y=(blanksubs*24), colour = plotID))+  
  geom_boxplot()+
  geom_point(aes(colour = plotID), size = 5,
             position = position_jitterdodge())+  
  geom_hline(yintercept = 0)+
  theme_classic()+
  theme(legend.position = "none",  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        axis.text.x = element_text(angle = 45, hjust = 1))+  
  scale_y_reverse(name = bquote('Slope of remaining fraction of coccoliths ('*day^-1*')'))+  
  scale_x_discrete(name = "Target cell concentration and lith:cell ratio",labels = c("1000; Lith:Cell (100:1)",
                                                                                     "10000; Lith:Cell (100:1)",
                                                                                     "100000; Lith:Cell (100:1)",  
                                                                                     "1000; Lith:Cell (1000:1)"))+  
  scale_colour_manual(values = c("100-1000" = "#0072B2",
                                 "100-10000" = "#0072B2",
                                 "100-100000" = "#0072B2",
                                 "1000-1000" ="#0072B2"))

hist(ratioloss_lith_slopes$blanksubs)
shapiro.test(ratioloss_lith_slopes$blanksubs)
kruskal.test(blanksubs~plotID, data = ratioloss_lith_slopes)
dunnTest(blanksubs~plotID, data = ratioloss_lith_slopes)

ratioloss_lith_slopes$blanksubs_day <- ratioloss_lith_slopes$blanksubs * 24
ads_stats <- summarySE(ratioloss_lith_slopes, measurevar = "blanksubs_day", groupvars = c("plotID"))


#multiplying the lith cell ratio by the derived removal rate###########################################################################
#this removes everything except for the t0
only_t0 <- filter(ads_data_dilution, hour < 23)
#removes the blank
only_t0 <- filter(only_t0, treatment !="blank")

#analyzing the growth rates over 24hrs for all of the samples
only_growth_comp <- filter(ads_data_dilution, hour < 25)

only_growth_comp <- filter(only_growth_comp, treatment !="blank")

#removes the unused 1000:1 treatment
only_growth_comp <- filter(only_growth_comp, number.of.liths.per.cell !="1000")

only_growth_comp$ID <- as.factor(paste(only_growth_comp$number.of.liths.per.cell, only_growth_comp$replicate, only_t0$target.density, only_growth_comp$treatment, sep="-"))
only_growth_comp$plotID <- as.factor(paste(only_growth_comp$number.of.liths.per.cell, only_growth_comp$target.density, sep="-"))

only_growth_comp$cell.concentration <- as.numeric(as.character(only_growth_comp$cell.concentration))

#plotting the cell densities
#rearranges the order or the treatments so the control appears first
only_growth_comp %>% #using a pipe here
  arrange(cell.concentration) %>%
  mutate(plotID = factor(plotID, levels=c("100-1000", "100-10000", "100-100000"))) %>%
  
  ggplot(aes(x=hour, y=lith.concentration, colour = plotID, shape = experiment.name))+  #with the pipe you don't need df
  geom_point(size = 3)+
  theme_bw()+theme(legend.text = element_text(lineheight = 2),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),  
                   axis.text.x = element_blank(),  
                   axis.title.x = element_blank())+  
  theme_classic()+  
  facet_wrap(~plotID, scales = 'free')

#running the linear model
only_growth_slopes <- only_growth_comp %>% group_by(experiment.name, replicate, plotID, liths) %>%
  do(model = lm(log(cell.concentration) ~ hour, data = .)) %>%
  mutate(coef=coef(model)["hour"])

#converts it into per day
only_growth_slopes$growth_day <- only_growth_slopes$coef *24

#plotting the 24hr growth rates
ggplot(only_growth_slopes, aes(x= plotID, y=growth_day, colour = liths))+  
  geom_boxplot()+  
  geom_point(aes(colour = liths, shape = experiment.name), size = 5,
             position = position_jitterdodge())+  
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_discrete(name = "Dilution experimental treatments", labels = c("100-1000" = "1000 (100:1)",  
                                                                         "100-10000" = "10000 (100:1)",  
                                                                         "100-100000" = "100000 (100:1)"))
scale_y_continuous(name = bquote('Specific Growth Rate ('*mu ~ day^-1*')'))+
  scale_colour_manual(values = c("100-1000" = "#0072B2",
                                 "100-10000" = "#0072B2",
                                 "100-100000" = "#0072B2"))

#removes the other experiment based on poor growth rates
only_t0 <-filter(only_t0, experiment.name  != "Lith Adsorption Dilution Series")

#making an new ID to calculate the removal rates
only_t0$ID <- as.factor(paste(only_t0$number.of.liths.per.cell, only_t0$replicate, only_t0$target.density, only_t0$treatment, sep="-"))

#making a new ID for plotting purposes
only_t0$plotID <- as.factor(paste(only_t0$number.of.liths.per.cell, only_t0$target.density, sep="-"))
only_t0$rephour <- as.factor(paste(only_t0$experiment.name, only_t0$target.density, sep="-"))

#removes the unused 1000:1 treatment
only_t0 <- filter(only_t0, plotID !="1000-1000")

#converting the concentrations into numeric values
only_t0$lith.concentration <- as.numeric(as.character(only_t0$lith.concentration))
only_t0$cell.concentration <- as.numeric(as.character(only_t0$cell.concentration))


require(data.table)
#combines the necessary data grames
setDT(only_t0)[setDT(ratioloss_lith_slopes), blanksubs_day := i.blanksubs_day, on=c("ID")]
only_t0$blanksubs_day_abs <- abs(only_t0$blanksubs_day)

#calculating the number of liths per cell at 0hr
only_t0$liths_per_cell <- only_t0$lith.concentration/only_t0$cell.concentration

#multiplying by the slopes derived from the ratio of removal for the lith concentrations for each replicate
only_t0$removal_rate <- only_t0$liths_per_cell*only_t0$blanksubs_day_abs

#plotting the removal rates
plot1 <- ggplot(only_t0,aes(x=plotID,y=removal_rate, colour = plotID))+  
  geom_boxplot()+  
  geom_point(aes(colour = plotID, shape = experiment.name), size = 3, position = position_jitterdodge())+  
  theme_classic()+theme(legend.text = element_text(lineheight = 2),  
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())
  # scale_x_discrete(name = "Dilution experimental treatments",  
  #                  labels = c("100-1000" = "1000 (100:1)",  
  #                             "100-10000" = "10000 (100:1)",  
  #                             "100-100000" = "100000 (100:1)"))+  
  # scale_y_continuous(name = bquote('Number of liths removed'*~cell^-1 ~day^-1*''), limits = c(0,8), breaks = c(0,2,4,6,8))+  
  # scale_colour_manual(values = c("100-1000" = "#0072B2",
  #                                "100-10000" = "#0072B2",
  #                                "100-100000" = "#0072B2"))
summaryBy(cell.concentration ~ rephour, data = only_t0, 
          FUN = list(mean, max, min, median, sd))

#here's the encounter rate calculatinos
#looking at encounter rates
beta_DS = 8.990647e-06 #the calculated beta for differential settling for naked cells and liths
#note that the beta_DS value was updated 08_25_2022 to reflect the change in coccolith sinking speed
#we used 0.205 m d-1 based on the Zhang et al., 2018 paper
beta_BM = 8.7734970414853e-7 #the calculated beta for brownian motion for naked cells and liths

beta_BM_DS = beta_DS + beta_BM

only_t0$lith_encounter_rate <- (only_t0$lith.concentration * beta_BM_DS)

plot2 <- ggplot(only_t0,aes(x=plotID,y=log10(lith_encounter_rate), colour = plotID))+  
  geom_boxplot()+  
  geom_point(aes(colour = plotID), size = 3, position = position_jitterdodge())+  
  theme_classic()+theme(legend.text = element_text(lineheight = 2),  
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())+
  scale_x_discrete(name = "Dilution experimental treatments",  
                   labels = c("100-1000" = "1000 (100:1)",  
                              "100-10000" = "10000 (100:1)",  
                              "100-100000" = "100000 (100:1)"))+  
  scale_y_continuous(name = bquote('Number of liths encountered'*~ml^-1 ~day^-1*'; log10'), limits = c(0,2))+  
  scale_colour_manual(values = c("100-1000" = "#0072B2",
                                 "100-10000" = "#0072B2",
                                 "100-100000" = "#0072B2"))

enc_stats <- summarySE(only_t0, measurevar = "lith_encounter_rate", groupvars = c("plotID"))

#this calculates the adsorption efficiency
#this was the old way, might not be correct based on one of the reviewers comments
only_t0$adsorption_efficiency <- (only_t0$removal_rate / only_t0$lith_encounter_rate)*100

#tried to change the calculation based on Heidi's suggestions
# only_t0$adsorption_rate <- (only_t0$blanksubs_day_abs/40)
# only_t0$encounters <- (only_t0$cell.concentration*only_t0$lith.concentration)*beta_BM_DS
# 
# only_t0$adsorption_efficiency <- only_t0$adsorption_rate/only_t0$encounters



plot3 <-  only_t0 %>% filter(adsorption_efficiency<100) %>% ggplot(aes(x=plotID,y=adsorption_efficiency, colour = plotID))+  
  geom_boxplot()+  
  geom_point(aes(colour = plotID), size = 3, position = position_jitterdodge())+  
  theme_classic()+theme(legend.text = element_text(lineheight = 2),  
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())+
  scale_x_discrete(name = "Dilution experimental treatments", labels = c("100-1000" = "1000 (100:1)",  
                                                                         "100-10000" = "10000 (100:1)",  
                                                                         "100-100000" = "100000 (100:1)"))+  
  scale_y_continuous(name = bquote('Adsorption efficiency %'*~cell^-1*''), limits = c(0,80), breaks = c(0,20,40,60,80))+  #units were corrected to % cell-1
  scale_colour_manual(values = c("100-1000" = "#0072B2",
                                 "100-10000" = "#0072B2",
                                 "100-100000" = "#0072B2"))

ads_eff_stats <- summarySE(only_t0, measurevar = "adsorption_efficiency", groupvars = c("plotID"))

ggarrange(plot1, plot2, plot3, align = 'hv', ncol = 2, nrow = 2, common.legend = FALSE, legend = 'none')

only_t0$labels <- rep (c("three", "four", "five"), each = 3, length = 9)
ggplot(only_t0, aes(x = removal_rate)) +
  geom_histogram(aes(color = labels, fill = labels), 
                 position = "identity", bins = 4, alpha = 0.4) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#0072B2")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#0072B2"))
hist(only_t0$removal_rate)
shapiro.test(only_t0$removal_rate) #data is greater than 0.05, normally distributed

aov_removal_rates <- aov(removal_rate~labels, data = only_t0)  

tukey_removal_rates <- TukeyHSD(aov_removal_rates, "labels", ordered = TRUE)
tukey_removal_rates <- HSD.test(aov_removal_rates, trt = 'labels')

#just in case I want to export the csv file
# only_t0 <- only_t0[,-(1:29),drop=FALSE]
# only_t0 <- only_t0[,-(19:20),drop=FALSE]
# 
# write.csv(only_t0, "figure_4_data.csv")
#viroliths########################################################################################################
#virolith experiment for Figure 4
#plotting bar graph for viroliths
lithbath <- subset(adsorption, experiment.name == "ViroLiths" | experiment.name == "viroliths 2", select = c(1:31))
lithbath <- dplyr::filter(lithbath, treatment != "type A")

#plotting the viral lith experiment
#summarizing the cells concentrations
lithbath$cell.concentration <- as.numeric(as.character(lithbath$cell.concentration))
lithbathcells <- summarySE(lithbath,measurevar = "cell.concentration",  
                           groupvars = c("experiment.name", "hour","treatment","liths"))
#summarizing the virus concentrations
lithbath$viral.concentration <- as.numeric(as.character(lithbath$viral.concentration))
lithbathvirus <- summarySE(lithbath,measurevar = "viral.concentration",  
                           groupvars = c("experiment.name", "hour","treatment","liths")) %>%  
  dplyr::filter(liths != "cells+liths")

#plotting percent calcified
lithbath$hour <- as.numeric(lithbath$hour)
virolith1 <- lithbath %>% #using a pipe here
  filter(hour < 49) %>%  
  filter(experiment.name != "viroliths 2") %>%
  arrange(percent.calcified) %>%
  mutate(treatment = factor(treatment, levels=c("control ",  
                                        "374 liths 50:1",  
                                        "374 viroliths"))) %>%  
  ggplot(aes(x=as.factor(hour),y=percent.calcified, colour=treatment))+  
  geom_boxplot()+
  geom_point(aes(colour = treatment),
             size = 4,
             position = position_jitterdodge())+
  theme_bw()+
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_classic()+
  scale_colour_manual("Treatment",  
                      labels = c("Control","CCMP374 + Liths","CCMP374 + VirusLiths"),  
                      values = c("control " = "#E69F00",
                                 "374 viroliths" = "#56B4E9",
                                 "374 liths 50:1" = "#0072B2"))+
  scale_y_continuous(name = "Percent calcified cells (Side scatter; SSC)", limits = c(0,75), breaks = c(0,25,50,75))+
  scale_x_discrete(name = "Hour", limits = c("0","24","48","72","96"))



#plotting the growth curve 
# lithbathcells %>% #using a pipe here
#   arrange(cell.concentration) %>%  
#   filter(experiment.name != "viroliths 2") %>%
#   mutate(liths = factor(liths, levels=c("control",  
#                                         "374 liths 50:1",  
#                                         "374 viroliths"))) %>%
#   ggplot(aes(x=hour, y=cell.concentration, colour=liths))+  
#   geom_point(size=5)+theme_bw()+  
#   theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
#         panel.grid.major = element_blank(),  
#         panel.grid.minor = element_blank())+  
#   theme_classic()+
#   geom_errorbar(aes(ymin=cell.concentration-se, ymax=cell.concentration+se), width= 5, size=0.75)+  
#   geom_line(aes(colour=liths,group=liths))+
#   scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96))+  
#   scale_y_continuous(name = bquote('Cell concentration ('*cells ~ mL^-1 ~ x10^6*')'),  
#                      breaks = c(0,1000000,2000000,3000000,4000000),  
#                      labels = c(0,1,2,3,4))+
#   scale_colour_manual("Treatment",
#                       values = c("control" = "#E69F00",
#                                  "374 viroliths" = "#56B4E9",
#                                  "374 liths 50:1" = "#0072B2"),
#                       labels = c("Control","CCMP374 + Liths","CCMP374 + VirusLiths"))

#plotting the growth curve as a boxplot so you can see the data points
virolith2 <- lithbath %>% #using a pipe here
  arrange(cell.concentration) %>%  
  filter(experiment.name != "viroliths 2") %>%
  mutate(treatment = factor(treatment, levels=c("control ",  
                                        "374 liths 50:1",  
                                        "374 viroliths"))) %>%  
  ggplot(aes(x=as.factor(hour),y=cell.concentration, colour=treatment))+  
  geom_boxplot()+
  geom_point(aes(colour = treatment),
             size = 4,
             position = position_jitterdodge())+  
  theme_bw()+
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_classic()+
  scale_colour_manual("Treatment",
                      values = c("control " = "#E69F00",
                                 "374 viroliths" = "#56B4E9",
                                 "374 liths 50:1" = "#0072B2"),
                      labels = c("Control","CCMP374 + Liths","CCMP374 + VirusLiths"))+
  scale_x_discrete(name = "Hour")+  
  scale_y_continuous(name = bquote('Cell concentration ('*cells ~ mL^-1 ~ x10^6*')'),  
                     breaks = c(0,1000000,2000000,3000000,4000000),  
                     labels = c(0,1,2,3,4))+  
  stat_summary(fun = median,
               geom = "line", size = 1,
               aes(group = treatment), position = position_jitterdodge())


#plotting virus data
# lithbathvirus %>% #using a pipe here
#   arrange(viral.concentration) %>%
#   mutate(liths = factor(liths, levels=c("control","viroliths"))) %>%
#   ggplot(aes(x=hour, y=viral.concentration, colour=liths))+  
#   geom_point(size=5)+theme_bw()+  
#   theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
#         panel.grid.major = element_blank(),  
#         panel.grid.minor = element_blank())+  
#   theme_classic()+
#   geom_errorbar(aes(ymin=viral.concentration-se, ymax=viral.concentration+se), width= 5, size=0.75)+  
#   geom_line(aes(colour=liths,group=liths))+
#   scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96))+  
#   scale_y_continuous(name = bquote('Viral concentration ('*virus ~ mL^-1 ~ x10^8*')'),  
#                      breaks = c(0,25000000,50000000,75000000,100000000),  
#                      labels = c(0,0.25,0.50,0.75,1))+
#   scale_colour_manual("Treatment",
#                       values = c("control" = "#E69F00",
#                                  "viroliths" = "#56B4E9"),
#                       labels = c("Control","CCMP374 + VirusLiths"))

#plotting as a boxplot
virolith3 <- lithbath %>% #using a pipe here
  filter(experiment.name != "viroliths 2") %>%  
  filter(treatment != "374 liths 50:1") %>%
  arrange(viral.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control ",
                                                "374 viroliths"))) %>%
  ggplot(aes(x=as.factor(hour),y=viral.concentration, colour=treatment))+  
  geom_boxplot()+
  geom_point(aes(colour = treatment),
             size = 4,
             position = position_jitterdodge())+  
  theme_bw()+
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_classic()+
  scale_colour_manual("Treatment",
                       values = c("control " = "#E69F00",
                                  "374 viroliths" = "#56B4E9"),
                       labels = c("Control","CCMP374 + VirusLiths"))+
  scale_x_discrete(name = "Hour")+  
  scale_y_continuous(name = bquote('Viral concentration ('*virus ~ mL^-1 ~ x10^8*')'),  
                     breaks = c(0,25000000,50000000,75000000,100000000),  
                     labels = c(0,0.25,0.50,0.75,1))+
  stat_summary(fun = median,
               geom = "line", size = 1,
               aes(group = treatment), position = position_jitterdodge())

ggarrange(virolith1, virolith2, virolith3, align = 'hv', label.x = 0, 
          nrow = 2, ncol = 2, hjust = -0.35, vjust = 26, 
          common.legend = TRUE, legend = "bottom")


#the next part is for the supplementalary information and has been incorporated into Extended Data Fig. 13
#for this experiment we compared the infection dynamics of viroliths with virolith free supernatant containing only desorbed free viruses

#this also includes the total viruses across the different steps involved with generating the viroliths but gave the ability to assess how many viruses were attached per coccolith
#need the virus concentrations
virus_wash <- read.csv("virus counts.csv")

#plotting the total viruses
virolith4 <- virus_wash %>%  
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
  scale_colour_manual(values = c("lysate" = "#000000",
                                 "3" = "#000000",
                                 "4" = "#000000",
                                 "5" = "#56B4E9",  
                                 "6" = "#009E73"))+
  scale_y_continuous(name = "Number of viruses; log10")+
  scale_x_discrete(name = "Virolith preperation step",  
                   labels = c("Lysate","Virolith\nresuspension","Virolith\nsupernatant",  
                              "Washed virolith\nresuspension", "Washed virolith\nsupernatant"),  
                   guide = guide_axis(angle = 45))


#plotting
virolith5 <- lithbath %>% #using a pipe here  
  filter(experiment.name == "viroliths 2") %>%  
  filter(hour != "120") %>%
  arrange(percent.calcified) %>%
  mutate(liths = factor(liths, levels=c("control",
                                        "cells+liths",
                                        "viroliths",  
                                        "supernatant"))) %>%
  ggplot(aes(x=as.factor(hour),y=percent.calcified, colour=liths))+  
  geom_boxplot()+
  geom_point(aes(colour = liths),
             size = 4,
             position = position_jitterdodge())+
  theme_bw()+
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_classic()+
  scale_colour_manual("Treatment",
                      labels = c("Control","Naked cells + coccoliths","Naked cells + viroliths", "Naked cells + free viruses"),
                      values = c("control" = "#E69F00",
                                 "viroliths" = "#56B4E9",
                                 "cells+liths" = "#0072B2",  
                                 "supernatant" = "#009E73"))+
  scale_y_continuous(name = "Percent calcified cells (Side scatter; SSC)", limits = c(0,30), breaks = c(0,5,10,15,20,25,30))+
  scale_x_discrete(name = "Hour")
  

#plotting the growth curve 
virolith6 <- lithbathcells %>% #using a pipe here  
  filter(experiment.name == "viroliths 2") %>%  
  filter(hour != "120") %>%
  arrange(cell.concentration) %>%
  mutate(liths = factor(liths, levels=c("control",
                                        "cells+liths",
                                        "viroliths",  
                                        "supernatant"))) %>%
  ggplot(aes(x=hour, y=cell.concentration, colour=liths))+
  geom_point(size=5)+theme_bw()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+
  geom_errorbar(aes(ymin=cell.concentration-se, ymax=cell.concentration+se), width= 5, size=0.75)+  
  geom_line(aes(colour=liths,group=liths))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96))+  
  scale_y_continuous(name = bquote('Cell concentration ('*cells ~ mL^-1 ~ x10^5*')'),  
                     breaks = c(0,300000,600000,900000),  
                     labels = c(0,3,6,9))+
  scale_colour_manual("Treatment",
                      labels = c("Control","Naked cells + coccoliths","Naked cells + viroliths", "Naked cells + free viruses"),
                      values = c("control" = "#E69F00",
                                 "viroliths" = "#56B4E9",
                                 "cells+liths" = "#0072B2",  
                                 "supernatant" = "#009E73"))
#plotting as a boxplot
lithbath %>% #using a pipe here
  filter(experiment.name == "viroliths 2") %>%  
  filter(hour < 120) %>%  
  arrange(cell.concentration) %>%
  mutate(liths = factor(liths, levels=c("control",
                                        "cells+liths",
                                        "viroliths",
                                        "supernatant"))) %>%
  ggplot(aes(x=as.factor(hour),y=cell.concentration, colour=liths))+  
  geom_boxplot()+
  geom_point(aes(colour = liths),
             size = 4,
             position = position_jitterdodge())+  
  theme_bw()+
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_classic()+
  scale_x_discrete(name = "Hour")+  
  scale_y_continuous(name = bquote('Cell concentration ('*cells ~ mL^-1 ~ x10^5*')'))+
  scale_colour_manual("Treatment",
                      # labels = c("Control","Naked cells + coccoliths","Naked cells + viroliths", "Naked cells + free viruses"),
                      values = c("control" = "#E69F00",
                                 "cells+liths" = "#0072B2",  
                                 "viroliths" = "#56B4E9",
                                 "supernatant" = "#009E73"))+
  stat_summary(fun = median,
               geom = "line", size = 1,
               aes(group = liths), position = position_jitterdodge())


#plotting virus data
virolith7 <- lithbathvirus %>% #using a pipe here  
  filter(experiment.name == "viroliths 2") %>%  
  filter(hour != "120") %>%
  arrange(viral.concentration) %>%
  mutate(liths = factor(liths, levels=c("control",
                                        "viroliths",  
                                        "supernatant"))) %>%
  ggplot(aes(x=hour, y=viral.concentration, colour=liths))+  
  geom_point(size=5)+theme_bw()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+
  geom_errorbar(aes(ymin=viral.concentration-se, ymax=viral.concentration+se), width= 5, size=0.75)+  
  geom_line(aes(colour=liths,group=liths))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96))+  
  scale_y_continuous(name = bquote('Viral concentration ('*virus ~ mL^-1 ~ x10^7*')'),  
                     breaks = c(0,10000000,20000000,30000000,40000000),  
                     labels = c(0,1,2,3,4))+
  scale_colour_manual("Treatment",
                    labels = c("Control","Naked cells + viroliths", "Naked cells + free viruses"),
                    values = c("control" = "#E69F00",
                               "viroliths" = "#56B4E9",
                               "supernatant" = "#009E73"))
#plotting as a boxplot
lithbath %>% #using a pipe here
  filter(experiment.name == "viroliths 2") %>%  
  filter(hour < 120) %>%  
  filter(liths != "cells+liths") %>%
arrange(viral.concentration) %>%
  mutate(liths = factor(liths, levels=c("control",
                                        "viroliths",  
                                        "supernatant"))) %>%
  ggplot(aes(x=as.factor(hour),y=viral.concentration, colour=liths))+  
  geom_boxplot()+
  geom_point(aes(colour = liths),
             size = 4,
             position = position_jitterdodge())+  
  theme_bw()+
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_classic()+
  scale_y_continuous(name = bquote('Viral concentration ('*virus ~ mL^-1 ~ x10^7*')'),  
                     breaks = c(0,10000000,20000000,30000000,40000000,50000000,60000000,70000000),  
                     labels = c(0,1,2,3,4,5,6,7))+
  scale_colour_manual("Treatment",
                      labels = c("Control","Naked cells + viroliths", "Naked cells + free viruses"),
                      values = c("control" = "#E69F00",
                                 "viroliths" = "#56B4E9",
                                 "supernatant" = "#009E73"))+
  scale_x_discrete(name = "Hour")+  
  stat_summary(fun = median,
               geom = "line", size = 1,
               aes(group = liths), position = position_jitterdodge())


ggarrange(virolith4, virolith5, virolith6, virolith7, fig13_cell_gam_output, fig13_virus_gam_output, align = 'hv', label.x = 0, 
          nrow = 3, ncol = 2, hjust = -0.35, vjust = 26, 
          common.legend = FALSE, legend = "none")


#infection with coccoliths###########################################################################################################################
#this was for sorted coccoliths
lith_infection$hour <- lith_infection$day*24
lith_infection %>% ggplot(aes(x=hour, y=log10(concentration), colour= treatment))+ geom_boxplot()+ geom_point(aes(colour = treatment), size = 3, position = position_jitterdodge())
same_MOI <- summarySE(lith_infection,measurevar = "concentration",  
                      groupvars = c("hour","treatment")) 
sorting <- same_MOI %>%  
  filter(treatment != "control") %>%  
  filter(hour <120) %>%  
  ggplot(aes(x=hour, y=concentration, colour= treatment, shape = treatment))+
  geom_point(size=5)+theme_bw()+   
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.ticks.length=unit(.25, "cm"),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.line = element_line(colour = 'black', size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymin=concentration-se, ymax=concentration+se), width= 5, size=0.75)+ 
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120))+
  scale_y_continuous(name = bquote('Cell concentration ('*cells ~ mL^-1*')'), limits = c(0,1200000), breaks=c(0,300000,600000,900000,1200000))+  
  scale_shape_manual(values = c(17,17))+
  scale_colour_manual("Treatment",
                      # labels = c("Naked cells + virus","Naked cells + coccoliths and virus"),
                      values = c("virus" = "#009E73",
                                 "liths and virus" = "#56B4E9"))

same_MOI_96h <- filter(lith_infection, hour == "96")
same_MOI_96h <- filter(same_MOI_96h, treatment != "control")
shapiro.test(same_MOI_96h$concentration) #pvalue is 0.1561, the data is normal will do a t-test
t.test(concentration ~ treatment, data = same_MOI_96h, var.equal = TRUE)
#at 96 hour viroliths and virus alone are significantly different with a pvalue of 0.005133

lithspike_infection %>%  
  filter(treatment != "cells and liths") %>% filter(treatment != "control") %>%  
  ggplot(aes(x=time, y=log10(cell.concentration), colour= treatment))+
  geom_point(size=5)+theme_bw()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+ 
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96))+
  scale_y_continuous(name = bquote('Cell concentration ('*cells ~ mL^-1*')'))+
  scale_colour_manual("Treatment",
                      labels = c("Naked cells + virus","Naked cells + coccoliths and virus"),
                      values = c("cells and virus" = "#E69F00",
                                 "cells liths virus" = "#56B4E9"))

infection <- lithspikestress %>% filter(treatment != "cells and liths") %>% filter(treatment != "control") %>% ggplot(aes(x=as.factor(time), y=log10(cell.concentration), colour= treatment))+ geom_boxplot()+ geom_point(aes(colour = treatment), size = 3, position = position_jitterdodge())

lithspike_infection <- summarySE(lithspikestress,measurevar = "cell.concentration",  
                                 groupvars = c("time","treatment")) 

lithspike_infection$time <- as.numeric(as.character(lithspikestress_infection$time))
infection <- lithspike_infection %>%  
  filter(treatment != "cells and liths") %>% filter(treatment != "control") %>%  
  ggplot(aes(x=time, y=log10(cell.concentration), colour= treatment))+
  geom_point(size=5)+theme_bw()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+ 
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96))+
  scale_y_continuous(name = bquote('Cell concentration ('*cells ~ mL^-1*')'))+
  scale_colour_manual("Treatment",
                      labels = c("Naked cells + virus","Naked cells + coccoliths and virus"),
                      values = c("cells and virus" = "#E69F00",
                                 "cells liths virus" = "#56B4E9"))

protection <- adsorption %>% subset(treatment == "cells+virus" | treatment == "lithcells+virus")%>% 
  ggplot(aes(x=as.factor(hour), y=log10(as.numeric(as.character(cell.concentration))), colour= treatment))+ geom_boxplot()+facet_wrap(~experiment.name)

lithspike_protection <- adsorption %>% subset(treatment == "cells+virus" | treatment == "lithcells+virus")%>% filter(experiment.name == "Lith Protection 2")
lithspike_protection$cell.concentration <- as.numeric(as.character(lithspike_protection$cell.concentration))
liths_protect <- summarySE(lithspike_protection,measurevar = "cell.concentration",  
                           groupvars = c("hour", "treatment")) 

protection <- liths_protect %>%  
  ggplot(aes(x=hour, y=log10(cell.concentration), colour= treatment))+
  geom_point(size=5)+theme_bw()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+ 
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,148,168))+
  scale_y_continuous(name = bquote('Cell concentration ('*cells ~ mL^-1*')'))+
  scale_colour_manual("Treatment",
                      labels = c("Naked cells + virus","Naked cells + coccoliths and virus"),
                      values = c("cells+virus" = "#E69F00",
                                 "lithcells+virus" = "#56B4E9"))

ggarrange(protection, infection, nrow = 1, ncol = 2, align = "h", legend = FALSE)

#comparing the previous experiment with coccoliths isolated via density centrifugation with coccolith:cell ratio of 50:1 and MOI of 5
virolith_p3 <- subset(adsorption, experiment.name == "viroliths 3" | experiment.name == "viroliths 4")
virolith_p3$cell.concentration <- as.numeric(as.character(virolith_p3$cell.concentration))

virolith_p3_cells <- summarySE(virolith_p3,measurevar = "cell.concentration",  
          groupvars = c("hour","treatment","moi")) 
virolith_p3_cells$moi <- as.character(as.factor(virolith_p3_cells$moi))

virolith_p3_liths <- summarySE(virolith_p3,measurevar = "lith.concentration",  
                               groupvars = c("hour","treatment","moi")) 
virolith_p3_liths$moi <- as.character(as.factor(virolith_p3_liths$moi))
virolith_p3_virus <- summarySE(virolith_p3,measurevar = "viral.concentration",  
                               groupvars = c("hour","treatment","moi"))
virolith_p3_virus$moi <- as.character(as.factor(virolith_p3_virus$moi))

cells <- virolith_p3_cells %>% filter(treatment != "control") %>% filter(treatment != "374 liths") %>%  
  filter(hour < 169) %>%
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control", "374 liths",
                                        "374 viroliths 5",
                                        "374 viroliths 0.1",
                                        "374 viroliths 0.01",
                                        "virus 5",
                                        "virus 0.1",
                                        "virus 0.01"))) %>%
  ggplot(aes(x=hour, y=cell.concentration, colour= treatment, shape = moi))+
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.ticks.length=unit(.25, "cm"),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.line = element_line(colour = 'black', size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymin=cell.concentration-se, ymax=cell.concentration+se), width= 5, size=0.75)+ 
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192))+
  scale_y_continuous(name = bquote('Cell concentration ('*cells ~ mL^-1*'; log10)'))+
  scale_colour_manual("Treatment",
                      # labels = c("Naked cells + coccoliths and virus","Naked cells + virus"),
                      values = c("control" = "#E69F00", "374 liths" = "#0072B2", "374 viroliths 5" = "#56B4E9",
                                 "374 viroliths 0.1" = "#56B4E9",
                                 "374 viroliths 0.01" = "#56B4E9",
                                 "virus 5" = "#009E73",  
                                 "virus 0.1" = "#009E73",
                                 "virus 0.01" = "#009E73"))

lith_infection$cell.concentration <- lith_infection$concentration
sorted_cells_96 <- filter(lith_infection, hour >90 & hour <119)
sorted_cells_96 <- filter(sorted_cells_96, treatment != "control")

sorted_cells_96$experiment.name <- rep("sorted_liths", length=6)
sorted_cells_96$moi <- rep("0.1", length=6) 


cells_96 <- subset(virolith_p3, treatment == "virus 5" | treatment == "374 viroliths 5")
cells_96 <- filter(cells_96, hour >90)
cells_168 <- virolith_p3 %>% filter(hour == "168") %>% filter(moi != "0")

cells1 <- merge(cells_96, cells_168, all = TRUE)

cells2 <- merge(cells1, sorted_cells_96, all = TRUE)


virolith_cells <- cells2 %>%
  arrange(cell.concentration) %>% 
  mutate(treatment = factor(treatment, levels=c("virus 5",
                                                "374 viroliths 5",
                                                "virus",  
                                                "liths and virus",
                                                "virus 0.1",
                                                "374 viroliths 0.1",
                                                "virus 0.01",
                                                "374 viroliths 0.01"))) %>%
  ggplot(aes(x=as.factor(treatment),y=log10(cell.concentration), colour=treatment))+  
  geom_boxplot()+
  geom_point(aes(colour = treatment),
             size = 4,
             position = position_jitterdodge())+theme_bw()+  
  # theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank())+
  theme_classic()+
  scale_y_continuous(name = bquote('Cell concentration ('*cells ~ mL^-1*'; log10)'))+
  scale_colour_manual("Treatment",
                      # labels = c("Naked cells + coccoliths and virus","Naked cells + virus"),
                      values = c( "374 viroliths 5"  = "#56B4E9",
                                  "virus 5" = "#009E73",
                                  "liths and virus" = "#56B4E9",
                                  "virus" = "#009E73",
                                 "374 viroliths 0.1" = "#56B4E9",
                                 "374 viroliths 0.01" = "#56B4E9",
                                 "virus 0.1" = "#009E73",
                                 "virus 0.01" = "#009E73"), guide = "none")



virolith_p3 %>% #using a pipe here
  filter(hour < 168) %>%
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control", "374 liths",
                                                "374 viroliths 5",
                                                "374 viroliths 0.1",
                                                "374 viroliths 0.01",
                                                "virus 5",
                                                "virus 0.1",
                                                "virus 0.01"))) %>%
  ggplot(aes(x=as.factor(hour),y=as.numeric(cell.concentration), colour=treatment, shape = as.factor(moi)))+  
  geom_boxplot()+
  geom_point(aes(colour = liths),
             size = 4,
             position = position_jitterdodge())+  
  theme_bw()+
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_classic()+
  scale_y_continuous(name = bquote('Cell concentration ('*cells ~ mL^-1*'; log10)'))+
  scale_colour_manual("Treatment",
                      # labels = c("Naked cells + coccoliths and virus","Naked cells + virus"),
                      values = c("control" = "#E69F00", "374 liths" = "#0072B2", "374 viroliths 5" = "#56B4E9",
                                 "374 viroliths 0.1" = "#56B4E9",
                                 "374 viroliths 0.01" = "#56B4E9",
                                 "virus 5" = "#009E73",  
                                 "virus 0.1" = "#009E73",
                                 "virus 0.01" = "#009E73"))
  scale_x_discrete(name = "Hour")+  
  stat_summary(fun = median,
               geom = "line", size = 1,
               aes(group = treatment), position = position_jitterdodge())


cells_5 <- subset(virolith_p3, treatment == "374 viroliths 5" | treatment == "virus 5")
cells_5 <- filter(cells_5, hour > 73 & hour < 97)
shapiro.test(cells_5$cell.concentration) #pvalue is greater than 0.05 (0.994), implying the data is normal and I'll use a t-test
t.test(cell.concentration ~ treatment, data = cells_5, var.equal = TRUE)
#viroliths and virus alone are very slightly significantly different at MOI 5- pvalue 0.04282

cells_0.1 <- subset(cells_168, treatment == "374 viroliths 0.1" | treatment == "virus 0.1")
cells_0.1 <- filter(cells_0.1, hour < 169)
shapiro.test(cells_0.1$cell.concentration) #pvalue is greater than 0.05 (0.272), implying the data is normal and I'll use a t-test
t.test(cell.concentration ~ treatment, data = cells_0.1, var.equal = TRUE)
#viroliths and virus alone are significantly different at 168 h (pvalue 0.00362)

cells_0.01 <- subset(cells_168, treatment == "374 viroliths 0.01" | treatment == "virus 0.01")
cells_0.01 <- filter(cells_0.01, hour < 169)
shapiro.test(cells_0.01$cell.concentration) #pvalue is slightly less than 0.05, and I'll use a Mann Whitney
wilcox.test(cell.concentration ~ treatment, data = cells_0.01, paired = FALSE) #not paired data so this is a Mann-Whitney for independent samples
#viroliths and virus alone are not significantly different at a ratio of 0.01- this is due to one replicate in the virus alone being different than the rest

liths <- virolith_p3_liths %>% filter(lith.concentration > 100000) %>% filter(hour < 169) %>% filter(treatment != "374 liths") %>%
  ggplot(aes(x=hour, y=lith.concentration, color = treatment, shape = moi))+
  geom_point(size=5)+theme_bw()+  
  # theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
  #       panel.grid.major = element_blank(),  
  #       panel.grid.minor = element_blank())+  
  theme_classic()+
  geom_errorbar(aes(ymin=lith.concentration-se, ymax=lith.concentration+se), width= 5, size=0.75)+ 
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168,192))+
  scale_y_continuous(name = bquote('Coccolith concentration ('*cells ~ mL^-1*')'))+
scale_colour_manual("Treatment", values = c("374 viroliths 5" = "#56B4E9",
                    "374 viroliths 0.1" = "#56B4E9",
                    "374 viroliths 0.01" = "#56B4E9"))

virus <- virolith_p3_virus %>% filter(treatment != "control") %>% filter(treatment != "374 liths") %>% filter(hour < 169) %>%
  arrange(viral.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control", "374 liths",
                                                "374 viroliths 5",
                                                "374 viroliths 0.1",
                                                "374 viroliths 0.01",
                                                "virus 5",
                                                "virus 0.1",
                                                "virus 0.01"))) %>%
  ggplot(aes(x=hour, y=viral.concentration, colour= treatment, shape= moi))+
  geom_point(size=5)+theme_bw()+  
  # theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
  #       panel.grid.major = element_blank(),  
  #       panel.grid.minor = element_blank())+  
  theme_classic()+
  geom_errorbar(aes(ymin=viral.concentration-se, ymax=viral.concentration+se), width= 5, size=0.75)+ 
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168,192))+
  scale_y_continuous(name = bquote('Virus concentration ('*cells ~ mL^-1*')'))+
scale_colour_manual("Treatment",
                    # labels = c("Naked cells + coccoliths and virus","Naked cells + virus"),
                    values = c("control" = "#E69F00", "374 liths" = "#0072B2", "374 viroliths 5" = "#56B4E9",
                               "374 viroliths 0.1" = "#56B4E9",
                               "374 viroliths 0.01" = "#56B4E9",
                               "virus 5" = "#009E73",  
                               "virus 0.1" = "#009E73",
                               "virus 0.01" = "#009E73"))

ggarrange(sorting, cells, nrow = 2, ncol = 2, align = "hv", legend = FALSE)

#enzymes###################################################################################################################
enzyme_data <- subset(adsorption, experiment.name == "Enzymes")

enzyme_data %>% filter(treatment != "calcein") %>% filter(hour < 49) %>%  
  arrange(percent.calcified) %>%
  mutate(treatment = factor(treatment, levels=c("control",  
                                        "374 liths",  
                                        "bleached 374 liths",  
                                        "374 amylase liths",  
                                        "374 proteinase liths",  
                                        "374 lipase liths"))) %>%
  ggplot(aes(x=as.factor(hour), y=percent.calcified, colour= treatment))+ geom_boxplot()+
  geom_point(aes(colour = treatment),  
             size = 4, 
             position = position_jitterdodge())+  
  theme_bw()+  
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+ 
  scale_color_manual(name = "Treatment", labels = c("control" = "Naked cells", "374 liths" = "Naked cells + coccoliths", "bleached 374 liths" = "Naked cells + oxidized coccoliths", "374 amylase liths" = "Naked cells + amylase coccoliths", "374 proteinase liths" = "Naked cells + proteinase coccoliths", "374 lipase liths" = "Naked cells + lipase coccoliths"), values = c("control" = "#88CCEE", "374 liths" = "#CC6677", "bleached 374 liths" = "#DDCC77", "374 amylase liths" = "#117733", "374 proteinase liths" = "#AA4499", "374 lipase liths" = "#44AA99"))+
  scale_y_continuous(name = "Percent calcified cells", limits = c(0,100), breaks = c(0,25,50,75,100))+
  scale_x_discrete(name = "Hour")

