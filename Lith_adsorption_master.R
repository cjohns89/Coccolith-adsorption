setwd("/Users/cjohns/Dropbox/Lith Adsorption//Spreadsheets")
adsorption <- read.csv("Lith_Adsorption_Mastersheet_20191011.csv")
lithspikestress <- read.csv("Lith Spike Infection II 20181308.csv")
measurements <- read.csv("percent_calcified_measurements.csv")
library(dplyr)
library(Rmisc)
library(dplyr)
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
#percent calcified via lith adsorption######################################################################################################################################################################################
#filtering everything that I've given 374 liths to
growth <- subset(adsorption, liths == "control" | liths == "cells+liths", select = c(1:21))
ggplot(growth, aes(x=hour, y=percent.calcified, colour = liths, shape = host.strain))+  
  geom_point(size = 5)+  
  scale_y_continuous(name = "Percent calcified")+
  scale_x_continuous(name = bquote('Cell concentration ('*cells ~ mL^-1*'); log10'))+
  labs(color = "Treatment")+
  scale_colour_manual(labels = c("Non-calcified cells\n+coccoliths", "Control"),
                      values = c("control" = "#E69F00",
                                 "cells+liths" = "#0072B2"))+  
  scale_shape_discrete(name = "Strain", labels = c("CCMP1516", "CCMP374", "Phaeocystis globosa"))


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


#uses percent calcified as detected by side scatter (SSC)
ssc <- ggplot(just374,aes(x=as.factor(hour), y=percent.calcified, colour= liths))+ geom_boxplot()+
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

just374_sum <- summarySE(just374, measurevar = "percent.calcified", groupvars = c("hour","liths"))

            
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
#making new IDs for the growth rates
just374$cell.concentration <- as.numeric(as.character((just374$cell.concentration))) #forcing a factor to a numeric can warp the actual numbers fo you need to change it to a character afterwards
just374$ID <- as.factor(paste(just374$liths, just374$replicate, just374$experiment.name, sep="-"))
just374$rephour <- as.factor(paste(just374$replicate, just374$hour, sep="-"))


#calc growth rates and transforming to lognormal
just374 <- data.table (just374, key= c("ID") )
just374 [, ratiolossln:= (log((cell.concentration/(cell.concentration[match("0", hour)])))/(hour - hour[match("0", hour)])), by= c("ID")]

ggplot(just374, aes(x=cell.concentration, y=log10(ratiolossln), colour = liths))+  
  geom_point(size = 5)

just374$ratiolossln <- just374$ratiolossln *24
is.na(just374) <- sapply(just374, is.infinite)



#plotting the specific growth rates after summarizing each replicate and finding the slope
just374$statID <- as.factor(paste(just374$liths, just374$replicate, just374$experiment.name, just374$treatment,sep="-"))

slopes <- just374 %>% group_by(experiment.name, statID, liths) %>%
  do(model = lm(log(cell.concentration) ~ hour, data = .)) %>%
  mutate(coef=coef(model)["hour"])

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

ggplot(slopes, aes(x=growth_day, fill=liths)) +
  geom_density(alpha=0.4)
ggplot(slopes, aes(x=liths, y = growth_day)) +
  geom_violin()

#summarizing the growth rates
growthrates <- group_by(just374, liths) %>%
  summarise(
    count = n(),
    mean = mean(ratiolossln, na.rm = TRUE),
    sd = sd(ratiolossln, na.rm = TRUE))

#stats
ggplot(slopes, aes(growth_day)) +
  geom_histogram(fill = "white", color = "grey30")+facet_wrap(~liths)
library(lawstat)
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

sum_stats_growth <- summarySE(slopes, measurevar = "growth_day",  
                                   groupvars = c("liths"))
#other hosts#######################################################################################################################################################################################
#this shows any host that wasn't CCMP374
otherhosts <- subset(growth, host.strain == "1516" | host.strain == "phaeocystis", select = c(1:21))
otherhosts <- filter(otherhosts, hour < 100)

ggplot(otherhosts,aes(x=as.factor(hour), y=percent.calcified, colour= liths))+ geom_boxplot()+
  geom_point(aes(colour = liths),  
             size = 4, 
             position = position_jitterdodge())+  
  facet_rep_wrap(~host.strain, scales = "fixed", repeat.tick.labels = FALSE)+
  theme_bw()+  
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+  
  scale_color_manual(name = "Treatment", labels = c("control" = "Naked", "cells+liths" = "Cells + Liths"), values = c("control" = "#E69F00", "cells+liths" = "#0072B2"))+  
  scale_y_continuous(name = "Percent calcified cells", limits = c(0,100), breaks = c(0,25,50,75,100))+  
  scale_x_discrete(name = "Hour")

sum_stats_other_hosts <- summarySE(otherhosts, measurevar = "percent.calcified",  
                                            groupvars = c("experiment.name","hour","liths"))
 
#bleached liths######################################################################################################################################################################################
#plotting the data for just the bleached experiments
bleachedliths <- subset(adsorption, experiment.name == "Organic Removal and Calcein" | experiment.name == "UNCW_beast", select = c(1:21))
bleachedliths <- subset(bleachedliths, treatment == "control" | treatment == "374 liths" | treatment == "calcein"  | treatment == "bleached 374 liths", select = c(1:21))

bleachedliths_24hrs <- filter(bleachedliths, hour < 40)

bleachedliths_24hrs %>%
  arrange(percent.calcified) %>%
  mutate(liths = factor(liths, levels=c("control", "cells+liths", "bleach"))) %>%
ggplot(aes(x=as.factor(hour), y=percent.calcified, colour= liths))+ geom_boxplot()+
  geom_point(aes(colour = liths, shape = experiment.name),  
             size = 4, 
             position = position_jitterdodge())+  
  theme_bw()+  
  theme(legend.position = c(0.1,0.88),legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+  
  scale_color_manual(name = "Treatment", labels = c("control" = "Naked", "cells+liths" = "Cells + Liths", "bleach" = "Cells + Oxidized Liths"), values = c("control" = "#E69F00", "cells+liths" = "#0072B2", "bleach" = "#D55E00"))+  
  scale_y_continuous(name = "Percent calcified cells", limits = c(0,100), breaks = c(0,25,50,75))+  
  scale_x_discrete(name = "Hour")+  
  scale_shape_discrete(name = "Experiment", labels = c("Organic Removal and Calcein" = "Oxy Liths 1", "UNCW_beast" = "Oxy Liths 2"))

hist(bleachedliths$percent.calcified)
ggplot(bleachedliths, aes(percent.calcified)) +
  geom_histogram(fill = "white", color = "grey30")+facet_wrap(~hour)


t0_bleach <- subset(bleachedliths, hour == "0", select = c(1:21))
ggplot(t0_bleach, aes(percent.calcified)) +
  geom_histogram(fill = "white", color = "grey30", binwidth = 0.5)
t24_bleach <- subset(bleachedliths, hour == "24", select = c(1:21))
ggplot(t24_bleach, aes(percent.calcified)) +
  geom_histogram(fill = "white", color = "grey30", binwidth = 0.5)

krusk_t0 <- kruskal.test(percent.calcified~liths, data = t0_bleach)
krusk_t24 <- kruskal.test(percent.calcified~liths, data = t24_bleach)

DT_t0 = dunnTest(percent.calcified ~ liths,
              data=t0_bleach,
              method="holm")

DT_t0 = DT_t0$res
DT_t0
cldList(comparison = DT_t0$Comparison,
        p.value    = DT_t0$P.adj,
        threshold  = 0.05)
DT_t0letters <- cldList(P.adj ~ Comparison, data = DT_t0, threshold = 0.05)


DT_t24 = dunnTest(percent.calcified ~ liths,
                 data=t24_bleach,
                 method="holm")
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
growth <- growthcurve %>% #using a pipe here
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
ros_data <- filter(ros_data, strain == "CCMP374")

#removing all of the NA ros values
ros_data <- ros_data[!is.na(ros_data$normalized.ros),]
ros_data <- filter(ros_data, time < 100)
combined_ros_data <- merge(ros_data, stress,   
                           by.x = c("date", "time" , "strain" , "normalized.ros"),  
                           by.y = c("date", "time" , "strain" , "normalized.ros"),  
                           all.x = TRUE,  
                           all.y= TRUE)
#looking for normality
ggplot(ros_data, aes(normalized.ros)) +
  geom_histogram(fill = "white", color = "grey30")

CI_ros <- wilcox.test(ros_data$normalized.ros,
            alternative="two.sided",
            correct=TRUE,
            conf.int=TRUE,
            conf.level=0.95)


medianros <- median(ros_data$normalized.ros)
#95% confidence interval

upperros <- log10(447.5)
lowerros <- log10(288.5)
medianros <- log10(305) 

#looking where stressed infected cells might fall
stressed_cells <- read.csv("Lith Spike Infection II 20181308.csv")
stressed_cells <- subset(stressed_cells, treatment == "cells liths virus" | treatment == "cells and virus", select = c(1:23))
stressed_cells <- filter(stressed_cells, time > 50)
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
# combined_sytox_data <- merge(sytox_data, stress,   
#                            by.x = c("date", "time" , "strain" , "normalized.sytox"),  
#                            by.y = c("date", "time" , "strain" , "normalized.sytox"),  
#                            all.x = TRUE,  
#                            all.y= TRUE)

ggplot(sytox_data, aes(normalized.sytox)) +
  geom_histogram(fill = "white", color = "grey30") #need to use the median, not normal

#median
mediansytox <- median(sytox_data$normalized.sytox)
#95% confidence interval
CI_sytox <- wilcox.test(sytox_data$normalized.sytox,
                  alternative="two.sided",
                  correct=TRUE,
                  conf.int=TRUE,
                  conf.level=0.95)
uppersytox <- 5.850005
lowersytox <- 3.999997

stress$deadcell_ml <- ((stress$normalized.sytox/100)*stress$cell.concentration)

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
# combined_daf_data <- merge(daf_data, stress,   
#                            by.x = c("date", "time" , "strain" , "normalized.daf"),  
#                            by.y = c("date", "time" , "strain" , "normalized.daf"),  
#                            all.x = TRUE,  
#                            all.y= TRUE)


#median
mediandaf <- log10(median(daf_data$normalized.daf)) #685.2319 is the mean of the values
#95% confidence interval
CI_daf <- wilcox.test(daf_data$normalized.daf,
                      alternative="two.sided",
                      correct=TRUE,
                      conf.int=TRUE,
                      conf.level=0.95)
upperdaf <- log10(6974)
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
#for this to run you need to run all of the code for the PI datasets
#plotting just the photophysiology data
pmax <- pidata_pmax %>%    
  filter(hour < 48) %>%
  arrange(pmax) %>%  
  mutate(culture = factor(culture, levels = c("Naked","Cells + Liths","Calcified"))) %>%  
  ggplot(aes(x=hour, y=pmax, color=culture))+  
  geom_boxplot()+  
  geom_point(aes(colour=culture), size = 4, position = position_jitterdodge())+
  theme_classic2()+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  scale_x_discrete(name= "Hour")+  
  scale_y_continuous(name = bquote(atop('Max photosynthetic rate', '(Pmax; '*electrons ~ s^-1 ~ PSII^-1*')')))+  
  scale_color_manual(name = "Type", values = c("Naked" = "#E69F00", "Cells + Liths" = "#0072B2", "Calcified" = "#CC79A7"))+  
  stat_compare_means(comparisons = )

hist(pidata_pmax$pmax)
shapiro.test(pidata_pmax$pmax)#the data is normally distributed, going to do an anova for each group

t0 <- filter(pidata_pmax, hour < 24)
t24 <- filter(pidata_pmax, hour >23 & hour < 48)

aov_t0 <- aov(pmax~culture, data = t0)
aov_t24 <- aov(pmax~culture, data = t24)

Tukey_t0 <- TukeyHSD(aov_t0, "culture", ordered = TRUE)

Tukey_t24 <- TukeyHSD(aov_t24, "culture", ordered = TRUE)
library("agricolae") #need to install questionr package prior, otherwise it will fail
library(multcomp)
library(multcompView)
t0_letters <- multcompLetters(extract_p(Tukey_t0$culture))
t24_letters <- multcompLetters(extract_p(Tukey_t24$culture))

#summary stats for pax
pmax_sum <- summarySE(pidata_pmax, measurevar = "pmax", groupvars = c("culture","hour"))


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
# morpho_48hr <- filter(morpho_48hr, experiment.name != "G_Oceanica Lith Adsorption")

morpho_48hr %>% #using a pipe here
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
hist(morpho_nocontrol$percent.calcified)

KW_morpho <- kruskal.test(morpho_48hr$percent.calcified, morpho_48hr$treatment, method = "bonferroni")
KW_morpho <- kruskal.test(morpho_48hr$percent.calcified, morpho_48hr$treatment, method = "holm")

DT_morpho <- dunnTest(percent.calcified ~ treatment, data = morpho_48hr)

DT_morpho_letters = DT_morpho$res

DT_morpho_letters <- cldList(P.adj ~ Comparison, data = DT_morpho_letters, threshold = 0.05)
DT_morpho_letters

morpho_sum <- summarySE(morphotypes, measurevar = "percent.calcified", groupvars = c("hour","date","liths","treatment"))


ggarrange(organics, morpho, align = 'hv', ncol = 2, nrow = 1)

#removal rates######################################################################################################################################################################################
#looking at the lith removal for each experiment###########
#removing the NA's 
lith_ads <- adsorption[!is.na(adsorption$lith.concentration), ]
lith_ads <- filter(lith_ads, lith.concentration != 0)
lith_ads <- filter(lith_ads, experiment.name != "1000yr old liths")
lith_ads <- filter(lith_ads, hour < 97)
lith_ads <- filter(lith_ads, host.strain != "1516")
lith_ads <- filter(lith_ads, experiment.name != "Organic Removal and Calcein")
lith_ads <- filter(lith_ads, experiment.name != "374 growth and adsorption")
lith_ads <- filter(lith_ads, treatment != "calcified")
lith_ads <- filter(lith_ads, number.of.liths.per.cell != "50:1 G. Oceanica; 50:1 Ehux")

# #making new IDs for the lith removal rates
# lith_ads$lith.concentration <- as.numeric(as.character((lith_ads$lith.concentration))) #forcing a factor to a numeric can warp the actual numbers fo you need to change it to a character afterwards
# lith_ads$ID <- as.factor(paste(lith_ads$liths, lith_ads$replicate, lith_ads$experiment.name, sep="-"))
# lith_ads$rephour <- as.factor(paste(lith_ads$replicate, lith_ads$hour, sep="-"))
# lith_ads$day <- lith_ads$hour / 24
#  
# # 
# #plotting the apparent growth rates after summarizing each replicate and finding the slope
# lith_ads$statID <- as.factor(paste(lith_ads$liths, lith_ads$replicate, lith_ads$experiment.name, lith_ads$treatment,sep="-"))
# 
# lith_slopes <- lith_ads %>% group_by(statID, treatment, experiment.name, liths, replicate, number.of.liths.per.cell) %>%
#   do(model = lm(log(lith.concentration) ~ hour, data = .)) %>%
#   mutate(coef=coef(model)["hour"])
# 
# lith_slopes$rephour <- as.factor(paste(lith_slopes$experiment.name, lith_slopes$replicate, sep="-"))
# 
# lith_slopes <- as.data.table(lith_slopes)
# 
# lith_slopes$removal_per_day <- lith_slopes$coef * 24
# 
# lith_slopes [, blanksubs:= (removal_per_day-(removal_per_day[match ("blank", treatment)])), by=c("rephour")]
# 
# lith_slopes <- filter(lith_slopes, blanksubs != "0")
# # lith_slopes %>%
# #   arrange(blanksubs) %>%
# #   mutate(liths = factor(liths, levels=c("control", "cells+liths"))) %>%
# #   
#   ggplot(lith_slopes, aes(x= experiment.name, y=blanksubs, colour = experiment.name))+  
#   geom_boxplot()+
#   geom_point(aes(colour = experiment.name), size = 5,
#              position = position_jitterdodge())+  
#   theme_classic()+
#   theme(legend.position = "none",  
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),  
#         axis.text.x = element_text(angle = 45, hjust = 1))+  
#     scale_y_reverse()
#   # scale_x_discrete(name = "Treatment", labels = c("Control","Non-calcified cells\n+coccoliths"))+
#   # scale_fill_manual(values = c("#E69F00","#0072B2"))+
#   # scale_colour_manual(labels = c("Control",
#   #                                "Non-calcified cells\n+coccoliths"),
#   #                     values = c("control" = "#E69F00",
#   #                                "cells+liths" = "#0072B2"))
  
#option two looking at the ratio of loss
#subsets the lith dilution series experiment and the others into their own datasheets
ads_data <- subset(lith_ads, experiment.name == "Phaeocystis" | experiment.name == "lith spike" |  
                      experiment.name == "G_Oceanica Lith Adsorption" | experiment.name == "Photo_physiology"
                     , select = 1:31)
ads_data_dilution <- subset(lith_ads, experiment.name == "Lith Adsorption Dilution Series" | experiment.name == "calcein_dilution_series", select = 1:31)
  
#we are goingn to do the slopes for the lith dilution series first 
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
ratioloss_lith_slopes <- as.data.table(filter(ratioloss_lith_slopes, treatment != "blank"))
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

#now for the other removal data
#make new IDs
ads_data$ID <- as.factor(paste(ads_data$replicate, ads_data$treatment, ads_data$experiment.name, sep="-"))
ads_data$rephour <- as.factor(paste(ads_data$experiment.name, ads_data$replicate, ads_data$hour, sep="-"))
  
  
#calc ratio loss and transforming to lognormal
ads_data<- data.table (ads_data, key= c("ID") )
ads_data [, ratiolossln:= log((lith.concentration/(lith.concentration[match("0", hour)]))), by= c("ID")]

  
#subtracting the blanks
#first fun the linear regressions
ratioloss_lith_slopes_ads <- ads_data %>% group_by(experiment.name, ID, treatment, liths, replicate, number.of.liths.per.cell) %>%
    do(model = lm(ratiolossln ~ hour, data = .)) %>%
    mutate(coef=coef(model)["hour"])
#make a new ID to average the blanks and treatments
ratioloss_lith_slopes_ads$avgID <- as.factor(paste(ratioloss_lith_slopes_ads$experiment.name, ratioloss_lith_slopes_ads$treatment, sep="-"))
#convert to a data table
ratioloss_lith_slopes_ads <- data.table(ratioloss_lith_slopes_ads, key = c("avgID"))
#takes the mean
ratioloss_lith_slopes_ads [, blanksavg:=   mean(coef), by = "avgID"]
#new ID to match up the blank averages with the correct treatment
ratioloss_lith_slopes_ads$rephour <- as.factor(paste(ratioloss_lith_slopes_ads$experiment.name, ratioloss_lith_slopes_ads$replicate, ratioloss_lith_slopes_ads$hour, sep="-"))
#subtracts the average slope of each blank  
ratioloss_lith_slopes_ads [, blanksubs:= (coef-(blanksavg[match ("blank", treatment)])), by=c("rephour")]
#removes the blank for plotting  
ratioloss_lith_slopes_ads <- as.data.table(filter(ratioloss_lith_slopes_ads, treatment != "blank"))

#plots as boxplot
ggplot(ratioloss_lith_slopes_ads, aes(x= experiment.name, y=blanksubs, colour = experiment.name))+  
    geom_boxplot()+
    geom_point(aes(colour = experiment.name), size = 5,
               position = position_jitterdodge())+  
    theme_classic()+
    theme(legend.position = "none",  
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),  
          axis.text.x = element_text(angle = 45, hjust = 1))+  
  scale_y_reverse(name = bquote('Slope of remaining fraction of coccoliths ('*hour^-1*')'), limits = c(0.01,-0.017))+  
  scale_x_discrete(name = "Experiment Name")
  # scale_colour_manual(name =  "Treatment", labels = c("1000; Lith:Cell Ratio (100:1)",
  #                                                     "10000; Lith:Cell Ratio (100:1)",
  #                                                     "100000; Lith:Cell Ratio (100:1)",  
  #                                                     "1000; Lith:Cell Ratio (1000:1)"),
  #                     values = c("100-1000" = "#333333",
  #                                "100-10000" = "#333333",
  #                                "100-100000" = "#333333",
  #                                "1000-1000" ="#333333"))



#binding the two data frames together to plot everything  
full_ads_data <- rbind.fill(ratioloss_lith_slopes,ratioloss_lith_slopes_ads)
#plots all data as boxplot
ggplot(full_ads_data, aes(x= experiment.name, y=blanksubs, colour = experiment.name))+  
    geom_boxplot()+
    geom_point(aes(colour = experiment.name), size = 5,
               position = position_jitterdodge())+  
  geom_hline(yintercept = 0)+
    theme_classic()+
    theme(legend.position = "none",  
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),  
          axis.text.x = element_text(angle = 45, hjust = 1))+  
  scale_y_reverse(name = bquote('Slope of remaining fraction of coccoliths ('*hour^-1*')'), limits = c(0.01,-0.017)) 
  scale_x_discrete(name = "Experiment Name", labels = c("G_Oceanica Lith Adsorption" = "G. Oceanica",
                                                        "Lith Adsorption Dilution Series" = "Lith Dilution Series",
                                                        "lith spike" = "Lith Spike",
                                                        "Phaeocystis" ="P. Globosa",
                                                        "Photo_physiology" ="Photophysiology",
                                                        "calcein_dilution_series" = "Calcein Dilution Series"))+
scale_colour_manual(values = c("G_Oceanica Lith Adsorption" = "#0072B2",
                               "Lith Adsorption Dilution Series" = "#0072B2",
                               "lith spike" = "#0072B2",
                               "Phaeocystis" ="#0072B2",
                               "Photo_physiology" ="#0072B2",
                               "calcein_dilution_series" = "#0072B2"))

#removing all the other experiments that weren't just 374
full_ads_data <- filter(full_ads_data, experiment.name != "G_Oceanica Lith Adsorption")
full_ads_data <- filter(full_ads_data, experiment.name != "Phaeocystis")
full_ads_data <- filter(full_ads_data, experiment.name != "Lith Adsorption Dilution Series")
full_ads_data <- filter(full_ads_data, experiment.name != "calcein_dilution_series")

all <- ggplot(full_ads_data, aes(x= experiment.name, y=blanksubs, colour = experiment.name))+  
  geom_boxplot()+
  geom_point(aes(colour = experiment.name), size = 5,
             position = position_jitterdodge())+  
  geom_hline(yintercept = 0)+
  theme_classic()+
  theme(legend.position = "none",  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        axis.text.x = element_text(angle = 45, hjust = 1))+  
  scale_y_reverse(name = bquote('Slope of remaining fraction of coccoliths ('*hour^-1*')'), limits = c(0.01,-0.017))+ 
  scale_x_discrete(name = "Experiment Name", labels = c("lith spike" = "Lith Spike",
                                                        "Photo_physiology" ="Photophysiology"))+
scale_colour_manual(values = c("lith spike" = "#000000",
                               "Photo_physiology" ="#418CFC"))

ggplot(full_ads_data, aes(x= liths, y=blanksubs))+  
  geom_boxplot()+
  geom_point(aes(colour = experiment.name), size = 5,
             position = position_jitterdodge())+  
  geom_hline(yintercept = 0)+
  theme_classic()+
  theme(legend.position = "none",  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        axis.text.x = element_text(angle = 45, hjust = 1))+  
  scale_y_reverse(name = bquote('Slope of remaining fraction of coccoliths ('*hour^-1*')'), limits = c(0.01,-0.017))

ggarrange(all, dilution, align = 'hv', ncol = 2, nrow = 1)
#encounter rates######################################################################################################################################################################################
#looking at encounter rates
beta_DS = 9.981893e-06 #the calculated encounter rate for naked cells and liths

ads_data$cell_encounter_rate <- (as.numeric(as.character(ads_data$cell.concentration)) * beta_DS)
ads_data$lith_encounter_rate <- (ads_data$lith.concentration * beta_DS)

# lithdilution$calcified_cell_conc <- (lithdilution$cell.concentration * lithdilution$percent.calcified)
# lithdilution$num_liths_cal_cell <- (lithdilution$lith.concentration / lithdilution$calcified_cell_conc)

ads_data <- filter(ads_data, treatment != "blank")
ads_data <- filter(ads_data, experiment.name != "Phaeocystis")
ads_data <- filter(ads_data, experiment.name != "G_Oceanica Lith Adsorption")
#plotting the lit encounter rates
all_enc <- ggplot(ads_data,aes(x=as.factor(hour),y=log10(lith_encounter_rate), colour = experiment.name))+  
  geom_boxplot()+  
  geom_point(aes(colour = experiment.name), size = 3, position = position_jitterdodge())+  
  theme_classic()+theme(legend.text = element_text(lineheight = 2),  
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  scale_x_discrete(name = "Hours")+  
  scale_y_continuous(name = bquote('Number of liths encountered; log10 '*~day^-1*''), limits = c(0,2.5))+  
  scale_colour_manual(values = c("lith spike" = "#000000",
                                 "Photo_physiology" ="#418CFC"))

# ggplot(lithdilution, 
#        aes(y=as.factor(lith.concentration), x=as.factor(cell.concentration))) +
#   geom_tile(aes(fill = lith_encounter_rate), colour = "grey50") +
#   labs (title = "host:lith encounters", fill="encounters per day") + 
#   xlab("Host density") + ylab("Lith density") + theme(legend.position = "bottom", legend.direction = "horizontal", legend.key.height=unit(1,"line"), legend.key.width = unit (2, "line"))

ads_data_dilution$cell_encounter_rate <- (as.numeric(as.character(ads_data_dilution$cell.concentration)) * beta_DS)
ads_data_dilution$lith_encounter_rate <- (ads_data_dilution$lith.concentration * beta_DS)

# lithdilution$calcified_cell_conc <- (lithdilution$cell.concentration * lithdilution$percent.calcified)
# lithdilution$num_liths_cal_cell <- (lithdilution$lith.concentration / lithdilution$calcified_cell_conc)

ads_data_dilution <- filter(ads_data_dilution, treatment != "blank")
#plotting the lith encounter rates
ads_data_dilution$plotID <- as.factor(paste(ads_data_dilution$number.of.liths.per.cell, ads_data_dilution$target.density, sep="-"))
ads_data_dilution <- filter(ads_data_dilution, hour < 24)
dilution_enc <- ggplot(ads_data_dilution,aes(x=plotID,y=log10(lith_encounter_rate), colour = plotID))+  
  geom_boxplot()+  
  geom_point(aes(colour = plotID), size = 3, position = position_jitterdodge())+  
  theme_classic()+theme(legend.text = element_text(lineheight = 2),  
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  scale_x_discrete(name = "Hours")+  
  scale_y_continuous(name = bquote('Number of liths encountered'*~cell^-1 ~day^-1*'; log10'), limits = c(0,2.5))+  
  scale_colour_manual(values = c("100-1000" = "#0072B2",
                                 "100-10000" = "#0072B2",
                                 "100-100000" = "#0072B2",
                                 "1000-1000" ="#0072B2"))

#encounter rate summary stats
enc_stats <- summarySE(ads_data_dilution, measurevar = "lith_encounter_rate", groupvars = c("plotID"))

encounter_rate_dilution_series <- summarySE(ads_data_dilution, measurevar = "lith_encounter_rate", groupvars = c("plotID", "hour"))

dilution_and_model <- left_join(ads_data_dilution, modata_liths, by = "plotID")
dilution_and_model <- melt(dilution_and_model, id.vars = c("experiment.name","date","hour","plotID","number.of.liths.per.cell",  
                                                        "cell.concentration","treatment"),  
                                       measure = c("lith_encounter_rate","E_BM_DS"),  #when I used measure.var it took all the columns, measure worked beter
                                       value.name = "enc",
                                       variable.name = "encounter_rate")

dilution_and_model$encounter_rate <- as.factor(dilution_and_model$encounter_rate)
dilution_enc <- ggplot(dilution_and_model,aes(x=plotID,y=log10(enc), colour = encounter_rate))+  
  geom_boxplot()+  
  geom_point(aes(colour = encounter_rate), size = 3, position = position_jitterdodge())+  
  theme_classic()+theme(legend.text = element_text(lineheight = 2),  
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())+
  scale_x_discrete(name = "Initial cell concentration\n(Lith:cell ratio)")+
  scale_y_continuous(name = bquote('Number of liths encountered'*~cell^-1 ~day^-1*'; log10'))+  
  scale_color_manual(name = "Type", labels = c("Experimental", "Model"),  
                       values = c("lith_encounter_rate" = "#111111","E_BM_DS" = "#777777"))


#arranging the plots
ggarrange(dilution_enc, dilution, align = 'h', ncol = 2, nrow = 2, common.legend = TRUE, legend = 'bottom')

removal_rates <- plot_grid(all,
          dilution, all_enc, dilution_enc,
          align = c('h','v'),
          label_x = 0.2,
          ncol = 2) -> new_p1

plot_grid(removal_rates,
          theoretical_enc,
          #labels = c('Fig A','Fig B'),
          label_x = 0.2,
          nrow = 2,  
          scale = c(1, 1, .5))

#viroliths########################################################################################################
#plotting bar graph for viroliths
lithbath <- dplyr::filter(adsorption, experiment.name == "ViroLiths")
lithbath <- dplyr::filter(lithbath, treatment != "type A")

pclithbath<- summarySE(lithbath, measurevar = "percent.calcified", groupvars = c("treatment","hour"))

#plotting
lithbath$hour <- as.numeric(lithbath$hour)
virolith1 <- lithbath %>% #using a pipe here
  arrange(percent.calcified) %>%
  mutate(liths = factor(liths, levels=c("control",  
                                                "cells+liths",  
                                                "viroliths"))) %>%  
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
                      labels = c("Control","CCMP374 + Liths","CCMP374 + VirusLiths"),  
                      values = c("control" = "#E69F00",
                                 "viroliths" = "#56B4E9",
                                 "cells+liths" = "#0072B2"))+
  scale_y_continuous(name = "Percent calcified cells (Side scatter; SSC)", limits = c(0,75), breaks = c(0,25,50,75))+
  scale_x_discrete(name = "Hour")
  


#plotting the viral lith experiment
lithbath$cell.concentration <- as.numeric(lithbath$cell.concentration)
lithbathcells <- summarySE(lithbath,measurevar = "cell.concentration",  
                           groupvars = c("hour","treatment","liths"))

lithbath$viral.concentration <- as.numeric(lithbath$viral.concentration)
lithbathvirus <- summarySE(lithbath,measurevar = "viral.concentration",  
                           groupvars = c("hour","treatment","liths")) %>%  
  dplyr::filter(liths != "cells+liths")



#plotting the growth curve 
virolith2 <- lithbathcells %>% #using a pipe here
  arrange(cell.concentration) %>%
  mutate(liths = factor(liths, levels=c("control",  
                                                "cells+liths",  
                                                "viroliths"))) %>%
  ggplot(aes(x=hour, y=cell.concentration, colour=liths))+  
  geom_point(size=5)+theme_bw()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+
  geom_errorbar(aes(ymin=cell.concentration-se, ymax=cell.concentration+se), width= 5, size=0.75)+  
  geom_line(aes(colour=liths,group=liths))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96))+  
  scale_y_continuous(name = bquote('Cell concentration ('*cells ~ mL^-1 ~ x10^6*')'),  
                     breaks = c(0,1000000,2000000,3000000,4000000),  
                     labels = c(0,1,2,3,4))+
  scale_colour_manual("Treatment",
                      values = c("control" = "#E69F00",
                                 "viroliths" = "#56B4E9",
                                 "cells+liths" = "#0072B2"),
                      labels = c("Control","CCMP374 + Liths","CCMP374 + VirusLiths"))

#plotting virus data
virolith3 <- lithbathvirus %>% #using a pipe here
  arrange(viral.concentration) %>%
  mutate(liths = factor(liths, levels=c("control","viroliths"))) %>%
  ggplot(aes(x=hour, y=viral.concentration, colour=liths))+  
  geom_point(size=5)+theme_bw()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+
  geom_errorbar(aes(ymin=viral.concentration-se, ymax=viral.concentration+se), width= 5, size=0.75)+  
  geom_line(aes(colour=liths,group=liths))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96))+  
  scale_y_continuous(name = bquote('Viral concentration ('*virus ~ mL^-1 ~ x10^8*')'),  
                     breaks = c(0,25000000,50000000,75000000,100000000),  
                     labels = c(0,0.25,0.50,0.75,1))+
  scale_colour_manual("Treatment",
                      values = c("control" = "#E69F00",
                                 "viroliths" = "#56B4E9"),
                      labels = c("Control","CCMP374 + VirusLiths"))

ggarrange(virolith1, virolith2, virolith3, align = 'hv', label.x = 0, 
          nrow = 2, ncol = 2, hjust = -0.35, vjust = 26, 
          common.legend = TRUE, legend = "bottom")

