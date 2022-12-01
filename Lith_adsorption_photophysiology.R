###read in file
library(tidyverse)
library(dplyr)
library(data.table)
library(phytotools)
library(Rmisc)
library(gridExtra)
library('cowplot')
library(gtable)
library(ggpubr)
library(grid)
library(ggthemes)
library(lemon)
setwd("/Users/Christopher Johns/Dropbox/Lith Adsorption/Spreadsheets/adsorptive_exchange")
phytophys <- read.csv("CJ191811.csv")

#removing unused columns 
phytophys <- subset(phytophys, select = -c(Alp4:Tau4)) 
phytophys <- subset(phytophys, select = -c(MTF.Alp1:TauAv3.1)) 

#adding necessary colunmns
#converting PAR to Photons/m2/sec
phytophys$photons <- ((phytophys$PAR/1000000)*6.02e23)

#converting sigma to m2
phytophys$sigma_m2 <- ((phytophys$Sigma*2.4)/1e20) #2.4 is a correction factor applied each time the FIRe is recallibrated, for the up to date callibration factor contact Max Gorbunov

#adding numerical identifieers to each row
phytophys$ID <- 1:nrow(phytophys)

#adds the sample number to each set
phytophys$SampleID <- rep(c("1","2","3","4","5","6","7","8","9"), each = 18, length=810)
phytophys$hour <- rep(c("0","24","48","72","96"), each = 162, length=810)
phytophys$Treatment <- rep(c("naked","cells + liths","calcified"), each = 54, length=810)

#needed to create a new idea in order to find the correct ratios for the ETR
phytophys$fake_ID <- as.factor(paste(phytophys$DATE, phytophys$SampleID, sep="-"))

#calculating the ratios
phytophys <- data.table (phytophys, key= c("fake_ID") )
phytophys [, photo_eff_ratio:= (Fv.Fm/(Fv.Fm[match("0", PAR)])), by= c("fake_ID")]

#calculating the ETR
phytophys [, ETR:= (photons*(sigma_m2[match("0", PAR)])*photo_eff_ratio), by= c("fake_ID")]
# phytophystrimmed <- filter(phytophys, PAR < 470)
# phytophystrimmed <- filter(phytophystrimmed, hour == "0"| hour == "96")

#calculating the NPQ
phytophys [, NPQ:= (((Fm[match("0", PAR)])-Fm)/(Fm[match("0", PAR)])), by= c("fake_ID")]

#Looking at the trimmed vs. untrimmed PI curves
phytophys%>%  
  arrange(ETR) %>%  
  mutate(Treatment = factor(Treatment, levels = c("naked","cells + liths","calcified"))) %>%
ggplot(aes(x = PAR, y = ETR, colour = Treatment))+  
  geom_point()+  
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  facet_rep_wrap(~hour, scales = "fixed", repeat.tick.labels = FALSE)+  
  theme_classic()+
  geom_smooth(method = 'loess', formula = 'y~x')+ #uses the default 95% confidence intervals
  scale_colour_manual(name = "Treatment", labels = c("naked" = "Naked", "cells + liths" = "Cells + Liths", "calcified" = "Calcified"), values = c("naked" = "#E69F00", "cells + liths" = "#0072B2", "calcified" = "#CC79A7"))+  
  scale_x_continuous(name = bquote('Photosynthetically active radiation (PAR; '*mu~mol ~ photons ~ m^-2 ~ s^-1*')'))+  
  scale_y_continuous(name = bquote('Electron transfer rate (ETR; '*electrons ~ s^-1 ~ PSII^-1*')'), expand = c(0,0))

#for the first 0-24h of the experiments trimming the datasets by removing two of the last light steps (i.e. leaving only 16) did not make a 
#huge difference in the trend of the data, the overall pmax values were just slightly lower and generally within the standard error of the measurements thmselves
#therefore I did not trim off the last two light steps, and I analyzed everything as is
#Looks at the non-photochemical quenching
phytophys %>%  
  arrange(NPQ) %>%  
  mutate(Treatment = factor(Treatment, levels = c("naked","cells + liths","calcified"))) %>%
  ggplot(aes(x = PAR, y = NPQ, colour = Treatment))+  
  geom_point()+  
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
                   panel.border = element_blank(), panel.background = element_blank())+
  facet_wrap(~hour)+  
  geom_smooth(method = 'loess', formula = 'y~x')+ 
  scale_colour_manual(name = "Treatment", labels = c("naked" = "Naked", "cells + liths" = "Cells + Liths", "calcified" = "Calcified"), values = c("naked" = "#E69F00", "cells + liths" = "#0072B2", "calcified" = "#CC79A7"))


pi_sigma <-
  phytophys %>% 
  group_by(SampleID) %>% 
  slice(1)

# ###converting to dataframe
# masterpidata <- as.data.frame(masterpi)
# ###creating unique ID
# masterpidata$ID <- as.factor(paste(masterpidata$culture, masterpidata$date, masterpidata$Rep, sep='_'))

###use to trim PI curves after decide which points to eliminate
masterpidata_t0 <- subset(phytophys, phytophys$hour %in% c('0'))

masterpidata_24hr <- subset(phytophys, phytophys$hour %in% c('24'))

masterpidata_48hr <- subset(phytophys, phytophys$hour %in% c('48'))

masterpidata_72hr <- subset(phytophys, phytophys$hour %in% c('72'))

masterpidata_96hr <- subset(phytophys, phytophys$hour %in% c('96'))


###plotting ETR
ggplot(data = masterpidata_t0, aes(x = PAR, y = ETR)) + geom_point() + facet_wrap(~SampleID)

ggplot(data = masterpidata_24hr, aes(x = PAR, y = ETR)) + geom_point() + facet_wrap(~SampleID) 

ggplot(data = masterpidata_48hr, aes(x = PAR, y = ETR)) + geom_point() + facet_wrap(~SampleID) 

ggplot(data = masterpidata_72hr, aes(x = PAR, y = ETR)) + geom_point() + facet_wrap(~SampleID) 

ggplot(data = masterpidata_96hr, aes(x = PAR, y = ETR)) + geom_point() + facet_wrap(~SampleID) 


###plotting NPQ
ggplot(data = masterpidata_t0, aes(x = PAR, y = NPQ)) + geom_point() + facet_wrap(~SampleID)

ggplot(data = masterpidata_24hr, aes(x = PAR, y = NPQ)) + geom_point() + facet_wrap(~SampleID) 

ggplot(data = masterpidata_48hr, aes(x = PAR, y = NPQ)) + geom_point() + facet_wrap(~SampleID) 

ggplot(data = masterpidata_72hr, aes(x = PAR, y = NPQ)) + geom_point() + facet_wrap(~SampleID) 

ggplot(data = masterpidata_96hr, aes(x = PAR, y = NPQ)) + geom_point() + facet_wrap(~SampleID) 


####Start here

####t0
names(masterpidata_t0) #id is unique to a given RLC
ID <- unique(masterpidata_t0$SampleID) #Hold unique ids
n <- length(ID) #5 unique pidata500N20191017_2trim
#Setup arrays and vectors to store data
#All pidata500N20191017_2trim in example have the same 11 PAR steps in the same order
alpha <- array(NA,c(n,4)) #4 is the number of columns 
ek <- array(NA,c(n,4))
ssr <- rep(NA,n)
residuals <- array(NA,c(n,18)) ##change to the number of steps, edit when trimming off data points to number of points left##
#Loop through individual pidata500N20191017_2trim
for (i in 1:n){
  #Get ith data
  PAR <- masterpidata_t0$PAR[masterpidata_t0$SampleID==ID[i]]
  ETR <- masterpidata_t0$ETR[masterpidata_t0$SampleID==ID[i]]
  #Call function
  myfit <- fitJP(PAR,ETR, normalize=FALSE)
  #Store data
  alpha[i,] <- myfit$alpha
  ek[i,] <- myfit$ek
  ssr[i] <- myfit$ssr
  residuals[i,] <- myfit$residuals
}

masterpidata_t0_res <- as.data.frame(rbind(alpha, ek))
names(masterpidata_t0_res) = c ("Estimate", "Std. Error", "t value", "p value")
masterpidata_t0_res$culture <- rep(c("Naked","Cells + Liths","Calcified"), each = 3, length=9) #change this to total reps
masterpidata_t0_res$parameter <- rep (c("alpha", "ek"), each=9) #change this to total reps
masterpidata_t0_res$hour <- rep (c("0"), length =1)
masterpidata_t0_res$ssr <- ssr #sum of the squared residuals
masterpidata_t0_res$ID <- unique(masterpidata_t0$SampleID) ###add to other codes
alpha <- masterpidata_t0_res %>% filter (parameter == "alpha")
ek<- masterpidata_t0_res  %>% filter (parameter == "ek")
masterpidata_t0_pmax <- as.data.frame (alpha$Estimate*ek$Estimate)
colnames(masterpidata_t0_pmax)[1] <- "pmax"
masterpidata_t0_pmax$culture <- alpha$culture
masterpidata_t0_pmax$ID <- unique(masterpidata_t0$SampleID) ###add to other codes
masterpidata_t0_pmax$hour <- rep (c("0"), length =1)


####24hr
names(masterpidata_24hr) #id is unique to a given RLC
ID <- unique(masterpidata_24hr$SampleID) #Hold unique ids
n <- length(ID) #5 unique pidata500N20191017_2trim
#Setup arrays and vectors to store data
#All pidata500N20191017_2trim in example have the same 11 PAR steps in the same order
alpha <- array(NA,c(n,4)) #4 is the number of columns 
ek <- array(NA,c(n,4))
ssr <- rep(NA,n)
residuals <- array(NA,c(n,18)) ##change to the number of steps, edit when trimming off data points to number of points left##
#Loop through individual pidata500N20191017_2trim
for (i in 1:n){
  #Get ith data
  PAR <- masterpidata_24hr$PAR[masterpidata_24hr$SampleID==ID[i]]
  ETR <- masterpidata_24hr$ETR[masterpidata_24hr$SampleID==ID[i]]
  #Call function
  myfit <- fitJP(PAR,ETR, normalize=FALSE)
  #Store data
  alpha[i,] <- myfit$alpha
  ek[i,] <- myfit$ek
  ssr[i] <- myfit$ssr
  residuals[i,] <- myfit$residuals
}

masterpidata_24hr_res <- as.data.frame(rbind(alpha, ek))
names(masterpidata_24hr_res) = c ("Estimate", "Std. Error", "t value", "p value")
masterpidata_24hr_res$culture <- rep(c("Naked","Cells + Liths","Calcified"), each = 3, length=9) #change this to total reps
masterpidata_24hr_res$parameter <- rep (c("alpha", "ek"), each=9) #change this to total reps
masterpidata_24hr_res$hour <- rep (c("24"), length =1)
masterpidata_24hr_res$ssr <- ssr #sum of the squared residuals
masterpidata_24hr_res$ID <- unique(masterpidata_24hr$SampleID) ###add to other codes
alpha <- masterpidata_24hr_res %>% filter (parameter == "alpha")
ek<- masterpidata_24hr_res  %>% filter (parameter == "ek")
masterpidata_24hr_pmax <- as.data.frame (alpha$Estimate*ek$Estimate)
colnames(masterpidata_24hr_pmax)[1] <- "pmax"
masterpidata_24hr_pmax$culture <- alpha$culture
masterpidata_24hr_pmax$ID <- unique(masterpidata_24hr$SampleID) ###add to other codes
masterpidata_24hr_pmax$hour <- rep (c("24"), length =1)

####48hr
names(masterpidata_48hr) #id is unique to a given RLC
ID <- unique(masterpidata_48hr$SampleID) #Hold unique ids
n <- length(ID) #5 unique pidata500N20191017_2trim
#Setup arrays and vectors to store data
#All pidata500N20191017_2trim in example have the same 11 PAR steps in the same order
alpha <- array(NA,c(n,4)) #4 is the number of columns 
ek <- array(NA,c(n,4))
ssr <- rep(NA,n)
residuals <- array(NA,c(n,18)) ##change to the number of steps, edit when trimming off data points to number of points left##
#Loop through individual pidata500N20191017_2trim
for (i in 1:n){
  #Get ith data
  PAR <- masterpidata_48hr$PAR[masterpidata_48hr$SampleID==ID[i]]
  ETR <- masterpidata_48hr$ETR[masterpidata_48hr$SampleID==ID[i]]
  #Call function
  myfit <- fitJP(PAR,ETR, normalize=FALSE)
  #Store data
  alpha[i,] <- myfit$alpha
  ek[i,] <- myfit$ek
  ssr[i] <- myfit$ssr
  residuals[i,] <- myfit$residuals
}

masterpidata_48hr_res <- as.data.frame(rbind(alpha, ek))
names(masterpidata_48hr_res) = c ("Estimate", "Std. Error", "t value", "p value")
masterpidata_48hr_res$culture <- rep(c("Naked","Cells + Liths","Calcified"), each = 3, length=9) #change this to total reps
masterpidata_48hr_res$parameter <- rep (c("alpha", "ek"), each=9) #change this to total reps
masterpidata_48hr_res$hour <- rep (c("48"), length =1)
masterpidata_48hr_res$ssr <- ssr #sum of the squared residuals
masterpidata_48hr_res$ID <- unique(masterpidata_48hr$SampleID) ###add to other codes
alpha <- masterpidata_48hr_res %>% filter (parameter == "alpha")
ek<- masterpidata_48hr_res  %>% filter (parameter == "ek")
masterpidata_48hr_pmax <- as.data.frame (alpha$Estimate*ek$Estimate)
colnames(masterpidata_48hr_pmax)[1] <- "pmax"
masterpidata_48hr_pmax$culture <- alpha$culture
masterpidata_48hr_pmax$ID <- unique(masterpidata_24hr$SampleID) ###add to other codes
masterpidata_48hr_pmax$hour <- rep (c("48"), length =1)

####72hr
names(masterpidata_72hr) #id is unique to a given RLC
ID <- unique(masterpidata_72hr$SampleID) #Hold unique ids
n <- length(ID) #5 unique pidata500N20191017_2trim
#Setup arrays and vectors to store data
#All pidata500N20191017_2trim in example have the same 11 PAR steps in the same order
alpha <- array(NA,c(n,4)) #4 is the number of columns 
ek <- array(NA,c(n,4))
ssr <- rep(NA,n)
residuals <- array(NA,c(n,18)) ##change to the number of steps, edit when trimming off data points to number of points left##
#Loop through individual pidata500N20191017_2trim
for (i in 1:n){
  #Get ith data
  PAR <- masterpidata_72hr$PAR[masterpidata_72hr$SampleID==ID[i]]
  ETR <- masterpidata_72hr$ETR[masterpidata_72hr$SampleID==ID[i]]
  #Call function
  myfit <- fitJP(PAR,ETR, normalize=FALSE)
  #Store data
  alpha[i,] <- myfit$alpha
  ek[i,] <- myfit$ek
  ssr[i] <- myfit$ssr
  residuals[i,] <- myfit$residuals
}

masterpidata_72hr_res <- as.data.frame(rbind(alpha, ek))
names(masterpidata_72hr_res) = c ("Estimate", "Std. Error", "t value", "p value")
masterpidata_72hr_res$culture <- rep(c("Naked","Cells + Liths","Calcified"), each = 3, length=9) #change this to total reps
masterpidata_72hr_res$parameter <- rep (c("alpha", "ek"), each=9) #change this to total reps
masterpidata_72hr_res$hour <- rep (c("72"), length =1)
masterpidata_72hr_res$ssr <- ssr #sum of the squared residuals
masterpidata_72hr_res$ID <- unique(masterpidata_72hr$SampleID) ###add to other codes
alpha <- masterpidata_72hr_res %>% filter (parameter == "alpha")
ek<- masterpidata_72hr_res  %>% filter (parameter == "ek")
masterpidata_72hr_pmax <- as.data.frame (alpha$Estimate*ek$Estimate)
colnames(masterpidata_72hr_pmax)[1] <- "pmax"
masterpidata_72hr_pmax$culture <- alpha$culture
masterpidata_72hr_pmax$ID <- unique(masterpidata_72hr$SampleID) ###add to other codes
masterpidata_72hr_pmax$hour <- rep (c("72"), length =1)

####96hr
names(masterpidata_96hr) #id is unique to a given RLC
ID <- unique(masterpidata_96hr$SampleID) #Hold unique ids
n <- length(ID) #5 unique pidata500N20191017_2trim
#Setup arrays and vectors to store data
#All pidata500N20191017_2trim in example have the same 11 PAR steps in the same order
alpha <- array(NA,c(n,4)) #4 is the number of columns 
ek <- array(NA,c(n,4))
ssr <- rep(NA,n)
residuals <- array(NA,c(n,18)) ##change to the number of steps, edit when trimming off data points to number of points left##
#Loop through individual pidata500N20191017_2trim
for (i in 1:n){
  #Get ith data
  PAR <- masterpidata_96hr$PAR[masterpidata_96hr$SampleID==ID[i]]
  ETR <- masterpidata_96hr$ETR[masterpidata_96hr$SampleID==ID[i]]
  #Call function
  myfit <- fitJP(PAR,ETR, normalize=FALSE)
  #Store data
  alpha[i,] <- myfit$alpha
  ek[i,] <- myfit$ek
  ssr[i] <- myfit$ssr
  residuals[i,] <- myfit$residuals
}

masterpidata_96hr_res <- as.data.frame(rbind(alpha, ek))
names(masterpidata_96hr_res) = c ("Estimate", "Std. Error", "t value", "p value")
masterpidata_96hr_res$culture <- rep(c("Naked","Cells + Liths","Calcified"), each = 3, length=9) #change this to total reps
masterpidata_96hr_res$parameter <- rep (c("alpha", "ek"), each=9) #change this to total reps
masterpidata_96hr_res$hour <- rep (c("96"), length =1)
masterpidata_96hr_res$ssr <- ssr #sum of the squared residuals
masterpidata_96hr_res$ID <- unique(masterpidata_96hr$SampleID) ###add to other codes
alpha <- masterpidata_96hr_res %>% filter (parameter == "alpha")
ek<- masterpidata_96hr_res  %>% filter (parameter == "ek")
masterpidata_96hr_pmax <- as.data.frame (alpha$Estimate*ek$Estimate)
colnames(masterpidata_96hr_pmax)[1] <- "pmax"
masterpidata_96hr_pmax$culture <- alpha$culture
masterpidata_96hr_pmax$ID <- unique(masterpidata_96hr$SampleID) ###add to other codes
masterpidata_96hr_pmax$hour <- rep (c("96"), length =1)




#####STOP HERE


####binding dataframes
#putting all the data together
pidata_res <- rbind (masterpidata_t0_res, masterpidata_24hr_res, masterpidata_48hr_res, masterpidata_72hr_res, masterpidata_96hr_res)
pidata_pmax <- rbind (masterpidata_t0_pmax, masterpidata_24hr_pmax, masterpidata_48hr_pmax, masterpidata_72hr_pmax, masterpidata_96hr_pmax)

#alpha the initial slopre of the curve
pidata_res_alpha <- filter(pidata_res, parameter == "alpha")
ps_data_res_stats_alpa <- summarySE(pidata_res_alpha, measurevar = "Estimate", groupvars = c("hour","culture"))

alpha <- ps_data_res_stats_alpa %>%  
  arrange(Estimate) %>%  
  mutate(culture = factor(culture, levels = c("Naked","Cells + Liths","Calcified"))) %>%  
  ggplot(aes(x=hour, y=Estimate, fill=culture))+  
  geom_bar(stat='identity',position="dodge")+theme_bw()+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=Estimate-se, ymax=Estimate+se), width=.5, size=0.5, position = position_dodge(width = 0.9))+  
  scale_x_discrete(name= "Hour")+  
  scale_y_continuous(name = "Alpha")+ 
  scale_fill_manual(name = "Type", values = c("Naked" = "#E69F00", "Cells + Liths" = "#0072B2", "Calcified" = "#CC79A7"))

alpha <- pidata_res_alpha %>%  
  filter(hour < 48) %>%
  arrange(Estimate) %>%  
  mutate(culture = factor(culture, levels = c("Naked","Cells + Liths","Calcified"))) %>%  
  ggplot(aes(x=hour, y=Estimate, color=culture))+  
  geom_boxplot()+  
  geom_point(aes(colour=culture), size = 4, position = position_jitterdodge())+
  theme_classic2()+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  scale_x_discrete(name= "Hour")+  
  scale_y_continuous(name = "alpha")+  
  scale_color_manual(name = "Type", values = c("Naked" = "#E69F00", "Cells + Liths" = "#0072B2", "Calcified" = "#CC79A7"))

hist(pidata_res_alpha$Estimate)
shapiro.test(pidata_res_alpha$Estimate)
alpha_t24 <- filter(pidata_res_alpha, hour >23 & hour < 48)

aov_alpha_t24 <- aov(Estimate~culture, data = alpha_t24)

Tukey_alpha_t24 <- TukeyHSD(aov_alpha_t24, "culture", ordered = TRUE)

#ek the minimum saturation irradiance
pidata_res_ek <- filter(pidata_res, parameter == "ek")
ps_data_res_stats_ek <- summarySE(pidata_res_ek, measurevar = "Estimate", groupvars = c("hour","culture"))


ps_data_res_stats_ek %>%  
  arrange(Estimate) %>%  
  mutate(culture = factor(culture, levels = c("Naked","Cells + Liths","Calcified"))) %>%  
  ggplot(aes(x=hour, y=Estimate, fill=culture))+  
  geom_bar(stat='identity',position="dodge")+theme_bw()+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=Estimate-se, ymax=Estimate+se), width=.5, size=0.5, position = position_dodge(width = 0.9))+  
  scale_x_discrete(name= "Hour")+  
  scale_y_continuous(name = "ek value")+  
  scale_fill_manual(name = "Type", values = c("Naked" = "#E69F00", "Cells + Liths" = "#0072B2", "Calcified" = "#CC79A7"))

ek <- pidata_res_ek %>%  
  filter(hour < 48) %>%
  arrange(Estimate) %>%  
  mutate(culture = factor(culture, levels = c("Naked","Cells + Liths","Calcified"))) %>%  
  ggplot(aes(x=hour, y=Estimate, color=culture))+  
  geom_boxplot()+  
  geom_point(aes(colour=culture), size = 4, position = position_jitterdodge())+
  theme_classic2()+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  scale_x_discrete(name= "Hour")+  
  scale_y_continuous(name = "Ek")+  
  scale_color_manual(name = "Type", values = c("Naked" = "#E69F00", "Cells + Liths" = "#0072B2", "Calcified" = "#CC79A7"))

hist(pidata_res_ek$Estimate)
shapiro.test(pidata_res_ek$Estimate)
ek_t24 <- filter(pidata_res_ek, hour >23 & hour < 48)

aov_ek_t24 <- aov(Estimate~culture, data = ek_t24)

Tukey_ek_t24 <- TukeyHSD(aov_ek_t24, "culture", ordered = TRUE)

#pmax maximum photosynthetic rate
pidata_pmax_stats <- summarySE(pidata_pmax, measurevar = "pmax", groupvars = c("hour","culture"))

pidata_pmax_stats %>%  
  arrange(pmax) %>%  
  mutate(culture = factor(culture, levels = c("Naked","Cells + Liths","Calcified"))) %>%  
  ggplot(aes(x=hour, y=pmax, fill=culture))+  
  geom_bar(stat='identity',position="dodge")+theme_bw()+theme_classic2()+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=pmax-se, ymax=pmax+se), width=.5, size=0.5, position = position_dodge(width = 0.9))+  
  scale_x_discrete(name= "Hour")+  
  scale_y_continuous(name = "Pmax")+  
  scale_fill_manual(name = "Type", values = c("Naked" = "#E69F00", "Cells + Liths" = "#0072B2", "Calcified" = "#CC79A7"))

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
  scale_y_continuous(name = "Pmax")+  
  scale_color_manual(name = "Type", values = c("Naked" = "#E69F00", "Cells + Liths" = "#0072B2", "Calcified" = "#CC79A7"))

hist(pidata_pmax$pmax)
shapiro.test(pidata_pmax$pmax)

pmax_24h <- filter(pidata_pmax, hour == "24")

ggplot(pmax_24h, aes(pmax)) +
  geom_histogram(fill = "white", color = "grey30")+facet_wrap(~culture)

aov_pmax_t24h <- aov(pmax~culture, data = pmax_24h)

Tukey_pmax_t24 <- TukeyHSD(aov_pmax_t24h, "culture", ordered = TRUE)

###finding sigma
t0_sigma <-
  masterpidata_t0 %>% 
  group_by(SampleID) %>% 
  slice(1)

hr24_sigma <-
  masterpidata_24hr %>% 
  group_by(SampleID) %>% 
  slice(1)

hr48_sigma <-
  masterpidata_48hr %>% 
  group_by(SampleID) %>% 
  slice(1)

hr72_sigma <-
  masterpidata_72hr %>% 
  group_by(SampleID) %>% 
  slice(1)

hr96_sigma <-
  masterpidata_96hr %>% 
  group_by(SampleID) %>% 
  slice(1)

sigma_master <- rbind(t0_sigma,hr24_sigma,hr48_sigma,hr72_sigma,hr96_sigma)

sigma_stats <- summarySE(sigma_master, measurevar = "Sigma", groupvars = c("hour","Treatment"))

#plotting sigma
sigma_stats %>%  
  arrange(Sigma) %>%  
  mutate(Treatment = factor(Treatment, levels = c("naked","cells + liths","calcified"))) %>%  
  ggplot(aes(x=hour, y=Sigma, fill=Treatment))+  
  geom_bar(stat='identity',position="dodge")+theme_bw()+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=Sigma-se, ymax=Sigma+se), width=.5, size=0.5, position = position_dodge(width = 0.9))+  
  scale_x_discrete(name= "Hour")+  
  scale_y_continuous(name = "Sigma")+  
  scale_fill_manual(name = "Treatment", labels = c("naked" = "Naked", "cells + liths" = "Cells + Liths", "calcified" = "Calcified"), values = c("naked" = "#E69F00", "cells + liths" = "#0072B2", "calcified" = "#CC79A7"))

sigma <- sigma_master %>%  
  filter(hour < 48) %>%
  arrange(Sigma) %>%  
  mutate(Treatment = factor(Treatment, levels = c("naked","cells + liths","calcified"))) %>%  
  ggplot(aes(x=hour, y=Sigma, color=Treatment))+  
  geom_boxplot()+  
  geom_point(aes(colour=Treatment), size = 4, position = position_jitterdodge())+
  theme_classic2()+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  scale_x_discrete(name= "Hour")+  
  scale_y_continuous(name = "Sigma")+  
  scale_color_manual(name = "Type", values = c("naked" = "#E69F00", "cells + liths" = "#0072B2", "calcified" = "#CC79A7"))

hist(sigma_master$Sigma)
shapiro.test(sigma_master$Sigma)
sigma_t24 <- filter(sigma_master, hour >23 & hour < 48)

aov_sigma_t24 <- aov(Sigma~Treatment, data = sigma_t24)

Tukey_sigma_t24 <- TukeyHSD(aov_sigma_t24, "Treatment", ordered = TRUE)

#ploting everything together
ggarrange(alpha, ek, pmax, sigma, align = 'hv',  
          nrow = 2, ncol = 2, common.legend = TRUE)



###############################################################################################################
#this code wasn't used in the final paper
setwd("/Users/Christopher Johns/Dropbox/Rscripts/CSV/20192312_calcified_naked_photophysiology")
calc_naked_phytophys <- read.csv("CJ201923.csv")


#removing unused columns 
calc_naked_phytophys <- subset(calc_naked_phytophys, select = -c(Alp4:Tau4)) 
calc_naked_phytophys <- subset(calc_naked_phytophys, select = -c(MTF.Alp1:TauAv3.1)) 

#adding necessary colunmns
#converting PAR to Photons/m2/sec
calc_naked_phytophys$photons <- ((calc_naked_phytophys$PAR/1000000)*6.02e23)

#converting sigma to m2
calc_naked_phytophys$sigma_m2 <- ((calc_naked_phytophys$Sigma*2)/1e20)

#adding numerical identifieers to each row
calc_naked_phytophys$ID <- 1:nrow(calc_naked_phytophys)

#adds the sample number to each set
calc_naked_phytophys <- calc_naked_phytophys %>% mutate(SampleID =  
                                                          ifelse(ID >=1 & ID <19, "1",  
                                                                 ifelse(ID >=19 & ID <37, "2",  
                                                                        ifelse(ID >=37 & ID <55, "3", ""))))#for some reason needed the extra "" at the end to make it stop failing

calc_naked_phytophys <- filter(calc_naked_phytophys, ID < 55) 



calc_naked_phytophys <- calc_naked_phytophys %>% mutate(Treatment =  
                                                          ifelse(SampleID >=1 & SampleID <2, "calcified",  
                                                                 ifelse(SampleID >=2 & SampleID <3, "calcified and naked",  
                                                                        ifelse(SampleID >=3 & SampleID <=4, "naked", ""))))


#needed to create a new idea in order to find the correct ratios for the ETR
calc_naked_phytophys$fake_ID <- as.factor(paste(calc_naked_phytophys$DATE, calc_naked_phytophys$SampleID, sep="-"))

#calculating the ratios
calc_naked_phytophys <- data.table (calc_naked_phytophys, key= c("fake_ID") )
calc_naked_phytophys [, photo_eff_ratio:= (Fv.Fm/(Fv.Fm[match("0", PAR)])), by= c("fake_ID")]

#calculating the ETR
calc_naked_phytophys [, ETR:= (photons*(sigma_m2[match("0", PAR)])*photo_eff_ratio), by= c("fake_ID")]

#calculating the NPQ
calc_naked_phytophys [, NPQ:= (((Fm[match("0", PAR)])-Fm)/(Fm[match("0", PAR)])), by= c("fake_ID")]


calc_naked_phytophystrimmed <- filter(calc_naked_phytophys, PAR < 502)
                                      

fig2 <- calc_naked_phytophystrimmed %>%  
  arrange(ETR) %>%  
  mutate(Treatment = factor(Treatment, levels = c("naked","calcified and naked","calcified"))) %>%
  ggplot(aes(x = PAR, y = ETR, colour = Treatment))+  
  geom_point()+  
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())+  
  theme_classic()+
  geom_smooth(method = 'loess', formula = 'y~x')+ #uses the default 95% confidence intervals
  scale_colour_manual(name = "Treatment", labels = c("calcified" = "Calcified", "calcified and naked" = "Calcified and Naked", "naked" = "Naked"), values = c("calcified" = "#CC79A7", "calcified and naked" = "#0072B2", "naked" = "#E69F00"))+  
  scale_x_continuous(name = bquote('Photosynthetically active radiation (PAR; '*mu~mol ~ photons ~ m^-2 ~ s^-1*')'))+  
  scale_y_continuous(name = bquote('Electron transfer rate (ETR; '*electrons ~ s^-1 ~ PSII^-1*')'), expand = c(0,0))

calc_naked_phytophys %>%  
  arrange(NPQ) %>%  
  mutate(Treatment = factor(Treatment, levels = c("calcified","calcified and naked","naked"))) %>%
  ggplot(aes(x = PAR, y = NPQ, colour = Treatment))+  
  geom_point()+  
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  facet_wrap(~hour)+  
  theme_classic()+
  geom_smooth(method = 'loess', formula = 'y~x')+ #uses the default 95% confidence intervals
  scale_colour_manual(name = "Treatment", labels = c("calcified" = "Calcified", "calcified and naked" = "Calcified and Naked", "naked" = "Naked"), values = c("calcified" = "#E69F00", "calcified and naked" = "#0072B2", "naked" = "#CC79A7"))+  
  scale_x_continuous(name = bquote('Photosynthetically active radiation (PAR; '*mu~mol ~ photons ~ m^-2 ~ s^-1*')'))+  
  scale_y_continuous(name = bquote('Electron transfer rate (ETR; '*electrons ~ s^-1 ~ PSII^-1*')'), expand = c(0,0))

ggarrange(fig1, fig2, align = 'hv',  
          nrow = 1, ncol = 2, common.legend = FALSE, legend = "bottom")

