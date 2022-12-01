## ------------------------------------------------------------------------
setwd("D:/R program")
source ("inspack.R")
#values needed 
library(ggplot2)
library(viridis)
library(dplyr)
library(ggpubr)

K= 1.38064852*(10)^-23 #m2 kg/ s2 K boltzmann constant
mu= 1.126*(10)^-3 #kg/m s dynamic viscosity in 18C
v= 1.099*(10)^-6 #m2/s kinematic viscosity in 18C
Reh_calc= 3E-6 #in m radius Ehux #based on Paasche 2002
Reh_naked= 2.5E-6 #in m radius Ehux #based on Paasche 2002
Reh_lith = 1.5E-6 #in m radius #based on Paasche 2002
Rehv= 90*(10)^-9 #in m radius virus, from Castberg et al 2002
Temp = 18+273.15 #temp in kelvin, here assuming 18C
Den_OcM = 1.05 #g/cm3 density organic cell matter
Den_CH2O= 1.025 #g/cm3 density seawater at 18C
Den_Li <- 2.6 #in g/cm3 as Jakob et al
Den_Nc <- 1.05 #based on Paasche 2002
Den_Cc <- 1.19 #based on Paasche 2002

CcNum <- 1*((10)^3)
NcNum <- 0.1*CcNum
ViNum <- (CcNum)*30
LiNum <- (CcNum)*50 #balch 1991 from cj's paper


## ------------------------------------------------------------------------
#Brownian motion (modata)
#1. make a data frame
require(utils)

modata <- expand.grid(group1= as.factor(c("Cc", "Nc", "Li", "Vi")), group2= as.factor(c("Cc", "Nc", "Li", "Vi")))

modata <- modata %>%
  mutate(rad1 = case_when(group1 == "Nc" ~ Reh_naked,
                         group1 == "Cc" ~ Reh_calc, 
                         group1 == "Li" ~ Reh_lith,
                         group1 == "Vi" ~ Rehv)) %>%
  mutate(rad2 = case_when(group2 == "Nc" ~ Reh_naked,
                          group2 == "Cc" ~ Reh_calc, 
                          group2 == "Li" ~ Reh_lith,
                          group2 == "Vi" ~ Rehv)) %>%
  mutate(count1 = case_when(group1 == "Nc" ~ NcNum,
                          group1 == "Cc" ~ CcNum, 
                          group1 == "Li" ~ LiNum,
                          group1 == "Vi" ~ ViNum))%>%
  mutate(count2 = case_when(group2 == "Nc" ~ NcNum,
                            group2 == "Cc" ~ CcNum, 
                            group2 == "Li" ~ LiNum,
                            group2 == "Vi" ~ ViNum))%>%
  mutate(den1 = case_when(group1 == "Nc" ~ Den_Nc,
                            group1 == "Cc" ~ Den_Cc, 
                            group1 == "Li" ~ Den_Li,
                            group1 == "Vi" ~ Den_CH2O))%>%
  mutate(den2 = case_when(group2 == "Nc" ~ Den_Nc,
                         group2 == "Cc" ~ Den_Cc, 
                         group2 == "Li" ~ Den_Li,
                         group2 == "Vi" ~ Den_CH2O))


#2. calculate beta (beta)
modata$beta_BM <-((2*(K*Temp*((modata$rad1+modata$rad2)^2)))/(3*mu*(modata$rad1*modata$rad2)))*86400*10^6 #cm3/day

#3. calculate encounters (E)
modata$E_BM <- modata$beta_BM*modata$count2

ggplot(modata, aes(group1, group2)) +
  geom_tile(aes(fill = log(E_BM)), colour = "grey50")+
  scale_fill_gradient(low = "white", high = "blue")

##arrange factors
modata$group1 <- factor(modata$group1, levels=c("Cc", "Nc", "Li", "Vi"))
modata$group2 <- factor(modata$group2, levels=c("Cc", "Nc", "Li", "Vi")) 

## ------------------------------------------------------------------------
##differential settling
#calc sinkvel
#these sinking velocities are theoretically calculated but the reviewer wants us to use a more realistic sinking velocity for coccoliths, they suggested 0.14 m d^-1 from Honjo, 1976, but we are going to use 0.204 from Zhang et al., 2018 instead 
modata$SinkVel1 <- ((2*((modata$rad1*100)^2)*(981)*(modata$den1-Den_CH2O))/(9*(mu*10)))*864 #meter per day
modata$SinkVel2 <- ((2*((modata$rad2*100)^2)*(981)*(modata$den2-Den_CH2O))/(9*(mu*10)))*864 #meter per day

#this replaces the theoretical sinking rate of 0.59 with 0.204 instead from Zhang et al., 2018
modata$SinkVel3 <- replace(modata$SinkVel1, modata$group1== "Li", 0.205)
modata$SinkVel4 <- replace(modata$SinkVel2, modata$group2== "Li", 0.205)

modata$beta_DS <- (pi*(((modata$rad1+modata$rad2)*100)^2)*(abs((modata$SinkVel3-modata$SinkVel4)/864)))*86400 #in encounters cm3/day




#beta will be the same for pairings in reverse

modata$E_DS <- modata$beta_DS*modata$count2 #note I'm using the Zhang sinking rate from above

resize.win(6,6)
ggplot(modata, aes(group1, group2)) +
  geom_tile(aes(fill = log10(E_DS)), colour = "grey50") +
  scale_fill_viridis(option="cividis") +
  ylim(rev(levels(modata$group2))) + labs (fill="log10 encounters per day") + 
  theme_Publication() + theme(axis.title = element_blank(), legend.position = "bottom", legend.direction = "horizontal", legend.key.height=unit(1,"line"), legend.key.width = unit (2, "line"))


## ------------------------------------------------------------------------
##turbulence
#disrate is cm2/s3

#make data frame

modataext <- modata[rep(seq_len(nrow(modata)), 7), ]
modataext$disrate <- rep_len (c (1 %o% (10)^(seq(-8,-2, 1))), length.out=112)

#the 4.2 coefficient differed from Heidi and was changed to 0.42, this coefficient is represented two different ways in the Kiorbe paper
#Kiorbe 1995- Planktivorous feeding in calm and turbulent environments, with emphasis on copepods
#4.2 was being used in the wrong way, it combined the two equations improperly
modataext$beta_turb <- (0.42*pi*((modataext$disrate/(v))^0.5)*((modataext$rad1+modataext$rad2)^3))*86400*10^6 

modataext$E_turb <- modataext$beta_turb*modataext$count2


## ------------------------------------------------------------------------
## add them all

modataext$beta_all <- modataext$beta_BM + modataext$beta_DS + modataext$beta_turb

ggplot (data=modataext %>% filter (group1 %in% c ("Nc")) %>% filter (!(group2 %in% c("Nc"))), aes(x=log10(disrate),y = log10(beta_all), color=group2)) + geom_line() + geom_point ()

#compute encounters

modataext$E_all <- modataext$beta_all*modataext$count2 #doublechecked with Heidi's model already
modataext$E_allmix <- modataext$beta_all*modataext$count1*modataext$count2
modataext$disratef <- as.factor (modataext$disrate)

ggplot (data=modataext %>% filter (group1 %in% c ("Nc")) %>% filter (!(group2 %in% c("Nc"))), aes(x=log10(disrate),y = log10(E_all), color=group2)) + geom_line() + geom_point ()

ggplot(modataext, aes(x=group2, y=group1)) +
  geom_tile(aes(fill = log10(E_all)), colour = "grey50") + geom_text(size=2.8, color="white", aes(label=formatC(E_all, format = "e", digits = 1))) +
  scale_fill_viridis(option="cividis") + 
  ylim(rev(levels(modata$group2))) +
  theme_Publication() +
  theme(axis.title = element_blank(), legend.position = "bottom", legend.title= element_blank(), legend.direction = "horizontal", legend.key.height=unit(1.5,"line"), legend.key.width = unit (2.5, "line")) + facet_grid(~disrate)

##for cj
modataext$disratecode <- factor (modataext$disrate, levels = c("1e-08", "0.001"),
                                 labels = c("calm", "stormy"))
resize.win (9,6)
#the encounter rates have been logged to better see the difference
enc <- ggplot(modataext %>% filter (disratecode %in% c("calm", "stormy")), aes(x=group2, y=group1)) +
  geom_tile(aes(fill = log10(E_all)))  + labs (fill= expression(log[10]~"encounters"~ "entity"^-1 ~ "day"^-1)) +
  scale_fill_gradient(low = "blue", high = "red") + 
  #geom_text(size=4, color="white", aes(label=formatC(E_all, format = "e", digits = 1))) +
  ylim(rev(levels(modata$group2))) +
  # theme_Publication2() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "bottom", legend.direction = "horizontal", legend.key.height=unit(1.5,"line"), legend.key.width = unit (2.5, "line")) + facet_grid(~disratecode)

#the encounter rates with the numnbers on each tile, has not been logged
ggplot(modataext %>% filter (disratecode %in% c("calm", "stormy")), aes(x=group2, y=group1)) +
  geom_tile(aes(fill = (E_all)))  + 
  scale_fill_viridis(option="cividis") + 
  geom_text(size=4, color="white", aes(label=formatC(E_all, digits = 3))) +
  ylim(rev(levels(modata$group2))) +
  # theme_Publication2() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "bottom", legend.direction = "horizontal", legend.key.height=unit(1.5,"line"), legend.key.width = unit (2.5, "line")) + facet_grid(~disratecode)

##days
modataext$days <- 1/(modataext$E_all)

time_to_enc <- ggplot(modataext %>% filter (disratecode %in% c("calm", "stormy")), aes(x=group2, y=group1)) +
  geom_tile(aes(fill = log10(days)))  + labs (fill= expression(log[10]~"days to encounter")) +
  scale_fill_gradient(low = "blue", high = "red") + 
  #geom_text(size=4, color="white", aes(label=formatC(E_all, format = "e", digits = 1))) +
  ylim(rev(levels(modata$group2))) +
  # theme_Publication2() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "bottom", legend.direction = "horizontal", legend.key.height=unit(1.5,"line"), legend.key.width = unit (2.5, "line")) + facet_grid(~disratecode)

ggarrange(enc,time_to_enc, nrow = 2, ncol = 1, legend = "none")

sum <- modataext %>% filter (disratecode %in% c("calm", "stormy"))

sum_im_dis_rate <- modataext %>% filter(group1 == "Li" & group2 == "Vi")

sum_im_dis_rate %>% ggplot(aes(x = disratef, y = E_all))+geom_point(size=5)

sum_min_frac_1 <- modataext %>% filter(group1 == "Cc" & group2 == "Li")
sum_min_frac_2 <- modataext %>% filter(group1 == "Cc" & group2 == "Vi")
sum_min_frac_3 <- modataext %>% filter(group1 == "Nc" & group2 == "Li")
sum_min_frac_4 <- modataext %>% filter(group1 == "Nc" & group2 == "Vi")

sum_min_frac_cal <- rbind(sum_min_frac_1,sum_min_frac_2)

sum_min_frac_cal <- data.table(sum_min_frac_cal, key= c("disratef"))
sum_min_frac_cal <- sum_min_frac_cal %>%
  mutate(ratio = E_all / lag(E_all, default = first(E_all)))

row_odd <- seq_len(nrow(sum_min_frac_cal)) %% 2              # Create row indicator
row_odd  

data_row_even_cal <- sum_min_frac_cal[row_odd == 0, ]             # Subset odd rows


sum_min_frac_nak <- rbind(sum_min_frac_3,sum_min_frac_4)
sum_min_frac_nak <- data.table(sum_min_frac_nak, key= c("disratef"))
sum_min_frac_nak <- sum_min_frac_nak %>%
  mutate(ratio = E_all / lag(E_all, default = first(E_all)))
row_odd <- seq_len(nrow(sum_min_frac_nak)) %% 2              # Create row indicator
row_odd  

data_row_even_nak <- sum_min_frac_nak[row_odd == 0, ]             # Subset odd rows


sum_im_dis_rate$ratio <- sum_im_dis_rate$E_all
sum_min_frac <- rbind(data_row_even_cal,data_row_even_nak,sum_im_dis_rate)
sum_min_frac$disratef <- as.numeric(as.character(sum_min_frac$disratef))

sum_min_frac %>% ggplot(aes(x = log10(disratef), y = (ratio*100), group = group1))+
  geom_line(aes(color = group1, group = group1), size = 1.5)+  
  theme_classic2()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  scale_color_manual(values = c("Nc" = "#E69F00", "Li" = "#000000", "Cc" = "#CC79A7"))+ 
  scale_y_continuous(name = "% population", limits = c(0,100))+
  scale_x_continuous(name = "Dissipation rate", limits = c(-6,-3))
  
#arrange both the encounter rates and time to encounters
ggarrange(enc, time_to_enc, nrow = 2, ncol = 1, align = 'hv', legend = TRUE )


#number of hours to encounter
sum$hours <- sum$days*24



#finding the percentage of viroliths under calm and storm conditions
#calcified cells encounter liths 
calm_liths_cal <- 2.338756e-01
stormy_liths_cal <- 1.585262e+01
#calcified cells encountering viruses
calm_virus_cal <- 4.513123e-01
stormy_virus_cal <- 3.485450e+00

fraction_viroliths_calcified_calm <- 1/(calm_liths_cal/calm_virus_cal)
fraction_viroliths_calcified_stormy <- 1/(stormy_liths_cal/stormy_virus_cal)
  
#naked cells encountering liths
calm_liths_nak <- 5.281986e-01
stormy_liths_nak <- 1.149774e+01
#naked cells encounter virus
calm_virus_nak <- 2.061095e-01
stormy_virus_nak <- 1.992840e+00

fraction_viroliths_naked_calm <- 1/(calm_liths_nak/calm_virus_nak)
fraction_viroliths_naked_stormy <- 1/(stormy_liths_nak/stormy_virus_nak)


#based on these data the formation rate is 15% of the coccoliths would be viroliths under calm conditions and 56% in stormy conditions
#calcified cells would require the formation rate to be 102% in calm and 21% in stormy
#naked cells would require the formation rate to be 56% in calm and 18% in stormy
#in this case with a slower sinking rate viroliths would only be the dominate mode in stormy conditions when the formation rate
#is higher than the minimum requirement


write.table (sum, "Postdoc-R/Exported Tables/sum_cjpaper_paasche.csv", sep=";", col.names=T, row.names=F)
