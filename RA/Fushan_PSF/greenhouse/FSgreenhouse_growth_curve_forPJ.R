#library----
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(plotrix)
library(RColorBrewer)

#function----
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#data preperation (all)----
MSS <- read.csv(file = '/Users/yupeitseng/Documents/RA/FS_PSF/green_house/FS_PSF greenhouse - MSS_height.csv', sep = ',', header = T) #import data
colnames(MSS)[9:19] <- 1:11
MSS[241:250, 2] <- rep("background", 10)
MSS[241:250, 3] <- rep("background", 10)

#data preperation (MACHZU)----
MSS_M <- MSS %>% filter(MSS[, 19] != "Dead") %>% filter(apical_wilt == 0) %>% filter(seedling_sp == "MACHZU")
MSS_M[, 9:19] <- sapply(MSS_M[, 9:19], as.numeric, na.rm = T)

dat_M <- data_summary(MSS_M, varname = "1", groupnames=c("soil_sp", "decay_cat", "soil_inoculum"))
colnames(dat_M)[4] <- "mean"
for (i in as.character(2:11)){
  dat <- data_summary(MSS_M, varname = i, groupnames=c("soil_sp", "decay_cat", "soil_inoculum"))
  colnames(dat)[4] <- "mean"
  dat_M <- rbind(dat_M, dat)
}
dat_M$time <- rep(seq(90, 230, 14), each = 13)

#in MACHZU soil (live)
dat_M_l <- dat_M %>% filter(soil_sp == "background" | soil_sp == "MACHZU") %>% filter(soil_inoculum == "background" | soil_inoculum == "live") %>% filter(decay_cat == "fresh" | decay_cat == "shortd")

#in MACHZU soil (sterile)
dat_M_s <- dat_M %>% filter(soil_sp == "background" | soil_sp == "MACHZU") %>% filter(soil_inoculum == "background" | soil_inoculum == "sterile") %>% filter(decay_cat == "fresh" | decay_cat == "shortd")

#plotting----
p_M_l <- ggplot(dat_M_l, aes(x=time, y=mean, group=decay_cat, color=decay_cat)) + 
  geom_line() +
  geom_point(size = 3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(0.05))+
  labs(x = "Days", y = "Height (cm)", colour = "Soil decay categories")+
  theme_classic()+
  scale_x_continuous(breaks = seq(90, 230, by = 14))+
  theme(text = element_text(size = 20)) 

p_M_s <- ggplot(dat_M_s, aes(x=time, y=mean, group=decay_cat, color=decay_cat)) + 
  geom_line() +
  geom_point(size = 3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(0.05))+
  labs(x = "Days", y = "Height (cm)", colour = "Soil decay categories")+
  theme_classic()+
  scale_x_continuous(breaks = seq(90, 230, by = 14))+
  theme(text = element_text(size = 20)) 
