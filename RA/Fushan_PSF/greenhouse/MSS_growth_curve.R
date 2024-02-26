
library(tidyverse)
library(ggplot2)
library(ggpubr)

MSS <- read.csv(file = '/Users/yupeitseng/Documents/RA/FS_PSF/green_house/FS_PSF greenhouse - MSS_height.csv', sep = ',', header = T)
colnames(MSS)[9:18] <- 1:10
MSS[241:250, 2] <- rep("background", 10)
MSS[241:250, 3] <- rep("background", 10)
MSS_E <- MSS %>% filter(MSS[, 18] != "Dead") %>% filter(apical_wilt == 0) %>% filter(seedling_sp == "ENGERO")
MSS_E[, 9:18] <- sapply(MSS_E[, 9:18], as.numeric)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#ENGERO----
dat_E <- data_summary(MSS_E, varname = "1", groupnames=c("soil_sp", "decay_cat", "soil_inoculum"))
colnames(dat_E)[4] <- "mean"
for (i in as.character(2:10)){
  dat <- data_summary(MSS_E, varname = i, groupnames=c("soil_sp", "decay_cat", "soil_inoculum"))
  colnames(dat)[4] <- "mean"
  dat_E <- rbind(dat_E, dat)
}
dat_E$time <- rep(1:10, each = 13)

dat_E_l <- dat_E %>% filter(soil_sp == "background" | soil_sp == "ENGERO") %>% filter(soil_inoculum == "background" | soil_inoculum == "live")

dat_E_s <- dat_E %>% filter(soil_sp == "background" | soil_sp == "ENGERO") %>% filter(soil_inoculum == "background" | soil_inoculum == "sterile")

dat_M_l <- dat_E %>% filter(soil_sp == "background" | soil_sp == "MACHZU") %>% filter(soil_inoculum == "background" | soil_inoculum == "live")

dat_M_s <- dat_E %>% filter(soil_sp == "background" | soil_sp == "MACHZU") %>% filter(soil_inoculum == "background" | soil_inoculum == "sterile")

p1<- ggplot(dat_E_l, aes(x=time, y=mean, group=decay_cat, color=decay_cat)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05))

p2<- ggplot(dat_E_s, aes(x=time, y=mean, group=decay_cat, color=decay_cat)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05))

p3<- ggplot(dat_M_l, aes(x=time, y=mean, group=decay_cat, color=decay_cat)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05))

p4<- ggplot(dat_M_s, aes(x=time, y=mean, group=decay_cat, color=decay_cat)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05))

figure <- ggarrange(p1, p2, p3, p4, labels = c("ENGERO_live", "ENGERO_sterile", "MACHZU_live", "MACHZU_sterile"),ncol = 2, nrow = 2)
annotate_figure(figure, top = text_grob("ENGERO", 
                                      color = "black", face = "bold", size = 14))


#MACHZU----

MSS_M <- MSS %>% filter(MSS[, 18] != "Dead") %>% filter(apical_wilt == 0) %>% filter(seedling_sp == "MACHZU")
MSS_M[, 9:18] <- sapply(MSS_M[, 9:18], as.numeric)

dat_M <- data_summary(MSS_M, varname = "1", groupnames=c("soil_sp", "decay_cat", "soil_inoculum"))
colnames(dat_M)[4] <- "mean"
for (i in as.character(2:10)){
  dat <- data_summary(MSS_M, varname = i, groupnames=c("soil_sp", "decay_cat", "soil_inoculum"))
  colnames(dat)[4] <- "mean"
  dat_M <- rbind(dat_M, dat)
}
dat_M$time <- rep(1:10, each = 13)

dat_E_l <- dat_M %>% filter(soil_sp == "background" | soil_sp == "ENGERO") %>% filter(soil_inoculum == "background" | soil_inoculum == "live")

dat_E_s <- dat_M %>% filter(soil_sp == "background" | soil_sp == "ENGERO") %>% filter(soil_inoculum == "background" | soil_inoculum == "sterile")

dat_M_l <- dat_M %>% filter(soil_sp == "background" | soil_sp == "MACHZU") %>% filter(soil_inoculum == "background" | soil_inoculum == "live")

dat_M_s <- dat_M %>% filter(soil_sp == "background" | soil_sp == "MACHZU") %>% filter(soil_inoculum == "background" | soil_inoculum == "sterile")

p1<- ggplot(dat_E_l, aes(x=time, y=mean, group=decay_cat, color=decay_cat)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05))

p2<- ggplot(dat_E_s, aes(x=time, y=mean, group=decay_cat, color=decay_cat)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05))

p3<- ggplot(dat_M_l, aes(x=time, y=mean, group=decay_cat, color=decay_cat)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05))

p4<- ggplot(dat_M_s, aes(x=time, y=mean, group=decay_cat, color=decay_cat)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05))

figure <- ggarrange(p1, p2, p3, p4, labels = c("ENGERO_live", "ENGERO_sterile", "MACHZU_live", "MACHZU_sterile"),ncol = 2, nrow = 2)
annotate_figure(figure, top = text_grob("MACHZU", 
                                        color = "black", face = "bold", size = 14))
