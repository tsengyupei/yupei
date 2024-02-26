install.packages('AICcmodavg')
library(rstatix)
library(ggpubr)
library(effects)
library(ggeffects)
library(lmtest)
library(AICcmodavg)

#Import data
RGR <- read.csv(file = '/Users/yupeitseng/Downloads/FS_PSF greenhouse - MSS_height.csv', sep = ',', header = T)
RGR$soil_sp <- as.character(RGR$soil_sp)
#MACHZU
RGR_M <- RGR[1:120, ]
RGR_M$soil_sp <- as.factor(RGR_M$soil_sp)
RGR_M %>% group_by(soil_sp, decay_cat, soil_inoculum) %>% summarize(mean_RGR = mean(RGR))
model_M <- lm(RGR ~ soil_sp + decay_cat + soil_inoculum, data = RGR_M)
model_M_2 <- lm(RGR ~ soil_sp * decay_cat * soil_inoculum, data = RGR_M)
lrtest(model_M, model_M_2)
aictab(cand.set = list(model_M, model_M_2))
summary(model_M)
car::Anova(model_M, type = "II")
# plot(allEffects(model_M), multiline=TRUE, ci.style="bars")

RGR_M$decay_cat <- factor(RGR_M$decay_cat, levels = c("fresh", "shortd","longd"))

M.plot <- ggeffect(model_M, terms = c("decay_cat","soil_sp","soil_inoculum")) %>% plot(ci = T, ci.style = "errorbar", connect.lines = T, line.size = 0.5,colors = c("#d8b365","#A3B18A"), show.title = F)

RGR_M_f <- RGR_M[which(RGR_M$soil_inoculum == 'live'), ]
model_M_f <- lm(RGR ~ soil_sp + decay_cat, data = RGR_M_f)
summary(model_M_f)
car::Anova(model_M_f)
car::Anova(model_M_f, type = "III")


RGR_M_s <- RGR_M[which(RGR_M$soil_inoculum == 'sterile'), ]
model_M_s <- lm(RGR ~ soil_sp + decay_cat, data = RGR_M_s)
summary(model_M_s)
car::Anova(model_M_s, type = "II")


##E all##
RGR_E <- RGR[121:240, ]
RGR_E <- RGR_E[-which(is.na(RGR_E$RGR)), ]
# RGR_E %>% group_by(soil_sp, decay_cat, soil_inoculum) %>% identify_outliers(RGR)
# RGR_E %>% group_by(soil_sp, decay_cat, soil_inoculum) %>% shapiro_test(RGR)
# shapiro_test()
RGR_E$decay_cat <- factor(RGR_E$decay_cat, levels = c("fresh", "shortd","longd"))

model_E <- lm(data = RGR_E, RGR ~ soil_sp + decay_cat + soil_inoculum)
model_E_2 <- lm(data = RGR_E, RGR ~ soil_sp * decay_cat * soil_inoculum)
lrtest(model_E, model_E_2)
summary(model_E)
# ggqqplot(residuals(model))

# RGR_E %>% ggqqplot("RGR") + facet_grid(soil_sp ~ decay_cat ~ soil_inoculum)

# RGR_E %>% levene_test(RGR ~ soil_sp * decay_cat * soil_inoculum)

RGR_E %>% group_by(soil_sp, decay_cat, soil_inoculum) %>% summarize(mean_RGR = mean(RGR))
model_E <- lm(RGR ~ soil_sp * decay_cat * soil_inoculum, data = RGR_E)
summary(model_E)
car::Anova(model_E, type = "III")

E.plot <- ggeffect(model_E, terms = c("decay_cat","soil_sp","soil_inoculum")) %>% plot(ci = T, ci.style = "errorbar", connect.lines = T, line.size = 0.5,colors = c("#d8b365","#A3B18A"), show.title = F)

#Check interaction
with(RGR_E, tapply(RGR, list(soil_sp, decay_cat), mean))
interaction.plot(x.factor = dat$trt, 
                 trace.factor = dat$gender,
                 response = dat$resp)



RGR_E_f <- RGR_E[which(RGR_E$soil_inoculum == 'live'), ]
model_E_f <- lm(RGR ~ soil_sp +  decay_cat, data = RGR_E_f)
summary(model_E_f)
car::Anova(model_E_f)
car::Anova(model_E_f, type = "III")

RGR_E_s <- RGR_E[which(RGR_E$soil_inoculum == 'sterile'), ]
model_E_s <- lm(RGR ~ soil_sp + decay_cat, data = RGR_E_s)
summary(model_E_s)
car::Anova(model_E_s, type = "III")
# p <- c()
# for (i in c(0:11)){
#   p <- c(p, t.test(RGR[(1+(i*20)):(10+(i*20)), 10], RGR[(11+(i*20)):(20+(i*20)), 10], paired=FALSE)$p.value)
# }



#ISS----
#Import data
ISS <- read.csv(file = '/Users/yupeitseng/Downloads/FS_PSF greenhouse - ISS_height.csv', sep = ',', header = T)
