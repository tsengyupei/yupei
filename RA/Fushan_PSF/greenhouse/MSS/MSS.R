library(patchwork)
library(PerformanceAnalytics)
library(ggplot2)
library(scatterpie)
library(gtools)
library(rpsychi)
library(tidyverse)
library(stringr)
library(rstatix)

#0. Read in the data and clean out empty rows. ------
biomass_FS <- read.csv("/Users/yupeitseng/Documents/RA/FS_PSF/green_house/MSS/FS_greenhouse_harvest - MSS.csv")
height <- read.csv('/Users/yupeitseng/Documents/RA/FS_PSF/green_house/FS_PSF greenhouse - MSS_height.csv')

height <- height %>% mutate_at(c(9:19), as.numeric) %>%  mutate(growth_rate = X5.29 - X12.26)

biomass_FS <- biomass_FS %>% select(c(1:7, 15:18)) %>% mutate(replicate = c(rep(1:10, 24), rep(1:5, 2)), above_mass = as.numeric(above_mass), below_mass = as.numeric(below_mass), total_mass = (above_mass+below_mass)) %>% select(-apical_wilt)

biomass_FS <- left_join(biomass_FS, height[, c(5, 19:21)])
colnames(biomass_FS)[13] <- 'height'
biomass_FS <- biomass_FS %>% mutate(root_shoot_ratio = below_mass/above_mass)
biomass_FS[which(biomass_FS$soil_sp == "bg"), c("decay_cat", "sterilization")] <- matrix(ncol = 2, nrow = 10, rep("bg", 20))
# write.csv(biomass_FS, "/Users/yupeitseng/Documents/RA/FS_PSF/R_markdown/perfomance_FS.csv")

#1. Data investigation----
#checking the variation of root shoot ratio within group----
# Box plot by group with error bars
ggplot(biomass_FS %>% filter(soil_sp != "bg"), aes(x = plant_sp, y = root_shoot_ratio, col = sterilization)) + 
  stat_boxplot(geom = "errorbar", # Error bars
               width = 0.25) +    # Bars width
  geom_boxplot()+
  scale_color_manual(values=c("#A7C084", "#B6A87C"), "Sterilization")+
  labs(x = "Plant species", y = "Root-shoot ratio") +
  theme(text = element_text(size = 25)) +
  facet_grid(decay_cat ~ soil_sp)
    

#histogram for data distribution----
biom_hist <- ggplot(biomass_FS) +
  geom_histogram(aes(x = total_mass)) +
  facet_grid(.~plant_sp + soil_sp + decay_cat)
biom_hist_log <- biom_hist + scale_x_log10() + ggtitle("logged x-axis")
biom_hist/biom_hist_log

chart.Correlation(biomass_FS[, c('above_mass', 'below_mass', 'total_mass', 'height', 'growth_rate')])

# Make a "Cleveland dotplot" type plot of the biomass values----
# Theme for plots
theme_gsk <- function() {
  theme_minimal()+
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.tag = element_text(face = "bold")
    ) 
}

cleveland_dotplot <- function(value = value, x_title = x_title)
  {biomass_for_cleveland <- biomass_FS %>% select("plant_sp", "soil_sp", "decay_cat", "sterilization", value) 
  colnames(biomass_for_cleveland)[5] <- 'performance'
  biomass_for_cleveland <-  biomass_for_cleveland %>% group_by(plant_sp, soil_sp, decay_cat, sterilization) %>%
  summarize(
    # Get the mean and SEM for each (focal*soil source) group:
    mean_bm = mean(performance, na.rm = T),
    se_bm = sd(performance, na.rm = T)/sqrt(n()),
    # Also get the median and IQRs
    median_bm = median(performance, na.rm = T),
    lqr = quantile(performance, 0.25, na.rm = T),
    uqr = quantile(performance, 0.75, na.rm = T),
    # Get outliers
    min_val = min(performance, na.rm = T),
    max_val = max(performance, na.rm = T),
    out_low = ifelse(min_val < lqr - 1.5*(uqr-lqr), min_val, NA),
    out_upp = ifelse(max_val > uqr + 1.5*(uqr-lqr), max_val, NA),
    n = n()) 

# add a column that sets the y-value of the points -- 
# needed, because I want to add the outlier points
biomass_for_cleveland$plot_y <- (c(seq(0.6, 1.1, by = .1),
                                      seq(1.6, 2.1, by = .1)+1,
                                      2.6+2,
                                      seq(3.6, 4.1, by = .1)+3,
                                      seq(4.6, 5.1, by = .1)+4,
                                      5.6+5))
biomass_for_cleveland <- biomass_for_cleveland %>% mutate(pointcol = paste(decay_cat, sterilization, sep = '_'), treatment = paste(plant_sp, soil_sp, sep = '_'))
# Make a new data frame that includes the outlier points
# Outliers here are defined as those points that are either
# lower than (LQR - 1.5*IQR), or higher than (UQR + 1.5*IQR)
biomass_raw <- biomass_FS %>% select("plant_sp", "soil_sp", "decay_cat", "sterilization", value) 
colnames(biomass_raw)[5] <- 'performance'
outliers <- left_join(biomass_raw, biomass_for_cleveland) %>%
  filter(performance < lqr - (1.5*(uqr-lqr)) | 
           performance > uqr + (1.5*(uqr-lqr))) %>%
  select(plant_sp, soil_sp, pointcol, performance, plot_y)


(biomass_dotplot <- ggplot(biomass_for_cleveland, 
                           aes(y = plot_y, x = median_bm, fill = pointcol, label = treatment, shape = pointcol)) + scale_x_log10() +
    # add median value
    geom_point(size = 4, stroke = .9) +
    # add dashed lines connecting outlier to median (if outlier exists)
    geom_errorbarh(aes(xmin = out_low, xmax = lqr),
                   height = 0, linetype = 3, size = .25) +
    geom_errorbarh(aes(xmin = uqr, xmax = out_upp),
                   height = 0, linetype = 3, size = .25) +
    # add outlying points (if they exist)
    geom_point(aes(y = plot_y, x = out_low), shape = 21, fill = "grey50", size = .5) +
    geom_point(aes(y = plot_y, x = out_upp), shape = 21, fill = "grey50", size = .5) +
    geom_point(data = outliers, aes(x = performance, y = plot_y),
               shape = 21, fill = "grey50", size = .5, inherit.aes = F) +
    
    # add solid line connecting median to lower and upper quantile
    geom_errorbarh(aes(xmin = lqr, xmax = uqr), height = .1) +
    
    # Set color and shape for all the points
    scale_fill_manual(values = c("cornsilk4", "cornsilk", "blue4", "cadetblue1","darkgreen", "darkseagreen", "darkseagreen1"),
                      name = "Inoculum\nsource") +
    scale_shape_manual(values = c(22, 21, 22, 21, 23, 22, 21),
                       name = "Inoculum\nsource") +
    # add thin line between each focal species
    geom_hline(yintercept = c(2, 4, 6, 8, 10), linetype = 1, size = 0.25) +
    # add species names along y-axis
    scale_y_continuous(breaks = c(1,3,5,7,9,11),
                       labels = (unique(biomass_for_cleveland$treatment))) +
    # misc. theme settings.
    ylab("Focal species_Soil species") +
    xlab(x_title) +
    theme_gsk() +
    theme(#legend.justification=c(1,0), legend.position=c(1,.15),
      legend.position = "top",
      # legend.background = element_rect(colour = "grey50"),
      legend.text = element_text(size = 8)) + 
    NULL)
}

cleveland_dotplot("above_mass", 'Aboveground biomass (mg)')
cleveland_dotplot("height", 'Height (mg)')

pdf('/Users/yupeitseng/Documents/RA/FS_PSF/green_house/MSS/MSS_cleveland_dotplot.pdf', height = 8, width = 12)

# for(i in c("above_mass", "below_mass", "total_mass", "height", "growth_rate"))
# {
#   print(cleveland_dotplot(i, i))
# }
# dev.off()

# Make a boxplot for raw data ----
biomass_FS_box <- biomass_FS %>%  filter(soil_sp != "bg") %>% select(-c("rosetting")) %>% mutate(soil_sp = sapply(soil_sp, FUN = function(x){ifelse(x == "E", "Soil_E", "Soil_M")}), plant_sp = sapply(plant_sp, FUN = function(x){ifelse(x == "E", "Plant_E", "Plant_M")}))

stat.test <- biomass_FS_box %>%
  group_by(plant_sp, soil_sp, decay_cat) %>%
  t_test(log_total_mass ~ sterilization) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

stat.test <- biomass_FS_box %>%
  group_by(plant_sp, soil_sp, sterilization) %>%
  t_test(log_total_mass ~ decay_cat) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

stat.test <- biomass_FS_box %>%
  group_by(plant_sp, decay_cat, sterilization) %>%
  t_test(total_mass ~ soil_sp) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

# stat.test <- biomass_FS_box %>%
#   group_by(plant_sp, decay_cat, sterilization) %>%
#   t_test(total_mass ~ soil_sp) %>%
#   adjust_pvalue(method = "bonferroni") %>%
#   add_significance()
# stat.test 

# Create box plots with significance levels
# Hide ns (non-significant)
plot <- ggboxplot(
  biomass_FS_box, x = "decay_cat", y = "total_mass", color = "sterilization",
  facet = c("plant_sp", "soil_sp")
) +
  stat_pvalue_manual(stat.test_1, hide.ns = TRUE)+
  stat_pvalue_manual(stat.test_2, hide.ns = TRUE)+
  scale_color_manual(values=c("#556B37", "#B6A87C"), "Experiment")+
  labs(x = "Decay category", y = "Total mass")+
  theme(text = element_text(size = 20))

biomass_FS_box_m <- biomass_FS_box[which(biomass_FS_box$plant_sp == "Plant_M"), ]
biomass_FS_box_e <- biomass_FS_box[which(biomass_FS_box$plant_sp == "Plant_E"), ]

# biomass_FS_box$log_total_mass <- (biomass_FS_box$total_mass)^(1/3)
# biomass_FS_box %>%
#   levene_test(log_total_mass~plant_sp*soil_sp*decay_cat*sterilization)
# plot(lm(total_mass~plant_sp*soil_sp*decay_cat*sterilization, data=biomass_FS_box))
# 
# biomass_FS_box_rm %>%
#   group_by(plant_sp,soil_sp, decay_cat, sterilization) %>%
#   shapiro_test(log_total_mass)

results_m<-aov(total_mass~soil_sp*decay_cat*sterilization, data=biomass_FS_box_m)
anova(results_m)
test <- biomass_FS_box_m %>% tukey_hsd(total_mass~soil_sp*decay_cat)
test
results_e<-aov(total_mass~soil_sp*decay_cat*sterilization, data=biomass_FS_box_e)
anova(results_e)

# Make a bar plot for raw data----
biomass_FS_st <- biomass_FS %>% select(c("soil_sp", "decay_cat", "sterilization", "plant_sp", "total_mass")) 
colnames(biomass_FS_st)[5] <- 'performance'
biomass_FS_st <- biomass_FS_st %>% group_by(plant_sp, soil_sp, decay_cat, sterilization) %>% summarise(se = sd(performance, na.rm = T)/sqrt(n()), mean = mean(performance, na.rm = T)) %>%mutate(treatment = paste(decay_cat, sterilization, sep = '_'), group = paste(plant_sp, soil_sp, sep = '_')) %>% filter(soil_sp != 'bg')

ggplot(biomass_FS_st, aes(x=treatment, y=mean, fill=decay_cat)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9))+
  labs(y= "Total mass", x = "Plant sp._Soil sp.", fill = "Decay catagories")+
  scale_fill_brewer(palette='Greens')+
  facet_wrap(~group)

#Using model fitting approach to invastigate microbial effect----
#with mean sterilization
biomass_FS_ms <- biomass_FS %>% select(c("soil_sp", "decay_cat", "sterilization", "plant_sp", "total_mass")) 
colnames(biomass_FS_ms)[5] <- 'performance' 
biomass_FS_ms <- biomass_FS_ms %>% left_join(biomass_FS_ms %>%  group_by(soil_sp, decay_cat, plant_sp, sterilization) %>% summarise(mean = mean(performance, na.rm = T)) %>% filter(sterilization == 'sterile'), by = c('soil_sp',  'decay_cat', 'plant_sp', 'sterilization')) 
biomass_FS_ms[which(is.na(biomass_FS_ms$mean)), 'mean'] <- biomass_FS_ms[which(is.na(biomass_FS_ms$mean)), 'performance']
E <- biomass_FS_ms %>% select(-performance) %>% rename(performance = mean)

#with real sterilization
E <- biomass_FS %>% select(c("soil_sp", "decay_cat", "sterilization", "plant_sp", "total_mass")) 
colnames(E)[5] <- 'performance' 

#later calculation
E <- E %>% mutate(performance = log(performance)) %>% 
  filter(soil_sp != 'bg') %>% 
  mutate(sterilization = ifelse(sterilization == 'live', 1, 0)) %>% 
  group_by(plant_sp, soil_sp, decay_cat) %>% 
  # group_by(plant_sp, soil_sp) %>%
  summarise(mean = summary(lm(performance~sterilization))$coefficients[2, 1], se = summary(lm(performance~sterilization))$coefficients[2, 2], p = summary(lm(performance~sterilization))$coefficients[2, 4], n = sum(!is.na(performance))) %>% 
  mutate(treatment = paste(plant_sp, soil_sp, sep = "_")) %>% 
  mutate(ano_x = ifelse(p<0.001, '***', ifelse(p<0.01, '**', ifelse(p<0.05, '*', NA))), ano_y = ifelse(p<0.05, mean+se, NA)) %>% 
  mutate(sd = se * sqrt(n))

E <- E %>% mutate(plant_sp = ifelse(plant_sp == "E", "Plant_E", "Plant_M"), soil_sp = ifelse(soil_sp == "E", "Soil_E", "Soil_M"))

ggbarplot(
  E, x = "decay_cat", y = "mean", fill = "decay_cat", col = NULL,
  facet = c("plant_sp", "soil_sp")
) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9))+
  labs(y= "Microbial effect", x = "Decay category", fill = "Decay catagory")+
  scale_fill_manual(values=c("#556B37", "#A7C084", "#CAD9B5"))+
  geom_text(aes(label=ano_x, group = decay_cat), col='black', position = position_dodge(width = 1), vjust = -0.5, size = 10)+
  theme(text = element_text(size = 20))

#Two way anova by mean
two_way_anova <- function(anova_dat = anova_dat){
  # anova_dat <- E %>% filter(plant_sp == "Plant_E")
  total_mean <- sum(anova_dat$mean*anova_dat$n)/sum(anova_dat$n)
  ss_soil <- sum(anova_dat %>% group_by(soil_sp) %>% summarise(mean = (sum(mean*n))/sum(n), n = sum(n)) %>% mutate(ss = ((mean-total_mean)^2)*n) %>% select("ss"))
  ss_decay <- sum(anova_dat %>% group_by(decay_cat) %>% summarise(mean = (sum(mean*n))/sum(n), n = sum(n)) %>% mutate(ss = ((mean-total_mean)^2)*n) %>% select("ss"))
  ss_between <- sum(((anova_dat$mean - total_mean)^2)*anova_dat$n)
  ss_soil_decay <- ss_between - ss_soil - ss_decay
  ss_within <- sum(as.data.frame(anova_dat %>% mutate(ss = ((sd)^2)*(n-1)) %>% select("ss"))$ss)
  df_within <- sum(anova_dat$n) - 5
  anova_within <- ss_within/df_within
  
  anova_tab <- matrix(ncol = 2, nrow = 3, c(ss_soil, ss_decay, ss_soil_decay, 1, 2, 2), byrow = F)
  colnames(anova_tab) <- c("ss", "df")
  anova_tab <- anova_tab %>% as.data.frame() %>% mutate(mean_squares = ss/df) %>% mutate(f = mean_squares/anova_within) %>% mutate(p = pf(f, df, (sum(anova_dat$n) - 5), lower.tail = FALSE)) %>% add_row(tibble_row(ss = ss_within, df = df_within, mean_squares = anova_within))
  rownames(anova_tab) <- c("soil", "decay", "soil:decay", "error")
  return(anova_tab)
}
two_way_anova(E %>% filter(plant_sp == "Plant_E"))
two_way_anova(E %>% filter(plant_sp == "Plant_M"))

## Bar plot
# ggplot(E, aes(x=treatment, y=mean, fill=decay_cat)) +
#   geom_bar(position=position_dodge(), stat="identity") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9))+
#   labs(y= "Microbial effect", x = "Plant sp._Soil sp.", fill = "Decay catagory")+
#   scale_fill_manual(values=c("#556B37", "#A7C084", "#CAD9B5"))+
#   geom_text(aes(label=ano_x, group = factor(decay_cat, levels = c("fresh", "shortd", "longd"))), col='black', position = position_dodge(width = 1), vjust = -0.5, size = 5)+
#   theme(text = element_text(size = 20))  

one_way_aov <- function(group = group){
  tab <- ind.oneway.second(m = filter(E, treatment == group)$mean,
                           sd = filter(E, treatment == group)$sd,
                           n= filter(E, treatment == group)$n) 
  # tab <- ind.oneway.second(m = filter(E, plant_sp == group)$mean,
  #                          sd = filter(E, plant_sp == group)$sd,
  #                          n= filter(E, plant_sp == group)$n) 
  (pf(tab$anova.table[1, 4], tab$anova.table[1, 2], tab$anova.table[2, 2], lower.tail = FALSE))
}

# for(i in unique(E$plant_sp)){
for(i in unique(E$treatment)){
  print(i)
  print(one_way_aov(i))
}

#Anova by raw data----
E_aov <- biomass_FS %>% select(c("soil_sp", "decay_cat", "sterilization", "plant_sp", "total_mass")) 
colnames(E_aov)[5] <- 'performance'
E_aov <- E_aov %>% mutate(performance = performance) %>% 
  filter(soil_sp != 'bg')

summary(lm(performance ~ soil_sp + decay_cat + sterilization + soil_sp*decay_cat + soil_sp*sterilization, data = as.data.frame(filter(E_aov, plant_sp == 'E'))))

summary(lm(performance ~ soil_sp + decay_cat + sterilization + soil_sp*decay_cat + soil_sp*sterilization, data = as.data.frame(filter(E_aov, plant_sp == 'M'))))

#02. Invasion growth rate----
calculate_IGR <- function(df) {
  df %>% 
    # In df, each row is one repilcate/rack, and each column
    # represents the growth of one species in one soil type.
    # e.g. "FE_ACWR" is the growth of FE in ACWR-cultivated soil.
    mutate(freshE_freshM_E = (E_M_fresh_live - E_M_fresh_sterile) - (M_M_fresh_live - M_M_fresh_sterile),
           freshE_freshM_M = (M_E_fresh_live - M_E_fresh_sterile) - (E_E_fresh_live - E_E_fresh_sterile),
           shortdE_shortdM_E = (E_M_shortd_live - E_M_shortd_sterile) - (M_M_shortd_live - M_M_shortd_sterile),
           shortdE_shortdM_M = (M_E_shortd_live - M_E_shortd_sterile) - (E_E_shortd_live - E_E_shortd_sterile),
           longdE_longdM_E = (E_M_longd_live - E_M_longd_sterile) - (M_M_longd_live - M_M_longd_sterile),
           longdE_longdM_M = (M_E_longd_live - M_E_longd_sterile) - (E_E_longd_live - E_E_longd_sterile),
           freshE_shortdM_E = (E_M_shortd_live - E_M_shortd_sterile) - (M_M_shortd_live - M_M_shortd_sterile),
           freshE_shortdM_M = (M_E_fresh_live - M_E_fresh_sterile) - (E_E_fresh_live - E_E_fresh_sterile),
           freshE_longdM_E = (E_M_longd_live - E_M_longd_sterile) - (M_M_longd_live - M_M_longd_sterile),
           freshE_longdM_M = (M_E_fresh_live - M_E_fresh_sterile) - (E_E_fresh_live - E_E_fresh_sterile),
           shortdE_longdM_E = (E_M_longd_live - E_M_longd_sterile) - (M_M_longd_live - M_M_longd_sterile),
           shortdE_longdM_M = (M_E_shortd_live - M_E_shortd_sterile) - (E_E_shortd_live - E_E_shortd_sterile),
           shortdE_freshM_E = (E_M_fresh_live - E_M_fresh_sterile) - (M_M_fresh_live - M_M_fresh_sterile),
           shortdE_freshM_M = (M_E_shortd_live - M_E_shortd_sterile) - (E_E_shortd_live - E_E_shortd_sterile),
           longdE_freshM_E = (E_M_fresh_live - E_M_fresh_sterile) - (M_M_fresh_live - M_M_fresh_sterile),
           longdE_freshM_M = (M_E_longd_live - M_E_longd_sterile) - (E_E_longd_live - E_E_longd_sterile),
           longdE_shortdM_E = (E_M_shortd_live - E_M_shortd_sterile) - (M_M_shortd_live - M_M_shortd_sterile),
           longdE_shortdM_M = (M_E_longd_live - M_E_longd_sterile) - (E_E_longd_live - E_E_longd_sterile)
    )%>%
    select(pair_id, freshE_freshM_E:longdE_shortdM_M) %>% 
    gather(pair_id, IGR, freshE_freshM_E:longdE_shortdM_M) %>% 
    mutate(plant_sp = sapply(strsplit(pair_id, "_"), function(x){x[3]}), treatment = paste(sapply(strsplit(pair_id, "_"), function(x){x[1]}), sapply(strsplit(pair_id, "_"), function(x){x[2]}), sep = '_'))
}

#pre-paired
IGR <- calculate_IGR(biomass_wide) %>% group_by(pair_id, plant_sp, treatment) %>% summarise(mean = mean(IGR, na.rm = T), se = sd(IGR, na.rm = T)/n()) 

#rnorm data generation----
biomass_wide_st <- cbind(apply(biomass_wide[, -1], 2, FUN = function(x){mean(x, na.rm = T)}), apply(biomass_wide[, -1], 2, FUN = function(x){sd(x, na.rm = T)}))
colnames(biomass_wide_st) <- c('mean', 'sd')
rownames(biomass_wide_st)

n_reps <- 1000
vec <- c()
for(i in rownames(biomass_wide_st)){
  set.seed(31415)
  vec <- c(vec, rnorm(n_reps, mean = biomass_wide_st[i, 1], 
                      sd = biomass_wide_st[i, 2]))
}
mat <- matrix(vec, ncol = 24, nrow = n_reps, byrow = F)
colnames(mat) <- rownames(biomass_wide_st)
mat <- as.tibble(mat) %>% mutate(pair_id = 1:n_reps)

IGR <- calculate_IGR(mat) %>% group_by(pair_id, plant_sp, treatment) %>% summarise(mean = mean(IGR, na.rm = T), sd = sd(IGR, na.rm = T))

#bootstrap data generation----
n_reps <- 1000
vec <- c()
for (i in colnames(biomass_wide)[-1]){
  set.seed(31415)
  for (j in 1:1000){
    vec <- c(vec, mean(sample((na.exclude(unlist(unname(as.vector(biomass_wide[, i]))))), 10, replace = T), na.rm = T))
  }
}
mat <- matrix(vec, ncol = 24, nrow = n_reps, byrow = F)
colnames(mat) <- colnames(biomass_wide)[-1]
mat <- as.tibble(mat) %>% mutate(pair_id = 1:n_reps)

IGR <- calculate_IGR(mat) %>% group_by(pair_id, plant_sp, treatment) %>% summarise(mean = mean(IGR, na.rm = T), sd = sd(IGR, na.rm = T))

# ggline(calculate_IGR(biomass_wide), x = "pair_id", y = "IGR", 
#        add = c("mean_se", "jitter"), 
#        # order = c("ctrl", "trt1", "trt2"),
#        ylab = "IGR", xlab = "Treatment")

#se
ggplot(IGR, aes(x=plant_sp, y=mean, fill = plant_sp)) +
  geom_errorbar(width=.1, aes(ymin=mean-se, ymax=mean+se)) +
  geom_errorbar(width=.1, aes(ymin=mean-se, ymax=mean+se), data=IGR) +
  geom_point(shape=21, size=3)+
  facet_wrap(~treatment)+
  scale_fill_manual(values=c("#556B37", "#B6A87C"), "Plant species")+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank())+
  xlab("Pair") +
  ylab("Invasion growth rate")+
  geom_hline(yintercept=0)+
  theme(text = element_text(size = 20))

#sd
ggplot(IGR, aes(x=plant_sp, y=mean, fill = plant_sp)) +
  geom_errorbar(width=.1, aes(ymin=mean-sd, ymax=mean+sd)) +
  geom_errorbar(width=.1, aes(ymin=mean-sd, ymax=mean+sd), data=IGR) +
  geom_point(shape=21, size=3)+
  facet_wrap(~treatment)+
  scale_fill_manual(values=c("#556B37", "#B6A87C"), "Plant species")+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank())+
  xlab("Pair") +
  ylab("Invasion growth rate")+
  geom_hline(yintercept=0)+
  theme(text = element_text(size = 20))

#03. Reshaping data for calculating Stabilization and FDs ----
# So first, let's make a column that is simply the log AGB
# of each pot, and then go from there. 
biomass_for_SF <- biomass_FS %>% filter(soil_sp != 'bg') %>% select("plant_sp", "soil_sp", "decay_cat", "sterilization", "pair_id", "total_mass") 
colnames(biomass_for_SF)[6] <- 'performance'
biomass_for_SF <- biomass_for_SF %>% mutate(log_performance = log10(performance))


# Now, we need to convert the "long format" biomass data frame  
# into a "wide format". Here's a function for that:
make_wide_biomass <- function(df) {
  # This takes the biomass dataframe, and converts into a wide
  # dataframe, such that each row is one replicate rack,
  # and each column represents the log AGB of one species in 
  # one soil background in that rack.
  df %>% mutate(pair = paste(plant_sp, soil_sp, decay_cat, sterilization, sep = "_")) %>%
    select(pair_id, pair, log_performance) %>%
    spread(pair, log_performance)
}
biomass_wide <- make_wide_biomass(biomass_for_SF)

# Calculate the Stabilization between each species pair ----
# Recall that following the definition of the m terms in Bever 1997,
# and the analysis of this model in Kandlikar 2019,
# stablization = -0.5*(log(m1A) + log(m2B) - log(m1B) - log(m2A))
# Here is a function that does this calculation 
# (recall that after the data reshaping above, 
# each column represents a given value of log(m1A))
calculate_stabilization <- function(df) {
  df %>% 
    # In df, each row is one repilcate/rack, and each column
    # represents the growth of one species in one soil type.
    # e.g. "FE_ACWR" is the growth of FE in ACWR-cultivated soil.
    mutate(freshE_freshM = -0.5*((E_E_fresh_live - E_E_fresh_sterile) - (E_M_fresh_live - E_M_fresh_sterile) - (M_E_fresh_live- M_E_fresh_sterile) + (M_M_fresh_live - M_M_fresh_sterile)),
           shortdE_shortdM = -0.5*((E_E_shortd_live - E_E_shortd_sterile) - (E_M_shortd_live - E_M_shortd_sterile) - (M_E_shortd_live - M_E_shortd_sterile) + (M_M_shortd_live - M_M_shortd_sterile)),
           longdE_longdM = -0.5*((E_E_longd_live - E_E_longd_sterile) - (E_M_longd_live - E_M_longd_sterile) - (M_E_longd_live - M_E_longd_sterile) + (M_M_longd_live - M_M_longd_sterile)),
           freshE_shortdM = -0.5*((E_E_fresh_live - E_E_fresh_sterile) - (E_M_shortd_live - E_M_shortd_sterile) - (M_E_fresh_live - M_E_fresh_sterile) + (M_M_shortd_live - M_M_shortd_sterile)),
           freshE_longdM = -0.5*((E_E_fresh_live - E_E_fresh_sterile) - (E_M_longd_live - E_M_longd_sterile) - (M_E_fresh_live - M_E_fresh_sterile) + (M_M_longd_live - M_M_longd_sterile)),
           shortdE_freshM = -0.5*((E_E_shortd_live - E_E_shortd_sterile) - (E_M_fresh_live - E_M_fresh_sterile) - (M_E_shortd_live - M_E_shortd_sterile) + (M_M_fresh_live - M_M_fresh_sterile)),
           shortdE_longdM = -0.5*((E_E_shortd_live - E_E_shortd_sterile) - (E_M_longd_live - E_M_longd_sterile) - (M_E_shortd_live - M_E_shortd_sterile) + (M_M_longd_live - M_M_longd_sterile)),
           longdE_freshM = -0.5*((E_E_longd_live - E_E_longd_sterile) - (E_M_fresh_live - E_M_fresh_sterile) - (M_E_longd_live - M_E_longd_sterile) + (M_M_fresh_live - M_M_fresh_sterile)),
           longdE_shortdM = -0.5*((E_E_longd_live - E_E_longd_sterile) - (E_M_shortd_live - E_M_shortd_sterile) - (M_E_longd_live - M_E_longd_sterile) + (M_M_shortd_live - M_M_shortd_sterile))
)%>%
    select(pair_id, freshE_freshM:longdE_shortdM) %>% 
    gather(pair_id, stabilization, freshE_freshM:longdE_shortdM)
}

# Calculating the Fitness difference between each species pair ----
# Similarly, we can now calculate the fitness difference between
# each pair. Recall that FD = 0.5*(log(m1A)+log(m1B)-log(m2A)-log(m2B))
# But recall that here, IT IS IMPORTANT THAT
# m1A = (m1_soilA - m1_fieldSoil)!

# The following function does this calculation:
calculate_fitdiffs <- function(df) {
  df %>% 
    mutate(freshE_freshM = 0.5*((E_E_fresh_live-E_E_fresh_sterile) + (E_M_fresh_live - E_M_fresh_sterile) - (M_E_fresh_live-M_E_fresh_sterile) - (M_M_fresh_live-M_M_fresh_sterile)),
           shortdE_shortdM = 0.5*((E_E_shortd_live-E_E_shortd_sterile) + (E_M_shortd_live-E_M_shortd_sterile) - (M_E_shortd_live-M_E_shortd_sterile) - (M_M_shortd_live-M_M_shortd_sterile)),
           longdE_longdM = 0.5*((E_E_longd_live-E_E_longd_sterile) + (E_M_longd_live-E_M_longd_sterile) - (M_E_longd_live-M_E_longd_sterile) - (M_M_longd_live-M_M_longd_sterile)),
           freshE_shortdM = 0.5*((E_E_fresh_live-E_E_fresh_sterile) + (E_M_shortd_live-E_M_shortd_sterile) - (M_E_fresh_live-M_E_fresh_sterile) - (M_M_shortd_live-M_M_shortd_sterile)),
           freshE_longdM = 0.5*((E_E_fresh_live-E_E_fresh_sterile) + (E_M_longd_live-E_M_longd_sterile) - (M_E_fresh_live-M_E_fresh_sterile) - (M_M_longd_live-M_M_longd_sterile)),
           shortdE_freshM = 0.5*((E_E_shortd_live-E_E_shortd_sterile) + (E_M_fresh_live-E_M_fresh_sterile) - (M_E_shortd_live-M_E_shortd_sterile) - (M_M_fresh_live-M_M_fresh_sterile)),
           shortdE_longdM = 0.5*((E_E_shortd_live-E_E_shortd_sterile) + (E_M_longd_live-E_M_longd_sterile) - (M_E_shortd_live-M_E_shortd_sterile) - (M_M_longd_live-M_M_longd_sterile)),
           longdE_freshM = 0.5*((E_E_longd_live-E_E_longd_sterile) + (E_M_fresh_live-E_M_fresh_sterile) - (M_E_longd_live-M_E_longd_sterile) - (M_M_fresh_live-M_M_fresh_sterile)),
           longdE_shortdM = 0.5*((E_E_longd_live-E_E_longd_sterile) + (E_M_shortd_live-E_M_shortd_sterile) - (M_E_longd_live-M_E_longd_sterile) - (M_M_shortd_live-M_M_shortd_sterile))
           ) %>%
    select(pair_id, freshE_freshM:longdE_shortdM) %>% 
    gather(pair_id, fitdiff_fld, freshE_freshM:longdE_shortdM)
}

#pre-paired----
fd_values <- calculate_fitdiffs(biomass_wide)
fd_values <- fd_values %>% filter(!(is.na(fitdiff_fld)))
stabilization_values <- calculate_stabilization(biomass_wide)
# Seven values are NA; we can omit these 
stabilization_values <- stabilization_values %>% 
  filter(!(is.na(stabilization)))
# Now, generate statistical summaries of SD and FD and plotting----
SD_FD <- function(stabilization_values = stabilization_values, fd_values = fd_values){
stabiliation_summary <- stabilization_values %>% group_by(pair_id) %>%
  summarize(mean_sd = mean(stabilization),
            sem_sd = sd(stabilization)/sqrt(n()),
            n_sd = n())

fitdiff_summary <- fd_values %>% group_by(pair_id) %>%
  summarize(mean_fd = mean(fitdiff_fld),
            sem_fd = sd(fitdiff_fld)/sqrt(n()),
            n_fd = n())

# Combine the two separate data frames.
sd_fd_summary <- left_join(stabiliation_summary, fitdiff_summary) %>% mutate(pair_id = lapply(lapply(strsplit(pair_id, '_'), FUN = function(x){(paste(substr(x, 1, 1), str_sub(x, -1, -1)))}), FUN = function(x){paste0(x, collapse = '_')}))

# Some of the FDs are negative, let's flip these to be positive
# and also flip the label so that the first species in the name
# is always the fitness superior.
sd_fd_summary <- sd_fd_summary %>% 
  # if mean_fd is < 0, the following command gets the absolute
  # value and also flips around the species code so that
  # the fitness superior is always the first species in the code
  mutate(pair_id = ifelse(mean_fd < 0,
                       paste0(str_extract(pair_id, "...$"),
                              "_",
                              str_extract(pair_id, "^...")), 
                       pair_id),
         mean_fd = abs(mean_fd))

# Make a column that indicates the net outcome (coex. or exclusion)
sd_fd_summary <- sd_fd_summary %>% mutate(
  stabilize = ifelse(mean_sd - 2*sem_sd > 0, "yes", "no"), 
  fitness = ifelse(mean_fd - 2*sem_fd > 0, "yes", "no"), 
  outcome = ifelse(mean_fd - 2*sem_fd >
                     mean_sd + 2*sem_sd, "exclusion", "neutral"),
  outcome2 = ifelse(mean_fd > mean_sd, "exclusion", "coexistence"),
  outcome2 = ifelse(mean_sd < 0, "exclusion or priority effect", outcome2)
)


# Plot SD vs. FD for pre-paired data----
sd_fd_xpyplot <- ggplot(sd_fd_summary, aes(x = mean_sd, y = mean_fd, label = pair_id, fill = outcome)) +
   # geom_ribbon(aes(x = seq(0, 1.5, length.out = 15),
   #                 ymin = rep(0, 15),
   #                 ymax = seq(0, 1.5, length.out = 15)),
   #             fill = "#FFF2CC")  +
   
   ylim(c(-.6, 1)) +
   xlim(c(-.5, 1)) +
   geom_abline(linetype = 2) +
   # geom_label() +
   geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
   geom_point(shape = 21, size = 3, stroke = 1) +
   scale_fill_manual(values = c("#A7C084", "white"),
                     labels = c("Lower bound of fitness difference\nestimate is greater than upper\nbound of stabilization estimate",
                                "Confidence intervals of fitness\ndifference and stabilization\nestimates overlap"),
                     name = "") + 
  geom_errorbar(aes(ymin = (mean_fd)-2*sem_fd,
                    ymax = (mean_fd)+2*sem_fd),
                size = .35) +
  geom_errorbarh(aes(xmin = mean_sd-2*sem_sd,
                     xmax = mean_sd+2*sem_sd),
                 size = .35) +
  ggrepel::geom_text_repel(segment.size = 0.025, color = "#556B37", size = 6) +
   xlab(bquote(atop("Microbially mediated stabilization",
                    -frac(1,2)~(m["1A"]-m["1B"]-m["2A"]+m["2B"]) ))) +
   ylab(bquote(atop("Microbially mediated fitness difference",
                    frac(1,2)~(m["1A"]+m["1B"]-m["2A"]-m["2B"])))) + 
   annotate("text", x = 0.8, y = .2, label = "coexistence", 
            vjust = 0, hjust = 1, size = 10, fontface = "bold.italic") +
   # annotate("text", x = 1.27, y = .15, label = "(assuming otherwise", 
   #          hjust = 1, size = 2.75, color = "grey45", fontface = "italic") + 
   # annotate("text", x = 1.27, y = .07, label = "equal competitors)", 
   #          hjust = 1, size = 2.75, color = "grey45", fontface = "italic") + 
   
   annotate("text", x = 0.8, y = 0.95, label = "exclusion", 
            vjust = 0, hjust = 1, size = 10, fontface = "bold.italic") + 
   # annotate("text", x = 1.27, y = 1.45, label = "(assuming otherwise", 
   #          hjust = 1, size = 2.75, color = "grey45", fontface = "italic") + 
   # annotate("text", x = 1.27, y = 1.37, label = "equal competitors)", 
   #          hjust = 1, size = 2.75, color = "grey45", fontface = "italic") + 
   
   theme_gsk() +
   theme(axis.line = element_line(size = 0),
         legend.justification=c(1,0), legend.position=c(1, 0.025),
         legend.background = element_rect(colour = "grey50"),
         legend.key.height = unit(1, 'cm'),
         legend.text = element_text(size = 7.5),
         legend.title = element_text(size = 0)) +
  theme(text = element_text(size = 20))
   NULL
  print(sd_fd_xpyplot)
  return(sd_fd_summary)
}
SD_FD(stabilization_values, fd_values)
#Sampling ----
#Data generation for rnorm and bootstrap----
#rnorm
biomass_wide_st <- cbind(apply(biomass_wide[, -1], 2, FUN = function(x){mean(x, na.rm = T)}), apply(biomass_wide[, -1], 2, FUN = function(x){sd(x, na.rm = T)}))
colnames(biomass_wide_st) <- c('mean', 'sd')
rownames(biomass_wide_st)

n_reps <- 1000
vec <- c()
for(i in rownames(biomass_wide_st)){
  set.seed(31415)
  vec <- c(vec, rnorm(n_reps, mean = biomass_wide_st[i, 1], 
        sd = biomass_wide_st[i, 2]))
}
mat <- matrix(vec, ncol = 24, nrow = n_reps, byrow = F)
colnames(mat) <- rownames(biomass_wide_st)
mat <- as.tibble(mat) %>% mutate(pair_id = 1:n_reps)

#bootstrap
n_reps <- 1000
vec <- c()
for (i in colnames(biomass_wide)[-1]){
  set.seed(31415)
  for (j in 1:1000){
    vec <- c(vec, mean(sample((na.exclude(unlist(unname(as.vector(biomass_wide[, i]))))), 10, replace = T), na.rm = T))
  }
}
mat <- matrix(vec, ncol = 24, nrow = n_reps, byrow = F)
colnames(mat) <- colnames(biomass_wide)[-1]
mat <- as.tibble(mat) %>% mutate(pair_id = 1:n_reps)

# Calculate SD & FD----
stabilization_values_norm <- calculate_stabilization(mat)
fd_values_norm <- calculate_fitdiffs(mat)

# Use newly defined functions to calculate SD and FD
SD_vec_by_species <- split(stabilization_values_norm, stabilization_values_norm$pair_id) %>% lapply(function(x) { x["pair_id"] <- NULL; x })
FD_vec_by_species <- split(fd_values_norm, fd_values_norm$pair_id) %>% lapply(function(x) { x["pair_id"] <- NULL; x })

SD_vec_by_run <- split(stabilization_values_norm, rownames(stabilization_values_norm)) %>% lapply(function(x) { x["pair_id"] <- NULL; x })
FD_vec_by_run <- split(fd_values_norm, rownames(fd_values_norm)) %>% lapply(function(x) { x["pair_id"] <- NULL; x })

compare <- function(FD_vec, SD_vec) {
  exclusion <- sum(abs(FD_vec) > abs(SD_vec))
  coexistence <- sum(abs(FD_vec) < SD_vec)
  priority <- nrow(FD_vec)-(exclusion+coexistence)
  return(c(n_reps_exclusion = exclusion, n_reps_coex = coexistence, 
           n_reps_pes = priority))
}

outcomes_by_pair <- t(mapply(compare, FD_vec_by_species, SD_vec_by_species))
outcomes_by_run <- t(mapply(compare, FD_vec_by_run, SD_vec_by_run))
# to_return <- cbind(to_return, outcomes_by_pair)

#Function setting for Chesson's plot with two-color ribbon----
MCT_plot_mat <- stabilization_values_norm %>% filter(pair_id %in% c("freshE_freshM", "shortdE_shortdM", "longdE_longdM")) %>% cbind(fd_values_norm %>% filter(pair_id %in% c("freshE_freshM", "shortdE_shortdM", "longdE_longdM")) %>% select(-"pair_id"))

#Pie plot----
fd_values_st <- fd_values_norm %>% group_by(pair_id) %>% summarise( mean_FD = mean(fitdiff_fld))
stabilization_values_st <- stabilization_values_norm %>% group_by(pair_id) %>% summarise(mean_SD = mean(stabilization))
live_ref_sims <- outcomes_by_pair %>% cbind(stabilization_values_st,  fd_values_st) %>% select(-c(4, 6))

#filter version
live_ref_sims <- live_ref_sims[c(1, 5, 9), ]
rownames(live_ref_sims) <- c("Alive E_Alive M", "Long decay E_Long Decay M", "Short decay E_Short decay M")

# Make a base plot that can be modified for each sub question  
base_coex_plot <- 
  ggplot() +
  geom_ribbon(aes(x = seq(0, 5, length.out = 100),
                  ymin = rep(0, length.out = 100),
                  ymax = seq(0, 5, length.out = 100)),
              fill = alpha("#CCEAC9", 1)) + #99D594 d9ead3ff
  geom_ribbon(aes(x = seq(0, -5, length.out = 100),
                  ymin = rep(0, length.out = 100),
                  ymax = seq(0, 5, length.out = 100)),
              fill = alpha("#FFFFBF", 1)) +
  annotate("segment", x = 0, y = 0, xend = 5, yend = 5, linetype = 2) +
  annotate("segment", x = 0, y = 0, xend = -5, yend = 5, linetype = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab(bquote(atop("",
                   -frac(1,2)~(m["1A"]-m["1B"]-m["2A"]+m["2B"]) ))) +
  ylab(bquote(atop("Fitness difference",
                   frac(1,2)~(m["1A"]+m["1B"]-m["2A"]-m["2B"])))) +
  annotate("segment", x = 0.2, xend = 1.25, y = -.12, yend = -.12, 
           colour = "grey25", size = 0.5, arrow = arrow(length = unit(0.075, "inches"))) + 
  annotate("segment", x = -0.2, xend = -1.25, y = -.12, yend = -.12, 
           colour = "grey25", size = 0.5, arrow = arrow(length = unit(0.075, "inches"))) +
  theme(axis.line = element_blank())

# define the data range of the plot
x_min <- with(live_ref_sims, min(mean_SD))-0.25
x_max <- with(live_ref_sims, max(mean_SD))+0.1
y_min <- with(live_ref_sims, min(abs(mean_FD)))-0.05
y_max <- with(live_ref_sims, max(abs(mean_FD)))+0.15

# adapt a previously defined base plot
base_coex_plot1 <- base_coex_plot
base_coex_plot1$layers[7:8] <- NULL
base_coex_plot1 <- base_coex_plot1 +   
  annotate("text", x = 0.1, y = -0.03, hjust = 0, vjust = 0, label = "Stabilization", color = "black", size = 4.8) +
  annotate("text", x = -0.01, y = -0.03, hjust = 1, vjust = 0, label = "Destabilization",  color = "black", size = 4.8) + 
  theme(axis.title.x = element_text(margin = margin(t = 0)))
# base_coex_plot1

# plotting
(scatterpie_FL <- base_coex_plot1 +
    scale_x_continuous(limits = c(x_min, x_max)) +
    scale_y_continuous(limits = c(y_min-0.05, y_max)) +
    annotate("segment", x = 0, y = 0, xend =  min(x_max, y_max),
             yend = min(x_max, y_max), 
             linetype = 2) +
    annotate("segment", x = 0, y = 0, xend =  max(x_min, -y_max),
             yend =  min(-x_min, y_max),
             linetype = 2) +
    geom_scatterpie(aes(x=mean_SD, y=abs(mean_FD), r = 0.015), 
                    data=live_ref_sims, size = 0.15,
                    cols=c("n_reps_exclusion", "n_reps_coex", "n_reps_pes")) +
    geom_point(aes(x=mean_SD, y=abs(mean_FD)), 
               data=live_ref_sims, size = 0.1)+
    geom_text(aes(x=mean_SD, y=abs(mean_FD), label = rownames(live_ref_sims)), data=live_ref_sims, size = 4)+
    scale_fill_manual(values = c( "#f1f4f9", "#48a3a5", "#f2c949"), labels =  c( "exclusion", "coexistence", "priority effect"), name = "") +
    # guides(fill = "none") +
    # labs(title = "(A) Species pairs with reference growth in live soil") + 
    coord_fixed() +
    theme(axis.line = element_line(size = 0),
          legend.position = "bottom",
          legend.background = element_rect(colour = "grey50", size = .3),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 8, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 10))
)
