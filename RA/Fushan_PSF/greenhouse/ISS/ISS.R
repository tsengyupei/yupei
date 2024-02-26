library(patchwork)
library(PerformanceAnalytics)
library(ggplot2)
library(scatterpie)
library(gtools)
library(rpsychi)
library(tidyverse)
library(stringr)
library(ggpubr)
library(rstatix)
library(dplyr)
library(ggnewscale)
library(ggrepel)
#0. Read in the data and clean out empty rows. ------
#ISS data----
ISS <- read.csv("/Users/yupeitseng/Documents/RA/FS_PSF/green_house/ISS/FS_greenhouse_harvest - ISS.csv")
ISS_height <- read.csv("/Users/yupeitseng/Documents/RA/FS_PSF/green_house/ISS/FS_PSF greenhouse - ISS_height.csv")

ISS_height <- ISS_height %>% filter(ISS_height[, 18] != "Dead") %>% filter(Wilt == 0) %>% mutate_at(c(8:18), as.numeric) %>%  mutate(growth_rate = X6.19 - X12.26)

ISS <- ISS %>% select(1:5, 13:14) %>% mutate(above_mass = as.numeric(above_mass), below_mass = as.numeric(below_mass), total_mass = (above_mass+below_mass))

ISS <- left_join(ISS, ISS_height[, c(1, 2, 18:20)])
colnames(ISS)[10] <- 'height'
ISS <- ISS %>% mutate(root_shoot_ratio = below_mass/above_mass, plant_sp = str_sub(plant_sp, 1, 1), soil_sp = str_sub(soil_sp, 1, 1))
colnames(ISS)[11] <- "apical_wilt"
ISS$pair_id <- rep(rep(letters[11:15], each = 2), 6)
ISS <- ISS[, c(1:5, 14, 9, 6:8, 10:13)]

#MSS data----
biomass_FS <- read.csv("/Users/yupeitseng/Documents/RA/FS_PSF/green_house/MSS/FS_greenhouse_harvest - MSS.csv")
height <- read.csv('/Users/yupeitseng/Documents/RA/FS_PSF/green_house/FS_PSF greenhouse - MSS_height.csv')

height <- height %>% mutate_at(c(9:19), as.numeric) %>%  mutate(growth_rate = X5.29 - X12.26)

biomass_FS <- biomass_FS %>% select(c(1:7, 15:18)) %>% mutate(replicate = c(rep(1:10, 24), rep(1:5, 2)), above_mass = as.numeric(above_mass), below_mass = as.numeric(below_mass), total_mass = (above_mass+below_mass)) %>% select(-apical_wilt)

biomass_FS <- left_join(biomass_FS, height[, c(5, 19:21)])
colnames(biomass_FS)[13] <- 'height'
biomass_FS <- biomass_FS %>% mutate(root_shoot_ratio = below_mass/above_mass)

#combine MSS and ISS data----
MSS <- biomass_FS %>% mutate(exp = rep("MSS", 250)) %>% filter(sterilization == "live" | soil_sp != "bg") %>% select(-c("sterilization", "rosetting", "replicate"))
ISS <- ISS %>% mutate(exp = rep("ISS", 60))
combine <- rbind(MSS, ISS)

# 01. Make a box plot for compare ISS and MSS----
response <- c("total_mass", "above_mass", "below_mass", "height", "growth_rate")
y_lab <- c("Total mass", "Above mass", "Below mass", "Height", "Growth")

pdf('/Users/yupeitseng/Documents/RA/FS_PSF/green_house/ISS/Compare_MSS_ISS.pdf', height = 8, width = 12)
for (i in response){
  for(j in y_lab){
    if (which(response == i) == which(y_lab == j)){
      
      combine_st <- combine %>% select(c("exp", "soil_sp", "decay_cat", "plant_sp", i)) 
      colnames(combine_st)[5] <- 'performance'
      combine_st$soil_sp <- sapply(combine_st$soil_sp, FUN = function(x){
        ifelse(x == "E", "Soil_E", "Soil_M")
      })
      combine_st$plant_sp <- sapply(combine_st$plant_sp, FUN = function(x){
        ifelse(x == "E", "Plant_E", "Plant_M")
      })
      
      stat.test <- combine_st %>%
        group_by(plant_sp, soil_sp, decay_cat) %>%
        t_test(performance ~ exp) %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance()
      stat.test 
      
      # Create box plots with significance levels
      # Hide ns (non-significant)
      stat.test <- stat.test %>% add_xy_position(x = "decay_cat")
      plot <- ggboxplot(
        combine_st, x = "decay_cat", y = "performance", color = "exp",
        facet = c("plant_sp", "soil_sp")
      ) +
        stat_pvalue_manual(stat.test, hide.ns = TRUE)+
        scale_color_manual(values=c("#A7C084", "#B6A87C"), "Experiment")+
        labs(x = "Decay category", y = j)+
        theme(text = element_text(size = 20))
      print(plot)
    }
  }
}
dev.off()

#02. Reshaping data for calculating Stabilization and FDs ----
# So first, let's make a column that is simply the log AGB
# of each pot, and then go from there. 
#MSS data
biomass_for_SF <- biomass_FS %>% filter(soil_sp != 'bg' & sterilization == "live") %>% select("plant_sp", "soil_sp", "decay_cat", "pair_id", "total_mass") 
colnames(biomass_for_SF)[5] <- 'performance'
biomass_for_SF <- biomass_for_SF %>% mutate(log_performance = log10(performance))

# ISS data
biomass_for_SF_ISS <- ISS %>% select("plant_sp", "soil_sp", "decay_cat", "pair_id", "total_mass") 
colnames(biomass_for_SF_ISS)[5] <- 'performance'
biomass_for_SF_ISS <- biomass_for_SF_ISS %>% mutate(log_performance = log10(performance))

# Now, we need to convert the "long format" biomass data frame  
# into a "wide format". Here's a function for that:
make_wide_biomass <- function(df) {
  # This takes the biomass dataframe, and converts into a wide
  # dataframe, such that each row is one replicate rack,
  # and each column represents the log AGB of one species in 
  # one soil background in that rack.
  biomass_for_SF %>% mutate(pair = paste(plant_sp, soil_sp, decay_cat, sep = "_")) %>%
    select(pair_id, pair, log_performance) %>%
    spread(key = pair, value = log_performance)
}
biomass_wide <- make_wide_biomass(biomass_for_SF)
biomass_wide_ISS <- make_wide_biomass(biomass_for_SF_ISS)
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
    mutate(freshE_freshM = -0.5*((E_E_fresh) - (E_M_fresh) - (M_E_fresh) + (M_M_fresh)),
           shortdE_shortdM = -0.5*((E_E_shortd) - (E_M_shortd) - (M_E_shortd) + (M_M_shortd)),
           longdE_longdM = -0.5*((E_E_longd) - (E_M_longd) - (M_E_longd) + (M_M_longd)),
           freshE_shortdM = -0.5*((E_E_fresh) - (E_M_shortd) - (M_E_fresh) + (M_M_shortd)),
           freshE_longdM = -0.5*((E_E_fresh) - (E_M_longd) - (M_E_fresh) + (M_M_longd)),
           shortdE_freshM = -0.5*((E_E_shortd) - (E_M_fresh) - (M_E_shortd) + (M_M_fresh)),
           shortdE_longdM = -0.5*((E_E_shortd) - (E_M_longd) - (M_E_shortd) + (M_M_longd)),
           longdE_freshM = -0.5*((E_E_longd) - (E_M_fresh) - (M_E_longd) + (M_M_fresh)),
           longdE_shortdM = -0.5*((E_E_longd) - (E_M_shortd) - (M_E_longd) + (M_M_shortd))
    )%>%
    select(pair_id, freshE_freshM:longdE_shortdM) %>% 
    gather(pair_id, stabilization, freshE_freshM:longdE_shortdM)
}

#Sampling ----
#Data generation for rnorm and bootstrap----
#rnorm
stab_norm <- function(biomass_wide){
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
  mat <- matrix(vec, ncol = 12, nrow = n_reps, byrow = F)
  colnames(mat) <- rownames(biomass_wide_st)
  mat <- as.tibble(mat) %>% mutate(pair_id = 1:n_reps)
  
  stabilization_values_norm <- calculate_stabilization(mat)
  return(stabilization_values_norm)
}
stab_norm_MSS <- stab_norm(biomass_wide  = biomass_wide)
stab_norm_ISS <- stab_norm(biomass_wide  = biomass_wide_ISS)
stab_norm_combine <- rbind(stab_norm_MSS, stab_norm_ISS)
stab_norm_combine$exp <- c(rep("MSS", 1, 9000), rep("ISS", 1, 9000))

stat.test <- stab_norm_combine %>%
  group_by(pair_id) %>%
  t_test(stabilization ~ exp) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

# Create box plots with significance levels
# Hide ns (non-significant)
stat.test <- stat.test %>% add_xy_position(x = "exp")
ggboxplot(
  stab_norm_combine, x = "exp", y = "stabilization", color = "exp",
  facet = c("pair_id")
) +
  stat_pvalue_manual(stat.test, hide.ns = TRUE)+
  scale_color_manual(values=c("#B6A87C", "#A7C084"), "Experiment")+
  labs(x = "Experiment", y = "Stabilization")+
  theme(text = element_text(size = 20))


#bootstrap
stab_boot <- function(biomass_wide){
  n_reps <- 1000
  vec <- c()
  for (i in colnames(biomass_wide)[-1]){
    set.seed(31415)
    for (j in 1:1000){
      vec <- c(vec, mean(sample((na.exclude(unlist(unname(as.vector(biomass_wide[, i]))))), 10, replace = T), na.rm = T))
    }
  }
  mat <- matrix(vec, ncol = 12, nrow = n_reps, byrow = F)
  colnames(mat) <- colnames(biomass_wide)[-1]
  mat <- as.tibble(mat) %>% mutate(pair_id = 1:n_reps)
  
  stabilization_values_norm <- calculate_stabilization(mat)
  return(stabilization_values_norm)
}

stab_boot_MSS <- stab_boot(biomass_wide  = biomass_wide)
stab_boot_ISS <- stab_boot(biomass_wide  = biomass_wide_ISS)
stab_boot_combine <- rbind(stab_boot_MSS, stab_boot_ISS)
stab_boot_combine$exp <- c(rep("MSS", 1, 9000), rep("ISS", 1, 9000))

stat.test <- stab_boot_combine %>%
  group_by(pair_id) %>%
  t_test(stabilization ~ exp) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

# Create box plots with significance levels
# Hide ns (non-significant)
stat.test <- stat.test %>% add_xy_position(x = "exp")
ggboxplot(
  stab_boot_combine, x = "exp", y = "stabilization", color = "exp",
  facet = c("pair_id")
) +
  stat_pvalue_manual(stat.test, hide.ns = TRUE)+
  scale_color_manual(values=c("#B6A87C", "#A7C084"), "Experiment")+
  labs(x = "Experiment", y = "Stabilization")+
  theme(text = element_text(size = 20))

#03. Checking the pH variation of ISS----
load("/Users/yupeitseng/Documents/RA/FS_PSF/data/fs.elev.RData")
soil <- read.csv("/Users/yupeitseng/Documents/RA/FS_PSF/data/fungi/field_soil_property_data.csv")
ISS_height <- read.csv("/Users/yupeitseng/Documents/RA/FS_PSF/green_house/ISS/FS_PSF greenhouse - ISS_height.csv")
ISS_height[which(ISS_height$soil_sample_tag == "171583"), 2] <- c("171563", "171563")
ISS_height[which(ISS_height$soil_sample_tag == "240230"), 2] <- c("243230", "243230")
ISS_pH <- left_join(ISS, ISS_height[, 1:2], by = "pot_id")
ISS_pH <- ISS_pH %>% select(-7) %>% rename(tag = soil_sample_tag.y)
ISS_pH$tag <- as.numeric(ISS_pH$tag)
ISS_pH <- left_join(ISS_pH, soil[, c(1, 4:12)], by = "tag")
ISS_pH$sp_decay <- paste(ISS_pH$soil_sp, ISS_pH$decay_cat, sep = "_")

#ISS pH
ggboxplot(
  ISS_pH, x = "decay_cat", y = "pH", color = "soil_sp"
  
) +
  # stat_pvalue_manual(stat.test, hide.ns = TRUE)+
  scale_color_manual(values=c("#A7C084", "#B6A87C"), "Species")+
  labs(x = "Decay category", y = "pH")+
  theme(text = element_text(size = 20))

#ISS
ggplot()+
  geom_contour(fs.elev$col, mapping = aes(x/20, y/20, z = elev), bins = 40, alpha = 0.5,colour = "black")+
  labs(x = "East-West distance (m)", y = "North-South distance (m)", shape = "Decay category")+
  geom_point(data = ISS_pH, aes(x = gx_new, y = gy_new, color = soil_sp))+
  scale_color_manual(values=c("#A7C084", "#B6A87C"), "Species")+
  new_scale_color() +
  geom_point(data = ISS_pH, aes(x = gx_new, y = gy_new, color = sp_decay, shape = decay_cat, size = 10))+
  geom_text_repel(data = ISS_pH, aes(x = gx_new, y = gy_new, label = ISS_pH$sample_ID), vjust = 0, nudge_y = 0.5)+
  scale_shape_manual(values=c(15, 16, 17))+
  scale_color_manual(values=c("#556B37", "#A7C084", "#CAD9B5", "#7D6241", "#B6A87C", "#D3C1AB"))+
  guides(color = FALSE, size = FALSE)+
  theme(text = element_text(size = 20)) 

#MSS pH
ggboxplot(
  filter(soil, for_soil_mix == "Yes"), x = "decay_cat", y = "pH", color = "sp", order = c("Alive", "2014_2019", "2009_2014")) +
  scale_x_discrete(labels = c("fresh", "shortd", "longd"))+
  # stat_pvalue_manual(stat.test, hide.ns = TRUE)+
  scale_color_manual(values=c("#A7C084", "#B6A87C"), "Species", labels = c("E", "M"))+
  labs(x = "Decay category", y = "pH")+
  theme(text = element_text(size = 20))

#MSS
soil$sp_decay <- paste(soil$sp, soil$decay_cat, sep = "_")
ggplot()+
  geom_contour(fs.elev$col, mapping = aes(x/20, y/20, z = elev), bins = 40, alpha = 0.5,colour = "black")+
  labs(x = "East-West distance (m)", y = "North-South distance (m)", shape = "Decay category")+
  geom_point(data = ISS_pH, aes(x = gx_new, y = gy_new, color = soil_sp))+
  scale_color_manual(values=c("#A7C084", "#B6A87C"), "Species")+
  new_scale_color() +
  geom_point(data =filter(soil, for_soil_mix == "Yes"), aes(x = gx_new, y = gy_new, color = sp_decay, shape = decay_cat, size = 10))+
  geom_text_repel(data = filter(soil, for_soil_mix == "Yes"), aes(x = gx_new, y = gy_new, label = filter(soil, for_soil_mix == "Yes")$sample_ID), vjust = 0, nudge_y = 0.5)+
  scale_shape_manual(values=c(15, 16, 17), labels = c("Long decay", "Short decay", "Alive"))+
  scale_color_manual(values=c("#556B37", "#A7C084", "#CAD9B5", "#7D6241", "#B6A87C", "#D3C1AB"))+
  guides(color = FALSE, size = FALSE)+
  theme(text = element_text(size = 20)) 
  