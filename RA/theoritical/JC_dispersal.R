#### Load library
library("deSolve")
library("tidyverse")
library(pracma)


#Dispersal----
#with times <- seq(0, 3000, by = 0.02) in the function!!!
Data_dispersal <-
  crossing(
    # d = c(0, 2, 5, 8, 10),
    d = seq(0, 5, 0.1),
    ga1 = 1.6,
    ga2 = 1, 
    gb1 = 1,
    gb2 = 2, 
    # alpha1 = c(1, 0.2, 0.6, 1.1, 1.5, 0), 
    # alpha2 = c(1, 0.5, 0.9, 1.4, 1.8, 0),
    #competition 
    alpha1 = c(1.1, 1.5, 1.3, 1.1, 1.1, 1.6, 1.7),
    alpha2 = c(1.4, 1.8, 1.6, 1.2, 1.3, 1.9, 2.0),
    phi1 = 10, 
    phi2 = 10, 
    c11 = 1.5, 
    c22 = 0.8, 
    c12 = 1, 
    c21 = 0.8, 
    sigma1a = -6, 
    sigma1b = -0.4, 
    sigma2a = -0.5, 
    sigma2b = -6,
    r1 = 0.1, r2 = 0.1
  ) %>%
  # mutate(Par = pmap(., list)) %>%
  # filter(sigma1a == sigma2b) %>% 
  mutate(competition = paste(alpha1, alpha2, sep = ",")) %>% 
  # filter(competition == "1,1" | competition == "0.2,0.5" | competition == "0.6,0.9" | competition == "1.1,1.4"| competition == "1.5,1.8"| competition == "0,0") %>% 
  filter(competition == "1.3,1.6" | competition == "1.1,1.2" | competition == "1.1,1.3" | competition == "1.1,1.4"| competition == "1.5,1.8" | competition == "1.6,1.9" | competition == "1.7,2.0") %>%
  select(-competition) %>% 
  mutate(Par = pmap(., c)) %>%
  mutate(Scenario = paste("dispersal = ", d, ", competition = ", paste(alpha1, alpha2, sep = ","))) %>%
  mutate(Sim = map(Par, SimulateODE)) %>%
  mutate(Result = map2(Scenario, Sim, bind_cols)) %>%
  select(d, alpha1, alpha2, Result) %>%
  unnest(Result) %>%
  rename("Scenario" = "...1") 
  # gather(key = Variable, value = Density, -c(Scenario, time))

# setwd("/home/po/Documents/Yupei/patch_feedback_model/data")
# save(Data_dispersal, file = "JC_data_dispersal.RData")

load("/home/po/Documents/Yupei/patch_feedback_model/data/JC_data_dispersal.RData")

png(filename="/home/po/Documents/Yupei/patch_feedback_model/plot/JC_dispersal_1500_3000.png", width = 1600, height = 1200)
Data_dispersal %>% 
  gather(key = Variable, value = Density, Ma1:P2) %>% 
  filter(d %in% seq(from = 0.5, to = 3.5, by = 0.1)) %>%
  filter(time %in% seq(from = 1500, to = 3000, by = 0.02)) %>%
  group_by(d, alpha1, alpha2) %>%
  ggplot() +
  geom_line(aes(x = time, 
                y = Density, 
                color = Variable))+
  theme_classic()+
  # geom_text(x = 1875, y = 0.2, aes(label = text), data = coexsist %>% filter(m_effect == "m_effect =  -0.32"))+
  theme(text = element_text(size = 15))+
  facet_grid(d~alpha1 + alpha2)
dev.off()



png(filename="/home/po/Documents/Yupei/patch_feedback_model/plot/JC_dispersal_bifurcation_PM.png", width = 1600, height = 1200)
Data_dispersal %>% 
  filter(time %in% seq(from = 1500, to = 3000, by = 0.02)) %>%
  group_by(d, alpha1, alpha2) %>% 
  summarise(across(Ma1:P2, list(min = min, max = max))) %>% 
  mutate(competition = paste(alpha1, alpha2, sep = ",")) %>% 
  pivot_longer(cols = Ma1_min:P2_max, names_to = "Variable", values_to = "Density") %>% 
  mutate(value = sub(".*_", "", Variable), variable = sub("\\_.*", "", Variable)) %>% 
  # filter(variable %in% c("P1", "P2")) %>%
  # filter(variable %in% c("Ma1", "Ma2", "Mb1", "Mb2")) %>%
  ggplot(aes(x = d, 
             y = Density)) +
  geom_line(aes(linetype = value, color = variable))+
  theme_classic()+
  theme(text = element_text(size = 20))+
  facet_wrap(~competition, ncol = 1, nrow = 6)
dev.off()

peak <- Data_dispersal %>% 
  filter(time %in% seq(from = 1500, to = 3000, by = 0.02)) %>%
  filter(d == 0.8 & alpha1 == 1.1 & alpha2 == 1.2) %>% 
  select(Mb1)
peak <- findpeaks(unlist(peak), npeaks = 1000)
peak[-c(1:2, (length(peak[,1])-1):length(peak[,1])), ]           
