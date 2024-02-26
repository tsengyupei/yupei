#### Load library
library("deSolve")
library("tidyverse")


#### Model specification
patch_feedback_model <- function(times, state, parms) {
  with(as.list(c(state, parms, times)), {
    dMa1_dt = Ma1*(ga1-((Ma1 + alpha1*Mb1)/((phi1 * P1) + 0.000000001))) - d*Ma1*P2 + d*Ma2*P1
    dMa2_dt = Ma2*(ga2-((Ma2 + alpha2*Mb2)/((phi2 * P2) + 0.000000001))) - d*Ma2*P1 + d*Ma1*P2
    dMb1_dt = Mb1*(gb1-((Mb1 + alpha2*Ma1)/((phi1 * P1) + 0.000000001))) - d*Mb1*P2 + d*Mb2*P1
    dMb2_dt = Mb2*(gb2-((Mb2 + alpha1*Ma2)/((phi2 * P2) + 0.000000001))) - d*Mb2*P1 + d*Mb1*P2
    
    dP1_dt = P1*(r1 - c11*P1 - c12*P2 + sigma1a*Ma1 + sigma1b*Mb1)
    dP2_dt = P2*(r2 - c22*P2 - c21*P1 + sigma2a*Ma2 + sigma2b*Mb2)
    
    return(list(c(dMa1_dt, dMa2_dt, dMb1_dt, dMb2_dt, dP1_dt, dP2_dt)))
  })
}


#### Prepare functions
SimulateODE <- function(parameters){
  times <- seq(0, 2000, by = 0.02)
  state <- c(Ma1 = 0.1, Ma2 = 0.1, Mb1 = 0.1, Mb2 = 0.1, P1 = 0.1, P2 = 0.1)
  # parms <- unlist(parameters)
  parms <- parameters
  pop_size <- ode(func = patch_feedback_model, times = times, y = state, parms = parms)
  return(pop_size)
}

#Jazan-Conell----
#### Create data using functional programming
Data <-
  crossing(
    d = c(0, 2, 5, 8, 10),
    ga1 = 1.6,
    ga2 = 1, 
    gb1 = 1,
    gb2 = 2, 
    alpha1 = c(1, 0.2, 0.6, 1.1, 1.5, 0),
    alpha2 = c(1, 0.5, 0.9, 1.4, 1.8, 0),
    phi1 = 0.1, 
    phi2 = 0.1, 
    c11 = 1, 
    c22 = 0.8, 
    c12 = 1, 
    c21 = 0.8, 
    sigma1a = c(-0.32, -6), 
    sigma1b = -0.4, 
    sigma2a = -0.5, 
    sigma2b = c(-0.32, -6),
    r1 = 0.1, r2 = 0.1
  ) %>%
  # mutate(Par = pmap(., list)) %>%
  filter(sigma1a == sigma2b) %>% 
  mutate(competition = paste(alpha1, alpha2, sep = ",")) %>% 
  filter(competition == "1,1" | competition == "0.2,0.5" | competition == "0.6,0.9" | competition == "1.1,1.4"| competition == "1.5,1.8"| competition == "0,0") %>%
  select(-competition) %>% 
  mutate(Par = pmap(., c)) %>%
  mutate(Scenario = paste("dispersal = ", d, ", competition = ", paste(alpha1, alpha2, sep = ","), ", m_effect = ", sigma1a)) %>%
  mutate(Sim = map(Par, SimulateODE)) %>%
  mutate(Result = map2(Scenario, Sim, bind_cols)) %>%
  select(d, alpha1, alpha2, sigma1a, Result) %>%
  unnest(Result) %>%
  rename("Scenario" = "...1")
save(Data, file = "/home/po/Documents/Yupei/patch_feedback_model/data/JC_data.RData")
load("/home/po/Documents/Yupei/patch_feedback_model/data/JC_data.RData")
#Coexsistence
coexsist <- Data %>% 
  filter(time %in% seq(from = 1500, to = 2000, by = 0.02)) %>% 
  gather(key = Variable, value = Density, Ma1:P2) %>% 
  group_by(d, alpha1, alpha2, sigma1a, Scenario, Variable) %>% 
  summarise(Mean_density = mean(Density)) %>% 
  pivot_wider(names_from = Variable, values_from = Mean_density) %>% 
  mutate(d = as.factor(d), alpha1 = as.factor(alpha1), alpha2 = as.factor(alpha2), sigma1a= as.factor(sigma1a))
coexsist[coexsist < 0.00000001] <- 0

text <- apply(coexsist[, 6:11], 1, function(x){names(which(x == 0))})
vec_text <- c()
for (i in 1:nrow(coexsist)) {
  vec_text <- c(vec_text, paste(unlist(text[[i]]), collapse = " "))
}
vec_text[which(vec_text == "")] <- "coexist"
coexsist <- coexsist %>% as.data.frame() %>% mutate(result = vec_text) 

#### Plot data
#all
Strong <- Data %>% 
  filter(sigma1a ==  -6)%>% 
  gather(key = Variable, value = Density, Ma1:P2) %>% 
  # filter(time %in% seq(from = 0, to = 10, by = 0.02)) %>%
  # filter(time %in% seq(from = 1500, to = 2000, by = 0.02)) %>%
  ggplot() +
  geom_line(aes(x = time, 
                y = Density, 
                color = Variable))+
  theme_classic()+
  geom_text(x = 1000, y = 0.1, aes(label = result), data = coexsist %>% filter(sigma1a  == -6))+
  theme(text = element_text(size = 15))+
  # facet_grid(factor(Dispersal, levels = c("dispersal =  0 ", "dispersal =  2 ", "dispersal =  5 ", "dispersal =  8 ", "dispersal =  10 "))~Competition)
  facet_grid(d~alpha1+alpha2)

png(filename="/home/po/Documents/Yupei/patch_feedback_model/plot/JC_strong.png", width = 1600, height = 1200)
Strong
dev.off()


Weak <- Data %>% 
  filter(sigma1a  == -0.32)%>% 
  gather(key = Variable, value = Density, Ma1:P2) %>% 
  # filter(time %in% seq(from = 0, to = 10, by = 0.02)) %>%
  # filter(time %in% seq(from = 1500, to = 2000, by = 0.02)) %>%
  ggplot() +
  geom_line(aes(x = time, 
                y = Density, 
                color = Variable))+
  theme_classic()+
  geom_text(x = 1000, y = 0.1, aes(label = result), data = coexsist %>% filter(sigma1a  == -0.32))+
  theme(text = element_text(size = 15))+
  # facet_grid(factor(Dispersal, levels = c("dispersal =  0 ", "dispersal =  2 ", "dispersal =  5 ", "dispersal =  8 ", "dispersal =  10 "))~Competition)
  facet_grid(d~alpha1+alpha2)

png(filename="/home/po/Documents/Yupei/patch_feedback_model/plot/JC_weak.png", width = 1600, height = 1200)
Weak
dev.off()

#1500-2000
Strong <- Data %>% 
  filter(sigma1a ==  -6)%>% 
  gather(key = Variable, value = Density, Ma1:P2) %>% 
  # filter(time %in% seq(from = 0, to = 10, by = 0.02)) %>%
  filter(time %in% seq(from = 1500, to = 2000, by = 0.02)) %>%
  ggplot() +
  geom_line(aes(x = time, 
                y = Density, 
                color = Variable))+
  theme_classic()+
  geom_text(x = 1750, y = 0.02, aes(label = result), data = coexsist %>% filter(sigma1a  == -6))+
  theme(text = element_text(size = 15))+
  # facet_grid(factor(Dispersal, levels = c("dispersal =  0 ", "dispersal =  2 ", "dispersal =  5 ", "dispersal =  8 ", "dispersal =  10 "))~Competition)
  facet_grid(d~alpha1+alpha2)

png(filename="/home/po/Documents/Yupei/patch_feedback_model/plot/JC_strong_1500_2000.png", width = 1600, height = 1200)
Strong
dev.off()


Weak <- Data %>% 
  filter(sigma1a  == -0.32)%>% 
  gather(key = Variable, value = Density, Ma1:P2) %>% 
  # filter(time %in% seq(from = 0, to = 10, by = 0.02)) %>%
  filter(time %in% seq(from = 1500, to = 2000, by = 0.02)) %>%
  ggplot() +
  geom_line(aes(x = time, 
                y = Density, 
                color = Variable))+
  theme_classic()+
  geom_text(x = 1750, y = 0.1, aes(label = result), data = coexsist %>% filter(sigma1a  == -0.32))+
  theme(text = element_text(size = 15))+
  # facet_grid(factor(Dispersal, levels = c("dispersal =  0 ", "dispersal =  2 ", "dispersal =  5 ", "dispersal =  8 ", "dispersal =  10 "))~Competition)
  facet_grid(d~alpha1+alpha2)

png(filename="/home/po/Documents/Yupei/patch_feedback_model/plot/JC_weak_1500_2000.png", width = 1600, height = 1200)
Weak
dev.off()

