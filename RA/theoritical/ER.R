#Enemy release----
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
  times <- seq(0, 10000, by = 0.1)
  state <- c(Ma1 = 0.1, Ma2 = 0.1, Mb1 = 0.1, Mb2 = 0.1, P1 = 0.1, P2 = 0.1)
  # parms <- unlist(parameters)
  parms <- parameters
  pop_size <- ode(func = patch_feedback_model, times = times, y = state, parms = parms)
  return(pop_size)
}

#### Create data using functional programming
Data <-
  crossing(
    d = c(0, 2, 5, 8, 10),
    ga1 = 0.5,
    ga2 = 0.1, 
    gb1 = 0.1,
    gb2 = 0.5, 
    alpha1 = c(0, 1, 0.2, 0.6, 1.1, 1.5), 
    alpha2 = c(0, 1, 0.5, 0.9, 1.4, 1.8),
    phi1 = 100000, 
    phi2 = 100000, 
    c11 = 1.5, 
    c22 = 1.5, 
    c12 = 1, 
    c21 = 1.5, 
    sigma1a = -0.5, 
    sigma1b = -0.5, 
    sigma2a = c(-2, 0), 
    sigma2b = c(-2, 0),
    r1 = 0.001, r2 = 0.001
  ) %>%
  # mutate(Par = pmap(., list)) %>%
  filter(sigma2a == sigma2b) %>% 
  mutate(competition = paste(alpha1, alpha2, sep = ",")) %>% 
  filter(competition == "0,0" | competition == "1,1" | competition == "0.2,0.5" | competition == "0.6,0.9" | competition == "1.1,1.4"| competition == "1.5,1.8") %>% 
  select(-competition) %>% 
  mutate(Par = pmap(., c)) %>%
  mutate(Scenario = paste("dispersal = ", d, ", competition = ", paste(alpha1, alpha2, sep = ","), ", m_effect = ", sigma2a)) %>%
  mutate(Sim = map(Par, SimulateODE)) %>%
  mutate(Result = map2(Scenario, Sim, bind_cols)) %>%
  select(d, sigma2a, alpha1, alpha2, Result) %>%
  unnest(Result) %>%
  rename("Scenario" = "...1") 
save(Data, file = "/Users/yupeitseng/Documents/RA/Theoritical/ER_data.RData")
load("/Users/yupeitseng/Documents/RA/Theoritical/ER_data.RData")

View(Data %>% filter(time %in% seq(9990, 10000, 0.1)) %>% filter(d == 2, sigma1b == 2, alpha1 == 1, alpha2 == 1)) 


#Coexsistence
coexsist <- Data %>% 
  filter(time %in% seq(from = 9500, to = 10000, by = 0.1)) %>% 
  gather(key = Variable, value = Density, Ma1:P2) %>% 
  group_by(d, alpha1, alpha2, sigma2a, Scenario, Variable) %>% 
  summarise(Mean_density = mean(Density)) %>% 
  pivot_wider(names_from = Variable, values_from = Mean_density) %>% 
  mutate(d = as.factor(d), alpha1 = as.factor(alpha1), alpha2 = as.factor(alpha2), sigma2a= as.factor(sigma2a))
coexsist[coexsist < 0.00000001] <- 0

text <- apply(coexsist[, 6:11], 1, function(x){names(which(x == 0))})
vec_text <- c()
for (i in 1:nrow(coexsist)) {
  vec_text <- c(vec_text, paste(unlist(text[[i]]), collapse = " "))
}
vec_text[which(vec_text == "")] <- "coexist"
coexsist <- coexsist %>% as.data.frame() %>% mutate(result = vec_text)

#plotting
png(filename="/Users/yupeitseng/Documents/RA/Theoritical/plot/ER_strong.png", width = 1600, height = 1200)
Data %>% 
  filter(sigma2a == 0) %>% 
  gather(key = Variable, value = Density, Ma1:P2) %>% 
  ggplot() +
  geom_line(aes(x = time, 
                y = Density, 
                color = Variable))+
  theme_classic()+
  geom_text(x = 5000, y = 3, aes(label = result), data = coexsist %>% filter(sigma2a == 0))+
  theme(text = element_text(size = 15))+
  facet_grid(d~alpha1 + alpha2)
dev.off()

png(filename="/Users/yupeitseng/Documents/RA/Theoritical/plot/ER_strong_9500_10000.png", width = 1600, height = 1200)
Data %>% 
  filter(sigma2a == 0) %>% 
  filter(time %in% seq(from = 9500, to = 10000, by = 0.1)) %>%
  gather(key = Variable, value = Density, Ma1:P2) %>% 
  ggplot() +
  geom_line(aes(x = time, 
                y = Density, 
                color = Variable))+
  theme_classic()+
  geom_text(x = 9750, y = 3, aes(label = result), data = coexsist %>% filter(sigma2a == 0))+
  theme(text = element_text(size = 15))+
  facet_grid(d~alpha1 + alpha2)
dev.off()

png(filename="/Users/yupeitseng/Documents/RA/Theoritical/plot/ER_weak.png", width = 1600, height = 1200)
Data %>% 
  filter(sigma2a == -2) %>% 
  gather(key = Variable, value = Density, Ma1:P2) %>% 
  ggplot() +
  geom_line(aes(x = time, 
                y = Density, 
                color = Variable))+
  theme_classic()+
  geom_text(x = 5000, y = 3, aes(label = result), data = coexsist %>% filter(sigma2a == -2))+
  theme(text = element_text(size = 15))+
  facet_grid(d~alpha1 + alpha2)
dev.off()

png(filename="/Users/yupeitseng/Documents/RA/Theoritical/plot/ER_weak_9500_10000.png", width = 1600, height = 1200)
Data %>% 
  filter(sigma2a == -2) %>% 
  filter(time %in% seq(from = 9500, to = 10000, by = 0.1)) %>%
  gather(key = Variable, value = Density, Ma1:P2) %>% 
  ggplot() +
  geom_line(aes(x = time, 
                y = Density, 
                color = Variable))+
  theme_classic()+
  geom_text(x = 9750, y = 0.0015, aes(label = result), data = coexsist %>% filter(sigma2a == -2))+
  theme(text = element_text(size = 15))+
  facet_grid(d~alpha1 + alpha2)
dev.off()