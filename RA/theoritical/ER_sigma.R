#Enemy release----
#### Load library
library("deSolve")
library("tidyverse")
library("ggplot2")


#### Model specification
patch_feedback_model <- function(times, state, parms) {
  with(as.list(c(state, parms, times)), {
    dMa1_dt = Ma1*(ga1-((Ma1 + alpha1*Mb1)/((phi1 * P1) + 0.000000001))) - d*Ma1*P2 + d*Ma2*P1
    dMa2_dt = Ma2*(ga2-((Ma2 + alpha2*Mb2)/((phi2 * P2) + 0.000000001))) - d*Ma2*P1 + d*Ma1*P2
    dMb1_dt = Mb1*(gb1-((Mb1 + alpha2*Ma1)/((phi1 * P1) + 0.000000001))) - d*Mb1*P2 + d*Mb2*P1
    dMb2_dt = Mb2*(gb2-((Mb2 + alpha1*Ma2)/((phi2 * P2) + 0.000000001))) - d*Mb2*P1 + d*Mb1*P2
    
    dP1_dt = P1*(r1 - c11*P1 - c12*P2 + sigma1*Ma1 + sigma1*Mb1)
    dP2_dt = P2*(r2 - c22*P2 - c21*P1 + sigma2*Ma2 + sigma2*Mb2)
    
    return(list(c(dMa1_dt, dMa2_dt, dMb1_dt, dMb2_dt, dP1_dt, dP2_dt)))
  })
}


#### Prepare functions
SimulateODE <- function(parameters){
  times <- seq(0, 1000000, by = 10)
  state <- c(Ma1 = 0.1, Ma2 = 0.1, Mb1 = 0.1, Mb2 = 0.1, P1 = 0.1, P2 = 0.1)
  # parms <- unlist(parameters)
  parms <- parameters
  pop_size <- ode(func = patch_feedback_model, times = times, y = state, parms = parms)
  return(pop_size)
}

#### Create data using functional programming
Data <-
  crossing(
    d = 1000,
    ga1 = 0.5,
    ga2 = 0.1, 
    gb1 = 0.1,
    gb2 = 0.5, 
    alpha1 = 0, 
    alpha2 = 0,
    phi1 = 1, 
    phi2 = 1, 
    c11 = 1.5, 
    c22 = 1.5, 
    c12 = 1, 
    c21 = 1.5, 
    sigma1 = seq(-0.02, -0.001, 0.001), 
    sigma2 = seq(-1000, 0, 100), 
    r1 = 0.01, r2 = 0.01
  ) %>%
  # mutate(Par = pmap(., list)) %>%
  mutate(Par = pmap(., c)) %>%
  mutate(Scenario = paste("sigma1 = ", sigma1, ", sigma2 = ", sigma2)) %>%
  mutate(Sim = map(Par, SimulateODE)) %>%
  mutate(Result = map2(Scenario, Sim, bind_cols)) %>%
  select(sigma1, sigma2, Result) %>%
  unnest(Result) %>%
  rename("Scenario" = "...1") 
  
save(Data, file = "/Users/yupeitseng/Documents/RA/Theoritical/ER_sigma_data.RData")
load("/Users/yupeitseng/Documents/RA/Theoritical/ER_data.RData")

Data %>% select() %>% pivot_wider(names_from = station, values_from = seen)


#Coexsistence
coexsist <- Data %>% 
  filter(time %in% seq(from = 99900, to = 100000, by = 10)) %>% 
  gather(key = Variable, value = Density, Ma1:P2) %>% 
  group_by(sigma1, sigma2, Scenario, Variable) %>% 
  summarise(Mean_density = mean(Density)) %>% 
  pivot_wider(names_from = Variable, values_from = Mean_density) %>% 
  mutate(sigma1 = as.factor(sigma1), sigma2 = as.factor(sigma2))
coexsist[coexsist < 0.0000001] <- 0

text <- apply(coexsist[, 4:9], 1, function(x){names(which(x == 0))})
vec_text <- c()
for (i in 1:nrow(coexsist)) {
  vec_text <- c(vec_text, paste(unlist(text[[i]]), collapse = " "))
}
vec_text[which(vec_text == "")] <- "coexist"
coexsist <- coexsist %>% as.data.frame() %>% mutate(result = vec_text)

coexsist %>% select(sigma1, sigma2, result) %>% 
  ggplot(aes(sigma1, sigma2, fill = result)) +
  geom_tile(aes(fill = result),colour = "white")

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