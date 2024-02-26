#### Load library
library("deSolve")
library("tidyverse")


#### Model specification
patch_feedback_model <- function(times, state, parms) {
  with(as.list(c(state, parms, times)), {
    dMa1_dt = Ma1*(ga1-((Ma1 + alpha1*Mb1)/((phi1 * P1)))) - d*Ma1 + d*Ma2
    dMa2_dt = Ma2*(ga2-((Ma2 + alpha2*Mb2)/((phi2 * P2)))) - d*Ma2 + d*Ma1
    dMb1_dt = Mb1*(gb1-((Mb1 + alpha2*Ma1)/((phi1 * P1)))) - d*Mb1 + d*Mb2
    dMb2_dt = Mb2*(gb2-((Mb2 + alpha1*Ma2)/((phi2 * P2) ))) - d*Mb2 + d*Mb1
    
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
    # d = c(0, 2, 5, 8, 10),
    d = c(1.6, 1.8, 2, 2.2, 2.4, 2.6),
    ga1 = 1.6,
    ga2 = 1, 
    gb1 = 1,
    gb2 = 2, 
    # alpha1 = c(1, 0.2, 0.6, 1.1, 1.5, 0), 
    # alpha2 = c(1, 0.5, 0.9, 1.4, 1.8, 0),
    #competition 
    alpha1 = c(1.1, 1.5, 1.3, 1.1, 1.1),
    alpha2 = c( 1.4, 1.8, 1.6, 1.2, 1.3),
    phi1 = 10, 
    phi2 = 10, 
    c11 = 1.5, 
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
  # filter(competition == "1,1" | competition == "0.2,0.5" | competition == "0.6,0.9" | competition == "1.1,1.4"| competition == "1.5,1.8"| competition == "0,0") %>% 
  filter(competition == "1.3,1.6" | competition == "1.1,1.2" | competition == "1.1,1.3" | competition == "1.1,1.4"| competition == "1.5,1.8") %>%
  select(-competition) %>% 
  mutate(Par = pmap(., c)) %>%
  mutate(Scenario = paste("dispersal = ", d, ", competition = ", paste(alpha1, alpha2, sep = ","), ", m_effect = ", sigma1a)) %>%
  mutate(Sim = map(Par, SimulateODE)) %>%
  mutate(Result = map2(Scenario, Sim, bind_cols)) %>%
  select(Result) %>%
  unnest(Result) %>%
  rename("Scenario" = "...1") %>%
  gather(key = Variable, value = Density, -c(Scenario, time))
str <- strsplit(Data$Scenario, ", ") 
Data <- Data %>% mutate(Dispersal = sapply(str, function(i)i[[1]]), Competition = sapply(str, function(i)i[[2]]), m_effect = sapply(str, function(i)i[[3]]))

#### Plot data
Strong <- Data %>% 
  filter(m_effect == "m_effect =  -6")%>% 
  # filter(time %in% seq(from = 0, to = 10, by = 0.02)) %>%
  # filter(time %in% seq(from = 1500, to = 2000, by = 0.02)) %>%
  ggplot() +
  geom_line(aes(x = time, 
                y = Density, 
                color = Variable))+
  theme_classic()+
  # geom_text(x = 1000, y = 0.04, aes(label = text), data = coexsist %>% filter(m_effect == "m_effect =  -6"))+
  theme(text = element_text(size = 15))+
  # facet_grid(factor(Dispersal, levels = c("dispersal =  0 ", "dispersal =  2 ", "dispersal =  5 ", "dispersal =  8 ", "dispersal =  10 "))~Competition)
  facet_grid(factor(Dispersal, levels = c("dispersal =  1.6 ", "dispersal =  1.8 ", "dispersal =  2 ", "dispersal =  2.2 ", "dispersal =  2.4 ", "dispersal =  2.6 "))~Competition)
  

Weak <- Data %>% 
  filter(m_effect == "m_effect =  -0.32")%>% 
  filter(time %in% seq(from = 0, to = 10, by = 0.02)) %>%
  # filter(time %in% seq(from = 1750, to = 2000, by = 0.02)) %>%
  ggplot() +
  geom_line(aes(x = time, 
                y = Density, 
                color = Variable))+
  theme_classic()+
  # geom_text(x = 1875, y = 0.2, aes(label = text), data = coexsist %>% filter(m_effect == "m_effect =  -0.32"))+
  theme(text = element_text(size = 15))+
  facet_grid(factor(Dispersal, levels = c("dispersal =  0 ", "dispersal =  2 ", "dispersal =  5 ", "dispersal =  8 ", "dispersal =  10 "))~Competition)

#check the result of specific scenerion within specific time period
Data %>% 
  filter(m_effect == "m_effect =  -6" & Competition == "competition =  1.1,1.4 " & Dispersal == "dispersal =  2 ")%>% 
  filter(time %in% seq(from = 1500, to = 1600, by = 0.02)) %>% 
  pivot_wider(names_from = Variable, values_from = Density) %>% View()

#Coexsistence
coexsist <- Data %>% 
  filter(time %in% seq(from = 1500, to = 2000, by = 0.02)) %>% 
  group_by(Scenario, Variable) %>% 
  summarise(Mean_density = mean(Density)) %>% 
  pivot_wider(names_from = Variable, values_from = Mean_density)
coexsist[coexsist < 0.00000001] <- 0
str <- strsplit(coexsist$Scenario, ", ") 
coexsist <- coexsist %>% as.data.frame() %>% mutate(Dispersal = sapply(str, function(i)i[[1]]), Competition = sapply(str, function(i)i[[2]]), m_effect = sapply(str, function(i)i[[3]]))

text <- apply(coexsist[, 2:7], 1, function(x){names(which(x == 0))})
vec_text <- c()
for (i in 1:nrow(coexsist)) {
  vec_text <- c(vec_text, paste(unlist(text[[i]]), collapse = " "))
}
vec_text[which(vec_text == "")] <- "coexist"
coexsist <- coexsist %>% mutate(text = vec_text) 


#Dispersal----
Data_dispersal <-
  crossing(
    # d = c(0, 2, 5, 8, 10),
    d = seq(1, 3, 0.1),
    ga1 = 1.6,
    ga2 = 1, 
    gb1 = 1,
    gb2 = 2, 
    # alpha1 = c(1, 0.2, 0.6, 1.1, 1.5, 0), 
    # alpha2 = c(1, 0.5, 0.9, 1.4, 1.8, 0),
    #competition 
    alpha1 = c(1.1, 1.5, 1.3, 1.1, 1.1),
    alpha2 = c(1.4, 1.8, 1.6, 1.2, 1.3),
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
  filter(sigma1a == sigma2b) %>% 
  mutate(competition = paste(alpha1, alpha2, sep = ",")) %>% 
  # filter(competition == "1,1" | competition == "0.2,0.5" | competition == "0.6,0.9" | competition == "1.1,1.4"| competition == "1.5,1.8"| competition == "0,0") %>% 
  filter(competition == "1.3,1.6" | competition == "1.1,1.2" | competition == "1.1,1.3" | competition == "1.1,1.4"| competition == "1.5,1.8") %>%
  select(-competition) %>% 
  mutate(Par = pmap(., c)) %>%
  mutate(Scenario = paste("dispersal = ", d, ", competition = ", paste(alpha1, alpha2, sep = ","))) %>%
  mutate(Sim = map(Par, SimulateODE)) %>%
  mutate(Result = map2(Scenario, Sim, bind_cols)) %>%
  select(d, alpha1, alpha2, Result) %>%
  unnest(Result) %>%
  rename("Scenario" = "...1") %>%
  gather(key = Variable, value = Density, -c(Scenario, time))
str_dispersal <- strsplit(Data_dispersal$Scenario, ", ") 
Data_dispersal <- Data_dispersal %>% mutate(Dispersal = sapply(str_dispersal, function(i) sub('.* =', '', i[[1]])), Competition = sapply(str_dispersal, function(i)sub('.* =', '', i[[2]])))
Data_dispersal %>% 
  filter(time %in% seq(from = 1500, to = 2000, by = 0.02)) %>%
  ggplot(aes(x = time, 
             y = Density, 
             color = Variable)) +
  geom_line()+
  theme_classic()+
  theme(text = element_text(size = 20))+
  facet_wrap(~Scenario, ncol = 4, nrow = 3)


#Mutual facilitation----
#### Create data using functional programming
Data <-
  crossing(
    d = c(0.08, 0.3, 0.5, 2, 5),
    ga1 = 0.15,
    ga2 = 0.07, 
    gb1 = 0.07,
    gb2 = 0.15, 
    alpha1 = c(1, 0.2, 0.6, 1.1, 1.5), 
    alpha2 = c(1, 0.5, 0.9, 1.4, 1.8),
    phi1 = 10, 
    phi2 = 10, 
    c11 = 1, 
    c22 = 1, 
    c12 = 0.5, 
    c21 = 1, 
    sigma1a = c(0, 1), 
    sigma1b = -0.4, 
    sigma2a = -0.5, 
    sigma2b = c(0, 1),
    r1 = 0.1, r2 = 0.1
  ) %>%
  # mutate(Par = pmap(., list)) %>%
  filter(sigma1a == sigma2b) %>% 
  mutate(competition = paste(alpha1, alpha2, sep = ",")) %>% 
  filter(competition == "1,1" | competition == "0.2,0.5" | competition == "0.6,0.9" | competition == "1.1,1.4"| competition == "1.5,1.8") %>% 
  select(-competition) %>% 
  mutate(Par = pmap(., c)) %>%
  mutate(Scenario = paste("dispersal = ", d, ", competition = ", paste(alpha1, alpha2, sep = ","), ", m_effect = ", sigma1a)) %>%
  mutate(Sim = map(Par, SimulateODE)) %>%
  mutate(Result = map2(Scenario, Sim, bind_cols)) %>%
  select(Result) %>%
  unnest(Result) %>%
  rename("Scenario" = "...1") %>%
  gather(key = Variable, value = Density, -c(Scenario, time))
str <- strsplit(Data$Scenario, ", ") 
Data <- Data %>% mutate(Dispersal = sapply(str, function(i)i[[1]]), Competition = sapply(str, function(i)i[[2]]), m_effect = sapply(str, function(i)i[[3]]))

#### Plot data
Strong <- Data %>% 
  filter(m_effect == "m_effect =  -6")%>% 
  # filter(time %in% seq(from = 0, to = 10, by = 0.02)) %>%
  filter(time %in% seq(from = 1500, to = 2000, by = 0.02)) %>%
  ggplot() +
  geom_line(aes(x = time, 
                y = Density, 
                color = Variable))+
  theme_classic()+
  geom_text(x = 1800, y = 0.075, aes(label = text), data = coexsist %>% filter(m_effect == "m_effect =  -6"))+
  theme(text = element_text(size = 15))+
  facet_grid(factor(Dispersal, levels = c("dispersal =  0 ", "dispersal =  2 ", "dispersal =  5 ", "dispersal =  8 ", "dispersal =  10 "))~Competition)

Weak <- Data %>% 
  filter(m_effect == "m_effect =  -0.32")%>% 
  filter(time %in% seq(from = 0, to = 10, by = 0.02)) %>%
  # filter(time %in% seq(from = 1750, to = 2000, by = 0.02)) %>%
  ggplot() +
  geom_line(aes(x = time, 
                y = Density, 
                color = Variable))+
  theme_classic()+
  # geom_text(x = 1875, y = 0.2, aes(label = text), data = coexsist %>% filter(m_effect == "m_effect =  -0.32"))+
  theme(text = element_text(size = 15))+
  facet_grid(factor(Dispersal, levels = c("dispersal =  0 ", "dispersal =  2 ", "dispersal =  5 ", "dispersal =  8 ", "dispersal =  10 "))~Competition)

#check the result of specific scenerion within specific time period
Data_neutral %>% 
  filter(m_effect == "m_effect =  -6" & Competition == "competition =  1,1 " & Dispersal == "dispersal =  2 ")%>% 
  filter(time %in% seq(from = 1500, to = 2000, by = 0.02)) %>% 
  pivot_wider(names_from = Variable, values_from = Density)

#Coexsistence
coexsist <- Data %>% 
  filter(time %in% seq(from = 1500, to = 2000, by = 0.02)) %>% 
  group_by(Scenario, Variable) %>% 
  summarise(Mean_density = mean(Density)) %>% 
  pivot_wider(names_from = Variable, values_from = Mean_density)
coexsist[coexsist < 0.00000001] <- 0
str <- strsplit(coexsist$Scenario, ", ") 
coexsist <- coexsist %>% as.data.frame() %>% mutate(Dispersal = sapply(str, function(i)i[[1]]), Competition = sapply(str, function(i)i[[2]]), m_effect = sapply(str, function(i)i[[3]]))

text <- apply(coexsist[, 2:7], 1, function(x){names(which(x == 0))})
vec_text <- c()
for (i in 1:nrow(coexsist)) {
  vec_text <- c(vec_text, paste(unlist(text[[i]]), collapse = " "))
}
vec_text[which(vec_text == "")] <- "coexist"
coexsist <- coexsist %>% mutate(text = vec_text) 