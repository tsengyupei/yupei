#Jazan-Conell----
#### Create data using functional programming
Data <-
  crossing(
    d = c(0, 5, 10, 13, 18, 20, 22, 28, 30, 32, 35, 40),
    ga1 = 1.6,
    ga2 = 1, 
    gb1 = 1,
    gb2 = 2, 
    alpha1 = c(1, 1.1, 1.1, 1.1, 1.1, 1.2, 1.2, 1.2, 1.2, 1.3, 1.3, 1.3, 1.4),
    alpha2 = c(1, 1.1, 1.2, 1.3, 1.4, 1.2, 1.3, 1.4, 1.5, 1.3, 1.4, 1.5, 1.7),
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
  filter(competition == "1,1" | competition == "1.1,1.1" | competition == "1.1,1.2" | competition == "1.1,1.3"| competition == "1.1,1.4"| competition == "1.2,1.2"| competition == "1.2,1.3"| competition == "1.2,1.4"| competition == "1.2,1.5"| competition == "1.3,1.3"| competition == "1.3,1.4"| competition == "1.3,1.5"| competition == "1.4,1.7") %>%
  select(-competition) %>% 
  mutate(Par = pmap(., c)) %>%
  mutate(Scenario = paste("dispersal = ", d, ", competition = ", paste(alpha1, alpha2, sep = ","), ", m_effect = ", sigma1a)) %>%
  mutate(Sim = map(Par, SimulateODE)) %>%
  mutate(Result = map2(Scenario, Sim, bind_cols)) %>%
  select(d, alpha1, alpha2, sigma1a, Result) %>%
  unnest(Result) %>%
  rename("Scenario" = "...1")
save(Data, file = "/Users/yupeitseng/Documents/RA/Theoritical/JC_full_dispersal_data.RData")
load("")
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