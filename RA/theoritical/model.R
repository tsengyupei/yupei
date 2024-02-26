library(deSolve)
library(tidyr)

ke <- function(times, state, parms) {
  with(as.list(c(state, parms, times)), {
    dMa1_dt = Ma1*(ga1-((Ma1)/((phi1 * P1) + 0.000000001))) 
    
    dMb2_dt = Mb2*(gb2-((Mb2)/((phi2 * P2) + 0.000000001))) 
    
    dP1_dt = P1*(r1 - c11*P1 - c12*P2 + sigma1a*Ma1 + sigma1b*Mb2)
    dP2_dt = P2*(r2 - c22*P2 - c21*P1 + sigma2a*Ma1 + sigma2b*Mb2)
    
    return(list(c(dMa1_dt, dMb2_dt, dP1_dt, dP2_dt)))
  })
}
parms <- c(ga1 = 0.5, gb2 = 0.5, phi1 = 0.1, phi2 = 0.1, c11 = 1.5, c22 = 1.5, c12 = 1, c21 = 1.5, sigma1a = -0.5, sigma1b = -0.5, sigma2a = 0, sigma2b = 0, r1 = 0.1, r2 = 0.1)
### Model specification
patch_feedback_model <- function(times, state, parms) {
  with(as.list(c(state, parms, times)), {
    dMa1_dt = Ma1*(ga1-((Ma1 + alpha1*Mb1)/((phi1 * P1) + 0.000000001))) - d*Ma1*P2 + d*Ma2*P1
    dMa2_dt = Ma2*(ga2-((Ma2 + alpha2*Mb2)/((phi2 * P2) + 0.000000001))) - d*Ma2*P1 + d*Ma1*P2
    dMb1_dt = Mb1*(gb1-((Mb1 + beta1*Ma1)/((phi1 * P1) + 0.000000001))) - d*Mb1*P2 + d*Mb2*P1
    dMb2_dt = Mb2*(gb2-((Mb2 + beta2*Ma2)/((phi2 * P2) + 0.000000001))) - d*Mb2*P1 + d*Mb1*P2
    
    dP1_dt = P1*(r1 - c11*P1 - c12*P2 + sigma1a*Ma1 + sigma1b*Mb1)
    dP2_dt = P2*(r2 - c22*P2 - c21*P1 + sigma2a*Ma2 + sigma2b*Mb2)
    
    return(list(c(dMa1_dt, dMa2_dt, dMb1_dt, dMb2_dt, dP1_dt, dP2_dt)))
  })
}

patch_feedback_model <- function(times, state, parms) {
  with(as.list(c(state, parms, times)), {
    dMa1_dt = ga1*Ma1*(1-((Ma1 + alpha1*Mb1)/((phi1 * P1) + 0.000000001))) - d*Ma1*P2 + d*Ma2*P1
    dMa2_dt = ga2*Ma2*(1-((Ma2 + alpha2*Mb2)/((phi2 * P2) + 0.000000001))) - d*Ma2*P1 + d*Ma1*P2
    dMb1_dt = gb1*Mb1*(1-((Mb1 + beta1*Ma1)/((phi1 * P1) + 0.000000001))) - d*Mb1*P2 + d*Mb2*P1
    dMb2_dt = gb2*Mb2*(1-((Mb2 + beta2*Ma2)/((phi2 * P2) + 0.000000001))) - d*Mb2*P1 + d*Mb1*P2
    
    dP1_dt = P1*(r1 - c11*P1 - c12*P2 + sigma1a*Ma1 + sigma1b*Mb1)
    dP2_dt = P2*(r2 - c22*P2 - c21*P1 + sigma2a*Ma2 + sigma2b*Mb2)
    
    return(list(c(dMa1_dt, dMa2_dt, dMb1_dt, dMb2_dt, dP1_dt, dP2_dt)))
  })
}

### Model parameters
#mutual need 50000 or 100000 steps, jc need 2000 steps

times <- seq(0, 1000000, by = 10)
times <- seq(0, 5000, by = 0.02)
state <- c(Ma1 = 0.1, Ma2 = 0.1, Mb1 = 0.1, Mb2 = 0.1, P1 = 0.1, P2 = 0.1)
state <- c(Ma1 = 0.1, Mb2 = 0.1, P1 = 0.1, P2 = 0.1)

#jazen-connell
parms <- c(ga1 = 1.6, ga2 = 1, gb1 = 1, gb2 = 2, alpha1 = 1.2, alpha2 = 1.5, beta1 = 1.5, beta2 = 1.2, d = 13, phi1 = 0.1, phi2 = 0.1, c11 = 1, c22 = 0.8, c12 = 1, c21 = 0.8, sigma1a = -6, sigma1b = -0.4, sigma2a = -0.5, sigma2b = -6, r1 = 0.1, r2 = 0.1)  #jazen-connell/high home soil negative microbial effect/no dispersal/inter>intra

#mutual facilitation
parms <- c(ga1 = 0.2, ga2 = 0.01, gb1 = 0.01, gb2 = 0.05, alpha1 = 0, alpha2 = 0, beta1 = 0, beta2 = 0, d = 8, phi1 = 1, phi2 = 1, c11 = 1, c22 = 1, c12 = 0.9, c21 = 1, sigma1a = 0.1, sigma1b = 3, sigma2a = 3, sigma2b = 0.1, r1 = 0.1, r2 = 0.1) 
#sigma= 0/3

#Enemy release
parms <- c(ga1 = 5, ga2 = 0.1, gb1 = 0.1, gb2 = 5, alpha1 = 0, alpha2 = 0, beta1 = 0, beta2 = 0, d = 1000, phi1 = 1, phi2 = 1, c11 = 1.5, c22 = 1.5, c12 = 1.4, c21 = 1.5, sigma1a = -0.001, sigma1b = -0.001, sigma2a = -1000, sigma2b = -1000, r1 = 0.01, r2 = 0.01)

#1000000

#Conditioning
parms <- c(ga1 = 0.3, ga2 = 0.01, gb1 = 0.01, gb2 = 2.5, alpha1 = 0, alpha2 = 0, beta1 = 0, beta2 = 0, d = 10, phi1 = 0.1, phi2 = 0.1, c11 = 1, c22 = 1, c12 = 1, c21 = 0.95, sigma1a = 1, sigma1b = 0.9, sigma2a = 0.1, sigma2b = 0.05, r1 = 0.1, r2 = 0.1)
0.025/ 0.01

prop <- ode(func = ke, times = times, y = state, parms = parms)
tail(prop)

prop <- ode(func = patch_feedback_model, times = times, y = state, parms = parms)
tail(prop)

head(prop)
prop[which(prop$time %in% seq(660, 800, by = 0.02)),]

# prop[49900:50000,]
prop[1600:2000,]
prop[9900:10000,]
prop[49900:50000,]
prop[249900:250000,]
prop[299900:300000,]
prop[19990:20000,]
### ODE solver
prop <- as.data.frame(prop)
# prop <- prop[75000:100000,]
prop <- prop[199900:200000,]

plot(prop$P1 ~ prop$time, type = 'l', ylim = c(min(prop[, 6:7]), max(prop[, 6:7])))
lines(prop$P2 ~ prop$time, col = 'grey')
plot(prop$Mb1 ~ prop$time, col = 'red', type = 'l', ylim = c(min(prop[, 2:5]), max(prop[, 2:5])))
lines(prop$Mb2 ~ prop$time, col = 'pink')
lines(prop$Ma1 ~ prop$time, col = 'darkblue')
lines(prop$Ma2 ~ prop$time, col = 'skyblue')
prop %>%
  as.data.frame() %>%
  pivot_longer(cols = -time, names_to = "species", values_to = "N") %>%
  ggplot(aes(x = time, y = N, color = species)) + 
  geom_line(size = 1.5) +
  theme_classic(base_size = 12) +
  labs(x = "Time", y = "Population size") +
  scale_x_continuous(limits = c(0, 100.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max(pop_size)*1.2), expand = c(0, 0)) +
  scale_color_brewer(name = NULL, palette = "Set1")
}
