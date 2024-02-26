#### Load library
library("deSolve")
library("tidyverse")
LV <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dN1_dt = N1 * (1 - alpha_11 * N1 - alpha_12 * N2)
    dN2_dt = N2 * (1 - alpha_22 * N2 - alpha_21 * N1)
    return(list(c(dN1_dt, dN2_dt)))
  })
}


#### Prepare functions
SimulateODE <- function(parameters, initial){
  times <- seq(0, 50, by = 0.1)
  state <- initial
  # parms <- unlist(parameters)
  parms <- parameters
  pop_size <- ode(func = LV, times = times, y = state, parms = parms)
  return(pop_size)
}


#### Create data using functional programming
Data <-
  crossing(
    alpha_11 = 0.01,
    alpha_22 = 0.01,
    alpha_12 = 0.003,
    alpha_21 = c(0.005, 0.015, 0.025, 0.035)
  ) %>%
  # mutate(Par = pmap(., list)) %>%
  mutate(Par = pmap(., c)) %>%
  mutate(Initial = map2(5, 10, function(a, b) c(N1 = a, N2 = b))) %>%
  mutate(Scenario = as.character(alpha_21)) %>%
  mutate(Sim = map2(Par, Initial, SimulateODE)) %>%
  mutate(Result = map2(Scenario, Sim, bind_cols)) %>%
  select(Result) %>%
  unnest(Result) %>%
  rename(Scenario = "...1") %>%
  gather(key = Variable, value = Density, -c(Scenario, time))


# #### Create data using classic approach with multiple objects
# times <- seq(0, 50, by = 0.1)
# state <- c(N1 = 5, N2 = 10)
# parms.1 <- c(alpha_11 = 0.01, alpha_21 = 0.005, alpha_22 = 0.01, alpha_12 = 0.003)
# parms.2 <- c(alpha_11 = 0.01, alpha_21 = 0.015, alpha_22 = 0.01, alpha_12 = 0.003)
# pop_size.1 <- as.data.frame(ode(func = LV, times = times, y = state, parms = parms.1))
# pop_size.2 <- as.data.frame(ode(func = LV, times = times, y = state, parms = parms.2))
# pop_size.1$Scenario <- "1"
# pop_size.2$Scenario <- "2"
# Data <- 
#   as.data.frame(rbind(pop_size.1, pop_size.2)) %>%
#   gather(key = Variable, value = Density, -c(Scenario, time))


#### Plot data: LV-competition under two scenarios
Data %>%
  ggplot(aes(x = time, 
             y = Density, 
             color = Scenario, 
             linetype = Variable)) +
  geom_line(linewidth = 1) + 
  theme_classic()
  