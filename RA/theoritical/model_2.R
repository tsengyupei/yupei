install.packages('deSolve')
library(deSolve)
library(tidyverse)

### Model specification
hierarchy_patch_feedback_model <- function(times, state, parms) {
  with(as.list(c(state, parms, times)), {
    
    
    dP1_dt = c1*P1*(1 - P1 - P2)*(Ra1/(Ra1 + Ra2)) -m1*P1 - c2*P1*P2
    dP2_dt = c2*P2*(1 - P1 - P2)*(Rb2/(Rb1 + Rb2)) -m2*P2 - c1*P1*P2
    
    dRa1_dt = ra1*Ra1*(1-((Ra1 + aab*Rb1)/(A * P1))) - da*Ra1 + da*Ra2
    dRa2_dt = ra2*Ra2*(1-((Ra2 + aab*Rb2)/(A * P2))) - da*Ra2 + da*Ra1
    dRb1_dt = rb1*Rb1*(1-((Rb1 + aba*Ra1)/(A * P1))) - db*Rb1 + db*Rb2
    dRb2_dt = rb2*Rb2*(1-((Rb2 + aba*Ra2)/(A * P2))) - db*Rb2 + db*Rb1
    return(list(c(dRa1_dt, dRa2_dt, dRb1_dt, dRb1_dt, dP1_dt, dP2_dt)))
    })
}

### Model parameters
times <- seq(0, 20, by = 0.1)
state <- c(P1 = 0.2, P2 = 0.3, Ra1 = 100, Ra2 = 100, Rb1 = 100, Rb2 = 100)
# parms <- c(ra1 = ra1, ra2 = ra2, rb1 = rb1, b2 = rb2, aab = aab, aba = aba, da = da, db = db, A = A, c1 = c1, c2 = c2, m1 = m1, m2 = m2)
parms <- c(ra1 = 0.5, ra2 = 0.1, rb1 = 0.1, rb2 = 0.5, aab = 0.005, aba = 0.005, da = 0.1, db = 0.1, A = 10000000, c1 = 0.5, c2 = 0.5, m1 = 0.1, m2 = 0.1)  

### ODE solver
patch_prop <- ode(func = hierarchy_patch_feedback_model, times = times, y = state, parms = parms, method = "rk4")


data(cars)
plot(dist~speed, data = cars)
