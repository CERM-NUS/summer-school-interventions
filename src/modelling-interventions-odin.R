## digraph SIR {
##     rankdir=LR; // horizontal layout
##     node [shape=rectangle, fontname=Arial];
##     S [label="S\n(Susceptible)"];
##     I [label="I\n(Infectious)"];
##     R [label="R\n(Removed)"];
## 
##     S -> I [label="βIS"];
##     I -> R [label="γI"];
## 
## }

## ------------------------------------------------------------------------------------------------------------------------------
#| lst-label: lst-sir-simulation
#| lst-cap: Implementation of the SIR model in odin
#| label: code-simulation
#| output: false


library(odin)
library(tidyverse)

# define the model
sir_ode <- odin::odin({
  # Derivatives
  deriv(S) <- -beta * S * I
  deriv(I) <- beta * S * I  - gamma * I
  deriv(R) <- gamma * I
  
  # Total population
  N <- S + I + R
  
  # Initial conditions
  initial(S) <- N_ini - I_ini
  initial(I) <- I_ini
  initial(R) <- 0
  
  # user defined values
  # user defined values can be given another value
  # while initializing
  N_ini <- user(5700000)
  I_ini <- user(1)
  
  # user defined parameters
  beta <- user(8.77e-8)
  gamma <- user(0.2)
})


## ------------------------------------------------------------------------------------------------------------------------------
#| lst-label: lst-sir-run
#| lst-cap: Running SIR model
#| label: fig-sir-output
#| fig-cap: "Simulation of the SIR model for 100 days"
#| warning: false


# initialize model
sir_mod <- sir_ode$new()
# how long to run
times <- seq(0,100)
# run the model
sir_out <- sir_mod$run(times)
# extracting data as tibble
df_sir_out <- as.data.frame(sir_out)
# plotting data
p_sir <- df_sir_out %>%
  ggplot(aes(x = t, color = name)) +
  geom_line(aes(y=S, color="S"), alpha = 0.95) +
  geom_line(aes(y=I, color="I"), alpha = 0.95) +
  geom_line(aes(y=R, color="R"), alpha = 0.95) +
  scale_colour_discrete(limits = c("S", "I", "R")) +
  xlab("Time") +
  ylab("Total number of Individuals\n") +
  labs(color="Compartment") +
  theme_classic(base_size = 20) +
  theme(legend.justification = 0, legend.position = c(0.65, 0.55)) 
p_sir


## digraph SEIR {
##     rankdir=LR; // horizontal layout
##     node [shape=rectangle, fontname=Arial];
##     S [label="S\n(Susceptible)"];
##     E [label="E\n(Latent\\Exposed)"];
##     I [label="I\n(Infectious)"];
##     R [label="R\n(Removed)"];
## 
##     S -> E [label="βIS"];
##     E -> I [label="λE"];
##     I -> R [label="γI"];
##     edge [dir=front, label="μ(S+E+I+R)"];
##     muS [style=invis, shape=none, height=0, width=0];
##     muS ->  S;
## 
##     edge [label="μR"];
##     { rank=same; R; deltaR [style=invis, shape=none, height=0, width=0]}
##     R -> deltaR ;
## 
##     edge [label="μI"];
##     { rank=same; I; deltaI [style=invis, shape=none, height=0, width=0]}
##     I -> deltaI;
## 
##      edge [label="μE"];
##     { rank=same; E; deltaE [style=invis, shape=none, height=0, width=0]}
##     E -> deltaE;
## 
##     edge [label="μS"];
##     { rank=same; S; deltaS [style=invis;shape=none, height=0, width=0]}
##     S -> deltaS ;
## 
## }
## 

## ------------------------------------------------------------------------------------------------------------------------------
#| lst-label: lst-seir-simulation
#| lst-cap: Implementation of the SEIR model in odin
#| label: seir-code-simulation
#| output: false

# define the model
seir_ode <- odin::odin({
  ## Derivatives
  deriv(S) <- -beta * S * I + mu*N  -mu*S
  deriv(E) <- beta * S * I - lambda * E - mu *E
  deriv(I) <- lambda * E  - gamma * I - mu * I
  deriv(R) <- gamma * I - mu * R
  
  # Total population
  N <- S + E + I + R
  
  # Initial conditions
  initial(S) <- N_ini - I_ini
  initial(I) <- I_ini
  initial(R) <- 0
  initial(E) <- 0
  
  # user defined values
  # user defined values can be given another value
  # while initializing
  N_ini <- user(5700000)
  I_ini <- user(1)
  
  # user defined parameters
  beta <- user(8.77e-8)
  gamma <- user(0.2)
  lambda <- user(0.5)
  mu <- user(3.3e-5)
})


## ------------------------------------------------------------------------------------------------------------------------------
#| lst-label: lst-seir-run
#| lst-cap: Run SEIR model
#| label: fig-seir-output
#| fig-cap: "Simulation of the SEIR model for 150 days"
#| warning: false


# initialize model
seir_mod <- seir_ode$new()
# how long to run
times <- seq(0,150)
# run the model
seir_out <- seir_mod$run(times)
# extracting data as tibble
df_seir_out <- as.data.frame(seir_out)
# plotting data
p_seir <- df_seir_out %>%
  ggplot(aes(x = t, color = name)) +
  geom_line(aes(y=S, color="S"), alpha = 0.95) +
  geom_line(aes(y=I, color="I"), alpha = 0.95) +
  geom_line(aes(y=R, color="R"), alpha = 0.95) +
  scale_colour_discrete(limits = c("S", "I", "R")) +
  xlab("Time") +
  ylab("Total number of Individuals\n") +
  labs(color="Compartment") +
  theme_classic(base_size = 20) +
  theme(legend.justification = 0, legend.position = c(0.05, 0.55)) 
p_seir


## ------------------------------------------------------------------------------------------------------------------------------
#| lst-label: lst-tv-simulation
#| lst-cap: Implementation of the SEIR model with a shock intervention at $t_1=50$
#| label: tv-code-simulation
#| output: false

# define the model
tv_ode <- odin::odin({
  ## Derivatives
  deriv(S) <- -beta * S * I + mu*N  -mu*S
  deriv(E) <- beta * S * I - lambda * E - mu *E
  deriv(I) <- lambda * E  - gamma * I - mu * I
  deriv(R) <- gamma * I - mu * R
  
  # Total population
  N <- S + E + I + R
  ## Time-varying reproduction number Rt
  beta <- if (t < change_time) beta_0 else beta_0 * reduction
  
  
  # Initial conditions
  initial(S) <- N_ini - I_ini
  initial(I) <- I_ini
  initial(R) <- 0
  initial(E) <- 0
  
  # user defined values
  # user defined values can be given another value
  # while initializing
  N_ini <- user(5700000)
  I_ini <- user(1)
  
  # user defined parameters
  beta_0 <- user(8.77e-8)
  gamma <- user(0.2)
  lambda <- user(0.5)
  mu <- user(3.3e-5)
  reduction <- user(0.8) # Reduction in R0 due to interventions
  change_time <- user(80) # Time of intervention
})


## ------------------------------------------------------------------------------------------------------------------------------
#| lst-label: lst-tv-run
#| lst-cap: Run SEIR model
#| label: fig-tv-output
#| fig-cap: "Simulation of the SEIR model with a shock intervention at $t_1=80$ aand reduction by 20%"
#| warning: false


# initialize model
tv_mod <- tv_ode$new()
# how long to run
times <- seq(0,150)
# run the model
tv_out <- tv_mod$run(times)
# extracting data as tibble
df_tv_out <- as.data.frame(tv_out)
# plotting data
p_tv <- df_tv_out %>%
  ggplot(aes(x = t, color = name)) +
  geom_line(aes(y=S, color="S"), alpha = 0.95) +
  geom_line(aes(y=I, color="I"), alpha = 0.95) +
  geom_line(aes(y=R, color="R"), alpha = 0.95) +
  geom_vline(xintercept = 80, linetype = "dashed", color = "grey", size = 1.5) +
  scale_colour_discrete(limits = c("S", "I", "R")) +
  xlab("Time") +
  ylab("Total number of Individuals\n") +
  labs(color="Compartment") +
  theme_classic(base_size = 20) +
  theme(legend.justification = 0, legend.position = c(0.05, 0.55)) 
p_tv

