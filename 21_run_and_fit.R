library(tidyverse)
library(rjags)
library(runjags) # Loading this first prevents warnings about JAGS version
library(posterior) # Useful for postprocessing posterior samples

source("01_simulation.R")
source("02_summarize_sim.R")
source("03_data_prep.R")

## Number of capture events
k <- 10
## Add the initial state and the final, trap-only event
K <- k + 2
## Potential population size
M <- 1500

##~~~~PRIORS AND CONSTRAINTS:
## Detection efficiency
p <- 0.02

## Entry rate; will sample
runsize <- 1857
## Proportion entering before each event (none in the final event; doesn't count
## residualized fish who all "enter" before the first event). Use k + 1 sequence
## because they are diff'd. This leaves k probabilities of entry and then
## gamma <- seq(-2, 2, length.out = k + 1) |>
##   pnorm() |>
##   diff()
## Entry rate from the data
gamma <- c(0.0830998826355685, 0.216758138094457, 0.0568802271749238,
           0.227777798875943, 0.0464788101688508, 0.0437682338469231,
           0.0432419946705794, 0.0461711834433321, 0.0430909529926357,
           0.0242107786896443)
## No fish enter the system in the zero'th or last periods
gamma <- c(0, gamma, 0)

## Proportion of hatchery fish in population
beta <- 0.51
## Sex ratio
sex_ratio <- 0.51
## Residualization ratios: males at 7.7% and females at 0.7%
resid_ratio <- c(0.077, 0.007)

## delta: probability of remaining in returning state
delta <- 0.8

## phi: probability of moving to trap
phi <- 0.65

pars <- list(gamma = gamma,
             beta = beta,
             sex_ratio = sex_ratio,
             resid_ratio = resid_ratio,
             delta = delta,
             phi = phi,
             p = p)

res <- sim_model(500, pars)
dat <- prep_data(res)
sim_check(dat)

params <- c("phi", "gamma", "p", "beta", "delta", "resid_ratio",
            "sex_ratio", "N_NOR_spawn", "N_HOR_spawn",
            "N_NOR_entered", "N_HOR_entered",
            "pHOS")

fit <- runjags::run.jags(model    = "phos-model.jags",
                         data     = dat$data,
                         inits    = dat$init,
                         monitor  = params,
                         thin     = 1,
                         n.chains = 1,
                         burnin   = 4000,
                         adapt    = 1000,
                         sample   = 2500,
                         ## method   = 'parallel'
                         )

post <- as_draws_rvars(fit$mcmc)
