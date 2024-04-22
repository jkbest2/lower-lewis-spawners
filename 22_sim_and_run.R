library(tidyverse)
library(posterior)
library(runjags)
library(rjags)

source("01_simulation.R")
source("02_summarize_sim.R")
source("03_data_prep.R")

fit16 <- read_rds("data/pHOS_mcmc_differentRE_04292018.rds")
post16 <- as_draws_rvars(fit16$mcmc)
med16 <- map(post16, median)

M <- 1857 # Total number of candidate fish

shared_pars <- prep_pars(
  k           = 10,                # Number of mark-recapture events; total time
                                   # tracked is k + 2
  gamma       = med16$gamma,       # Probability of potential -> available trant
  beta        = med16$beta,        # pHOS
  sex_ratio   = med16$sex_ratio,   # Ratio of males to females
  resid_ratio = med16$resid_ratio, # Sex-specific proportion of residuals
  delta       = med16$delta,       # Probability of available -> available transition
  phi         = med16$phi,         # Probability of available -> trap transition
  p           = NA                 # Will be replaced for each scenario
)

ps <- c(0.005, 0.01, 0.02, 0.03)
pfitnames <- c("p005_fit", "p01_fit", "p02_fit", "p03_fit")

if (FALSE) { # Re-simulate observations?
  set.seed(20240308)
  if (!dir.exists("fit")) dir.create("fit")

  ## Generate parameter lists that include the varying value of p
  sim_pars <- map(ps, \(p) assign_in(shared_pars, "p", p))

  ## Simulate a collection of individuals and trajectories once, then simulate
  ## different observations conditional on the trajectories. Detection probability
  ## is not used until making observations, which occurs within the `sim_model`
  ## function, so the `shared_pars` can be used until then. `sim_pars` is
  ## generated to make extracting initial values for parameters easier later on.
  indivs <- replicate(M, sim_indiv(shared_pars), simplify = FALSE)
  indiv_df <- purrr::list_transpose(indivs) |>
    tibble::as_tibble()
  zsfull <- purrr::map(indivs, sim_eco, pars = shared_pars)
  sims <- map(sim_pars, \(pars) sim_model(indivs = indivs,
                                          zsfull = zsfull,
                                          pars = pars))

  mod_data <- map2(sims, sim_pars, prep_data)
  write_rds(mod_data, "data/mod_data.rds")
} else {
  mod_data <- read_rds("data/mod_data.rds")
}


params <- c("phi", "gamma", "p", "beta", "delta", "resid_ratio",
            "sex_ratio", "N_NOR_spawn", "N_HOR_spawn",
            "N_NOR_entered", "N_HOR_entered",
            "pHOS")

samples_per_step <- 500
n_steps <- 5

todo_df <- tibble(
  p = ps,
  p_idx = seq_along(ps),
  pfitname = pfitnames) |>
  expand_grid(step_idx = seq_len(n_steps)) |>
  mutate(fn = file.path("fit", paste0(pfitname, step_idx, ".rds"))) |>
  filter(!file.exists(fn)) |>
  arrange(step_idx, p_idx)

for (fit_idx in seq_len(nrow(todo_df))) {
  p_idx <- todo_df$p_idx[fit_idx]
  pfn <- todo_df$fn[fit_idx]
  fit <- run.jags(
    model = "phos-model.jags",
    data = mod_data[[p_idx]]$data,
    inits = mod_data[[p_idx]]$init,
    monitor = params,
    thin = 10,
    n.chains = 4,
    burnin = 4000,
    adapt = 1000,
    sample = samples_per_step,
    method = "parallel")
  write_rds(fit, pfn)
  message("Completed ", fit_idx, " of ", nrow(todo_df))
}

## for (step_idx in 2:n_steps) {
##   pfn <- file.path("fit", paste0(pfitnames, step_idx, ".rds"))
##   for (p_idx in seq_along(ps)) {
##     fit <- run.jags(
##       model = "phos-model.jags",
##       data = mod_data[[p_idx]]$data,
##       inits = mod_data[[p_idx]]$init,
##       monitor = params,
##       thin = 10,
##       n.chains = 4,
##       burnin = 4000,
##       adapt = 1000,
##       sample = samples_per_step,
##       method = "parallel")
##     write_rds(fit, pfn[p_idx])
##   }
## }

## jags_fit <- list()
## for (p_idx in 1:length(ps)) {
##   jags_fit[p_idx] <- run.jags(
##     model = "phos-model.jags",
##     data = mod_data[[p_idx]]$data,
##     inits = mod_data[[p_idx]]$init,
##     monitor = params,
##     thin = 10,
##     n.chains = 4,
##     burnin = 4000,
##     adapt = 1000,
##     sample = samples_per_step,
##     method = "parallel")
##   write_rds(jags_fit[[p_idx]], file.path("fit", pfitnames[p_idx]))
## }
## for (step_idx in 2:n_steps) {
##   for (p_idx in seq_along(ps)) {
##     jags_fit2[[p_idx]] <- extend.jags(
##       runjags.object = jags_fit2[[p_idx]],
##       combine = TRUE,
##       adapt = 1000,
##       thin = 10,
##       burnin = 0,
##       sample = samples_per_step,
##       method = "parallel")
##     write_rds(jags_fit2[[p_idx]], file.path("fit", pfitnames[p_idx]))
##   }
## }

## jags_fit2 <- map(file.path("fit", pfitnames), read_rds)

## fit <- read_rds("fit/p005_fit.rds")
