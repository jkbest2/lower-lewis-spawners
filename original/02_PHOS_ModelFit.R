library(tidyverse)
library(magrittr)
library(rjags)
library(lubridate)
library(mcmcplots)
library(parallel)
library(doParallel)
library(random)

#Process Data
source("original/00_PHOS_DataPreparation.R")

#Load Model Code and Functions
source("original/01_PHOS_functions_and_JAGSmodel.R")

## Fix the missing 4 states after capture in the trap (and removal from the system)
fix_trapped_state <- function(input) {
  z <- input$z
  K <- ncol(z)
  for (r in seq_len(nrow(z))) {
    merwin <- which(z[r, ] == 4)
    if (length(merwin) > 0) {
      first_merwin <- min(merwin)
      z[r, first_merwin:K] <- 4
    }
  }
  input$z <- z
  input
}
pHOS.Input <- fix_trapped_state(pHOS.Input)

pHOS.inits <- function(){
  list(z      = pHOS.ms.init.z(pHOS.Input$z),
       sex    = pHOS.or.init(pHOS.Input$sex),
       resid  = pHOS.or.init(pHOS.Input$resid),
       origin = pHOS.or.init(pHOS.Input$origin))
  }

params <- c("phi", "gamma", "p", "beta", "delta", "resid_ratio", 
            "sex_ratio", "N_NOR_spawn", "N_HOR_spawn",
            "N_NOR_entered", "N_HOR_entered",
            "pHOS")


m.PHOS <- runjags::run.jags(model    = pHOS_Multistate,
                            data     = pHOS.Input, 
                            inits    = pHOS.inits,
                            monitor  = params,
                            thin     = 10, # 25,000 samples will be drawn but 2500 will be saved
                            n.chains = 4,
                            burnin   = 4000, 
                            adapt    = 1000, 
                            sample   = 2500, 
                            method   = 'parallel'
                            )

fit <- runjags::extend.jags(m.PHOS,
                            adapt = 1000,
                            thin = 10,
                            burnin = 0,
                            sample = 1000,
                            combine = FALSE)
saveRDS(fit, "pHOS_mcmc.rds")

fit <- runjags::extend.jags(fit,
                            adapt = 1000,
                            thin = 10,
                            burnin = 0,
                            sample = 500,
                            combine = TRUE)
saveRDS(fit, "pHOS_mcmc2.rds")


#Run diagnostic plots
## mcmcplot(coda::as.mcmc.list(m.PHOS))

library(posterior)
post <- posterior::as_draws_rvars(m.PHOS$mcmc)
post2 <- as_draws_rvars(fit$mcmc)
write_rds(post2, "data/post2016.rds")
post_rhat <- map(post2, rhat)
post_ess <- map(post2, ess_bulk)
post_median_mcse <- map(post2, mcse_median)
post_median <- map(post2, median)


obs_df <- tibble(
  n_capture = apply(pHOS.Input$y, 1, \(r) sum(r %in% c(2, 3))),
  captured = n_capture > 0,
  recaptured = n_capture > 1,
  trapped = apply(pHOS.Input$y, 1, \(r) any(r == 4))) |>
  summarize(
    n_cap = sum(captured),
    n_recap = sum(recaptured),
    n_trapped = sum(trapped)) |>
  mutate(n = 357 + n_trapped) # Estimate from the original memo
write_rds(obs_df, "data/obs16_summary.rds")
