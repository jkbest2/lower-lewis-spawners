library(tidyverse)
library(posterior)

source("01_simulation.R")
source("02_summarize_sim.R")
source("03_data_prep.R")

## Import the results of the fit to the 2016 data. Use the posterior medians as
## parameter values for the simulated data sets. Should give dynamics reasonably
## similar to the 2016 data set for comparison.
post16 <- read_rds("data/post2016.rds")
med16 <- map(post16, median)

## Set up the simulation parameters, with 10 capture events (so tracking 12
## total time steps).
sim_pars <- prep_pars(
  k           = 10,                # Number of mark-recapture events; total time
                                   # tracked is k + 2
  gamma       = med16$gamma,       # Probability of potential -> available trant
  beta        = med16$beta,        # pHOS
  sex_ratio   = med16$sex_ratio,   # Ratio of males to females
  resid_ratio = med16$resid_ratio, # Sex-specific proportion of residuals
  delta       = med16$delta,       # Probability of available -> available transition
  phi         = med16$phi,         # Probability of available -> trap transition
  p           = NA                 # Will be replaced for each simulation
)

##
sim_df <- expand_grid(
  rep = 1:20,
  p = c(0.005, 0.01, 0.015, 0.02, 0.025, 0.03),
  M = c(500, 1000, 1500, 2000, 2500)) |>
  mutate(pars = map(p, \(p) assign_in(sim_pars, "p", p)),
         sim = map2(M, pars, sim_model),
         handled = map(sim, sim_handled))

handled_df <- sim_df |>
  select(rep, p, M, handled) |>
  unnest(handled)

handle_summ <- handled_df |>
  filter(final_state != 1) |>
  summarize(n = sum(n),
            n_cap = sum(n_cap),
            n_recap = sum(n_recap),
            n_trapped = sum(n_trapped),
            .by = c(rep, p, M)) |>
  mutate(pct = paste0(100 * p, "%"))

obs16_df <- read_rds("data/obs16_summary.rds")

ggplot(handle_summ, aes(x = n, y = n_cap, color = pct)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE, method = lm) +
  geom_point(data = obs16_df, color = "black", size = 2) +
  annotate(geom = "label",
           x = 1030, y = 67,
           label = "2016 Observed") +
  labs(x = "Total run (lower- and upper-Lewis)",
       y = "Number of individuals captured",
       color = "Probability\n of capture")

## t_cap <- apply(pHOS.Input$y, 2, \(.) sum(. %in% c(2, 3), na.rm = TRUE))

## pwt <- (post16$p * t_cap[2:11] / sum(t_cap[2:11])) |>
##   rvar_sum() |> median()

ggplot(handle_summ, aes(x = n, y = n_recap, color = pct)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE, method = lm) +
  geom_point(data = obs16_df, color = "black", size = 2) +
  annotate(geom = "label",
           x = 1060, y = 3.2,
           label = "2016 Observed") +
  labs(x = "Total run (lower- and upper-Lewis)",
       y = "Number of individuals recaptured",
       color = "Probability\n of capture")
