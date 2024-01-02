library(tidyverse)

n_events <- 11
n_individuals <- 350

p_entry <- 0.1

df <- tibble(
  id = seq_len(n_individuals),
  sex = sample(c("F", "M"), n_individuals, replace = TRUE, prob = c(0.45, 0.55)),
  origin = sample(c("N", "H"), n_individuals, replace = TRUE, prob = c(0.49, 0.51)),
  migrant = rbinom(n_individuals, 1, 0.4),
  resid = rbinom(n_individuals, 1, 0.07),
  state = 0,
  entry_t = ifelse(resid, 1, sample(seq_len(n_events))),
  trap_t = ifelse(migrant, sample(seq_len(n_events)), Inf)
)

p_tangle = 0.1
df2 <- expand_grid(
  df,
  event = seq_len(n_events)
) |>
  mutate(
    avail = event >= entry_t & event < trap_t) |>
  mutate(
    n_avail = sum(avail), .by = event) |>
  mutate(
    tangle_obs = rbinom(n_individuals * n_events, 1, avail * p_tangle),
    trap_obs = event == trap_t)


df2_long <- df2 |>
  select(id, event, tangle_obs) |>
  pivot_wider(
  names_from = event,
  values_from = tangle_obs)

df2 |>
  summarize(tangle_obs = sum(tangle_obs),
            trap_obs = sum(trap_obs),
            .by = event)
