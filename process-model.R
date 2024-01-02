library(tidyverse)

arrival <- 100   # Number of fish arriving each day

par <- list(
  frac_lower = 0.5,     # Fraction lower vs upper Lewis destination
  lower_tts = 7,    # Expected time to spawn
  lower_res = 7,        # Expected residence time on redd
  lower_redd = 1 / 1.6, # Expected redds per individual
  ## Upper fish parameters
  upper_ttt = 7,        # Expected days to the Merwin Trap
  upper_spf = 0.1       # Fraction of time upper Lewis fish spend in spawning areas
)

riv_l <- rbinom(1, arrival, par$frac_lower)

init <- c(
  riv_l = riv_l,
  sp_l = 0,
  res_l = 0,
  redd_l = 0,
  riv_u = arrival - riv_l,
  trap_u = 0
)

sim_step <- function(state, par) {
  st1 <- state

  ## How many transition from river to spawning
  riv_to_sp_l <- rbinom(1, state["riv_l"], pexp(1, 1 / par$lower_tts))
  ## How many stay in the river
  riv_to_riv_l <- state["riv_l"] - riv_to_sp_l
  ## All spawners yesterday are residents today
  sp_to_res_l <- state["sp_l"]
  ## How many residents die
  res_to_dead_l <- rbinom(1, state["res_l"], pexp(1, 1 / par$lower_res))
  ## How many remain on the spawning grounds
  res_to_res_l <- state["res_l"] - res_to_dead_l
  ## How many redds are produced
  sp_to_redd_l <- rpois(1, riv_to_sp_l * par$lower_redd)

  st1["riv_l"] <- riv_to_riv_l
  st1["sp_l"] <- riv_to_sp_l
  st1["res_l"] <- sp_to_res_l + res_to_res_l
  st1["redd_l"] <- sp_to_redd_l

  riv_to_trap_u <- rbinom(1, state["riv_u"], pexp(1, 1 / par$upper_ttt))
  riv_to_riv_u <- state["riv_u"] - riv_to_trap_u

  st1["riv_u"] <- riv_to_riv_u
  st1["trap_u"] <- riv_to_trap_u

  st1
}

sim <- function(init, par) {
  out <- list(init)
  i <- 1
  while (any(out[[i]] != 0)) {
    out[[i + 1]] <- sim_step(out[[i]], par)
    i <- i + 1
  }
  out
}

st <- sim(init, par)

st_df <- st |>
  list_transpose() |>
  as_tibble() |>
  mutate(day = seq_along(riv_l))

st_df |>
  pivot_longer(cols = -day,
               names_to = "state",
               values_to = "count") |>
  mutate(dest = str_sub(state, -1, -1)) |>
  ggplot(aes(x = day, y = count, color = state)) +
  geom_line() +
  facet_wrap(~ dest, ncol = 1)

st_df |>
  mutate(in_river  = riv_l + sp_l + res_l + riv_u,
         cum_trap = cumsum(trap_u)) |>
  select(day, in_river, cum_trap) |>
  pivot_longer(-day) |>
  ggplot(aes(x = day, y = value, color = name)) +
  geom_line()
