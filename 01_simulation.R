## Functions that perform simulations of individual trajectories. To simulate a
## single trajectory:

## 1. Define parameters in a list
## 2. Simulate an individual using `sim_indiv`
## 3. Simulate the state process of that individual using `sim_eco`
## 4. Simulate the observation process using `sim_obs`
## 5. Get the observed states using `sim_obseco`
## 6. Get the availability using `sim_avail`

## These steps should be repeated for the number of individuals required. This
## can be done using `apply` or `purrr::map`-style functions, or with the
## provided `sim_model` function provided. Note that these trajectories will
## likely include individuals that never enter the system and a larger number
## that are never observed. This is consistent with the estimation model. The
## trajectories will also be a list of individual states etc., and must be
## postprocessed in order to be used in the JAGS model.
#
sim_indiv <- function(pars) {
  sex <- rbinom(1, 1, pars$sex_ratio)
  resid <- rbinom(1, 1, pars$resid_ratio[sex + 1])
  ## Not currently using origin for anything, so just make it 50/50
  origin <- rbinom(1, 1, 0.51)
  list(sex = sex, resid = resid, origin = origin)
}

sim_eco <- function(indiv, pars, k = 10) {
  gamma <- pars$gamma
  delta <- pars$delta
  phi <- pars$phi

  z <- rep(NA_integer_, k + 2)
  z[1] <- 1L # Always start in the available (1) state

  ## Second state can transition to returner (2) state
  z[2] <- sample(1:4, 1, FALSE, c(1 - gamma[2], gamma[2], 0, 0))

  ## Second through penultimate states are determined by transition
  ## probabilities
  for (t in 2:(k + 1)) {
    if (z[t - 1] == 1) {
      z[t] <- sample(1:4, 1, FALSE, c(1 - gamma[t], gamma[t], 0, 0))
    } else if (z[t - 1] == 2) {
      z[t] <- sample(1:4, 1, FALSE,
                     c(0, delta, (1 - delta) * (1 - phi), (1 - delta) * phi))
    } else if (z[t - 1] == 3) {
      z[t] <- 3
    } else if (z[t - 1] == 4) {
      z[t] <- 4
    }
  }

  ## Final transition moves all (available) fish in returner state to either
  ## spawner or trap state, though only trapped fish will be observed.
  if (z[k + 1] == 2) {
    z[k + 2] <- sample(1:4, 1, FALSE, c(0, 0, 1 - phi, phi))
  } else {
    z[k + 2] <- z[k + 1]
  }

  ## Return the state vector
  z
}

sim_obs <- function(z, pars) {
  p <- c(1, pars$p, pars$p, 1) # Probability of observation by state
  y <- rlang::rep_along(z, NA)

  for (t in seq_along(z)) {
    if (z[t] == 1) {
      y[t] <- 1
    } else if (z[t] == 2 && t < length(z)) {
      ## Include t < length(z) condition because states 2 and 3 cannot be
      ## observed in the final time period
      y[t] <- sample(1:4, 1, FALSE, c(1 - p[2], p[2], 0, 0))
    } else if (z[t] == 3 && t < length(z)) {
      y[t] <- sample(1:4, 1, FALSE, c(1 - p[3], 0, p[3], 0))
    } else if (z[t] == 4) {
      ## Observed in trap only after transition; not available after that
      if (z[t - 1] != 4) {
        y[t] <- 4
      } else {
        y[t] <- 1
      }
    } else {
      ## Default to no observation (covers final observation)
      y[t] <- 1
    }
  }
  y
}

sim_obseco <- function(z, y) {
  n <- length(z)

  zobs <- rlang::rep_along(z, NA)
  zobs[1] <- 1

  for (t in 2:n) {
    if (y[t] == 3) {
      ## No imputation if observed in state 3
      zobs[t] <- z[t]
    } else if (y[t] %in% c(2, 4)) {
      zobs[t] <- z[t]
      ## Check for previous observations so that the "returning" state can be
      ## imputed. Use the max because if there is more than one returning
      ## observation the intervening nonobserved states are already imputed
      prev_obs <- head(y, t - 1)
      if (2 %in% prev_obs) {
        prev_retobs <- max(which(prev_obs == 2))
        zobs[prev_retobs:(t - 1)] <- 2
      }
    } else if (zobs[t - 1] %in% c(3, 4)) {
      ## If previously *observed* as 3 or 4, continue that forward
      zobs[t] <- z[t]
    }
  }

  zobs
}

sim_avail <- function(z) {
  a <- rlang::rep_along(z, 1)
  n <- length(z)

  ## Only need to mark as unavailable if collected in the trap *before* the
  ## final event
  if (any(z[seq_len(n - 1)] == 4)) {
    trap_t <- min(which(z == 4))
    a[(trap_t + 1):length(z)] <- 0
  }
  stopifnot(length(a) == length(z))
  a
}

sim_model <- function(n, pars) {
  indivs <- replicate(n, sim_indiv(pars), simplify = FALSE)
  indiv_df <- purrr::list_transpose(indivs) |>
    tibble::as_tibble()
  zsfull <- purrr::map(indivs, sim_eco, pars = pars, k = 10)
  ys <- purrr::map(zsfull, sim_obs, pars = pars)
  zs <- purrr::map2(zsfull, ys, sim_obseco)
  as <- purrr::map(zsfull, sim_avail)

  list(z = sim_tomat(zs),
       zfull = sim_tomat(zsfull),
       y = sim_tomat(ys),
       avail = sim_tomat(as),
       origin = indiv_df$origin,
       sex = indiv_df$sex,
       resid = indiv_df$resid)
}
