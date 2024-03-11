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

prep_pars <- function(k = 10, gamma, beta, sex_ratio, resid_ratio, delta, phi,
                      p) {
  if (length(gamma) == 1)
    gamma <- rep(gamma, k)
  if (length(gamma) != k)
    stop("gamma must be length 1 or k")

  if (length(beta) != 1)
    stop("beta must be length 1")

  if (length(sex_ratio) != 1)
    stop("sex_ratio must be length 1")

  if (length(resid_ratio) == 1)
    resid_ratio <- rep(resid_ratio, 2)
  if (length(resid_ratio) != 2)
    stop("resid_ratio must be length 2")

  if (length(delta) == 1)
    delta <- matrix(rep(delta, (k - 1) * 4), ncol = 4)
  if (length(delta) == k - 1)
    delta <- matrix(rep(delta, 4), ncol = 4)
  if (nrow(delta) != k - 1 | ncol(delta) != 4)
    stop("delta must be length, 1, k - 1, or be a matrix with k - 1 rows and 4 columns")

  if (length(phi) == 1)
    phi <- matrix(rep(phi, (k) * 4), ncol = 4)
  if (length(phi) == k)
    phi <- matrix(rep(phi, 4), ncol = 4)
  if (nrow(phi) != k | ncol(phi) != 4)
    stop("phi must be length, 1, k, or be a matrix with k rows and 4 columns")

  ## Change this later for variable detection probabilities
  if (length(p) != 1)
    stop("p must be length 1")

  list(k = k,
       gamma = gamma,
       beta = beta,
       sex_ratio = sex_ratio,
       resid_ratio = resid_ratio,
       delta = delta,
       phi = phi,
       p = p)
}

sim_indiv <- function(pars) {
  sex <- rbinom(1, 1, pars$sex_ratio)
  resid <- rbinom(1, 1, pars$resid_ratio[sex + 1])
  ## Not currently using origin for anything, so just make it 50/50
  origin <- rbinom(1, 1, 0.51)
  group <- 1 + origin + 2 * resid
  list(sex = sex, resid = resid, origin = origin, group = group)
}

## gamma is probability of entering the system. This is the same for all groups,
## and is length k, representing transitions *to* idx + 1. So gamma[1]
## represents the probability of transitioning to available state at time 2.

## delta is probability of remaining in available state. It is length k - 1 for
## each group, as this transition from available -> available can only occur at
## transitions transition times 2 -> 3, 3 -> 4, ..., k -> k + 1. Note that times
## here are indexed by the time at the end of the transition.

## phi is the probability of entering the trap. It is length k for each group,
## as it can only occur after time 2.
sim_eco <- function(indiv, pars) {
  k <- pars$k
  gamma <- pars$gamma
  delta <- pars$delta[, indiv$group]
  phi <- pars$phi[, indiv$group]

  z <- rep(NA_integer_, k + 2)
  z[1] <- 1L # Always start in the available (1) state

  ## Second through penultimate states are determined by transition
  ## probabilities
  for (t in 2:(k + 1)) {
    if (z[t - 1] == 1) {
      z[t] <- sample(1:4, 1, FALSE, c(1 - gamma[t - 1], gamma[t - 1], 0, 0))
    } else if (z[t - 1] == 2) {
      ## FIXME Is the indexing correct here? Concerned about phi in particular!
      z[t] <- sample(1:4, 1, FALSE,
                     c(0, delta[t - 2], (1 - delta[t - 2]) * (1 - phi[t - 2]), (1 - delta[t - 2]) * phi[t - 2]))
    } else if (z[t - 1] == 3) {
      z[t] <- 3
    } else if (z[t - 1] == 4) {
      z[t] <- 4
    }
  }

  ## Final transition moves all (available) fish in returner state to either
  ## spawner or trap state, though only trapped fish will be observed.
  if (z[k + 1] == 2) {
    z[k + 2] <- sample(1:4, 1, FALSE, c(0, 0, 1 - phi[k], phi[k]))
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

sim_model <- function(n = NULL, indivs = NULL, zsfull = NULL, pars) {
  ## Check inputs for consistency
  if (is.null(n) && is.null(indivs)) {
    stop("Must provide either n or indivs")
  }
  if (!is.null(n) && !is.null(indivs) && length(indivs) != n) {
    stop("If both specified, length(indivs) must be n")
  }
  if (is.null(indivs) && !is.null(zsfull)) {
    stop("indivs must be provided if zsfull is provided")
   }
  if (is.null(n) & !is.null(indivs)) {
    n <- length(indivs)
  }

  ## If not provided, generate a list of individual fish
  if (is.null(indivs)) {
    indivs <- replicate(n, sim_indiv(pars), simplify = FALSE)
  }
  ## Summarize in to data frame so that identity vectors can be extracted later
  indiv_df <- purrr::list_transpose(indivs) |>
    tibble::as_tibble()
  ## Generate state vectors if not provided
  if (is.null(zsfull)) {
    zsfull <- purrr::map(indivs, sim_eco, pars = pars)
  }
  ## Generate observation, oberved states, and availability vectors
  ys <- purrr::map(zsfull, sim_obs, pars = pars)
  zs <- purrr::map2(zsfull, ys, sim_obseco)
  as <- purrr::map(zsfull, sim_avail)

  ## Combine into matrices for passing to JAGS
  list(z = sim_tomat(zs),
       zfull = sim_tomat(zsfull),
       y = sim_tomat(ys),
       avail = sim_tomat(as),
       origin = indiv_df$origin,
       sex = indiv_df$sex,
       resid = indiv_df$resid)
}
