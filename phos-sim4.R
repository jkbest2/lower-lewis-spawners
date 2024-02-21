library(tidyverse)

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
  y <- rep_along(z, NA)

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

  zobs <- rep_along(z, NA)
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
  a <- rep_along(z, 1)
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

sim_tomat <- function(l) {
  m <- Reduce(rbind, l)
  colnames(m) <- paste("Event", 0:(ncol(m) - 1), sep = "_")
  rownames(m) <- NULL
  m
}

sim_model <- function(n, pars) {
  indivs <- replicate(n, sim_indiv(pars), simplify = FALSE)
  indiv_df <- list_transpose(indivs) |>
    as_tibble()
  zsfull <- map(indivs, sim_eco, pars = pars, k = 10)
  ys <- map(zsfull, sim_obs, pars = pars)
  zs <- map2(zsfull, ys, sim_obseco)
  as <- map(zsfull, sim_avail)

  list(z = sim_tomat(zs),
       zfull = sim_tomat(zsfull),
       y = sim_tomat(ys),
       avail = sim_tomat(as),
       origin = indiv_df$origin,
       sex = indiv_df$sex,
       resid = indiv_df$resid)
}

plot_obs <- function(res) {
  par(mfrow = c(1, 2))
  ## Collected in the river each observation event
  barplot(colSums(res$y == 2 | res$y == 3, na.rm = TRUE), main = "In river")
  ## Collected in the trap
  barplot(colSums(res$y == 4, na.rm = TRUE), main = "Trap")
}

number_handled <- function(res) {
  y <- res$y
  ## Total number of fish handled in river
  mark <- apply(y, 1, \(r) 2 %in% r | 3 %in% r) |> sum()
  ## Total number of fish trapped
  trap <- apply(y, 1, \(r) 4 %in% r) |> sum()
  ## Total number of *recaptures*
  recap <- apply(y, 1, \(r) sum(r > 1) > 1) |> sum()

  list(marked = mark, trapped = trap, recaptured = recap)
}

  }

}

prep_data <- function(res) {
  obs <- apply(res$y, 1, \(r) any(r != 1))

  ## Split each observation type into observed and unobserved for data and inits
  ## respectively
  z_init <- res$zfull
  z_init[!is.na(res$z)] <- NA

  origin_obs <- res$origin
  origin_obs[!obs] <- NA
  origin_init <- res$origin
  origin_init[obs] <- NA

  sex_obs <- res$sex
  sex_obs[!obs] <- NA
  sex_init <- res$sex
  sex_init[obs] <- NA

  resid_obs <- res$resid
  resid_obs[!obs] <- NA
  resid_init <- res$resid
  resid_init[obs] <- NA

  data <- list(M = nrow(res$y),
               K = ncol(res$y),
               origin = as.integer(origin_obs),
               sex = as.integer(sex_obs),
               resid = as.integer(resid_obs),
               z = res$z,
               y = res$y,
               avail = res$avail)
  init <- function() {
    list(z = z_init,
         origin = origin_init,
         sex = sex_init,
         resid = resid_init)
    }

  list(data = data,
       init = init)
}

sim_check <- function(dat) {
  M <- dat$data$M
  K <- dat$data$K

  init <- dat$init()

  ## Check data dimensions
  stopifnot(
    all(dim(dat$data$y) == c(M, K)),
    all(dim(dat$data$z) == c(M, K)),
    all(dim(dat$data$avail) == c(M, K)),
    length(dat$data$origin) == M,
    length(dat$data$sex) == M,
    length(dat$data$resid) == M
    )
  ## Check initial value dimensions
  stopifnot(
    all(dim(init$z) == c(M, K)),
    length(init$origin) == M,
    length(init$sex) == M,
    length(init$resid) == M
  )
  ## Check that there is one NA for data and one for intial value
  stopifnot(
    all(xor(is.na(dat$data$z), is.na(init$z))),
    all(xor(is.na(dat$data$origin), is.na(init$origin))),
    all(xor(is.na(dat$data$sex), is.na(init$sex))),
    all(xor(is.na(dat$data$resid), is.na(init$resid)))
  )
  invisible(dat)
}

library(rjags)
library(runjags) # Loading this first prevents warnings about JAGS version
library(posterior) # Useful for postprocessing posterior samples

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
