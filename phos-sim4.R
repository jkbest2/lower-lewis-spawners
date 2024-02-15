library(tidyverse)

## Number of capture events
k <- 10
## Add the initial state and the final, trap-only event
K <- k + 2
## Potential population size
M <- 1200

##~~~~PRIORS AND CONSTRAINTS:
## Detection efficiency
p <- 0.02

## Entry rate; will sample
runsize <- 1200
## Proportion entering before each event (none in the final event; doesn't count
## residualized fish who all "enter" before the first event). Use k + 1 sequence
## because they are diff'd. This leaves k probabilities of entry and then
gamma <- seq(-2, 2, length.out = k + 1) |>
  pnorm() |>
  diff()
## No fish enter the system in the zero'th or last periods
gamma <- c(0, gamma, 0)

## uniform prior on proportion of hatchery fish in population
beta <- 0.51
## Sex ratio assumed 55%
sex_ratio <- 0.55
## Residualization ratios: males at 7.7% and females at 0.7%
resid_ratio <- c(0.077, 0.007)

## delta: probability of remaining in returning state
delta <- 0.2

## phi: probability of moving to trap
phi <- 0.66

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
  z[1] <- 1L # Always start in the "available" state

  ## Second state can transition to "river" (2) state, either as a residual or as a returner
  if (indiv$resid || rbinom(1, 1, gamma[2])) {
    z[2] <- 2L
  } else {
    z[2] <- 1L
  }

  ## Third through penultimate states are determined by transition probabilities
  for (t in 2:k) {
    if (z[t] == 1) {
      ## Transition to river wp gamma, otherwise remain available
      z[t + 1] <- sample(1:4, 1, FALSE, c(1 - gamma[t + 1], gamma[t + 1], 0, 0))
    } else if (z[t] == 2) {
      z[t + 1] <- sample(1:4, 1, FALSE, c(0, delta, 1 - delta - phi, phi))
    } else if (z[t] == 3) {
      z[t + 1] <- 3
    } else if (z[t] == 4) {
      z[t + 1] <- 4
    }
  }

  ## Final transition can only be to the trap and only if the fish is already a returner
  if (z[k + 1] == 2) {
    z[k + 2] <- sample(1:4, 1, FALSE, c(0, 1 - phi, 0, phi))
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
  if (any(z == 4)) {
    trap_t <- min(which(z == 4))
    a[(trap_t + 1):length(z)] <- 0
  }
  a
}

sim_tomat <- function(l) {
  m <- Reduce(rbind, l)
  colnames(m) <- paste("Event", 0:(ncol(m) - 1), sep = "_")
  rownames(m) <- NULL
  m
}

i <- sim_indiv(pars)
z <- sim_eco(i, pars, k = 10)
y <- sim_obs(z, pars)
zobs <- sim_obseco(z, y)
a <- sim_avail(z)


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

res <- sim_model(1200, pars)
plot_obs(res)
number_handled(res)

plot_obs <- function(res) {
  par(mfrow = c(1, 2))
  ## Collected in the river each observation event
  barplot(colSums(res$y == 2 | res$y == 3, na.rm = TRUE), main = "In river")
  ## Collected in the trap
  barplot(colSums(res$y == 4, na.rm = TRUE), main = "Trap")
}

## check_results <- function(res) {
##   nas <- lapply(res, anyNA)
##   if (any(lapply(res, anyNA)) {
##       stop("NA in results")
##   }
## }

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

pad_results <- function(res, add = 1000) {
  which_obs <- apply(res$y, 1, \(r) any(r > 1)) |> which()
  origin <- c(res$origin[which_obs],
              rep(NA, add))
  sex <- c(res$sex[which_obs],
           rep(NA, add))
  resid <- c(res$resid[which_obs],
             rep(NA, add))
  z <- rbind(res$z[which_obs, ],
             matrix(NA, nrow = add, ncol = ncol(res$z)))
  y <- rbind(res$y[which_obs, ],
             matrix(1, nrow = add, ncol = ncol(res$y)))
  avail <- rbind(res$avail[which_obs, ],
                 matrix(TRUE, nrow = add, ncol = ncol(y)))

  if (!all.equal(nrow(y), nrow(z), nrow(avail),
                 length(origin), length(sex),
                 length(resid))) {
    stop("Not all have same number of rows")
  }

  list(origin = origin,
       sex = sex,
       resid = resid,
       z = z,
       y = y,
       avail = avail)
}

prep_data <- function(res, add = 1000) {
  res_pad <- pad_results(res, add = add)

  zobs <- res_pad$z
  ## zobs[res_pad$y == 1] <- NA

  list(M = nrow(res_pad$y),
       K = ncol(res_pad$y),
       origin = as.integer(res_pad$origin),
       sex = as.integer(res_pad$sex),
       resid = as.integer(res_pad$resid),
       z = zobs,
       y = res_pad$y,
       avail = res_pad$y)

}


library(R2jags)

dat <- prep_data(res)

fit <- jags(data = dat,
            inits = NULL,
            parameters.to.save = NULL,
            model.file = "phos-model.jags",
            n.iter = 1000, n.chains = 1)
