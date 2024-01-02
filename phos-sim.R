## Number of capture events
K <- 10
M <- 1000
avail <- matrix(1, nrow = M, ncol = K)

##~~~~PRIORS AND CONSTRAINTS:
## prior on tangle-net detetion efficiency which we expect to be low
mean.p <- rnorm(1, -1.5, 1)
## prior on standard deviation for detection efficieny random effect
sigma.p <- abs(rnorm(1, 0, 0.5))
## JAGS operates on precision rather than variance
## tau.p <- pow(sigma.p, -2)

## mean entry probability
mean.gamma <- rt(1, 7) * sqrt(1 / 0.1)
## standard deviation for entry probability random effects (per event)
sigma.gamma <- abs(rnorm(1, 0, sqrt(1 / 0.5)))
## tau.gamma <- pow(sigma.gamma, -2)

beta_delta <- rep(NA, 4)
##intercept term for probability of becoming spawner
beta_delta[1] <- rt(1, 7) * sqrt(1 / 0.1)
## offset for natural origin
beta_delta[2] <- rt(1, 7) * sqrt(1 / 0.25)
## offest for residuals
beta_delta[3] <- rt(1, 7) * sqrt(1 / 0.25)
## interaction of origin and residual status
beta_delta[4] <- rt(1, 7) * sqrt(1 / 0.25)

## standard deviation for probability of becoming spawner random effect
sigma.delta <- abs(rnorm(1, 0, sqrt(1 / 0.5)))
## tau.delta <- pow(sigma.delta, -2)

beta_phi <- rep(NA, 4)
## intercept term for probability of moving to trap
beta_phi[1] <- rt(1, 7) * sqrt(1 / 0.1)
beta_phi[2] <- rt(1, 7) * sqrt(1 / 0.25)
beta_phi[3] <- rt(1, 7) * sqrt(1 / 0.25)
beta_phi[4] <- rt(1, 7) * sqrt(1 / 0.25)
## standard deviation for probability of moving to trap random effect
sigma.phi <- abs(rnorm(1, 0, sqrt(1 / 0.5)))
## tau.phi <- pow(sigma.phi, -2)

## uniform prior on proportion of hatchery fish in population
## Made more informative based on past results (close to 50/50)
beta <- rbeta(1, 10, 10)
## more informative prior on sex ratio since it should be close to 50%
sex_ratio <- rbeta(1, 5, 5)
## % of male population that residualize
resid_ratio <- rep(NA, 2)
## expect this to be low, this prior has a mean of 10%
resid_ratio[1] <- rbeta(1, 1, 9)
## % of female population that residualize
resid_ratio[2] <- rbeta(1, 1, 9)

## gamma: probability of entering from potential population
## Event 0 through Event K - 2
e.gamma <- rnorm(K - 2, 0, sigma.gamma)
gamma <- plogis(mean.gamma + e.gamma)

## delta: probability of remaining in returning state
## Event 1 through Event K - 2
e.delta <- matrix(NA, nrow = K - 3, ncol = 4)
delta <- matrix(NA, nrow = K - 3, ncol = 4)
for (t in 1:(K - 3)) {
  ## Draw different random effect for each category of fish
  ## (More flexible than single shared random effect)
  e.delta[t, 1] <- rnorm(1, 0, sigma.delta)
  e.delta[t, 2] <- rnorm(1, 0, sigma.delta)
  e.delta[t, 3] <- rnorm(1, 0, sigma.delta)
  e.delta[t, 4] <- rnorm(1, 0, sigma.delta)
  ## probability of transitioning to kelt for hatchery anadramous
  delta[t, 1] <- plogis(beta_delta[1] + e.delta[t, 1])
  ## probability of transitioning to spawner for natural anadramous
  delta[t, 2] <- plogis(beta_delta[1] + beta_delta[2] + e.delta[t, 2])
  ## probability of transitioning to spawner for hatchery residuals
  delta[t, 3] <- plogis(beta_delta[1] + beta_delta[3] + e.delta[t, 3])
  ## probability of transitioning to spawner for natural residuals
  delta[t, 4] <- plogis(beta_delta[1] + beta_delta[2] + beta_delta[3] +
                        beta_delta[4] + e.delta[t, 4])
}

## phi: probability of moving to trap
## Event 1 through Event K - 1
e.phi <- matrix(NA, nrow = K - 2, ncol = 4)
phi <- matrix(NA, nrow = K - 2, ncol = 4)
for (t in 1:(K - 2)) {
  e.phi[t, 1] <- rnorm(1, 0, sigma.phi)
  e.phi[t, 2] <- rnorm(1, 0, sigma.phi)
  e.phi[t, 3] <- rnorm(1, 0, sigma.phi)
  e.phi[t, 4] <- rnorm(1, 0, sigma.phi)
  ## Draw different random effect for each category of fish
  ## (More flexible than single shared random effect)
  phi[t, 1] <- plogis(beta_phi[1] + e.phi[t, 1])
  ## probability of transitioning to the trap for hatchery anadramous
  phi[t, 2] <- plogis(beta_phi[1] + beta_phi[2] + e.phi[t, 2])
  ## probability of transitioning to the trap for natural anadramous
  phi[t, 3] <- plogis(beta_phi[1] + beta_phi[3] + e.phi[t, 3])
  ## probability of transitioning to the trap for hatchery residuals
  phi[t, 4] <- plogis(beta_phi[1] + beta_phi[2] + beta_phi[3] +
                      beta_phi[4] + e.phi[t, 4])
  ## probability of transitioning to the trap for natural residuals
}

## p: detection probability
## Event 1 through Event K - 1
e.p <- rep(NA, K - 2)
p <- rep(NA, K - 2)
for (t in 1:(K - 2)) {
  e.p[t] <- rnorm(1, 0, sigma.p)
  ## Detection efficiency varies among capture events around a common mean
  p[t] <- plogis(mean.p + e.p[t])
}

## Individual covariates
origin <- rep(NA, M)
sex <- rep(NA, M)
resid <- rep(NA, M)
## Group index, determines which offsets are included in the model for each individual
g <- rep(NA, M)
## Transition matrices
ps <- array(NA, dim = c(K - 1, M, 4, 4))
## Observation matrices
po <- array(NA, dim = c(K - 1, M, 4, 4))
for (i in 1:M){
  origin[i] <- rbinom(1, 1, beta) ## 0: HOR, 1: NOR
  sex[i] <- rbinom(1, 1, sex_ratio) ## 0: M, 1: F
  resid[i] <- rbinom(1, 1, resid_ratio[1 + sex[i]]) ## 0: anadramous, 1: residual

  g[i] <- 1 + 1 * origin[i] + 2 * resid[i]

  ##~~ TRANSITION ~~#
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## Indices as follows:
  ## [time, individual, current state, next state]
  ## Initial entry to Lewis River, k = 1
  ## Event 0 -> Event 1
  ps[1, i, 1, 1] <- 1 - gamma[t]
  ps[1, i, 1, 2] <- gamma[t]
  ps[1, i, 1, 3] <- 0
  ps[1, i, 1, 4] <- 0

  ps[1, i, 2, 1] <- 0
  ps[1, i, 2, 2] <- 1
  ps[1, i, 2, 3] <- 0
  ps[1, i, 2, 4] <- 0

  ps[1, i, 3, 1] <- 0
  ps[1, i, 3, 2] <- 0
  ps[1, i, 3, 3] <- 1
  ps[1, i, 3, 4] <- 0

  ps[1, i, 4, 1] <- 0
  ps[1, i, 4, 2] <- 0
  ps[1, i, 4, 3] <- 0
  ps[1, i, 4, 4] <- 1

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ##~~ OBSERVATION ~~#
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## [time, individual, actual state, observed state]
  ## Event 1, k = 1
  po[1, i, 1, 1] <- 1
  po[1, i, 1, 2] <- 0
  po[1, i, 1, 3] <- 0
  po[1, i, 1, 4] <- 0

  po[1, i, 2, 1] <- 1 - p[t] * avail[i, 2]
  po[1, i, 2, 2] <- p[t] * avail[i, 2]
  po[1, i, 2, 3] <- 0
  po[1, i, 2, 4] <- 0

  po[1, i, 3, 1] <- 1 - p[t] * avail[i, 2]
  po[1, i, 3, 2] <- 0
  po[1, i, 3, 3] <- p[t] * avail[i, 2]
  po[1, i, 3, 4] <- 0

  po[1, i, 4, 1] <- 1 - avail[i, 2]
  po[1, i, 4, 2] <- 0
  po[1, i, 4, 3] <- 0
  po[1, i, 4, 4] <- avail[i, 2]

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  for (k in 2:(K - 2)) {
    ##~~ TRANSITION ~~#
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ## Event 1 -> Event 2, ..., Event K - 2 -> Event K - 1
    ps[k, i, 1, 1] <- 1 - gamma[k]
    ps[k, i, 1, 2] <- gamma[k]
    ps[k, i, 1, 3] <- 0
    ps[k, i, 1, 4] <- 0

    ps[k, i, 2, 1] <- 0
    ps[k, i, 2, 2] <- (delta[k - 1, g[i]] * avail[i, k]) + (1 - avail[i, k])
    ps[k, i, 2, 3] <- (1 - delta[k - 1, g[i]]) * (1 - phi[k - 1, g[i]]) * (avail[k, t])
    ps[k, i, 2, 4] <- (1 - delta[k - 1, g[i]]) * (phi[k - 1, g[i]]) * (avail[k, t])

    ps[k, i, 3, 1] <- 0
    ps[k, i, 3, 2] <- 0
    ps[k, i, 3, 3] <- 1
    ps[k, i, 3, 4] <- 0

    ps[k, i, 4, 1] <- 0
    ps[k, i, 4, 2] <- 0
    ps[k, i, 4, 3] <- 0
    ps[k, i, 4, 4] <- 1

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ##~~ OBSERVATION ~~#
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ## Event 2, ..., Event K - 1
    po[k, i, 1, 1] <- 1
    po[k, i, 1, 2] <- 0
    po[k, i, 1, 3] <- 0
    po[k, i, 1, 4] <- 0

    po[k, i, 2, 1] <- 1 - p[t] * avail[i, k + 1]
    po[k, i, 2, 2] <- p[t] * avail[i, k + 1]
    po[k, i, 2, 3] <- 0
    po[k, i, 2, 4] <- 0

    po[k, i, 3, 1] <- 1 - p[t] * avail[i, k + 1]
    po[k, i, 3, 2] <- 0
    po[k, i, 3, 3] <- p[t] * avail[i, k + 1]
    po[k, i, 3, 4] <- 0

    po[k, i, 4, 1] <- 1 - avail[i, k + 1]
    po[k, i, 4, 2] <- 0
    po[k, i, 4, 3] <- 0
    po[k, i, 4, 4] <- avail[i, k + 1]
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  }

  ##~~ TRANSITION ~~#
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## Event K - 1 -> Event K
  ## Turn all fish remaining on spawning grounds into spawners
  ps[K - 1, i, 1, 1] <- 1
  ps[K - 1, i, 1, 2] <- 0
  ps[K - 1, i, 1, 3] <- 0
  ps[K - 1, i, 1, 4] <- 0

  ps[K - 1, i, 2, 1] <- 0
  ps[K - 1, i, 2, 2] <- 1 - avail[i, K - 1]
  ps[K - 1, i, 2, 3] <- (1 - phi[K - 2, g[i]]) * avail[i, K - 1]
  ps[K - 1, i, 2, 4] <- phi[K - 2, g[i]] * avail[i, K - 1]

  ps[K - 1, i, 3, 1] <- 0
  ps[K - 1, i, 3, 2] <- 0
  ps[K - 1, i, 3, 3] <- 1
  ps[K - 1, i, 3, 4] <- 0

  ps[K - 1, i, 4, 1] <- 0
  ps[K - 1, i, 4, 2] <- 0
  ps[K - 1, i, 4, 3] <- 0
  ps[K - 1, i, 4, 4] <- 1

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ##~~ OBSERVATION ~~#
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## Event K
  ## Only trap observations after last tangle-net event.
  po[K - 1, i, 1, 1] <- 1
  po[K - 1, i, 1, 2] <- 0
  po[K - 1, i, 1, 3] <- 0
  po[K - 1, i, 1, 4] <- 0

  po[K - 1, i, 2, 1] <- 1
  po[K - 1, i, 2, 2] <- 0
  po[K - 1, i, 2, 3] <- 0
  po[K - 1, i, 2, 4] <- 0

  po[K - 1, i, 3, 1] <- 1
  po[K - 1, i, 3, 2] <- 0
  po[K - 1, i, 3, 3] <- 0
  po[K - 1, i, 3, 4] <- 0

  po[K - 1, i, 4, 1] <- 1 - avail[i, K]
  po[K - 1, i, 4, 2] <- 0
  po[K - 1, i, 4, 3] <- 0
  po[K - 1, i, 4, 4] <- avail[i, K]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
}

## rcat <- function(x, p) {
##   sample(x, prob = p)
## }

z <- matrix(NA, nrow = M, ncol = K)
y <- matrix(NA, nrow = M, ncol = K)
avail <- matrix(NA, nrow = M, ncol = K)
avail[, 1] <- TRUE
###~~~ LIKELIHOOD:
for (i in 1:M){
  if (resid[i] || rbinom(1, 1, gamma[1])) {
    z[i, 1] <- 2
    y[i, 1] <- rbinom()
  } else {
    z[i, 1] <- 1
    y[i, 1] <- 1
  }
  for (k in 2:K){
    z[i, k] <- sample(x = 1:4,
                      size = 1,
                      replace = FALSE,
                      prob = ps[k - 1, i, z[i, k - 1], ])
    ## z[i, k] ~ dcat(ps[k - 1, i, z[i, k - 1], ])
    y[i, k] <- sample(1:4,
                      size = 1,
                      replace = FALSE,
                      prob = po[k - 1, i, z[i, k], ])
    ## y[i, k] ~ dcat(po[k - 1, i, z[i, k], ])
  }
}

fix_trap_obs <- function(r) {
  idx <- which(r == 4)
  ## Replace the repeated trap observations with non-observations because the
  ## individual is not available after the first trap observation
  r[idx[-which.min(idx)]] <- 1
  r
}

for (r in seq_len(nrow(y)))


###~~~ DERIVED PARAMETERS:
## for (i in 1:M){
##   NOR_entered[i] <- (1 - equals(z[i, K], 1)) * equals(origin[i], 1)
##   ## If latent state is not equal to 1 at last event, and latent origin state
##   ## is equal to 1 count as natural origin fish returning to Lewis River

##   HOR_entered[i] <- (1 - equals(z[i, K], 1)) * equals(origin[i], 0)
##   ## If latent state is not equal to 1 at last event, and latent origin state
##   ## is equal to 0 count as hatchery origin fish returning to Lewis River

##   NOR_spawn[i] <- equals(z[i, K], 3) * equals(origin[i], 1)
##   ## If latent state is equal to 3 at last event, and latent origin state is
##   ## equal to 0 count as natural origin spawner

##   HOR_spawn[i] <- equals(z[i, K], 3) * equals(origin[i], 0)
##   ## If latent state is equal to 3 at last event, and latent origin state is
##   ## equal to 0 count as hatchery origin spawner
## }

## N_NOR_entered <- sum(NOR_entered[])
## N_HOR_entered <- sum(HOR_entered[])
## N_NOR_spawn <- sum(NOR_spawn[])
## N_HOR_spawn <- sum(HOR_spawn[])

## pHOS <- N_HOR_spawn/(N_HOR_spawn + N_NOR_spawn)

###---------------------------------------------

## Individual covariates
origin <- rep(NA, M)
sex <- rep(NA, M)
resid <- rep(NA, M)
## Availability
## Group index, determines which offsets are included in the model for each individual
g <- rep(NA, M)
## Transition matrices
## ps <- array(NA, dim = c(K - 1, M, 4, 4))
## Observation matrices
## po <- array(NA, dim = c(K - 1, M, 4, 4))
avail <- matrix(NA, nrow = M, ncol = K)
z <- matrix(NA, nrow = M, ncol = K)
y <- matrix(NA, nrow = M, ncol = K)
for (i in 1:M) {
  origin[i] <- rbinom(1, 1, beta) ## 0: HOR, 1: NOR
  sex[i] <- rbinom(1, 1, sex_ratio) ## 0: M, 1: F
  resid[i] <- rbinom(1, 1, resid_ratio[1 + sex[i]]) ## 0: anadramous, 1: residual

  g[i] <- 1 + 1 * origin[i] + 2 * resid[i]

  ## Initial state; available if residual or enters before observation event 1
  avail[i, 1] <- TRUE
  if (resid[i] || rbinom(1, 1, gamma[1])) {
    z[i, 1] <- 2
    y[i, 1] <- sample(1:4, 1, TRUE, c(1 - p[1], p[1], 0, 0))
  } else {
    z[i, 1] <- 1 # Still potential; not in the river yet
    y[i, 1] <- 1 # Not observed
  }

  for (k in 2:(K - 2)) {
    ## Update availability
    if (z[i, k - 1] == 4.1) {
      avail[i, k] <- FALSE
    } else {
      avail[i, k] <- TRUE
    }

    ## Update state
    if (z[i, k - 1] == 1) {
      ps <- c(1 - gamma[k], gamma[k], 0, 0)
      z[i, k] <- sample(1:4, 1, TRUE, ps)
    } else if (z[i, k - 1] == 2) {
      ps <- c(
        0,
        delta[k - 1, g[i]],
        (1 - delta[k - 1, g[i]]) * (1 - phi[k - 1, g[i]]),
        (1 - delta[k - 1, g[i]]) * (phi[k - 1, g[i]]))
      z[i, k] <- sample(1:4, 1, TRUE, ps)
    } else if (z[i, k - 1] == 3) {
      z[i, k] <- 3
    } else if (z[i, k - 1] >= 4) {
      z[i, k] <- 4.1
    }

    ## Update observations
    if (z[i, k] == 1) {
      y[i, k] <- 1
    } else if (z[i, k] == 2) {
      po <- c(1 - p[k], p[k], 0, 0)
      y[i, k] <- sample(1:4, 1, FALSE, po)
    } else if (z[i, k] == 3) {
      po <- c(1 - p[k], 0, p[k], 0)
      y[i, k] <- sample(1:4, 1, FALSE, po)
    } else if (z[i, k] == 4) {
      y[i, k] <- 4
    } else if (z[i, k] == 4.1) {
      ## If state is 4.1, already trapped, it will not be available and will not
      ## be observed
      y[i, k] <- 1
    }
  }
}

  ##~~ TRANSITION ~~#
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## Indices as follows:
  ## [time, individual, current state, next state]
  ## Event 0 -> Event 1
  ps[1, i, 1, 1] <- 1 - gamma[t]
  ps[1, i, 1, 2] <- gamma[t]
  ps[1, i, 1, 3] <- 0
  ps[1, i, 1, 4] <- 0

  ps[1, i, 2, 1] <- 0
  ps[1, i, 2, 2] <- 1
  ps[1, i, 2, 3] <- 0
  ps[1, i, 2, 4] <- 0

  ps[1, i, 3, 1] <- 0
  ps[1, i, 3, 2] <- 0
  ps[1, i, 3, 3] <- 1
  ps[1, i, 3, 4] <- 0

  ps[1, i, 4, 1] <- 0
  ps[1, i, 4, 2] <- 0
  ps[1, i, 4, 3] <- 0
  ps[1, i, 4, 4] <- 1

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ##~~ OBSERVATION ~~#
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## [time, individual, actual state, observed state]
  ## Event 1, k = 1
  po[1, i, 1, 1] <- 1
  po[1, i, 1, 2] <- 0
  po[1, i, 1, 3] <- 0
  po[1, i, 1, 4] <- 0

  po[1, i, 2, 1] <- 1 - p[t] * avail[i, 2]
  po[1, i, 2, 2] <- p[t] * avail[i, 2]
  po[1, i, 2, 3] <- 0
  po[1, i, 2, 4] <- 0

  po[1, i, 3, 1] <- 1 - p[t] * avail[i, 2]
  po[1, i, 3, 2] <- 0
  po[1, i, 3, 3] <- p[t] * avail[i, 2]
  po[1, i, 3, 4] <- 0

  po[1, i, 4, 1] <- 1 - avail[i, 2]
  po[1, i, 4, 2] <- 0
  po[1, i, 4, 3] <- 0
  po[1, i, 4, 4] <- avail[i, 2]

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  for (k in 2:(K - 2)) {
    ##~~ TRANSITION ~~#
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ## Event 1 -> Event 2, ..., Event K - 2 -> Event K - 1
    ps[k, i, 1, 1] <- 1 - gamma[k]
    ps[k, i, 1, 2] <- gamma[k]
    ps[k, i, 1, 3] <- 0
    ps[k, i, 1, 4] <- 0

    ps[k, i, 2, 1] <- 0
    ps[k, i, 2, 2] <- (delta[k - 1, g[i]] * avail[i, k]) + (1 - avail[i, k])
    ps[k, i, 2, 3] <- (1 - delta[k - 1, g[i]]) * (1 - phi[k - 1, g[i]]) * (avail[k, t])
    ps[k, i, 2, 4] <- (1 - delta[k - 1, g[i]]) * (phi[k - 1, g[i]]) * (avail[k, t])

    ps[k, i, 3, 1] <- 0
    ps[k, i, 3, 2] <- 0
    ps[k, i, 3, 3] <- 1
    ps[k, i, 3, 4] <- 0

    ps[k, i, 4, 1] <- 0
    ps[k, i, 4, 2] <- 0
    ps[k, i, 4, 3] <- 0
    ps[k, i, 4, 4] <- 1

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ##~~ OBSERVATION ~~#
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ## Event 2, ..., Event K - 1
    po[k, i, 1, 1] <- 1
    po[k, i, 1, 2] <- 0
    po[k, i, 1, 3] <- 0
    po[k, i, 1, 4] <- 0

    po[k, i, 2, 1] <- 1 - p[t] * avail[i, k + 1]
    po[k, i, 2, 2] <- p[t] * avail[i, k + 1]
    po[k, i, 2, 3] <- 0
    po[k, i, 2, 4] <- 0

    po[k, i, 3, 1] <- 1 - p[t] * avail[i, k + 1]
    po[k, i, 3, 2] <- 0
    po[k, i, 3, 3] <- p[t] * avail[i, k + 1]
    po[k, i, 3, 4] <- 0

    po[k, i, 4, 1] <- 1 - avail[i, k + 1]
    po[k, i, 4, 2] <- 0
    po[k, i, 4, 3] <- 0
    po[k, i, 4, 4] <- avail[i, k + 1]
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  }

  ##~~ TRANSITION ~~#
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## Event K - 1 -> Event K
  ## Turn all fish remaining on spawning grounds into spawners
  ps[K - 1, i, 1, 1] <- 1
  ps[K - 1, i, 1, 2] <- 0
  ps[K - 1, i, 1, 3] <- 0
  ps[K - 1, i, 1, 4] <- 0

  ps[K - 1, i, 2, 1] <- 0
  ps[K - 1, i, 2, 2] <- 1 - avail[i, K - 1]
  ps[K - 1, i, 2, 3] <- (1 - phi[K - 2, g[i]]) * avail[i, K - 1]
  ps[K - 1, i, 2, 4] <- phi[K - 2, g[i]] * avail[i, K - 1]

  ps[K - 1, i, 3, 1] <- 0
  ps[K - 1, i, 3, 2] <- 0
  ps[K - 1, i, 3, 3] <- 1
  ps[K - 1, i, 3, 4] <- 0

  ps[K - 1, i, 4, 1] <- 0
  ps[K - 1, i, 4, 2] <- 0
  ps[K - 1, i, 4, 3] <- 0
  ps[K - 1, i, 4, 4] <- 1

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ##~~ OBSERVATION ~~#
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## Event K
  ## Only trap observations after last tangle-net event.
  po[K - 1, i, 1, 1] <- 1
  po[K - 1, i, 1, 2] <- 0
  po[K - 1, i, 1, 3] <- 0
  po[K - 1, i, 1, 4] <- 0

  po[K - 1, i, 2, 1] <- 1
  po[K - 1, i, 2, 2] <- 0
  po[K - 1, i, 2, 3] <- 0
  po[K - 1, i, 2, 4] <- 0

  po[K - 1, i, 3, 1] <- 1
  po[K - 1, i, 3, 2] <- 0
  po[K - 1, i, 3, 3] <- 0
  po[K - 1, i, 3, 4] <- 0

  po[K - 1, i, 4, 1] <- 1 - avail[i, K]
  po[K - 1, i, 4, 2] <- 0
  po[K - 1, i, 4, 3] <- 0
  po[K - 1, i, 4, 4] <- avail[i, K]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
}

## rcat <- function(x, p) {
##   sample(x, prob = p)
## }

z <- matrix(NA, nrow = M, ncol = K)
y <- matrix(NA, nrow = M, ncol = K)
avail <- matrix(NA, nrow = M, ncol = K)
avail[, 1] <- TRUE
###~~~ LIKELIHOOD:
for (i in 1:M){
  if (resid[i] || rbinom(1, 1, gamma[1])) {
    z[i, 1] <- 2
    y[i, 1] <- rbinom()
  } else {
    z[i, 1] <- 1
    y[i, 1] <- 1
  }
  for (k in 2:K){
    ps <- c(1 - gamma[t], gamma[t], 0, 0)
    z[i, k] <- sample(x = 1:r)
    z[i, k] <- sample(x = 1:4,
                      size = 1,
                      replace = FALSE,
                      prob = ps[k - 1, i, z[i, k - 1], ])
    ## z[i, k] ~ dcat(ps[k - 1, i, z[i, k - 1], ])
    y[i, k] <- sample(1:4,
                      size = 1,
                      replace = FALSE,
                      prob = po[k - 1, i, z[i, k], ])
    ## y[i, k] ~ dcat(po[k - 1, i, z[i, k], ])
  }
}

  ps[1, i, 1, 1] <- 1 - gamma[t]
  ps[1, i, 1, 2] <- gamma[t]
  ps[1, i, 1, 3] <- 0
  ps[1, i, 1, 4] <- 0

  ps[1, i, 2, 1] <- 0
  ps[1, i, 2, 2] <- 1
  ps[1, i, 2, 3] <- 0
  ps[1, i, 2, 4] <- 0

  ps[1, i, 3, 1] <- 0
  ps[1, i, 3, 2] <- 0
  ps[1, i, 3, 3] <- 1
  ps[1, i, 3, 4] <- 0

  ps[1, i, 4, 1] <- 0
  ps[1, i, 4, 2] <- 0
  ps[1, i, 4, 3] <- 0
  ps[1, i, 4, 4] <- 1

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ##~~ OBSERVATION ~~#
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## [time, individual, actual state, observed state]
  ## Event 1, k = 1
  po[1, i, 1, 1] <- 1
  po[1, i, 1, 2] <- 0
  po[1, i, 1, 3] <- 0
  po[1, i, 1, 4] <- 0

  po[1, i, 2, 1] <- 1 - p[t] * avail[i, 2]
  po[1, i, 2, 2] <- p[t] * avail[i, 2]
  po[1, i, 2, 3] <- 0
  po[1, i, 2, 4] <- 0

  po[1, i, 3, 1] <- 1 - p[t] * avail[i, 2]
  po[1, i, 3, 2] <- 0
  po[1, i, 3, 3] <- p[t] * avail[i, 2]
  po[1, i, 3, 4] <- 0

  po[1, i, 4, 1] <- 1 - avail[i, 2]
  po[1, i, 4, 2] <- 0
  po[1, i, 4, 3] <- 0
  po[1, i, 4, 4] <- avail[i, 2]

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  for (k in 2:(K - 2)) {
    ##~~ TRANSITION ~~#
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ## Event 1 -> Event 2, ..., Event K - 2 -> Event K - 1
    ps[k, i, 1, 1] <- 1 - gamma[k]
    ps[k, i, 1, 2] <- gamma[k]
    ps[k, i, 1, 3] <- 0
    ps[k, i, 1, 4] <- 0

    ps[k, i, 2, 1] <- 0
    ps[k, i, 2, 2] <- (delta[k - 1, g[i]] * avail[i, k]) + (1 - avail[i, k])
    ps[k, i, 2, 3] <- (1 - delta[k - 1, g[i]]) * (1 - phi[k - 1, g[i]]) * (avail[k, t])
    ps[k, i, 2, 4] <- (1 - delta[k - 1, g[i]]) * (phi[k - 1, g[i]]) * (avail[k, t])

    ps[k, i, 3, 1] <- 0
    ps[k, i, 3, 2] <- 0
    ps[k, i, 3, 3] <- 1
    ps[k, i, 3, 4] <- 0

    ps[k, i, 4, 1] <- 0
    ps[k, i, 4, 2] <- 0
    ps[k, i, 4, 3] <- 0
    ps[k, i, 4, 4] <- 1

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ##~~ OBSERVATION ~~#
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ## Event 2, ..., Event K - 1
    po[k, i, 1, 1] <- 1
    po[k, i, 1, 2] <- 0
    po[k, i, 1, 3] <- 0
    po[k, i, 1, 4] <- 0

    po[k, i, 2, 1] <- 1 - p[t] * avail[i, k + 1]
    po[k, i, 2, 2] <- p[t] * avail[i, k + 1]
    po[k, i, 2, 3] <- 0
    po[k, i, 2, 4] <- 0

    po[k, i, 3, 1] <- 1 - p[t] * avail[i, k + 1]
    po[k, i, 3, 2] <- 0
    po[k, i, 3, 3] <- p[t] * avail[i, k + 1]
    po[k, i, 3, 4] <- 0

    po[k, i, 4, 1] <- 1 - avail[i, k + 1]
    po[k, i, 4, 2] <- 0
    po[k, i, 4, 3] <- 0
    po[k, i, 4, 4] <- avail[i, k + 1]
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  }

  ##~~ TRANSITION ~~#
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## Event K - 1 -> Event K
  ## Turn all fish remaining on spawning grounds into spawners
  ps[K - 1, i, 1, 1] <- 1
  ps[K - 1, i, 1, 2] <- 0
  ps[K - 1, i, 1, 3] <- 0
  ps[K - 1, i, 1, 4] <- 0

  ps[K - 1, i, 2, 1] <- 0
  ps[K - 1, i, 2, 2] <- 1 - avail[i, K - 1]
  ps[K - 1, i, 2, 3] <- (1 - phi[K - 2, g[i]]) * avail[i, K - 1]
  ps[K - 1, i, 2, 4] <- phi[K - 2, g[i]] * avail[i, K - 1]

  ps[K - 1, i, 3, 1] <- 0
  ps[K - 1, i, 3, 2] <- 0
  ps[K - 1, i, 3, 3] <- 1
  ps[K - 1, i, 3, 4] <- 0

  ps[K - 1, i, 4, 1] <- 0
  ps[K - 1, i, 4, 2] <- 0
  ps[K - 1, i, 4, 3] <- 0
  ps[K - 1, i, 4, 4] <- 1

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ##~~ OBSERVATION ~~#
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## Event K
  ## Only trap observations after last tangle-net event.
  po[K - 1, i, 1, 1] <- 1
  po[K - 1, i, 1, 2] <- 0
  po[K - 1, i, 1, 3] <- 0
  po[K - 1, i, 1, 4] <- 0

  po[K - 1, i, 2, 1] <- 1
  po[K - 1, i, 2, 2] <- 0
  po[K - 1, i, 2, 3] <- 0
  po[K - 1, i, 2, 4] <- 0

  po[K - 1, i, 3, 1] <- 1
  po[K - 1, i, 3, 2] <- 0
  po[K - 1, i, 3, 3] <- 0
  po[K - 1, i, 3, 4] <- 0

  po[K - 1, i, 4, 1] <- 1 - avail[i, K]
  po[K - 1, i, 4, 2] <- 0
  po[K - 1, i, 4, 3] <- 0
  po[K - 1, i, 4, 4] <- avail[i, K]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
}
