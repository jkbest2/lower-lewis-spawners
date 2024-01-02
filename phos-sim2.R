## Number of capture events
K <- 11
## Potential population size
M <- 1200

##~~~~PRIORS AND CONSTRAINTS:
## Detection efficiency
p <- 0.02

## Entry rate; will sample
runsize <- 1200
## Proportion entering before each event (none in the final event; doesn't count
## residualized fish who all "enter" before the first event)
gamma <- pnorm(seq(-2, 2, length.out = 11)) |>
  diff() |>
  c(0)

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
             phi = phi)

sim_model <- function(p, runsize, pars) {
  gamma <- pars$gamma
  beta <- pars$beta
  sex_ratio <- pars$sex_ratio
  resid_ratio <- pars$resid_ratio
  delta <- pars$delta
  phi <- pars$phi

###---------------------------------------------
  ## Individual covariates
  entry <- rmultinom(1, runsize, gamma)
  entry_k <- rep(seq_len(K), entry)
  origin <- rep(NA, runsize)
  sex <- rep(NA, runsize)
  resid <- rep(NA, runsize)
  ## Availability
  avail <- matrix(TRUE, nrow = runsize, ncol = K)
  z <- matrix(NA, nrow = runsize, ncol = K)
  y <- matrix(NA, nrow = runsize, ncol = K)
  for (i in seq_len(runsize)) {
    origin[i] <- rbinom(1, 1, beta) ## 0: HOR, 1: NOR
    sex[i] <- rbinom(1, 1, sex_ratio) ## 0: M, 1: F
    resid[i] <- rbinom(1, 1, resid_ratio[1 + sex[i]]) ## 0: anadramous, 1: residual

    g[i] <- 1 + 1 * origin[i] + 2 * resid[i]

    ## Initial state; available if residual or enters before observation event 1
    avail[i, 1] <- TRUE
    if (resid[i] || entry_k[i] == 1) {
      z[i, 1] <- 2
      y[i, 1] <- sample(1:4, 1, TRUE, c(1 - p, p, 0, 0))
    } else {
      z[i, 1] <- 1 # Still potential; not in the river yet
      y[i, 1] <- 1 # Not observed
    }

    for (k in 2:(K - 1)) {
      ## Update availability
      if (z[i, k - 1] == 4.1) {
        avail[i, k] <- FALSE
      } else {
        avail[i, k] <- TRUE
      }

      ## Update state
      if (entry_k[i] > k) {
        z[i, k] <- 1
      } else if (entry_k[i] == k) {
        z[i, k] <- 2
      } else if (z[i, k - 1] == 2) {
        ps <- c(
          0,
          delta,
          (1 - delta) * (1 - phi),
          (1 - delta) * phi)
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
        po <- c(1 - p, p, 0, 0)
        y[i, k] <- sample(1:4, 1, TRUE, po)
      } else if (z[i, k] == 3) {
        po <- c(1 - p, 0, p, 0)
        y[i, k] <- sample(1:4, 1, TRUE, po)
      } else if (z[i, k] == 4) {
        y[i, k] <- 4
      } else if (z[i, k] == 4.1) {
        ## If state is 4.1, already trapped, it will not be available and will not
        ## be observed
        y[i, k] <- 1
      }
    }

    ## Final transition period; all fish end up as either spawners (3) or in the
    ## trap (4)
    if (z[i, K - 1] == 1) {
      z[i, K] <- 1
    } else if (z[i, K - 1] == 2) {
      ps <- c(0, 0, 1 - phi, phi)
      z[i, K] <- sample(1:4, 1, TRUE, ps)
    } else if (z[i, K - 1] == 4) {
      z[i, K] <- 4.1
      avail[i, K] <- FALSE
    } else if (z[i, K - 1] > 2) {
      z[i, K] <- z[i, K - 1]
    }

    ## Final observation; only the trap is operating so those fish are the only
    ## observations
    if (z[i, K] != 4) {
      y[i, K] <- 1
    } else {
      y[i, K] <- 4
    }
  }

  zfl <- floor(z)

  list(z = zfl,
       y = y,
       avail = avail,
       origin = origin,
       sex = sex,
       resid = resid)
}

res <- sim_model(0.02, 1200, pars)
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
             matrix(NA, nrow = add, ncol = ncol(z)))
  y <- rbind(res$y[which_obs, ],
             matrix(1, nrow = add, ncol = ncol(y)))
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
