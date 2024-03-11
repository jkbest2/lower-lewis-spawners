## These functions reshape the simulated data set to conform to the requirements
## of the JAGS model, including combining trajectories into matrices with rows
## for each individual with `sim_tomat` and extracting the required data list
## and initial values from the simulated data. The initial values are provided
## as a function (that always returns the same values) so that multiple chains
## can be used easily. The `sim_check` function is also provided to ensure that
## all data and initial value components are the correct dimensions and have
## missing values where appropriate.
sim_tomat <- function(l) {
  m <- Reduce(rbind, l)
  colnames(m) <- paste("Event", 0:(ncol(m) - 1), sep = "_")
  rownames(m) <- NULL
  m
}

prep_data <- function(res, pars = NULL) {
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
  if (is.null(pars)) {
    init <- function() {
      list(z = z_init,
           origin = origin_init,
           sex = sex_init,
           resid = resid_init)
    }
  }
  else {
    init <- function() {
      list(z = z_init,
           origin = origin_init,
           sex = sex_init,
           resid = resid_init,
           mean.p = qlogis(pars$p),
           mean.gamma = mean(qlogis(pars$gamma)),
           beta_delta = apply(pars$delta, 2, \(c) mean(qlogis(c))),
           beta_phi = apply(pars$phi, 2, \(c) mean(qlogis(c))),
           beta = pars$beta,
           sex_ratio = pars$sex_ratio,
           resid_ratio = pars$resid_ratio)
    }
  }
  list(init = init, data = data)
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
