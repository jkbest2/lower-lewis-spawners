
###Model derived from Chapter 10 of Kery and Schaub "Bayesian Population Analysis Using WinBUGS"

###Function to set initial values for ecological state matrix for unobserved fish
###(Avoids incosistent parent nodes)

pHOS.ms.init.z <-function(ch){

  #return ch.init if all ch are 1 (totally unobserved fish)
  # ~ 70% should stay unentered
  all_ones_init <- function(ch.row){
    ch.init.row = rep(NA, length(ch.row))

    if(runif(1) < 0.7){
     ch.init.row = rep(1, length(ch.row))
    } else {
      if (runif(1) < 0.5) {
        ch.init.row = c(1, rep(2, length(ch.row) - 2), 3)
      } else {
        ripe = sample(2:length(ch.row), 1)
        ch.init.row[1] = 1
        ch.init.row[2:(ripe-1)] = 2
        ch.init.row[3:length(ch.row)] = 3
      }
    }
    return(ch.init.row)
  }
  
  #return ch.init if fish is captured in returning state, but not
  #as a kelt. because we don't have availability information which
  #would limit transitions, just assume fish remains in remaining state at all time
  #periods but the last
  no_kelt_init  <- function(ch.row){
    ch.init.row = ch.row
    first = min(which(as.numeric(ch[i,]) == 2))
    
    for(i in first:(length(ch.row)-1)){
      ch.init.row[i] = ifelse(is.na(ch.row[i]), 2, ch.row[i])
    }
    
    if (first > 2) {
      enter = sample(2:(first-1), 1)
      ch.init.row[2:enter] = 1
      ch.init.row[enter:(first-1)] = 2
    }
      
    return(ch.init.row)
  }
  
  #return ch.init if fish is captured as a kelt
  kelt_init     <- function(ch.row){
    ch.init.row = ch.row
    first_kelt = min(which(as.numeric(ch.row) == 3))
    if (any(na.omit(ch.row) == 2)){
      first = min(which(as.numeric(ch.row) == 2))
        for(i in first:(first_kelt-1)){
          ch.init.row[i] = ifelse(is.na(ch.row[i]), 2, ch.row[i])
        }
      if (first > 2) {
        enter = sample(2:(first-1), 1)
        ch.init.row[2:enter] = 1
        ch.init.row[enter:(first-1)] = 2
        }
    } else {
      enter = sample(2:(first_kelt-1), 1)
      ch.init.row[2:enter] = 1
      ch.init.row[enter:(first_kelt-1)] = 2
    }
    return(ch.init.row)
    
  }

  ch.init <- matrix(rep(NA, length(ch)), nrow=dim(ch)[1])

  for(i in 1:dim(ch)[1]){
    # first check if there are any unknown states
    if( any(is.na(ch[i,])) ){
      # second check if fish is recaptured at all
      if( all(na.omit(ch[i,]) == 1) ) {
        ch.init[i,]  = all_ones_init(ch[i,])
        #third check if fish is ever identified as a kelt
      } else if ( any(na.omit(ch[i,2:(dim(ch)[2]-1)])  == 3) ) {
        ch.init[i,]  = kelt_init(ch[i,])
      } else{
        ch.init[i,]  = no_kelt_init(ch[i,])  
      }
    } else {
      ch.init[i,] = ch[i,]
    }
    for (j in 1:dim(ch)[2]){
      if( !is.na(ch[i,j]) ) ch.init[i,j] = NA 
    }
  }
  return(ch.init)
}

##Function to set initial values for origin of unobserved fish

pHOS.or.init <- function(or){
  or.init <- rep(NA, length(or))
  for(i in 1:length(or)){
    if(is.na(or[i])){
      or.init[i] <- rbinom(size=1, p = .5, n = 1)
    }
  }
  return(or.init)
}

##Multistate model stored as text string

pHOS_Multistate = "
model {
  #-----------------------------------------------
  #~~~~DATA:
  # M: number of fish in augmented data
  # K: number of capture events
  # origin: [M] vector taking values 0 or 1
  # sex: [M] vector taking values 0 or 1
  # resid: [M] vector taking values 0 or 1
  # z: [M x K] matrix taking values NA, 1, 2, 3, or 4
  # y: [M x K] matrix taking values 1, 2, 3, or 4
  # avail: [M x K] matrix taking values 0 or 1

  #-----------------------------------------------
  #####STATES:
  # Ecological States: 1 - Potential Population
  #                    2 - Returning Fish
  #                    3 - Kelt/Spawner
  #                    4 - Trap
  
  # Observed States:   1 - Not Observed
  #                    2 - Tangle Net trapped
  #                    3 - Tangle Net trapped Kelt
  #                    4 - Trapped
  
  #####PARAMETERS:
  #beta: origin ratio
  #sex_ratio: male-female ratio
  #resid_ratio: residual ratio
  #gamma: entry probability
  #delta: probability of becoming a spawner
  #phi: probability of moving to trap
  #p: tangle-net capture probability 
  
  ####PRIORS AND CONSTRAINTS:
  
  mean.p  ~ dnorm(-1.5, 1)        #prior on tangle net detetion efficiency which we expect to be low
  sigma.p ~ dnorm(0, 0.5) T(0,) #prior on standard deviation for detection efficieny random effect (per event)
  tau.p   <- pow(sigma.p, -2)
  
  mean.gamma  ~ dt(0, .1, 7)
  sigma.gamma ~ dnorm(0, 0.5) T(0,) #standard deviation for probability of survival random effect (release group)
  tau.gamma   <- pow(sigma.gamma, -2)
  
  beta_delta[1] ~ dt(0, .1, 7)
  beta_delta[2] ~ dt(0, .25, 7)
  beta_delta[3] ~ dt(0, .25, 7)
  beta_delta[4] ~ dt(0, .25, 7)
  
  sigma.delta ~ dnorm(0, 0.5) T(0,) #standard deviation for probability of survival random effect (release group)
  tau.delta   <- pow(sigma.delta, -2)
  
  beta_phi[1] ~ dt(0, .1, 7)
  beta_phi[2] ~ dt(0, .25, 7)
  beta_phi[3] ~ dt(0, .25, 7)
  beta_phi[4] ~ dt(0, .25, 7)
  sigma.phi   ~ dnorm(0, .5) T(0,) #standard deviation for probability of survival random effect (release group)
  tau.phi     <- pow(sigma.phi, -2)
  
  beta           ~ dbeta(1, 1)
  sex_ratio      ~ dbeta(5, 5)
  resid_ratio[1] ~ dbeta(1, 9) # % of male population that residualize
  resid_ratio[2] ~ dbeta(1, 9) # % of female population that residualize
  
  ## gamma: probability of entering from potential population
  ## Event 0 through Event K - 2 
  for (t in 1:(K - 2)) {
    e.gamma[t] ~ dnorm(0, tau.gamma)
  
    logit(gamma[t]) <- mean.gamma + e.gamma[t]
  }
  
  ## delta: probability of remaining in returning state
  ## Event 1 through Event K - 2
  for (t in 1:(K - 3)) {
    e.delta[1, t] ~ dnorm(0, tau.delta)
    e.delta[2, t] ~ dnorm(0, tau.delta)
    e.delta[3, t] ~ dnorm(0, tau.delta)
    e.delta[4, t] ~ dnorm(0, tau.delta)
  
    logit(delta[t, 1]) <- beta_delta[1] + e.delta[1, t] 
    # different probability of transitioning to kelt for hatchery origin anadramous fish
    logit(delta[t, 2]) <- beta_delta[1] + beta_delta[2] + e.delta[2, t] 
    # different probability of transitioning to spawner for natural origin anadramous fish
    logit(delta[t, 3]) <- beta_delta[1] + beta_delta[3] + e.delta[3, t] 
    # different probability of transitioning to spawner for hatchery origin residual fish
    logit(delta[t, 4]) <- beta_delta[1] + beta_delta[2] + beta_delta[3] + beta_delta[4] + e.delta[4, t] 
    # different probability of transitioning to spawner for natural origin residual fish
  }
  
  #phi: probability of moving to trap
  ## Event 1 through Event K - 1
  for (t in 1:(K - 2)) {
    e.phi[1, t] ~ dnorm(0, tau.phi)
    e.phi[2, t] ~ dnorm(0, tau.phi)
    e.phi[3, t] ~ dnorm(0, tau.phi)
    e.phi[4, t] ~ dnorm(0, tau.phi)
  
    logit(phi[t, 1]) <- beta_phi[1] + e.phi[1, t]  
    # probability of transitioning to the trap for hatchery origin anadramous fish
    logit(phi[t, 2]) <- beta_phi[1] + beta_phi[2] + e.phi[2, t]  
    # probability of transitioning to the trap for natural origin anadramous fish
    logit(phi[t, 3]) <- beta_phi[1] + beta_phi[3] + e.phi[3, t] 
    # different probability of transitioning to the trap for hatchery origin residual fish
    logit(phi[t, 4]) <- beta_phi[1] + beta_phi[2] + beta_phi[3] + beta_phi[4] + e.phi[4, t]  
    #probability of transitioning to the trap for natural origin residual fish
  }
  
  #p: detection probability
  #Event 1 through Event K - 1
  for (t in 1:(K - 2)) {
    e.p[t] ~ dnorm(0, tau.p)
  
    logit(p[t]) <- mean.p + e.p[t] #this means that detection efficiency is the same for each tangle net recapture
  }
  
  
  ####TRANSITION & OBSERVATION MATRICES
  
  for (i in 1:M){
    origin[i] ~ dbern(beta) # 0 - HOR, 1 - NOR
    sex[i]    ~ dbern(sex_ratio)  # 0 - M, 1 - F
    resid[i]  ~ dbern(resid_ratio[1 + sex[i]]) # 0 - anadramous, 1 - residual
  
    g[i] <- 1 + 1 * origin[i] + 2 * resid[i]  
  
    ## Initial entry to Lewis River (all fish start as potential population)
    ##TRANSITION ## Event 0 -> Event 1
    ###########################
    #[time, individual, current state, next state]
  
    ps[1, i, 1, 1] <- 1 - gamma[1]
    ps[1, i, 1, 2] <- gamma[1]
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
    ##########################
  
    ##OBSERVATION## Event 1
    ###########################
    #[time, individual, actual state, observed state]
  
    po[1, i, 1, 1] <- 1
    po[1, i, 1, 2] <- 0
    po[1, i, 1, 3] <- 0
    po[1, i, 1, 4] <- 0
  
    po[1, i, 2, 1] <- 1 - p[1] * avail[i, 2]
    po[1, i, 2, 2] <- p[1] * avail[i, 2]
    po[1, i, 2, 3] <- 0
    po[1, i, 2, 4] <- 0
  
    po[1, i, 3, 1] <- 1 - p[1] * avail[i, 2]
    po[1, i, 3, 2] <- 0
    po[1, i, 3, 3] <- p[1] * avail[i, 2]
    po[1, i, 3, 4] <- 0
  
    po[1, i, 4, 1] <- 1 - avail[i, 2]
    po[1, i, 4, 2] <- 0
    po[1, i, 4, 3] <- 0
    po[1, i, 4, 4] <- avail[i, 2]
    ###########################
  
    for (k in 2:(K - 2)) { 
      ##TRANSITION ## Event 1 -> Event 2, ..., Event K - 2 -> Event K - 1
      ###########################
      ps[k, i, 1, 1] <- 1 - gamma[k]
      ps[k, i, 1, 2] <- gamma[k]
      ps[k, i, 1, 3] <- 0
      ps[k, i, 1, 4] <- 0
  
      ps[k, i, 2, 1] <- 0
      ps[k, i, 2, 2] <- (delta[k - 1, g[i]] * avail[i, k]) + (1 - avail[i, k])
      ps[k, i, 2, 3] <- (1 - delta[k - 1, g[i]]) * (1 - phi[k - 1, g[i]]) * (avail[i, k])
      ps[k, i, 2, 4] <- (1 - delta[k - 1, g[i]]) * (phi[k - 1, g[i]]) * (avail[i, k])
  
      ps[k, i, 3, 1] <- 0
      ps[k, i, 3, 2] <- 0
      ps[k, i, 3, 3] <- 1
      ps[k, i, 3, 4] <- 0
  
      ps[k, i, 4, 1] <- 0
      ps[k, i, 4, 2] <- 0
      ps[k, i, 4, 3] <- 0
      ps[k, i, 4, 4] <- 1
      ###########################
  
      ##OBSERVATION## Event 2, ..., Event K - 1
      ###########################
      po[k, i, 1, 1] <- 1
      po[k, i, 1, 2] <- 0
      po[k, i, 1, 3] <- 0
      po[k, i, 1, 4] <- 0
  
      po[k, i, 2, 1] <- 1 - p[k] * avail[i, k + 1]
      po[k, i, 2, 2] <- p[k] * avail[i, k + 1]
      po[k, i, 2, 3] <- 0
      po[k, i, 2, 4] <- 0
  
      po[k, i, 3, 1] <- 1 - p[k] * avail[i, k + 1]
      po[k, i, 3, 2] <- 0
      po[k, i, 3, 3] <- p[k] * avail[i, k + 1]
      po[k, i, 3, 4] <- 0
  
      po[k, i, 4, 1] <- 1 - avail[i, k + 1]
      po[k, i, 4, 2] <- 0
      po[k, i, 4, 3] <- 0
      po[k, i, 4, 4] <- avail[i, k + 1]
      ###########################
    }
  
    ##TRANSITION ## Event K - 1 -> Event K
    ###########################
    # Turn all fish remaining on spawning grounds into spawners
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
    ###########################
  
    ##OBSERVATION## Event K
    ###########################
    # No tangle-netting, only trap observations after last tangle-netting event.
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
    ###########################
  }
  
  ####LIKELIHOOD:
  for (i in 1:M){
    for (k in 2:K){
      z[i, k] ~ dcat(ps[k - 1, i, z[i, k - 1], ])
      y[i, k] ~ dcat(po[k - 1, i, z[i, k], ])
    }
  }
  
  ###DERIVED PARAMETERS:
  for (i in 1:M){
    NOR_entered[i] <- (1 - equals(z[i,K], 1)) * equals(origin[i], 1)
    HOR_entered[i] <- (1 - equals(z[i,K], 1)) * equals(origin[i], 0)
    NOR_spawn[i]   <- equals(z[i,K], 3) * equals(origin[i], 1)
    HOR_spawn[i]   <- equals(z[i,K], 3) * equals(origin[i], 0)
  }
  
  N_NOR_entered <- sum(NOR_entered[])
  N_HOR_entered <- sum(HOR_entered[])
  
  N_NOR_spawn <- sum(NOR_spawn[])
  N_HOR_spawn <- sum(HOR_spawn[])
  
  pHOS <- N_HOR_spawn/(N_HOR_spawn + N_NOR_spawn)
}
"

## Function to generate posterior predictive check
## Operates on a per row basis
post_pred_pHOS <- function(N, draws_row){
  params <- draws_row %>% 
    gather(Parameter, Value)
  
  sex_ratio <- params %>% 
    filter(Parameter == "sex_ratio") %>% 
    pull()
  
  beta <- params %>% 
    filter(Parameter == "beta") %>% 
    pull()
  
  resid_ratio <- params %>% 
    filter(grepl("resid_ratio", Parameter)) %>% 
    pull()
  
  gamma <-  params %>% 
    filter(grepl("gamma", Parameter)) %>% 
    pull()
  
  p <-  params %>% 
    filter(Parameter %in% c("p[1]", "p[2]", "p[3]", "p[4]", "p[5]", 
                            "p[6]", "p[7]", "p[8]", "p[9]", "p[10]")) %>% 
    pull()
  
  delta <- matrix(ncol = 4, params %>% 
                    filter(grepl("delta", Parameter)) %>% 
                    pull())
  
  phi <- matrix(ncol = 4, params %>% 
                  filter(grepl("phi", Parameter)) %>% 
                  pull())
  
  trans_array <- array(NA, c(10, 4, 3))
  
  trans_array[, , 1] <- rbind(delta, c(0,0,0,0)) # 2 -> 2
  trans_array[, , 2] <- rbind((1 - delta), c(1, 1, 1,1)) * (1 - phi) # 2 -> 3
  trans_array[, , 3] <- rbind((1 - delta), c(1, 1, 1,1)) * (phi) # 2 -> 4
  
  eco_sim <- data_frame(origin_state = rbinom(n = N, size = 1, prob = beta), #0 - HOR, 1 - NOR,
                        sex_state    = rbinom(n = N, size = 1, prob = sex_ratio), #0 - M, 1 - F,
                        resid_state  = if_else(sex_state == 0,  
                                               rbinom(N, size = 1, prob = resid_ratio[1]),
                                               rbinom(N, size = 1, prob = resid_ratio[2])
                        ) ) %>% 
    mutate(state_trait = 1 + 1 * origin_state + 2 * resid_state ,
           Event_01    =  1 + rbinom(N, size = 1, prob = gamma[1]),
           Event_02    =  
             case_when(
               Event_01 == 1                    ~ Event_01 + rbinom(N, size = 1, prob = gamma[2]),
               Event_01 == 2 & state_trait == 1 ~ Event_01 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[1, 1, ]),
               Event_01 == 2 & state_trait == 2 ~ Event_01 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[1, 2, ]),
               Event_01 == 2 & state_trait == 3 ~ Event_01 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[1, 3, ]),
               Event_01 == 2 & state_trait == 4 ~ Event_01 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[1, 4, ])),
           Event_03    =  
             case_when(
               Event_02 == 1                    ~ Event_02 + rbinom(N, size = 1, prob = gamma[3]),
               Event_02 == 2 & state_trait == 1 ~ Event_02 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[2, 1, ]),
               Event_02 == 2 & state_trait == 2 ~ Event_02 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[2, 2, ]),
               Event_02 == 2 & state_trait == 3 ~ Event_02 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[2, 3, ]),
               Event_02 == 2 & state_trait == 4 ~ Event_02 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[2, 4, ]),
               Event_02 == 3                    ~ Event_02,
               Event_02 == 4                    ~ Event_02),
           
           Event_04    =  
             case_when(
               Event_03 == 1                    ~ Event_03 + rbinom(N, size = 1, prob = gamma[4]),
               Event_03 == 2 & state_trait == 1 ~ Event_03 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[3, 1, ]),
               Event_03 == 2 & state_trait == 2 ~ Event_03 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[3, 2, ]),
               Event_03 == 2 & state_trait == 3 ~ Event_03 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[3, 3, ]),
               Event_03 == 2 & state_trait == 4 ~ Event_03 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[3, 4, ]),
               Event_03 == 3                    ~ Event_03,
               Event_03 == 4                    ~ Event_03),
           
           Event_05    =  
             case_when(
               Event_04 == 1                    ~ Event_04 + rbinom(N, size = 1, prob = gamma[5]),
               Event_04 == 2 & state_trait == 1 ~ Event_04 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[4, 1, ]),
               Event_04 == 2 & state_trait == 2 ~ Event_04 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[4, 2, ]),
               Event_04 == 2 & state_trait == 3 ~ Event_04 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[4, 3, ]),
               Event_04 == 2 & state_trait == 4 ~ Event_04 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[4, 4, ]),
               Event_04 == 3                    ~ Event_04,
               Event_04 == 4                    ~ Event_04),
           
           Event_06    =  
             case_when(
               Event_05 == 1                    ~ Event_05 + rbinom(N, size = 1, prob = gamma[6]),
               Event_05 == 2 & state_trait == 1 ~ Event_05 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[5, 1, ]),
               Event_05 == 2 & state_trait == 2 ~ Event_05 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[5, 2, ]),
               Event_05 == 2 & state_trait == 3 ~ Event_05 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[5, 3, ]),
               Event_05 == 2 & state_trait == 4 ~ Event_05 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[5, 4, ]),
               Event_05 == 3                    ~ Event_05,
               Event_05 == 4                    ~ Event_05),
           
           Event_07    =  
             case_when(
               Event_06 == 1                    ~ Event_06 + rbinom(N, size = 1, prob = gamma[7]),
               Event_06 == 2 & state_trait == 1 ~ Event_06 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[6, 1, ]),
               Event_06 == 2 & state_trait == 2 ~ Event_06 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[6, 2, ]),
               Event_06 == 2 & state_trait == 3 ~ Event_06 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[6, 3, ]),
               Event_06 == 2 & state_trait == 4 ~ Event_06 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[6, 4, ]),
               Event_06 == 3                    ~ Event_06,
               Event_06 == 4                    ~ Event_06),
           
           Event_08    =  
             case_when(
               Event_07 == 1                    ~ Event_07 + rbinom(N, size = 1, prob = gamma[8]),
               Event_07 == 2 & state_trait == 1 ~ Event_07 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[7, 1, ]),
               Event_07 == 2 & state_trait == 2 ~ Event_07 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[7, 2, ]),
               Event_07 == 2 & state_trait == 3 ~ Event_07 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[7, 3, ]),
               Event_07 == 2 & state_trait == 4 ~ Event_07 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[7, 4, ]),
               Event_07 == 3                    ~ Event_07,
               Event_07 == 4                    ~ Event_07),
           
           Event_09    =  
             case_when(
               Event_08 == 1                    ~ Event_08 + rbinom(N, size = 1, prob = gamma[9]),
               Event_08 == 2 & state_trait == 1 ~ Event_08 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[8, 1, ]),
               Event_08 == 2 & state_trait == 2 ~ Event_08 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[8, 2, ]),
               Event_08 == 2 & state_trait == 3 ~ Event_08 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[8, 3, ]),
               Event_08 == 2 & state_trait == 4 ~ Event_08 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[8, 4, ]),
               Event_08 == 3                    ~ Event_08,
               Event_08 == 4                    ~ Event_08),
           
           Event_10    =  
             case_when(
               Event_09 == 1                    ~ Event_09 + rbinom(N, size = 1, prob = gamma[10]),
               Event_09 == 2 & state_trait == 1 ~ Event_09 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[9, 1, ]),
               Event_09 == 2 & state_trait == 2 ~ Event_09 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[9, 2, ]),
               Event_09 == 2 & state_trait == 3 ~ Event_09 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[9, 3, ]),
               Event_09 == 2 & state_trait == 4 ~ Event_09 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[9, 4, ]),
               Event_09 == 3                    ~ Event_09,
               Event_09 == 4                    ~ Event_09),
           
           Event_11    =  
             case_when(
               Event_10 == 1                    ~ Event_10,
               Event_10 == 2 & state_trait == 1 ~ Event_09 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[10, 1, ]),
               Event_10 == 2 & state_trait == 2 ~ Event_09 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[10, 2, ]),
               Event_10 == 2 & state_trait == 3 ~ Event_09 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[10, 3, ]),
               Event_10 == 2 & state_trait == 4 ~ Event_09 + sample(c(0, 1, 2), size = N, replace = T, 
                                                                    prob = trans_array[10, 4, ]),
               Event_10 == 3                    ~ Event_10,
               Event_10 == 4                    ~ Event_10)) %>% 
    mutate(Event_11 = if_else(Event_10 == 4, 0, Event_11),
           Event_10 = if_else(Event_09 == 4, 0, Event_10),
           Event_09 = if_else(Event_08 == 4, 0, Event_09),
           Event_08 = if_else(Event_07 == 4, 0, Event_08),
           Event_07 = if_else(Event_06 == 4, 0, Event_07),
           Event_06 = if_else(Event_05 == 4, 0, Event_06),
           Event_05 = if_else(Event_04 == 4, 0, Event_05),
           Event_04 = if_else(Event_03 == 4, 0, Event_04),
           Event_03 = if_else(Event_02 == 4, 0, Event_03),
           Event_02 = if_else(Event_01 == 4, 0, Event_02))
  
  obs_sim <- eco_sim %>% 
    mutate(Event_01 = case_when(Event_01 %in% c(0,1) ~ 0,
                                Event_01 %in% c(2,3) ~ Event_01 * rbinom(N, 1, p[1]),
                                Event_01 == 4        ~ 4),
           Event_02 = case_when(Event_02 %in% c(0,1) ~ 0,
                                Event_02 %in% c(2,3) ~ Event_02 * rbinom(N, 1, p[2]),
                                Event_02 == 4        ~ 4),
           Event_03 = case_when(Event_03 %in% c(0,1) ~0,
                                Event_03 %in% c(2,3) ~ Event_03 * rbinom(N, 1, p[3]),
                                Event_03 == 4        ~ 4),
           Event_04 = case_when(Event_04 %in% c(0,1) ~ 0,
                                Event_04 %in% c(2,3) ~ Event_04 * rbinom(N, 1, p[4]),
                                Event_04 == 4        ~ 4),
           Event_05 = case_when(Event_05 %in% c(0,1) ~ 0,
                                Event_05 %in% c(2,3) ~ Event_05 * rbinom(N, 1, p[5]),
                                Event_05 == 4        ~ 4),
           Event_06 = case_when(Event_06 %in% c(0,1) ~ 0,
                                Event_06 %in% c(2,3) ~ Event_06 * rbinom(N, 1, p[6]),
                                Event_06 == 4        ~ 4),
           Event_07 = case_when(Event_07 %in% c(0,1) ~ 0,
                                Event_07 %in% c(2,3) ~ Event_07 * rbinom(N, 1, p[7]),
                                Event_07 == 4        ~ 4),
           Event_08 = case_when(Event_08 %in% c(0,1) ~ 0,
                                Event_08 %in% c(2,3) ~ Event_08 * rbinom(N, 1, p[8]),
                                Event_08 == 4        ~ 4),
           Event_09 = case_when(Event_09 %in% c(0,1) ~ 0,
                                Event_09 %in% c(2,3) ~ Event_09 * rbinom(N, 1, p[9]),
                                Event_09 == 4        ~ 4),
           Event_10 = case_when(Event_10 %in% c(0,1) ~ 0,
                                Event_10 %in% c(2,3) ~ Event_10 * rbinom(N, 1, p[10]),
                                Event_10 == 4        ~ 4),
           Event_11 = if_else(Event_11 == 4, 4, 0))
  
  return(obs_sim %>% 
           gather(-origin_state, -sex_state, -resid_state, -state_trait,
                  key = event, value = obs_state) %>% 
           group_by(event, origin_state, sex_state, resid_state) %>% 
           summarise(n_return = sum(obs_state == 2),
                     n_spawn  = sum(obs_state == 3),
                     n_trap   = sum(obs_state == 4)) %>%  
           ungroup())
}

