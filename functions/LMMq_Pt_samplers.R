
# x sampler ----
nimSamplerXmove <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    S <- dim(model$latentObservation)[2]
    n <- length(model$latentIndex)
    D <- control$D
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    
    # ----------------------------------------------- Intro ----
    
    latentObsProp <- model$latentObservation
    latentIdxProp <- model$latentIndex
    tauProp <- model$tau
    
    log_MH_ratio <- 0
    skip <- TRUE
    
    nu0 <- nimNumeric(S, value = 0.0)
    nu2 <- nimNumeric(S, value = 0.0)
    
    # ----
    for(nb_modif in 1:D){
    
    
    choice <- runif(1, 0, 1)
    # ----------------------------------------------- add error ----
    if(choice <= 0.5){
      # What time is available for adding an error
      available_fully <- getAvailable0(latentObsProp, latentIdxProp)
      if(any(available_fully)){
        skip <- FALSE
        # sample t, nu2, nu0 and nu1 ----
        t <- nimSample(available_fully)
        existing <- latentIdxProp > 0
        has0 <- latentObsProp[, t] == 0
        which_nu0 <- nimNumeric()
        which_nu0 <- which((has0 * existing) > 0)
        nu0i <- nimSample(which_nu0)
        nu0 <- latentObsProp[nu0i, ]
        all_nu1i <- nimNumeric()
        all_nu1i <- which(latentIdxProp == 1+3^(t-1))
        als <- 1 - pnorm(model$a[1] * tauProp[all_nu1i, t] + model$b[1])
        r_nu1i <- rcat(1, als/sum(als))
        nu1i <- all_nu1i[r_nu1i]
        nu2 <- nu0
        nu2[t] <- 2
        
        # update nu0 as nu2 ----
        latentIdxProp[nu0i] <- nimGetLatentIndex(nu2)
        latentObsProp[nu0i, t] <- 2
        tauProp[nu0i, t] <- tauProp[nu1i, t]
        # put nu1 as unexisting ----
        latentIdxProp[nu1i] <- 0
        latentObsProp[nu1i, t] <- 0
        tauProp[nu1i, t] <- 0
        
        # q(x|xprop) / q(xprop|x) ----
        av2 <- getAvailable2(latentObsProp)
        all_nu2i <- nimNumeric()
        all_nu2i <- which(latentObsProp[, t] == 2)
        als2 <- pnorm(model$a[1] * tauProp[all_nu2i, t] + model$b[1])
        # [X | X'] / [X' | X]
        log_qx <- log(1/length(av2)) + log(1/sum(als2))
        log_qxprop <- log(1/length(available_fully)) + 
                      log(1/sum(als)) + log(1/length(which_nu0))
        log_qratio <- log_qx - log_qxprop
        # [Y | X', Z', N', ...] / [Y | X, Z, N, ...]
        N <- sum(model$latentIndex > 0)
        log_Lratio <- - log(N) - sum(log(1-model$capture))
        # MH ratio
        log_MH_ratio <- log_MH_ratio + log_Lratio + log_qratio
      }
    }
    
    # ----------------------------------------------- remove error ----
    else {
      
      # What time is available for adding an error
      available_2 <- getAvailable2(latentObsProp)
      if(any(available_2)){
        skip <- FALSE
        # sample t, nu0, nu1 and nu2 ----
        t <- nimSample(available_2)
        all_nu2i <- nimNumeric()
        all_nu2i <- which(latentObsProp[, t] == 2)
        als <- pnorm(model$a[1] * tauProp[all_nu2i, t] + model$b[1])
        r_nu2i <- rcat(1, als/sum(als))
        nu2i <- all_nu2i[r_nu2i]
        nu2 <- latentObsProp[nu2i, ]
        nu0 <- nu2
        nu0[t] <- 0
        nu1i <- which(latentIdxProp == 0)[1]
        
        # put nu1 as seen ----
        latentIdxProp[nu1i] <- 1+3^(t-1)
        latentObsProp[nu1i, t] <- 1
        tauProp[nu1i, t] <- tauProp[nu2i, t]
        # update nu2 as nu0 ----
        latentIdxProp[nu2i] <- nimGetLatentIndex(nu0)
        latentObsProp[nu2i, t] <- 0
        tauProp[nu2i, t] <- 0
        
        # q(x|xprop) / q(xprop|x) ----
        av0 <- getAvailable0(latentObsProp, latentIdxProp)
        existing <- latentIdxProp > 0
        has0 <- latentObsProp[, t] == 0
        which_nu0 <- nimNumeric()
        which_nu0 <- which((has0 * existing) > 0)
        possib0 <- length(which_nu0)
        all_nu1i <- nimNumeric()
        all_nu1i <- which(latentIdxProp == 1+3^(t-1))
        als1 <- 1 - pnorm(model$a[1] * tauProp[all_nu1i, t] + model$b[1])
        # [X | X'] / [X' | X]
        log_qx <- log(1/length(av0)) + 
          log(1/sum(als1)) + log(1/length(which_nu0))
        log_qxprop <- log(1/length(available_2)) + log(1/sum(als))
        log_qratio <- log_qx - log_qxprop
        # [Y | X', Z', N', ...] / [Y | X, Z, N, ...]
        Nprop <- sum(latentIdxProp > 0)
        log_Lratio <- log(Nprop) + sum(log(1-model$capture))
        # MH ratio
        log_MH_ratio <- log_MH_ratio + log_Lratio + log_qratio
      }
    }
    
    }
    # ------------------------------ metrolis ratio and acceptance ----
    
    if(!skip){
      
      u <- runif(1, 0, 1)
      if(u < exp(log_MH_ratio)){
        model[["N"]] <<- sum(latentIdxProp >  0)
        model[["D"]] <<- sum(latentIdxProp == 1)
        model[["latentObservation"]] <<- latentObsProp
        model[["latentIndex"]] <<- latentIdxProp
        model[["tau"]] <<- tauProp
        
        model[["nbErr"]] <<- sum(latentObsProp == 2)
        
        copy(from = model, to = mvSaved, row = 1,
             nodes = c(calcNodes, "N", "D", "tau", "latentIndex", "nbErr"), logProb = T)
      }
    }
  },
  methods = list(
    getAvailable0 = function(latObsProp = double(2), 
                             latIdxProp = double(1)){
      existing <- latIdxProp > 0
      available_t <- logical(S)
      available_nu0 <- logical(S)
      for(t in 1:S){
        if(any(latIdxProp == 1+3^(t-1))) 
          available_t[t] <- TRUE
        if(any( latObsProp[existing, t] == 0 )) 
          available_nu0[t] <- TRUE
      }
      available_fully <- nimNumeric()
      available_fully <- which(available_t*available_nu0 == 1)
      return(available_fully)
      returnType(double(1))
    },
    getAvailable2 = function(latentObsProp = double(2)){
      available_t <- logical(S)
      for(t in 1:S){
        if(any( latentObsProp[, t] == 2 )) 
          available_t[t] <- TRUE
      }
      available_fully <-nimNumeric()
      available_fully <- which(available_t)
      return(available_fully)
      returnType(double(1))
    },
    reset = function () {})
)



# Sampler unseen indiv ----

nimSamplerX0 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    D <- if(!is.null(control$D))         control$D         else 10
    S <- dim(model$latentObservation)[2]
  },
  run = function() {
    
    x0move <- nimSample(c(-D:(-1), 1:D))
    if(model$D[1] + x0move >= 0){
      
      N <- sum(model$latentIndex > 0)
      Nprop <- N+x0move
      Dprop <- model$D[1] + x0move
      lprobHist <- sum(log(1-model$capture))
      # nb unseen binomial  N'! / n0'!(N'-n0')! pSeen^(N'-n0') pUnseen^n0'
      #                     N!  / n0!(N-n0)!    pSeen^(N-n0)   pUnseen^n0
      #               = (N'! / N) (n0! / n0'!) pUnseen^(n0'-n0)
      log_MH_ratio <- lfactorial(Nprop) - lfactorial(N) +
        lfactorial(model$D[1]) - lfactorial(Dprop) +
        x0move * lprobHist
      
      u <- runif(1, 0, 1)
      if(u < exp(log_MH_ratio)){
        latIdxProp <- model$latentIndex
        idx <- nimNumeric()
        nbMove <- abs(x0move)
        if(x0move > 0){
          idx <- which(latIdxProp == 0)[1:nbMove]
          latIdxProp[idx] <- 1
        }
        else{
          idx <- which(latIdxProp == 1)[1:nbMove]
          latIdxProp[idx] <- 0
        }
        
        model[["latentIndex"]] <<- latIdxProp
        model[["D"]] <<- Dprop
        model[["N"]] <<- Nprop
        
        copy(from = model, to = mvSaved, row = 1,
             nodes = c("latentIndex", "N", "D"), logProb = T)
      }
    }
  },
  methods = list(
    reset = function () {} )
)

# Pt sampler ----
nimUpdatePt <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    S <- dim(model$latentObservation)[2]
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    
    apt <- nimNumeric(S, value = 1)
    bpt <- nimNumeric(S, value = 1)
    
    lat_obs_tmp <- model$latentObservation[which(model$latentIndex > 0), ]
    
    for(t in 1:S){
      apt[t] <- sum(lat_obs_tmp[, t] > 0)  # nb de captures total (> 0 dans les histo latents)
      bpt[t] <- sum(lat_obs_tmp[, t] == 0) # nb de non-capture total (0 dans les histo latents)
    }
    
    new_p <- nimNumeric(S)
    for(t in 1:S){
      new_p[t] <- rbeta(1, apt[t], bpt[t])
    }
    
    
    model[[target]] <<- new_p
    copy(from = model, to = mvSaved, row = 1,
         nodes = calcNodes, logProb = T)
  },
  methods = list( reset = function () {} )
)

# a sampler ----
nimUpdatea <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    
    n <- nrow(model$tau)
    S <- ncol(model$tau)
    nb_identification <- sum(model$latentObservation > 0)
    
    prior_mua <- control$prior_mua
    prior_sda <- control$prior_sda
    sda <- sqrt( ( 1/(prior_sda^2) + sum(model$tau^2) )^(-1) )
  },
  run = function() {
    
    tau_u <- get_tau_u()
    
    mua <- sda^2 * (prior_mua/(prior_sda^2) + 
                      sum( tau_u[, 1] * (tau_u[, 2] - model$b[1]) )
    )
    
    new_a <- rnorm(1, mua, sda)
    
    model[[target]] <<- new_a
    copy(from = model, to = mvSaved, row = 1,
         nodes = calcNodes, logProb = T)
  },
  methods = list( 
    get_tau_u = function(){
      res <- nimMatrix(nrow = nb_identification, ncol = 2, value = 0)
      idx <- which(model$latentIndex > 1)
      nb_id <- 1
      for(i in 1:length(idx)){
        I <- idx[i]
        for(t in 1:S){
          if(model$latentObservation[I, t] > 0){
            res[nb_id, 1] <- model$tau[I, t]
            id_u <- model$latentObservation[I, t]
            post_mean_u <- model$a[1] * model$tau[I, t] + model$b[1]
            if(id_u == 1)
              res[nb_id, 2] <- nimRnormTrunc(1, post_mean_u, lower = 0)[1]
            else
              res[nb_id, 2] <- nimRnormTrunc(1, post_mean_u, upper = 0)[1]
            nb_id <- nb_id + 1
          }
        }
      }
      return(res)
      returnType(double(2))
    },
    reset = function () {} )
)

# b sampler ----
nimUpdateb <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    
    n <- nrow(model$tau)
    S <- ncol(model$tau)
    nb_identification <- sum(model$latentObservation > 0)
    
    prior_sdb <- control$prior_sdb
    prior_mub <- control$prior_mub
    sdb <- sqrt( prior_sdb^2 / (nb_identification*prior_sdb^2 + 1) )
  },
  run = function() {
    
    tau_u <- get_tau_u()
    
    mub <- sdb^2 * (prior_mub/(prior_sdb^2) + sum(tau_u[,2] - model$a[1] *  tau_u[,1]))
    
    new_b <- rnorm(1, mub, sdb)
    new_alpha <- pnorm(model$a[1] * tau_u[,1] + new_b)
    
    model[[target]] <<- new_b
    model[["alpha_mean"]] <<- mean(new_alpha)
    copy(from = model, to = mvSaved, row = 1,
         nodes = c(calcNodes, "alpha_mean"), logProb = T)
  },
  methods = list( 
    get_tau_u = function(){
      res <- nimMatrix(nrow = nb_identification, ncol = 2, value = 0)
      idx <- which(model$latentIndex > 1)
      nb_id <- 1
      for(i in 1:length(idx)){
        I <- idx[i]
        for(t in 1:S){
          if(model$latentObservation[I, t] > 0){
            res[nb_id, 1] <- model$tau[I, t]
            id_u <- model$latentObservation[I, t]
            post_mean_u <- model$a[1] * model$tau[I, t] + model$b[1]
            if(id_u == 1)
              res[nb_id, 2] <- nimRnormTrunc(1, post_mean_u, lower = 0)[1]
            else
              res[nb_id, 2] <- nimRnormTrunc(1, post_mean_u, upper = 0)[1]
            nb_id <- nb_id + 1
          }
        }
      }
      return(res)
      returnType(double(2))
    },
    reset = function () {} )
)

# ab sampler ----

nimUpdateBab <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    
    S <- ncol(model$tau)
    nb_identification <- sum(model$latentObservation > 0)
    
    priorMuB <- control$priorMuB
    priorSdB <- control$priorSdB
    
    priorSdBinv <- inverse(priorSdB)
    W <- matrix(rep(1, 2*nb_identification), ncol = 2)
    W[,1] <- model$tau[which(model$tau > 0)]
    tWW <- t(W) %*% W
    
    SdB <- inverse(priorSdBinv + tWW) 
    chSdB <- chol(SdB)
  },
  run = function() {
    
    tau_1_u <- align_tau_1_u()
    MuBtmp <- SdB %*% (priorSdBinv %*% priorMuB + t(tau_1_u[, 1:2]) %*% (tau_1_u[, 3]))
    MuB <- c(MuBtmp[1,1], MuBtmp[2,1])
    ab <- rmnorm_chol(n = 1, mean = MuB, cholesky = chSdB, prec_param = FALSE)
    
    new_alpha <- pnorm(ab[1] * tau_1_u[1:nb_identification, 1] + ab[2])
    
    model[["a"]] <<- ab[1]
    model[["b"]] <<- ab[2]
    model[["alpha_mean"]] <<- mean(new_alpha)
    copy(from = model, to = mvSaved, row = 1,
         nodes = c(calcNodes, "alpha_mean"), logProb = T)
  },
  methods = list( 
    get_u = function(tau = double(0), u = double(0)){
      utilde <- 0
      if(u > 0){
        post_mean_u <- model$a[1] * tau + model$b[1]
        if(u == 1)
          utilde <- nimRnormTrunc(1, post_mean_u, lower = 0)[1]
        else
          utilde <- nimRnormTrunc(1, post_mean_u, upper = 0)[1]
      }
      return(utilde)
      returnType(double(0))
    },
    align_tau_1_u = function(){
      tau_1_u <- nimMatrix(nrow = nb_identification, ncol = 3, value = 0)
      tau_1_u[, 2] <- 1
      
      idx <- which(model$latentIndex > 1)
      nb_id <- 1
      for(i in 1:length(idx)){
        I <- idx[i]
        for(t in 1:S){
          if(model$latentObservation[I, t] > 0){
            tau_1_u[nb_id, 1] <- model$tau[I, t]
            tau_1_u[nb_id, 3] <- get_u(model$tau[I, t], 
                                       model$latentObservation[I, t])
          }
        }
      }
      
      return(tau_1_u)
      returnType(double(2))
    },
    reset = function () {} )
)