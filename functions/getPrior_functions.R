estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

nimRnormTrunc <- nimbleFunction(
  run = function(n = double(0), 
                 mu=double(0, default = 0),
                 sd=double(0, default = 1), 
                 lower=double(0, default=-Inf), 
                 upper=double(0, default=Inf)){
    
    # inverse transformation method
    
    # Generate 10,000 uniform distributed sample from [0, 1]
    q_simulation = runif(n, min=0, max=1)
    
    # Put in back to inverse distribution
    q_temp_simulation <- (pnorm(lower, mean = mu, sd=sd ) + 
                            q_simulation * (pnorm(upper, mean = mu, sd=sd) - 
                                              pnorm(lower, mean = mu, sd=sd)))
    
    res <- qnorm(q_temp_simulation, 
                 mean = mu, sd=sd)
    
    return(res)
    returnType(double(1))
  }
)


nimGetUtilde <- nimbleFunction(
  run = function(tau = double(1), u = double(1), a = double(0), b = double(0)){
    utilde <- nimNumeric(length(tau))
    for(i in 1:length(tau)){
      if(u[i] == 1) 
        utilde[i] <- nimRnormTrunc(1, a * tau[i] + b, lower = 0)[1]
      else{
        utilde[i] <- nimRnormTrunc(1, a * tau[i] + b, upper = 0)[1]
      }
    }
    return(utilde)
    returnType(double(1))
  }
)


nimRsample <- nimbleRcall(function(x = double(1), 
                                   size = double(0), 
                                   replace = logical(0),
                                   prob = double(1)){}, 
                          double(1), "sample")


chainningGen1 <- nimbleFunction(
  setup = function(tau, prior_sda) {
    dta_observed <- tau
    dta_observed[which(tau > 0)] <- 1
    tau_arrangedCapt <- t(tau)[which(t(tau) > 0)]
    tau_oneCapt <- apply(tau[which(apply(dta_observed, 1, sum) == 1), ], 1, sum)
    tau_oneCapt_index <- as.double(which(tau_arrangedCapt %in% tau_oneCapt))
    
    nbCapture <- sum(tau > 0)
    nbIndivMax <- sum(apply(dta_observed, 1, sum) > 0)
    nbUniqCapt <- sum(apply(dta_observed, 1, sum) == 1)
    if(nbUniqCapt/sum(dta_observed) > 0.3){
      nbMaxErr <- floor(sum(dta_observed)*0.3)
      nbIndivMin <- nbIndivMax - nbMaxErr
    }
    else{
      nbIndivMin <- sum(apply(dta_observed, 1, sum) > 1)
    }
    
    prior_mub <- (qnorm(1/nbIndivMax) + qnorm(1/nbIndivMin))/2
    prior_mua <- (qnorm(1 - round(nbUniqCapt/2)/sum(dta_observed)) - prior_mub) /
      mean(tau_arrangedCapt)
    prior_sdb <- abs(qnorm(1/nbIndivMax) - qnorm(1/nbIndivMin))*2
    sda <- sqrt( 1 / (sum(tau^2) + 1/(prior_sda^2)) )
    sdb <- sqrt( prior_sdb^2 / (nbCapture*prior_sdb^2 + 1) )
  },
  run <- function(nbIter = double(0), 
                  nbErr = double(0),
                  nthin = double(0),
                  nbBurn = double(0)){
    prior_mub <- qnorm(1/(nbIndivMax-nbErr+1) )
    a_next <- prior_mua
    b_next <- rnorm(1, prior_mub, sdb)
    samples <- nimMatrix(value = 0.0, nrow=nbIter, ncol = 2)
    
    for(iter in 1:(nbIter+nbBurn)){
      for(t in 1:nthin){
        # sample u
        alphaTmp <- pnorm(a_next * tau_oneCapt + b_next)
        p <- (1-alphaTmp)/sum(1-alphaTmp)
        u <- rep(1, nbCapture)
        samp <- nimRsample(tau_oneCapt_index,
                           nbErr,
                           replace = FALSE,
                           prob = p)
        u[samp] <- 0
        
        # update a
        u_tilde <- nimGetUtilde(tau_arrangedCapt, u, a_next, b_next)
        mua <- sda^2 * (prior_mua/(prior_sda^2) + sum(tau_arrangedCapt * (u_tilde - b_next)))
        a_next <- rnorm(1, mua, sda)
        
        # update b
        b_next <- rnorm(1, prior_mub, sdb)
        # u_tilde <- nimGetUtilde(tau_arrangedCapt, u, a_next, b_next)
        # mub <- sdb^2 * (prior_mub/(prior_sdb^2) + sum(u_tilde - a_next * tau_arrangedCapt))
        # b_next <- rnorm(1, mub, sdb)
      }
      if(iter > nbBurn) samples[iter-nbBurn, ] <- c(a_next, b_next)
    }
    return(samples)
    returnType(double(2))
  }
)

nimEstimPriorGen <- nimbleFunction(
  setup = function(tau, id, prior_mua, prior_sda, prior_mub, prior_sdb) {
    
    nb_id <- length(tau)
    
    a <- prior_mua
    b <- prior_mub
    
    sda <- sqrt( ( 1/(prior_sda^2) + sum(tau^2) )^(-1) )
    sdb <- sqrt( prior_sdb^2 / (nb_id*prior_sdb^2 + 1) )
  },
  run = function(niter = double(0), nburn = double(0)) {
    
    res <- nimMatrix(ncol = 2, nrow = niter+nburn)
    
    for(it in 1:(niter+nburn)){
      
      utilde <- nimGetUtilde(tau, id, a, b)
      mua <- sda^2 * (prior_mua/(prior_sda^2) + sum( tau * (utilde - b) ) )
      a <<- rnorm(1, mua, sda)
      
      utilde <- nimGetUtilde(tau, id, a, b)
      mub <- sdb^2 * (prior_mub/(prior_sdb^2) + sum(utilde - a *  tau))
      b <<- rnorm(1, mub, sdb)
      
      res[it, ] <- c(a, b)
    }
    
    return(res[(nburn+1):(nburn+niter), ])
    returnType(double(2))
  },
  methods = list(
    reset = function () {} )
)
