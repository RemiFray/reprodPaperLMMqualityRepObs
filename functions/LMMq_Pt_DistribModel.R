
# ------- Distribution ------ ----

dLMMqPt <- nimbleFunction(
  run = function(x = double(2), 
                 tau = double(2),
                 capture = double(1),
                 latentIndex = double(1),
                 a = double(0), 
                 b = double(0),
                 S = double(0), 
                 log = integer(0, default = 1) 
  ) {
    
    indexs <- which(latentIndex > 0)
    
    logProbData <- 0
    
    for(i in 1:length(indexs)){
      I <- indexs[i]
      hist <- x[I,]
      probHist <- 1
      for(t in 1:S){
        if(hist[t] == 0)
          probHist <- probHist * (1-capture[t])
        else if(hist[t] == 1) {
          alpha <- pnorm(a * tau[I, t] + b)
          probHist <- probHist * capture[t] * alpha
        }
        else{
          alpha <- pnorm(a * tau[I, t] + b)
          probHist <- probHist * capture[t] * (1-alpha)
        }
      }
      logProbData <- logProbData + log(probHist) 
    }
    
    if(log) return(logProbData)
    return(exp(logProbData))
    returnType(double(0))
  })

rLMMqPt <- nimbleFunction(
  run = function(n = integer(),
                 tau = double(2), capture = double(1), 
                 latentIndex = double(1),
                 a = double(0), b = double(0), S = double(0)) {
    returnType(double(2))
    return(tau) ## dummy behavior
  })

deregisterDistributions("dLMMqPt")
registerDistributions(
  list(dLMMqPt = list(BUGSdist = "dLMMqPt(tau, capture, latentIndex,
                        a, b, S)",
                            types =c("value = double(2)", 
                                     "tau = double(2)", 
                                     "capture = double(1)", 
                                     "latentIndex = double(1)"
                                     ))))


# ------- Model Code ------ ----

LMMqPt_Code <- nimbleCode({
  
  # priors
  for(t in 1:S){
    capture[t] ~ dbeta(1.0, 1.0)
  }
  
  a ~ dnorm(mean = prior_mua, sd = prior_sda)
  b ~ dnorm(mean = prior_mub, sd = prior_sdb)

  alpha_mean <- alphaMean
  alpha_sd <- alphaSd
  
  
  N ~ dunif(1, 10e4)
  D <- d
  
  nbErr <- NbErr
  
  # Likelihood, x
  latentObservation[1:n, 1:S] ~ dLMMqPt(tau[1:n, 1:S], capture[1:S], 
                                        latentIndex[1:n],
                                        a, b, S)
  
})

