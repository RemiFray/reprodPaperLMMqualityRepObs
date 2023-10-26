library(nimble)
library(tidyverse)

# ----- Functions, samplers, distrib and model ----- ----

source("./functions/LMM_Pt_DistribModel.R")
source("./functions/LMM_Pt_functions.R")
source("./functions/LMM_Pt_samplers.R")

# ----- Parameters ----- ----

modelType <- "LMMPt"
S <- 7

alpha_test <- c(0.8, 0.9, 0.95)
p_test <- c(0.4, 0.3, 0.2, 0.1)
N <- 500

priorType <- "uninformative"
priorAlpha <- c(1, 1)

bias_plan <- data.frame(a = rep(alpha_test, each = length(p_test)*100),
                        p = rep(rep(p_test, each = 100), length(alpha_test)),
                        iter = rep(rep(1:100, length(p_test)), length(alpha_test)))

burnin <- 100000
nthin <- 100
niter <- 10000*nthin+burnin


# ------- File management ------ ----

tmpDir <- paste0("tmp/", modelType, "_S", S, "_a", alpha, "_noPrior")
if(!file.exists(tmpDir)) dir.create(tmpDir, recursive = TRUE)

resDir <- paste0("./results/noRep/", modelType, "_S",S, "_", priorType)
if(!file.exists(resDir)) dir.create(resDir, recursive = TRUE)


# ------- Big loop ------ ----

for(line in 1:nrow(bias_plan)){
  
  # ------- Parameters ------ ----  
  
  alpha <- bias_plan$a[line]
  capture <- rep(bias_plan$p[line], S)
  
  iter <- bias_plan$iter[line]
  
  # ------- Data ------ ----
  
  load(file = paste0("./data/LMMq_Simul_S",S, "/LMMqSimul_S", S, "_N", N, "_a", alpha, 
                     "_p", capture[1], "_iter", iter, ".Rdata"))
  
  dta_observed <- res_observed$observed
  summaryDta <- getSummaryCMR(dta_observed)
  
  # MATRIX OF ORDERED OBSERVED HISTORIES
  observation <- as.data.frame(t(sapply(summaryDta$history, 
                                        function(h) eval(parse(text=h)))))
  rownames(observation) <- NULL
  observation$index <- apply(observation, 1, nimGetLatentIndex)
  observation <- observation %>% 
    arrange(index) %>% select(-index) %>% as.matrix()
  colnames(observation) <- NULL
  
  # MATRIX OF HISTORIES AND SPACE FOR OTHERS
  n <- nrow(observation)
  M <- min(sum(summaryDta$Freq[-1]), 3^S-n)
  latentObservation <- matrix(nrow=n+M, ncol = S, data=0)
  latentObservation[1:n,] <- observation
  latentIndex <- numeric(n+M)
  latentIndex[1:n] <- summaryDta$index
  
  xInit <- numeric(n+M)
  xInit[1:n] <- summaryDta$Freq
  
  # MATRIX A (just to check if xInit if valid)
  A <- getA(latentObservation, observation)
  
  # INITIAL X TO START MCMC
  set.seed(1248)
  xInit2 <- xInit
  latObs2 <- latentObservation
  latIdx2 <- latentIndex
  for(j in 1:40){
    tmp <- addError(xInit2, latObs2, latIdx2, S, n, M, n+M)
    xInit2 <- tmp$x
    latObs2 <- tmp$latObs
    latIdx2 <- tmp$latIdx
    
    A2 <- getA(latObs2, observation)
    if(!(all(xInit2 >= 0) & all(summaryDta$Freq[-1] == A2 %*% xInit2)))
      print(j)
  }
  cat(all(summaryDta$Freq[-1] == A %*% xInit), "  ")
  A2 <- getA(latObs2, observation)
  cat(all(xInit2 >= 0) & all(summaryDta$Freq[-1] == A2 %*% xInit2), "  ")
  print(sum(xInit-xInit2))
  
  # ------- Constructing nimble model ------ ----
  
  
  LMMPtConsts <- list(S = S, nbLatentObs=n+M)
  
  # Initialisation
  LMMPtInits <- function(i) list(capture = rep(0.5, S),
                                 alpha = 0.5, 
                                 N = c(sum(xInit),sum(xInit2))[[i]],
                                 x = list(xInit, xInit2)[[i]],
                                 latentObservation = list(latentObservation, latObs2)[[i]],
                                 latentIndex = list(latentIndex, latIdx2)[[i]])
  
  LMMPtModel <- nimbleModel(code = LMMPtCode, name = "Pt", 
                            constants = LMMPtConsts,
                            data = list(), inits = LMMPtInits(1))
  
  # Compilation modèle OK
  CLMMPt <- compileNimble(LMMPtModel, 
                          dirName = tmpDir,
                          showCompilerOutput = FALSE)
  
  
  # ------- MCMC configuration ------ ----
  
  # Configuration of MCMC
  LMMPtConf <- configureMCMC(LMMPtModel,
                             monitors = c("N", "D", "capture", "alpha"), 
                             thin = 1,
                             # enableWAIC = T,
                             print = TRUE)
  
  LMMPtConf$removeSampler("x", "capture", "alpha", "N")
  LMMPtConf$addSampler(target = c("alpha"), type = nimUpdateA, 
                       control = list(prior = priorAlpha))
  LMMPtConf$addSampler(target = paste0("capture[1:",S,"]"),
                       type = nimUpdatePt)
  LMMPtConf$addSampler(target = c("x"), type = nimSamplerXmove,
                       control = list(D = 1, latentIndex=as.double(latentIndex),
                                      n = n, M = M))
  LMMPtConf$addSampler(target = c("x"), type = nimSamplerX0,
                       control = list(D = 5))
  LMMPtConf$printSamplers()
  
  # Compilation MCMC
  LMMPtMCMC <- buildMCMC(LMMPtConf)
  
  
  CLMMPtMCMC <- compileNimble(LMMPtMCMC, project = LMMPtModel,
                              dirName = tmpDir,
                              showCompilerOutput = F)
  
  # ------- run MCMC ------ ----
  
  inits <- list(LMMPtInits(1), LMMPtInits(2))
  
  system.time(
    samples <- runMCMC(CLMMPtMCMC, niter = niter, nburnin = burnin, 
                       thin = nthin, nchains = 2,inits = inits, setSeed = c(777, 1234))
  )
  
  save(samples,
       file = paste0(resDir, 
                     "/LMM_S", S, "_", priorType, "_N", N, "_a", alpha, 
                     "_p", capture[1], "_iter", iter, ".Rdata"))
  
}



