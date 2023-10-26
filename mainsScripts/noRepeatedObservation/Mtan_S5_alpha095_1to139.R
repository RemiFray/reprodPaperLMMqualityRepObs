library(nimble)
library(tidyverse)

# ----- Functions, samplers, distrib and model ----- ----

source("./functions/LMMq_Pt_DistribModel.R")
source("./functions/LMMq_Pt_functions.R")
source("./functions/LMMq_Pt_samplers.R")
source("./functions/getPrior_functions.R")

# ----- Parameters ----- ----

modelType <- "LMMqPtab"
S <- 5

alpha <- 0.95
p_test <- c(0.4, 0.3, 0.2, 0.1)
N <- 500

priorType <- "uninformative"

bias_plan <- data.frame(p = rep(p_test, each = 100),
                        iter = rep(1:100, length(p_test)))

burnin <- 40000
nthin <- 100
niter <- 10000*nthin+burnin


# ------- File management ------ ----

tmpDir <- paste0("tmp/", modelType, "_S", S, "_a", alpha, "_noPrior")
if(!file.exists(tmpDir)) dir.create(tmpDir, recursive = TRUE)

resDir <- paste0("./results/noRep/", modelType, "_S",S, "_", priorType)
if(!file.exists(resDir)) dir.create(resDir, recursive = TRUE)


# ------- Big loop ------ ----

for(line in 1:139){
  
  # ------- Parameters ------ ----  

  capture <- rep(bias_plan$p[line], S)
  
  iter <- bias_plan$iter[line]
  
  # ------- Data ------ ----
  
  load(file = paste0("./data/LMMq_Simul_S",S, "/LMMqSimul_S", S, 
                     "_N", N, "_a", alpha, "_p", capture[1], 
                     "_iter", iter, ".Rdata"))
  
  dta_observed <- res_observed$observed
  tau <- res_observed$tau
  
  # PRIORS PROBIT
  source("./functions/getPrior_main.R")
  priorMuB <- matrix(c(res_observed$prior_mua, 
                       res_observed$prior_mub), ncol = 1)
  priorSdB <- matrix(c(res_observed$prior_sda, 0, 0, res_observed$prior_sdb), 
                     ncol = 2)
  
  # MATRIX OF OBSERVED HISTORIES WITH SPACE FOR OTHER INDIV
  observation <- dta_observed[which(apply(dta_observed, 1, sum) > 0),]
  nb_obs <- nrow(observation)
  n <- 4*N
  tau <- tau[which(apply(dta_observed, 1, sum) > 0),]
  tauTmp <- tau
  tau <- matrix(nrow = n, ncol = S, data=0)
  tau[1:nb_obs, ] <- tauTmp
  latentObservation <- matrix(nrow = n, ncol = S, data=0)
  latentObservation[1:nb_obs,] <- observation
  latentIndex <- numeric(n)
  latentIndex[1:nb_obs] <- apply(observation, 1, nimGetLatentIndex)
  unseen_init <- sample(5:50, 1)
  latentIndex[which(latentIndex == 0)[1:unseen_init]] <- 1
  rm(tauTmp, nb_obs)
  
  # DIFFERENT INITIALIZATION 
  set.seed(1248)
  latObs2 <- latentObservation
  latIdx2 <- latentIndex
  tau2 <- tau
  nbErr <- 0
  for(j in 1:40){
    tmp <- addError(tau2, latObs2, latIdx2, S)
    latObs2 <- tmp$latObs
    latIdx2 <- tmp$latIdx
    tau2 <- tmp$tau
    nbErr <- nbErr + tmp$err
  }
  
  nb_err_init2 <- sum(latentIndex>0) - sum(latIdx2>0)
  cat("Error added to second init: ", nb_err_init2, "\n")
  
  
  # ------- Constructing nimble model ------ ----
  
  
  LMMqPt_Consts <- list(S = S, n = n, 
                        prior_mua = priorMuB[1,1], prior_mub = priorMuB[2,1],  
                        prior_sda = priorSdB[1,1], prior_sdb = priorSdB[2,2], 
                        nb_identification = sum(res_observed$observed))
  
  # Initialisation
  a_init <- priorMuB[1,1]
  b_init <- c(qnorm(1/sum(apply(res_observed$observed, 1, sum) > 0)),
       qnorm(1/(sum(apply(res_observed$observed, 1, sum) > 0)-nb_err_init2)))
  LMMqPt_Inits <- function(i) list(capture = rep(0.5, S),
                                   a = a_init, b = b_init[i],
                                   alphaMean = 0.5, alphaSd = 0.1,
                                   N = c(sum(latentIndex>0), sum(latIdx2>0))[i],
                                   d = c(unseen_init, sum(latIdx2 == 1))[i],
                                   NbErr = c(0, nb_err_init2)[i],
                                   tau = list(tau, tau2)[[i]],
                                   latentObservation = list(latentObservation, 
                                                            latObs2)[[i]],
                                   latentIndex = list(latentIndex, latIdx2)[[i]]
  )
  
  LMMqPt_Model <- nimbleModel(code = LMMqPt_Code, name = "LMMqPt", 
                                constants = LMMqPt_Consts,
                                data = list(), inits = LMMqPt_Inits(1))
  LMMqPt_Model$getLogProb("latentObservation")
  
  # Compilation modèle OK
  CLMMqPAPt_Model <- compileNimble(LMMqPt_Model, 
                                   dirName = tmpDir,
                                   showCompilerOutput = FALSE)
  
  
  
  # ------- MCMC configuration ------ ----
  
  # Configuration of MCMC
  LMMqPt_Conf <- configureMCMC(LMMqPt_Model,
                               monitors = c("N", "D", "nbErr", "capture", 
                                            "alpha_mean", "a", "b"), 
                               monitors2 = c("latentIndex"),
                               thin = 1, thin2 = 5*nthin, 
                               print = FALSE)
  
  LMMqPt_Conf$removeSampler("capture", "a", "b", "N")
  LMMqPt_Conf$addSampler(target = c("a"), type = nimUpdatea,
                         control = list(prior_mua = priorMuB[1,1],
                                        prior_sda = priorSdB[1,1]
                         ))
  LMMqPt_Conf$addSampler(target = c("b"), type = nimUpdateb,
                         control = list(prior_mub = priorMuB[2,1],
                                        prior_sdb = priorSdB[2,2]
                         ))
  LMMqPt_Conf$addSampler(target = paste0("capture[1:",S,"]"),
                         type = nimUpdatePt)
  LMMqPt_Conf$addSampler(target = c("latentObservation"), type = nimSamplerXmove,
                         control = list(n = n, D = 1))
  LMMqPt_Conf$addSampler(target = c("latentObservation"), type = nimSamplerX0,
                         control = list(D = 5))
  LMMqPt_Conf$printSamplers()
  
  LMMqPt_MCMC <- buildMCMC(LMMqPt_Conf)
  
  CLMMqPt_MCMC <- compileNimble(LMMqPt_MCMC, project = LMMqPt_Model,
                                  dirName = tmpDir,
                                  showCompilerOutput = F, resetFunctions = T)
  
  
  # ------- run MCMC ------ ----
  
  inits <- list(LMMqPt_Inits(1), LMMqPt_Inits(2))
  
  system.time(
    samples <- runMCMC(CLMMqPt_MCMC, niter = niter, nburnin = burnin, 
                       thin = nthin, nchains = 2,inits = inits, 
                       setSeed = c(777, 1234))
  )
  
  save(samples,
       file = paste0(resDir,  
                     "/LMM_S", S, "_", priorType, "_N", N, "_a", alpha, 
                     "_p", capture[1], "_iter", iter, ".Rdata"))
  
  rm(LMMqPt_Model, LMMqPt_MCMC, CLMMqPt_Model, CLMMqPt_MCMC)
  
  unlink(paste0(tmpDir, "/*"))
}



