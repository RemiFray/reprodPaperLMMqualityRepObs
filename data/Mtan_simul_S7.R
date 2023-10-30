library(nimble)
library(tidyverse)
library(EnvStats)

# ----- Functions, samplers, distrib and model ----- ----

source("./functions/LMMq_Pt_functions.R")

# ----- Parameters ----- ----

S <- 7
n_rate <- 1.5

tau_mean <- 0.6
tau_sd <- 0.1
trunc_min <- 0.4
trunc_max <- 1
tau_median <- tau_mean +
  qnorm( ( pnorm( (trunc_min - tau_mean)/tau_sd ) +
             pnorm( (trunc_max - tau_mean)/tau_sd )) / 2 ) * tau_sd

# should lead to approx alpha mean of 0.8, 0.9 and 0.95
# (whith b = qnorm(1/1001) = -3.09)
alphaMean <- c(0.8, 0.9, 0.95)
alpha_median <- read.csv("data/alphaMedianToUse.csv")
alpha_min <- 0.6


p_test <- c(0.4, 0.3, 0.2, 0.1)
N_test <- c(500)

bias_plan <- data.frame(N = rep(N_test, each = length(p_test)*length(alphaMean)),
                        aMean = rep(rep(alphaMean, each = length(p_test)), length(N_test)),
                        a50 = numeric(length(N_test)*length(alphaMean)*length(p_test)),
                        p = rep(rep(p_test, length(alphaMean)), length(N_test)))
for(i in 1:nrow(bias_plan)){
  l <- alpha_median %>% 
    filter(S == S, p == bias_plan$p[i], N == bias_plan$N[i])
  if(bias_plan$aMean[i] == 0.8) bias_plan$a50[i] <- l[1, 6]
  else if(bias_plan$aMean[i] == 0.9) bias_plan$a50[i] <- l[1, 7]
  else bias_plan$a50[i] <- l[1, 8]
}

dtaDir <- paste0("./data/LMMq_Simul_S",S)
if(!file.exists(dtaDir)) dir.create(dtaDir, recursive = TRUE)


# ------- Big loop ------ ----


for(line in 1:nrow(bias_plan)){
  set.seed(1234)
  
  N <- bias_plan$N[line]
  capture <- rep(bias_plan$p[line], S)
  a_mean <- bias_plan$aMean[line]
  alpha_50 <- bias_plan$a50[line]
  a <- (qnorm(alpha_50) - qnorm(alpha_min)) / (tau_median - trunc_min)
  
  for(iter in 1:100){
    
    # ------- Preparing data ------ ----
    
    capture_data <- simulCapture(N, S, capture)
    
    b <- qnorm(1/(sum( apply(capture_data, 1, sum) > 0 )+1))
    
    tau_values = rnormTrunc(sum(capture_data), tau_mean, tau_sd, 
                            min = trunc_min, max = trunc_max)
    tau <- capture_data
    tau[which(tau > 0)] <- tau_values
    
    res_observed <- simulIdentification(capture_data, a, b, tau)
    res_observed$true_capture <- capture_data
    res_observed$a <- a
    res_observed$b <- b
    
    
    save(res_observed,
         file = paste("./data/LMMq_Simul_S",S, "/LMMqSimul_S", S, "_N", N, "_a", a_mean, 
                      "_p", capture[1], "_iter", iter, ".Rdata", sep = ""))
  }
}

