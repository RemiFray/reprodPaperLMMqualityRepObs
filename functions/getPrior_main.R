# script to be sourced after 
#   - loading functions/getPrior.R
#   - simulating data to get a list "res_observed"

prior_sda <- 5

chainning1 <- chainningGen1(res_observed$tau, prior_sda)
cChainning1 <- suppressMessages(compileNimble(chainning1))

nbUniqCapt <- sum(apply(res_observed$observed, 1, sum) == 1)
if(nbUniqCapt/sum(res_observed$observed) > 0.3){
  nbMaxErr <- floor(sum(res_observed$observed)*0.3)
} else
  nbMaxErr <- nbUniqCapt

nbIter <- 10
nthin <- 1
# 
# res <- matrix(nrow = nbMaxErr*nbIter, ncol = 2,
#               dimnames = list(NULL, c("a", "b")))
# cat("|-----------|-----------|\n")
# c <- 0.04
# for(j in c(1:nbMaxErr)){
# res[((j-1)*nbIter+1):(j*nbIter), ] <-
#   cChainning1$run(nbIter = nbIter,
#                   nbErr = j,
#                   nthin = nthin,
#                   nbBurn = 200)
# if(j/nbMaxErr > c){ cat("-") ; c <- c + 0.04}
# }
# cat("\n")
# res_observed[["prior_mua"]] <- estimate_mode(res[, 1])
# res_observed[["prior_mub"]] <- mean(res[, 2])
# res_observed[["prior_sda"]] <- sd(res[, 1])
# res_observed[["prior_sdb"]] <- sd(res[, 2])
# 
# res_observed[["prior_mina"]] <- min(res[, 1])
# res_observed[["prior_maxa"]] <- max(res[, 1])
# 
# res_observed[["prior_minb"]] <- min(res[, 2])
# res_observed[["prior_maxb"]] <- max(res[, 2])


# min (1 error)
restmp <- cChainning1$run(nbIter = nbIter,
                          nbErr = 1,
                          nthin = nthin,
                          nbBurn = 200)
res_observed[["prior_maxa"]] <- max(restmp[, 1])
res_observed[["prior_minb"]] <- min(restmp[, 2])
# max (nbErrMax)
restmp <- cChainning1$run(nbIter = nbIter,
                          nbErr = nbMaxErr,
                          nthin = nthin,
                          nbBurn = 200)
res_observed[["prior_mina"]] <- min(restmp[, 1])
res_observed[["prior_maxb"]] <- max(restmp[, 2])

res_observed[["prior_mua"]] <- mean(c(res_observed$prior_mina,
                                      res_observed$prior_maxa))
res_observed[["prior_mub"]] <- mean(c(res_observed$prior_minb,
                                      res_observed$prior_maxb))
res_observed[["prior_sda"]] <- (res_observed$prior_maxa -
                                  res_observed$prior_mina)/4
res_observed[["prior_sdb"]] <- (res_observed$prior_maxb -
                                  res_observed$prior_minb)/4
