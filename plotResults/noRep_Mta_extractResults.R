# -------- Packages and functions -------- ----
library(dplyr)
library(MCMCvis)

wideSummary <- function(tmp){
  # if("D" %in% rownames(tmp)) idx <- c(2:4,6,8)
  # else idx <- c(1:3,5,7)
  
  rnames <- gsub("\\[|, |\\]", "", rownames(tmp), perl=TRUE)
  rnames <- gsub("capture", "p", rnames, perl=TRUE)
  rnames <- gsub("alpha", "a", rnames, perl=TRUE)
  cnames <- gsub("%", "", colnames(tmp), perl=TRUE)
  allCols <- as.character(sapply(rnames, function(x) paste(x, cnames, sep="_")))

  df <- as.matrix(tmp) %>% 
    t() %>% as.numeric() %>% 
    matrix(nrow=1) %>% as.data.frame()
  colnames(df) <- allCols
  df
}


# -------- Parameters -------- ----

S <- c(5, 7, 9)
alpha_test <- c(0.8, 0.9, 0.95)
p_test <- c(0.1, 0.2, 0.3, 0.4)
N <- 500
model <- "LMMPt"
prior <- c("uninformative")


# -------- Loop for extraction -------- ----

# 6 : params
# D, N, alpha, psi (9), capture(S*3), delta(3) *7
# -> (1+1+1+9+s*3+3)*5 = (15+3*s)*7
# => 216 colonnes pour S=9

for(s in S){
  cat(s, "\n")
  testSummaries <- as.data.frame(matrix(ncol = 6+(15+3*s)*7, nrow=0))
  for(pr in prior){
    cat("\t", pr, "\n")
      for(p in p_test){
        cat("\t\t\t", p, "\n")
        for(a in alpha_test){
          directory <- paste0("./results/norep/", model, "_S", s, "_", pr)
          # if(!file.exists(directory)) cat(directory, "\n")
          files <- list.files(directory)
          for(iter in 1:100){
            fileExp <- paste0("N", N, "_a", a, "_p", p, "_iter", iter, ".Rdata")
            filePath <- files[which(grepl(fileExp, files))]
            load(file.path(directory, filePath))
            
            tmp <- data.frame(Prior=pr, S=s, N=N, p=p, 
                       alpha=a, iter=iter) %>% 
              cbind(wideSummary(MCMCsummary(samples)))
            testSummaries <- rbind(testSummaries, tmp)
            
          }
        }
      }
    }
  
  write.csv(testSummaries,
            file=paste0("./results/noRep/", model, "_S", s, "_Results.csv"), 
            row.names = F)
  
}


# -------- Merge tables -------- ----

resS5 <- read.csv(paste0("./results/noRep/", model, "_S5_Results.csv"))
resS7 <- read.csv(paste0("./results/noRep/", model, "_S7_Results.csv"))
resS9 <- read.csv(paste0("./results/noRep/", model, "_S9_Results.csv"))

allres <-  resS5 %>% 
  rbind(resS7[,!grepl("p[6-9]", colnames(resS7))]) %>% 
  rbind(resS9[,!grepl("p[6-9]", colnames(resS9))])

write.csv(allres, file=paste0("./results/noRep/", model, "_AllResults.csv"), 
          row.names = F)
  
  