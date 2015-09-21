##' Testsuite
##' 
##' Tests for the VineCopula package
##' 
##' @author Dr. Ulf Schepsmeier
##' 

## Main function

library(VineCopula)

source("../tests/testRun.r")
source("../tests/testCheck.r")

# BiCopPar2Tau
results_BiCopPar2Tau <- testRunBiCopPar("BiCopPar2Tau")
check_BiCopPar2Tau <- testCheck(results_BiCopPar2Tau)
if(!all(check_BiCopPar2Tau)){
  print(check_BiCopPar2Tau)
} else {
  rm(results_BiCopPar2Tau)
  gc()
}

# BiCopPar2Beta
results_BiCopPar2Beta <- testRunBiCopPar("BiCopPar2Beta")
check_BiCopPar2Beta <- testCheck(results_BiCopPar2Beta)
if(!all(check_BiCopPar2Beta)){
  print(check_BiCopPar2Beta)
} else {
  rm(results_BiCopPar2Beta)
  gc()
}


