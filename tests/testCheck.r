##' Testsuite - Check
##' 
##' Run several tests for the BiCop-functions of the VineCopula-package
##' 
##' @author Dr. Ulf Schepsmeier
##' @param results list of results returned from testRun*

testCheck <- function(results){
  ## length of results
  n <- length(results)
  
  check <- rep(TRUE, n)
  
  for(i in 1:n){
    ## Check 1: is.na
    if(any(is.na(results[[i]]))) check[i] <- FALSE
    ## Check 2: is.nan
    if(any(is.nan(results[[i]]))) check[i] <- FALSE
    ## Check 3: is.infinite
    if(any(is.infinite(results[[i]]))) check[i] <- FALSE
    ## Check 4: in range
    if(names(results)[i] %in% c(1:10,13,14,16:20,104,114,204,214)){
      if(any( results[[i]] < 0 || results[[i]] > 1 ) ) check[i] <- FALSE
    } else {
      if(any( results[[i]] > 0 || results[[i]] < -1 ) ) check[i] <- FALSE
    }
    ## check for jumps
    ## TODO
  }
  
  return(check)
}