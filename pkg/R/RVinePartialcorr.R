# functions for partial correlations

# specific partial correlation from a covariance or correlation matrix
# 'given' is vector indices for the given variables
# j,k are indices for the conditioning variables
partcor=function(S,given,j,k) {
  S11=S[given,given]
  jk=c(j,k)
  S12=S[given,jk]
  S21=S[jk,given]
  S22=S[jk,jk]
  if(length(given)>1) { tem=solve(S11,S12); Om212=S21%*%tem }
  else { tem=S12/S11; Om212=outer(S21,tem) }
  om11=1-Om212[1,1]
  om22=1-Om212[2,2]
  om12=S[j,k]-Om212[1,2]
  om12/sqrt(om11*om22)
}

#============================================================

# correlations to partial correlations and vice
# versa for general R-vine with RVineMatrix

# param: cor, correlation matrix
# param: RVM, RVineMatrix defining the strucutre of the RVine
# result: RVineMatrix with transformed partial correlations
cor2pcor=function(corMat, RVM) {
  d <- nrow(corMat)
  stopifnot(d == nrow(RVM$Matrix))
  if(d<=2) 
    return(corMat)
  
  pp <- matrix(0, d, d)
  A  <- matrix(0, d, d)
  
  # rotate towards notation in Kurowicka and Joe (2011), p. 9
  for (i in 1:d) {
    A[i,] <- rev(RVM$Matrix[d-i+1,])
  }
  
  # following algorithm is credited to Harry Joe
  for(j in 2:d) { # j <- 2
    pp[1,j] <- corMat[A[1,j],j]
  }
  
  # tree 2
  for(j in 3:d) {
    a1 <- A[1,j]
    a2 <- A[2,j] 
    pp[2,j] <- (corMat[j,a2]-corMat[j,a1]*corMat[a1,a2])/sqrt((1-corMat[j,a1]^2)*(1-corMat[a1,a2]^2))
  }
  
  # remaining trees
  for(ell in 3:(d-1)) { 
    for(j in (ell+1):d) {
      given <- A[1:(ell-1),j]
      pp[ell,j] <- partcor(corMat,given,A[ell,j],j)  # assuming A[j,j]=j
    }
  }
  
  # re-rotate towards VineCopula notation
  pc <- pp
  for (i in 1:d) {
    pc[i,] <- rev(pp[d-i+1,])
  }
  
  RVM$par <- pc
  return(RVM)
}


# corMat <- matrix(c(1.00, 0.17, 0.15, 0.14, 0.13,
#                    0.17, 1.00, 0.30, 0.28, 0.05,
#                    0.15, 0.30, 1.00, 0.17, 0.05,
#                    0.14, 0.28, 0.17, 1.00, 0.04,
#                    0.13, 0.05, 0.05, 0.04, 1.00),5,5)
# 
# 
#  
# RVM <- vineCopula(5L)@RVM
# str(RVM)
# 
# cor2pcor.cvine(corMat)
# 
# newRVM <- cor2pcor(corMat, RVM)
# newRVM$family <- matrix(1,5,5)
# 
# # classic
# round(cor(rmvnorm(10000,rep(0,5),corMat))-corMat,2)
# # vine
# round(cor(qnorm(rCopula(10000, vineCopula(newRVM))))-corMat,2)
# 
# pcor2cor(cor2pcor(corMat, RVM))
# round(pcor2cor(RVM=newRVM)-corMat,2)
# cor2pcor(corMat, RVM)$par
# 

# library(gstat)
# data(sic2004)
# coordinates(sic.val) <- ~x+y
# hist(variogramLine(vgm(0.8,"Mat",150000,0.2),covariance=T,dist_vector=spDists(sic.val[,])))
# 
# 
# tauFromCor <- function(x) {
#   tau(normalCopula(x))
# }
# 
# plot(1:200*2500,sapply(variogramLine(vgm(0.8,"Mat",150000,0.2),covariance=T,dist_vector=1:200*2500)[,2],
#       tauFromCor),type="l",ylim=c(0,1))
# 
# points(1:200*2500,variogramLine(vgm(0.8,"Mat",150000,0.2),covariance=T,dist_vector=1:200*2500)[,2],
#        type="l", col="red")



# generate correlation matrix based on partial correlations of R-vine
# with vine array A that has 1:d on diagonal;

pcor2cor <- function(RVM) {
  d=nrow(RVM$Matrix)
  
  A <- matrix(0,d,d)
  # rotate towards notation in Kurowicka and Joe (2011), p. 9
  for (i in 1:d) {
    A[i,] <- rev(RVM$Matrix[d-i+1,])
  }
  
  pc <- matrix(0,d,d)
  for (i in 1:d) {
    pc[i,] <- rev(RVM$par[d-i+1,])
  }
  
  if(d<=2) { 
    corMat <- matrix(c(1,pc[1,2],pc[1,2],1))
    return(corMat) 
  }
  
  corMat <- matrix(0,d,d)
  diag(corMat) <- 1
  for(j in 2:d) { 
    a1 <- A[1,j]
    corMat[a1,j] <- pc[1,j]
    corMat[j,a1] <- pc[1,j]
  }
  
  # tree 2
  for(j in 3:d) { 
    a1 <- A[1,j]
    a2 <- A[2,j] 
    corMat[j,a2] <- corMat[j,a1]*corMat[a1,a2]+pc[2,j]*sqrt((1-corMat[j,a1]^2)*(1-corMat[a1,a2]^2))
    corMat[a2,j] <- corMat[j,a2]
  }
  
  # remaining trees
  for(ell in 3:(d-1)) { 
    for(j in (ell+1):d) { 
      given <- A[1:(ell-1),j]
      S11 <- corMat[given,given]
      anew <- A[ell,j]
      jk <- c(anew,j)
      S12 <- corMat[given,jk]
      S21 <- corMat[jk,given]
      S22 <- corMat[jk,jk]
      tem <- solve(S11,S12)
      Om212 <- S21%*%tem
      om11 <- 1-Om212[1,1]
      om22 <- 1-Om212[2,2]
      tem12 <- pc[ell,j]*sqrt(om11*om22)
      corMat[anew,j] <- tem12+Om212[1,2]
      corMat[j,anew] <- corMat[anew,j]
    }
  }
  corMat
}