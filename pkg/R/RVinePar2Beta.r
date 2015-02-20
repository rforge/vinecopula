RVinePar2Beta <- function(RVM) {
    
    if (is(RVM)[1] != "RVineMatrix") 
        stop("'RVM' has to be an RVineMatrix object.")
    
    taus <- RVM$par
    n <- dim(RVM)
    
    for (i in 2:n) {
        for (j in 1:(i - 1)) {
            taus[i, j] <- BiCopPar2Beta(RVM$family[i, j], 
                                        RVM$par[i, j],
                                        RVM$par2[i, j])
        }
    }
    
    return(taus)
}
