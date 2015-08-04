BiCopSim <- function(N, family, par, par2 = 0, obj = NULL) {
    ## extract family and parameters if BiCop object is provided
    if (missing(family))
        family <- NA
    if (missing(par))
        par <- NA
    # for short hand usage extract obj from family
    if (class(family) == "BiCop")
        obj <- family
    if (!is.null(obj)) {
        stopifnot(class(obj) == "BiCop")
        family <- obj$family
        par <- obj$par
        par2 <- obj$par2
    }
    
    ## sanity checks for family and parameters
    BiCopCheck(family, par, par2)
    
    tmp <- .C("pcc", 
              as.integer(N), 
              as.integer(2), 
              as.integer(family), 
              as.integer(1), 
              as.double(par),
              as.double(par2), 
              as.double(rep(0, N * 2)),
              PACKAGE = "VineCopula")[[7]]
    
    ## return results
    U <- matrix(tmp, ncol = 2)
    U
}
