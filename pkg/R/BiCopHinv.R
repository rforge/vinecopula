BiCopHinv <- function(u1, u2, family, par, par2 = 0, obj = NULL) {
    ## sanity checks for u1, u2
    if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
        stop("u1 and/or u2 are not set or have length zero.")
    if (length(u1) != length(u2)) 
        stop("Lengths of 'u1' and 'u2' do not match.")
    if (any(c(u1, u2) > 1) || any(c(u1, u2) < 0)) 
        stop("Data has be in the interval [0,1].")
    
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
    
    ## check for family/parameter consistency
    BiCopCheck(family, par, par2)
    
    ## hinv(u2 | u1)
    hinv1 <- .C("Hinv1",
                as.integer(family),
                as.integer(length(u2)),
                as.double(u2),
                as.double(u1),
                as.double(par),
                as.double(par2),
                as.double(rep(0, length(u1))),
                package = "VineCopula")[[7]]
    
    ## hinv(u2 | u1)
    hinv2 <- .C("Hinv2",
                as.integer(family),
                as.integer(length(u2)),
                as.double(u1),
                as.double(u2),
                as.double(par),
                as.double(par2),
                as.double(rep(0, length(u1))),
                package = "VineCopula")[[7]]
    
    
    ## return results
    list(hinv1 = hinv1, hinv2 = hinv2)
}

