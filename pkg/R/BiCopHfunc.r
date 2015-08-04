###### BiCopHfunc 
# Input: 
# u1,u2 copula data 
# family copula family 
# par copula parameter 
# par2 copula parameter 2 
# 
# Output:
# hfunc1 h-function h(u1,u2) 
# hfunc2 h-function h(u2,u1) 

BiCopHfunc <- function(u1, u2, family, par, par2 = 0, obj = NULL) {
    ## sanity checks for u1, u2
    if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
        stop("u1 and/or u2 are not set or have length zero.")
    if (length(u1) != length(u2)) 
        stop("Lengths of 'u1' and 'u2' do not match.")
    if (any(u1 > 1) || any(u1 < 0)) 
        stop("Data has be in the interval [0,1].")
    if (any(u2 > 1) || any(u2 < 0)) 
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

    ## sanity checks for family and parameters
    BiCopCheck(family, par, par2)
    
    ## h(u2|u1)
    hfunc1 <- .C("Hfunc1", 
                 as.integer(family),
                 as.integer(length(u1)), 
                 as.double(u2), 
                 as.double(u1), 
                 as.double(par),
                 as.double(par2), 
                 as.double(rep(0, length(u1))),
                 PACKAGE = "VineCopula")[[7]]
    ## h(u1|u2)
    hfunc2 <- .C("Hfunc2", 
                 as.integer(family),
                 as.integer(length(u1)), 
                 as.double(u1),
                 as.double(u2), 
                 as.double(par),
                 as.double(par2), 
                 as.double(rep(0, length(u1))),
                 PACKAGE = "VineCopula")[[7]]
    
    ## return results
    list(hfunc1 = hfunc1, hfunc2 = hfunc2)
}
