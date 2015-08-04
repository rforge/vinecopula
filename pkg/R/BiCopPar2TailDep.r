BiCopPar2TailDep <- function(family, par, par2 = 0, obj = NULL) {
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
    if (is.na(family) || is.na(par)) 
        stop("Provide either 'family' and 'par' or 'obj'")
    if (!(family %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                        13, 14, 16, 17, 18, 19, 20,
                        23, 24, 26, 27, 28, 29, 30, 33, 34, 36, 37, 38, 39, 40, 
                        41, 42, 51, 52, 61, 62, 71, 72,
                        104, 114, 124, 134, 204, 214, 224, 234))) 
        stop("Copula family not implemented.")
    if (c(2, 7, 8, 9, 10,
          17, 18, 19, 20, 
          27, 28, 29, 30, 
          37, 38, 39, 40,
          42, 52, 62, 72,
          104, 114, 124, 134, 
          204, 214, 224, 234) %in% family && par2 == 0) 
        stop("For t-, BB1, BB6, BB7, BB8 and Tawn copulas, 'par2' must be set.")
    if (c(1, 3, 4, 5, 6, 11, 13, 14, 16, 23, 24, 26, 33, 34, 36, 41, 51, 61, 71) %in% 
            family && length(par) < 1) 
        stop("'par' not set.")
    BiCopCheck(family, par, par2)
    
    if (family == 0 | family == 1 | family == 5 | family %in% c(23, 24, 26, 27, 28, 29,
                                                                30, 33, 34, 36, 37, 38, 39,
                                                                40, 124, 134, 224, 234)) {
        lower <- 0
        upper <- 0
    } else if (family == 2) {
        lower <- 2 * pt((-sqrt(par2 + 1) * sqrt((1 - par)/(1 + par))), df = par2 + 
                            1)
        upper <- lower
    } else if (family == 3) {
        lower <- 2^(-1/par)
        upper <- 0
    } else if (family == 4 | family == 6) {
        lower <- 0
        upper <- 2 - 2^(1/par)
    } else if (family == 7) {
        lower <- 2^(-1/(par * par2))
        upper <- 2 - 2^(1/par2)
    } else if (family == 8) {
        lower <- 0
        upper <- 2 - 2^(1/(par * par2))
    } else if (family == 9) {
        lower <- 2^(-1/par2)
        upper <- 2 - 2^(1/par)
    } else if (family == 10) {
        lower <- 0
        if (par2 == 1) 
            upper <- 2 - 2^(1/par) else upper <- 0
    } else if (family == 13) {
        lower <- 0
        upper <- 2^(-1/par)
    } else if (family == 14 | family == 16) {
        lower <- 2 - 2^(1/par)
        upper <- 0
    } else if (family == 17) {
        lower <- 2 - 2^(1/par2)
        upper <- 2^(-1/par * par2)
    } else if (family == 18) {
        lower <- 2 - 2^(1/(par * par2))
        upper <- 0
    } else if (family == 19) {
        lower <- 2 - 2^(1/par)
        upper <- 2^(-1/par2)
    } else if (family == 20) {
        if (par2 == 1) 
            lower <- 2 - 2^(1/par) else lower <- 0
        
        upper <- 0
    } else if (family == 104) {
        par3 <- 1
        upper <- par2 + par3 - 2 * ((0.5 * par2)^par + (0.5 * par3)^par)^(1/par)
        lower <- 0
    } else if (family == 114) {
        par3 <- 1
        lower <- par2 + par3 - 2 * ((0.5 * par2)^par + (0.5 * par3)^par)^(1/par)
        upper <- 0
    } else if (family == 204) {
        par3 <- par2
        par2 <- 1
        upper <- par2 + par3 - 2 * ((0.5 * par2)^par + (0.5 * par3)^par)^(1/par)
        lower <- 0
    } else if (family == 214) {
        par3 <- par2
        par2 <- 1
        lower <- par2 + par3 - 2 * ((0.5 * par2)^par + (0.5 * par3)^par)^(1/par)
        upper <- 0
    }
    
    return(list(lower = lower, upper = upper))
}