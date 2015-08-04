BiCopPar2Beta <- function(family, par, par2 = 0, obj = NULL) {
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
    if (family == 2) 
        stop("The CDF of the t-copula is not implemented.")
    if (!(family %in% c(0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 16, 17, 18, 19, 20, 
                        23, 24, 26, 27, 28, 29, 30, 33, 34, 36, 37, 38, 39, 40, 41,
                        51, 61, 71, 104, 114, 124, 134, 204, 214, 224, 234))) 
        stop("Copula family not implemented.")
    if (family %in% c(7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40, 
                      104, 114, 124, 134, 204, 214, 224, 234) && par2 == 0) 
        stop("For BB1, BB6, BB7, BB8 and Tawn copulas, 'par2' must be set.")
    if (family %in% c(1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36, 41, 51, 
                      61, 71) && length(par) < 1) 
        stop("'par' not set.")
    BiCopCheck(family, par, par2)
    
    ## calculate beta
    4 * BiCopCDF(0.5, 0.5, family, par, par2) - 1
}