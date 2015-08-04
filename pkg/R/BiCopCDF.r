BiCopCDF <- function(u1, u2, family, par, par2 = 0, obj = NULL) {
    ## sanity checks for u1, u2
    if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
        stop("u1 and/or u2 are not set or have length zero.")
    if (any(u1 > 1) || any(u1 < 0)) 
        stop("Data has be in the interval [0,1].")
    if (any(u2 > 1) || any(u2 < 0)) 
        stop("Data has be in the interval [0,1].")
    if (length(u1) != length(u2)) 
        stop("Lengths of 'u1' and 'u2' do not match.")
    
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
    if (is.na(family) | is.na(par))
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
    
    res <- rep(NA, length(u1))
    
    ## CDFs for the different families
    if (family == 0) {
        res <- u1 * u2
    } else if (family == 1) {
        cdf <- function(u, v) pmvnorm(upper = c(qnorm(u), qnorm(v)), 
                                      corr = matrix(c(1,   par, par, 1), 2, 2))
        res <- mapply(cdf, u1, u2, SIMPLIFY = TRUE)
        # }else if(family == 2){ par2=round(par2) cdf = function(u,v)
        # pmvt(upper=c(qt(u,df=par2),qt(v,df=par2)), corr=matrix(c(1,par,par,1),2,2),
        # df=par2) res = mapply(cdf, u1, u2, SIMPLIFY=TRUE)
    } else if (family %in% c(3:10, 41)) {
        res <- .C("archCDF",
                  as.double(u1), 
                  as.double(u2), 
                  as.integer(length(u1)), 
                  as.double(c(par, par2)), 
                  as.integer(family), 
                  as.double(rep(0, length(u1))), 
                  PACKAGE = "VineCopula")[[6]]
    } else if (family %in% c(13, 14, 16:20, 51)) {
        res <- u1 + u2 - 1 + .C("archCDF",
                                as.double(1 - u1),
                                as.double(1 - u2), 
                                as.integer(length(u1)),
                                as.double(c(par, par2)), 
                                as.integer(family - 10),
                                as.double(rep(0, length(u1))),
                                PACKAGE = "VineCopula")[[6]]
    } else if (family %in% c(23, 24, 26:30, 61)) {
        res <- u2 - .C("archCDF", 
                       as.double(1 - u1),
                       as.double(u2), 
                       as.integer(length(u1)), 
                       as.double(c(-par, -par2)),
                       as.integer(family - 20),
                       as.double(rep(0, length(u1))),
                       PACKAGE = "VineCopula")[[6]]
    } else if (family %in% c(33, 34, 36:40, 71)) {
        res <- u1 - .C("archCDF",
                       as.double(u1), 
                       as.double(1 - u2),
                       as.integer(length(u1)), 
                       as.double(c(-par, -par2)),
                       as.integer(family - 30),
                       as.double(rep(0, length(u1))), 
                       PACKAGE = "VineCopula")[[6]]
    } else if (family %in% c(104, 114, 124, 134, 204, 214, 224, 234)) {
        # maybe replace by C-Code auxiliary functions ###
        ta <- function(t, par, par2, par3) {
            (par2 * t)^par + (par3 * (1 - t))^par
        }
        ######## Pickands A
        A <- function(t, par, par2, par3) {
            (1 - par3) * (1 - t) + (1 - par2) * t + ta(t, par, par2, par3)^(1/par)
        }
        
        w <- function(u, v) {
            log(v)/log(u * v)
        }
        C <- function(u, v, par, par2, par3) {
            (u * v)^A(w(u, v), par, par2, par3)
        }
        
        if (family == 104) {
            par3 <- 1
            res <- C(u1, u2, par, par2, par3)
        } else if (family == 114) {
            par3 <- 1
            res <- u1 + u2 - 1 + C(1 - u1, 1 - u2, par, par2, par3)
        } else if (family == 124) {
            par3 <- 1
            res <- u2 - C(1 - u1, u2, -par, par2, par3)
        } else if (family == 134) {
            par3 <- 1
            res <- u1 - C(u1, 1 - u2, -par, par2, par3)
        } else if (family == 204) {
            par3 <- par2
            par2 <- 1
            res <- C(u1, u2, par, par2, par3)
        } else if (family == 214) {
            par3 <- par2
            par2 <- 1
            res <- u1 + u2 - 1 + C(1 - u1, 1 - u2, par, par2, par3)
        } else if (family == 224) {
            par3 <- par2
            par2 <- 1
            res <- u2 - C(1 - u1, u2, -par, par2, par3)
        } else if (family == 234) {
            par3 <- par2
            par2 <- 1
            res <- u1 - C(u1, 1 - u2, -par, par2, par3)
        }
    }
    
    ## return results
    res
}
