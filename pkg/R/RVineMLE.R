RVineMLE <- function(data, RVM, start = RVM$par, start2 = RVM$par2, maxit = 200, max.df = 30, 
                     max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6), BB7 = c(5, 6),  BB8 = c(6, 1)),
                     grad = FALSE, hessian = FALSE, se = FALSE, ...) {
  
    ## sanity checks
    if (is(RVM)[1] != "RVineMatrix") 
        stop("'RVM' has to be an RVineMatrix object.")
    if (maxit <= 0) 
        stop("'maxit' has to be greater than zero.")
    if (max.df <= 2) 
        stop("The upper bound for the degrees of freedom parameter has to be larger than 2.")
    if (!is.list(max.BB)) 
        stop("'max.BB' has to be a list.")
    if (max.BB$BB1[1] < 0.001) 
        stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
    if (max.BB$BB1[2] < 1.001) 
        stop("The upper bound for the second parameter of the BB1 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB6[1] < 1.001) 
        stop("The upper bound for the first parameter of the BB6 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB6[2] < 1.001) 
        stop("The upper bound for the second parameter of the BB6 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB7[1] < 1.001) 
        stop("The upper bound for the first parameter of the BB7 copula should be greater than 1.001 (lower bound for estimation).")
    if (max.BB$BB7[2] < 0.001) 
        stop("The upper bound for the second parameter of the BB7 copula should be greater than 0.001 (lower bound for estimation).")
    if (max.BB$BB8[1] < 1.001) 
        stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
    if (max.BB$BB8[2] < 0.001 || max.BB$BB8[2] > 1) 
        stop("The upper bound for the second parameter of the BB1 copula should be in the interval [0,1].")
    
    ## sanity checks for start parameters
    Matrix <- RVM$Matrix
    if (!all(start %in% c(0, NA))) {
        for (i in 2:dim(Matrix)[1]) {
            for (j in 1:(i - 1)) {
                if ((RVM$family[i, j] == 1 || RVM$family[i, j] == 2) && abs(start[i, j]) >= 1) 
                    stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
                if (RVM$family[i, j] == 2 && start2[i, j] <= 1) 
                    stop("The degrees of freedom parameter of the t-copula has to be larger than 1.")
                if ((RVM$family[i, j] == 3 || RVM$family[i, j] == 13) && start[i, j] <= 0) 
                    stop("The parameter of the Clayton copula has to be positive.")
                if ((RVM$family[i, j] == 4 || RVM$family[i, j] == 14) && start[i, j] < 1) 
                    stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
                if ((RVM$family[i, j] == 6 || RVM$family[i, j] == 16) && start[i, j] <= 1) 
                    stop("The parameter of the Joe copula has to be in the interval (1,oo).")
                if (RVM$family[i, j] == 5 && start[i, j] == 0) 
                    stop("The parameter of the Frank copula has to be unequal to 0.")
                if ((RVM$family[i, j] == 7 || RVM$family[i, j] == 17) && start[i, j] <= 0) 
                    stop("The first parameter of the BB1 copula has to be positive.")
                if ((RVM$family[i, j] == 7 || RVM$family[i, j] == 17) && start2[i, j] < 1) 
                    stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
                if ((RVM$family[i, j] == 8 || RVM$family[i, j] == 18) && start[i, j] <= 0) 
                    stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
                if ((RVM$family[i, j] == 8 || RVM$family[i, j] == 18) && start2[i, j] < 1) 
                    stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
                if ((RVM$family[i, j] == 9 || RVM$family[i, j] == 19) && start[i, j] < 1) 
                    stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
                if ((RVM$family[i, j] == 9 || RVM$family[i, j] == 19) && start2[i, j] <= 0) 
                    stop("The second parameter of the BB7 copula has to be positive.")
                if ((RVM$family[i, j] == 10 || RVM$family[i, j] == 20) && start[i, j] < 1) 
                    stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
                if ((RVM$family[i, j] == 10 || RVM$family[i, j] == 20) && (start2[i, j] <= 0 || start2[i, j] > 1)) 
                    stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
                if ((RVM$family[i, j] == 23 || RVM$family[i, j] == 33) && start[i, j] >= 0) 
                    stop("The parameter of the rotated Clayton copula has to be negative.")
                if ((RVM$family[i, j] == 24 || RVM$family[i, j] == 34) && start[i, j] > -1) 
                    stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
                if ((RVM$family[i, j] == 26 || RVM$family[i, j] == 36) && start[i, j] >= -1) 
                    stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
                if ((RVM$family[i, j] == 27 || RVM$family[i, j] == 37) && start[i, j] >= 0) 
                    stop("The first parameter of the rotated BB1 copula has to be negative.")
                if ((RVM$family[i, j] == 27 || RVM$family[i, j] == 37) && start2[i, j] > -1) 
                    stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
                if ((RVM$family[i, j] == 28 || RVM$family[i, j] == 38) && start[i, j] >= 0) 
                    stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
                if ((RVM$family[i, j] == 28 || RVM$family[i, j] == 38) && start2[i, j] > -1) 
                    stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
                if ((RVM$family[i, j] == 29 || RVM$family[i, j] == 39) && start[i, j] > -1) 
                    stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
                if ((RVM$family[i, j] == 29 || RVM$family[i, j] == 39) && start2[i, j] >= 0) 
                    stop("The second parameter of the rotated BB7 copula has to be negative.")
                if ((RVM$family[i, j] == 30 || RVM$family[i, j] == 40) && start[i, j] > -1) 
                    stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
                if ((RVM$family[i, j] == 30 || RVM$family[i, j] == 40) && (start2[i, j] >= 0 || start2[i, j] < (-1))) 
                    stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
                if ((RVM$family[i, j] == 104 || RVM$family[i, j] == 114 || RVM$family[i, j] == 204 || RVM$family[i, j] == 214) && start[i, j] < 1) 
                    stop("Please choose 'par' of the Tawn copula in [1,oo).")
                if ((RVM$family[i, j] == 104 || RVM$family[i, j] == 114 || RVM$family[i, j] == 204 || RVM$family[i, j] == 214) && (start2[i, j] < 0 ||  start2[i, j] > 1)) 
                    stop("Please choose 'par2' of the Tawn copula in [0,1].")
                if ((RVM$family[i, j] == 124 || RVM$family[i, j] == 134 || RVM$family[i, j] == 224 || RVM$family[i, j] == 234) && start[i, j] > -1) 
                    stop("Please choose 'par' of the Tawn copula in (-oo,-1].")
                if ((RVM$family[i, j] == 124 || RVM$family[i, j] == 134 || RVM$family[i, j] == 224 || RVM$family[i, j] == 234) && (start2[i, j] < 0 ||  start2[i, j] > 1)) 
                    stop("Please choose 'par2' of the Tawn copula in [0,1].")
                
                if (grad == TRUE) {
                    if (RVM$family[i, j] %in% c(7:10, 17:20, 27:30, 37:40, 104, 114, 124, 134, 204, 214, 224, 234)) {
                        message("The combination 'grad=TRUE' and a copula family of the vector (7:10,17:20,27:30,37:40,104,114,124,134,204,214,224,234) is not possible. The algorithm will continue with 'grad=FALSE'.")
                        grad <- FALSE
                    }
                }
            }
        }
    }
    
    ## sanity checks for input data
    data <- as.matrix(data)
    if (any(data > 1) || any(data < 0)) 
        stop("Data has be in the interval [0,1].")
    d <- dim(data)[2]
    T <- dim(data)[1]
    n <- d
    N <- T
    if (n != dim(RVM)) 
        stop("Dimensions of 'data' and 'RVM' do not match.")
    if (T < 2) 
        stop("Number of observations has to be at least 2.")
    
    ## normalization of R-vine matrix
    o <- diag(RVM$Matrix)
    oldRVM <- RVM
    RVM <- normalizeRVineMatrix(RVM)
    data <- data[, o[length(o):1]]
    
    
    n <- dim(RVM)
    N <- dim(data)[1]
    
    ## sequential estimation of start parameters if not provided
    if (all(start == 0)) {
        est_start <- RVineSeqEst(data, RVM, max.df = max.df, max.BB = max.BB)
        start <- est_start$RVM$par
        start2 <- est_start$RVM$par2
    }
    
    ## Position of parameters in R-vine matrix
    posParams <- (RVM$family > 0)
    posParams2 <- (RVM$family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20,
                                     27, 28, 29, 30, 37, 38, 39, 40,
                                     104, 114, 124, 134, 204, 214, 224, 234))
    posParams[is.na(posParams)] <- FALSE
    posParams2[is.na(posParams2)] <- FALSE
    
    ## number of parameters
    nParams <- sum(posParams, na.rm = TRUE)
    nParams2 <- sum(posParams2, na.rm = TRUE)
    
    ## vectors of start parameters and corresponding pair-copula families
    startpar <- double(nParams + nParams2)
    Copula.Types <- RVM$family[posParams]
    startpar[1:nParams] <- start[posParams]
    if (nParams2 > 0) {
        startpar[(nParams + 1):(nParams + nParams2)] <- start2[posParams2]
    }
    
    ## lower and upper bounds
    lb <- double(nParams + nParams2)
    ub <- double(nParams + nParams2)
    for (i in 1:nParams) {
        if (Copula.Types[i] %in% c(1, 2)) {
            # Normal
            lb[i] <- -0.98
            ub[i] <- 0.98
        } else if (Copula.Types[i] %in% c(3, 13)) {
            # clayton
            lb[i] <- 1e-04
            ub[i] <- Inf
        } else if (Copula.Types[i] %in% c(23, 33)) {
            # rotated clayton
            lb[i] <- -Inf
            ub[i] <- -1e-04
        } else if (Copula.Types[i] %in% c(4, 14)) {
            # gumbel
            lb[i] <- 1.0001
            ub[i] <- Inf
        } else if (Copula.Types[i] %in% c(24, 34)) {
            # rotated gumbel
            lb[i] <- -Inf
            ub[i] <- -1.0001
        } else if (Copula.Types[i] == 5) {
            # frank
            lb[i] <- -Inf
            ub[i] <- Inf
        } else if (Copula.Types[i] %in% c(6, 16)) {
            # joe
            lb[i] <- 1.0001
            ub[i] <- Inf
        } else if (Copula.Types[i] %in% c(26, 36)) {
            # rotated joe
            lb[i] <- -Inf
            ub[i] <- -1.0001
        } else if (Copula.Types[i] %in% c(7, 17)) {
            # bb1
            lb[i] <- 0.001
            ub[i] <- max.BB$BB1[1]
        } else if (Copula.Types[i] %in% c(8, 18)) {
            # bb6
            lb[i] <- 1.001
            ub[i] <- max.BB$BB6[1]
        } else if (Copula.Types[i] %in% c(9, 19)) {
            # bb7
            lb[i] <- 1.001
            ub[i] <- max.BB$BB7[1]
        } else if (Copula.Types[i] %in% c(10, 20)) {
            # bb8
            lb[i] <- 1.001
            ub[i] <- max.BB$BB8[1]
        } else if (Copula.Types[i] %in% c(27, 37)) {
            # rotated bb1
            lb[i] <- -max.BB$BB1[1]
            ub[i] <- -0.001
        } else if (Copula.Types[i] %in% c(28, 38)) {
            # rotated bb6
            lb[i] <- -max.BB$BB6[1]
            ub[i] <- -1.001
        } else if (Copula.Types[i] %in% c(29, 39)) {
            # rotated bb7
            lb[i] <- -max.BB$BB7[1]
            ub[i] <- -1.001
        } else if (Copula.Types[i] %in% c(30, 40)) {
            # rotated bb8
            lb[i] <- -max.BB$BB8[1]
            ub[i] <- -1.001
        } else if (Copula.Types[i] %in% c(43, 44)) {
            # double Clayton, Gumbel
            lb[i] <- -0.9999
            ub[i] <- 0.9999
        } else if (Copula.Types[i] %in% c(104, 114, 204, 214)) {
            # Tawn
            lb[i] <- 1.001
            ub[i] <- 20
        } else if (Copula.Types[i] %in% c(124, 134, 224, 234)) {
            # Tawn
            lb[i] <- -20
            ub[i] <- -1.001
        }
    }
    
    if (nParams2 > 0) {
        todo <- which(Copula.Types %in% c(2, 7:10, 17:20, 27:30, 37:40, 104, 114, 124, 134, 204, 214, 224, 234))
        
        for (i in 1:nParams2) {
            if (Copula.Types[todo[i]] == 2) {
                # t
                lb[nParams + i] <- 2.001
                ub[nParams + i] <- max.df
            } else if (Copula.Types[todo[i]] %in% c(7, 17)) {
                # bb1
                lb[nParams + i] <- 1.001
                ub[nParams + i] <- max.BB$BB1[2]
            } else if (Copula.Types[todo[i]] %in% c(8, 18)) {
                # bb6
                lb[nParams + i] <- 1.001
                ub[nParams + i] <- max.BB$BB6[2]
            } else if (Copula.Types[todo[i]] %in% c(9, 19)) {
                # bb7
                lb[nParams + i] <- 0.001
                ub[nParams + i] <- max.BB$BB7[2]
            } else if (Copula.Types[todo[i]] %in% c(10, 20)) {
                # bb8
                lb[nParams + i] <- 0.001
                ub[nParams + i] <- max.BB$BB1[2]
            } else if (Copula.Types[todo[i]] %in% c(27, 37)) {
                # rotated bb1
                lb[nParams + i] <- -max.BB$BB1[2]
                ub[nParams + i] <- -1.001
            } else if (Copula.Types[todo[i]] %in% c(28, 38)) {
                # rotated bb6
                lb[nParams + i] <- -max.BB$BB6[2]
                ub[nParams + i] <- -1.001
            } else if (Copula.Types[todo[i]] %in% c(29, 39)) {
                # rotated bb7
                lb[nParams + i] <- -max.BB$BB7[2]
                ub[nParams + i] <- -1.001
            } else if (Copula.Types[todo[i]] %in% c(30, 40)) {
                # rotated bb8
                lb[nParams + i] <- -max.BB$BB8[2]
                ub[nParams + i] <- -0.001
            } else if (Copula.Types[todo[i]] %in% c(104, 114, 124, 134, 204, 214, 224, 234)) {
                # Tawn
                lb[nParams + i] <- 0.001
                ub[nParams + i] <- 0.99
            }
        }
    }
    
    ## log-likelihood function to be maximized
    optim_LL <- function(parm, data, posParams, posParams2, Copula.Types, start, start2, RVM, calcupdate = NA) {
        
        nParams <- sum(posParams, na.rm = TRUE)
        nParams2 <- sum(posParams2, na.rm = TRUE)
        
        matrixParams <- start
        matrixParams2 <- start2
        
        matrixParams[posParams] <- parm[1:nParams]
        
        if (nParams2 > 0) {
            matrixParams2[posParams2] <- parm[(nParams + 1):(nParams + nParams2)]
        }
        
        ll <- RVineLogLik(data, RVM, par = matrixParams, par2 = matrixParams2)
        
        
        if (is.finite(ll$loglik)) {
            return(ll$loglik)
        } else {
            if (is.na(ll$loglik)) {
                message(parm)
            }
            message(ll$loglik)
            return(-10^305)
        }
    }
    
    ## gradient
    ableitung <- function(parm, data, posParams, posParams2, Copula.Types, start, start2, RVM, calcupdate) {
        nParams <- sum(posParams, na.rm = TRUE)
        nParams2 <- sum(posParams2, na.rm = TRUE)
        
        matrixParams <- start
        matrixParams2 <- start2
        
        matrixParams[posParams] <- parm[1:nParams]
        if (nParams2 > 0) {
            matrixParams2[posParams2] <- parm[(nParams + 1):(nParams + nParams2)]
        }
        
        grad <- RVineGrad(data = data, 
                          RVM = RVM, 
                          par = matrixParams,
                          par2 = matrixParams2, 
                          posParams = posParams)$gradient
        return(grad)
    }
    
    ## default values for parscale (see optim)
    pscale <- numeric()
    for (i in 1:nParams) {
        pscale[i] <- ifelse(Copula.Types[i] %in% c(1, 2, 43, 44), 0.01, 1)
    }
    pscale2 <- numeric()
    if (nParams2 > 0) {
        for (i in 1:nParams2) {
            pscale2[i] <- ifelse(Copula.Types[i] %in% c(104, 114, 124, 134, 204, 214, 224, 234), 0.05, 1)
        }
        pscale <- c(pscale, pscale2)
    }
    
    ## (default) values for control parameters of optim
    ctrl <- list(fnscale = -1, 
                 maxit = maxit,
                 trace = 1,
                 parscale = pscale, 
                 factr = 1e+08)
    ctrl <- modifyList(ctrl, list(...))
    
    ## optimization
    if (all(Copula.Types %in% c(0, 1, 2, 3:6, 13, 14, 16, 23, 24, 26, 33, 34, 36, 43, 44)) && grad == TRUE) {
        
        if (hessian == TRUE || se == TRUE) {
            out1 <- optim(par = startpar, 
                          fn = optim_LL, 
                          gr = ableitung,
                          data = data,
                          posParams = posParams,
                          posParams2 = posParams2, Copula.Types = Copula.Types, 
                          start = start,
                          start2 = start2, 
                          RVM = RVM,
                          method = "L-BFGS-B",
                          control = ctrl,
                          lower = lb, 
                          upper = ub,
                          hessian = TRUE)
        } else {
            out1 <- optim(par = startpar,
                          fn = optim_LL,
                          gr = ableitung,
                          data = data, 
                          posParams = posParams, 
                          posParams2 = posParams2, 
                          Copula.Types = Copula.Types, 
                          start = start, 
                          start2 = start2, 
                          RVM = RVM,
                          method = "L-BFGS-B",
                          control = ctrl, 
                          lower = lb,
                          upper = ub)
        }
    } else {
        if (hessian == TRUE || se == TRUE) {
            out1 <- optim(par = startpar, 
                          fn = optim_LL, 
                          data = data,
                          posParams = posParams,
                          posParams2 = posParams2,
                          Copula.Types = Copula.Types,
                          start = start, 
                          start2 = start2, 
                          RVM = RVM, 
                          calcupdate = NA, 
                          method = "L-BFGS-B",
                          control = ctrl, 
                          lower = lb, 
                          upper = ub,
                          hessian = TRUE)
        } else {
            out1 <- optim(par = startpar,
                          fn = optim_LL,
                          data = data, 
                          posParams = posParams,
                          posParams2 = posParams2,
                          Copula.Types = Copula.Types,
                          start = start, 
                          start2 = start2,
                          RVM = RVM, 
                          calcupdate = NA,
                          method = "L-BFGS-B", 
                          control = ctrl,
                          lower = lb,
                          upper = ub)
        }
    }
    
    ## list for final output
    out <- list()

    out$value <- out1$value
    out$convergence <- out1$convergence
    out$message <- out1$message
    out$counts <- out1$counts
    
    if (hessian == TRUE) 
        out$hessian <- out1$hessian
    
    if (se == TRUE) 
        out1$se <- sqrt((diag(solve(-out1$hessian))))
    
    out$RVM <- oldRVM
    
    kk <- 1
    for (ll in 1:nParams) {
        out1$par[ll] <- out1$par[ll]
        if (Copula.Types[ll] %in% c(2, 7:10, 17:20, 27:30, 37:40,
                                    104, 114, 124, 134, 204, 214, 224, 234)) {
            out1$par[nParams + kk] <- out1$par[nParams + kk]
            kk <- kk + 1
        }
    }
    
    out$RVM$par[posParams] <- out1$par[1:nParams]
    out$RVM$par2[posParams2] <- out1$par[(nParams + 1):(nParams + nParams2)]
    
    if (se == TRUE) {
        out$se <- matrix(0, d, d)
        out$se2 <- matrix(0, d, d)
        out$se[posParams] <- out1$se[1:nParams]
        out$se2[posParams2] <- out1$se[(nParams + 1):(nParams + nParams2)]
    }
    
    return(out)
}
