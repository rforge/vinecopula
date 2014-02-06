####################
##  vine copulas  ##
####################

validVineCopula = function(object) {
  dim <- object@dimension
  if( dim <= 2)
    return("Number of dimension too small (>2).")
  if(length(object@copulas)!=(dim*(dim-1)/2))
    return("Number of provided copulas does not match given dimension.")
  if(!any(unlist(lapply(object@copulas,function(x) is(x,"copula")))))
    return("Not all provided copulas in your list are indeed copulas.")
  return (TRUE)
}

setOldClass("RVineMatrix")

setClass("vineCopula",
         representation = representation(copulas="list", dimension="integer", 
                                         RVM="RVineMatrix"),
         prototype = prototype(RVM=structure(list(),class="RVineMatrix")),
         validity = validVineCopula,
         contains = list("copula")
)

# constructor
vineCopula <- function (RVM, type="CVine") { # RVM <- 4L
  if (is.integer(RVM)) {# assuming a dimension
    stopifnot(type %in% c("CVine","DVine"))
    if (type=="CVine")
      RVM <- C2RVine(1:RVM,rep(0,RVM*(RVM-1)/2),rep(0,RVM*(RVM-1)/2))
    if (type=="DVine")
      RVM <- D2RVine(1:RVM,rep(0,RVM*(RVM-1)/2),rep(0,RVM*(RVM-1)/2))
  }
  
  stopifnot(class(RVM)=="RVineMatrix") 
  
  ltr <- lower.tri(RVM$Matrix)
  copDef <- cbind(RVM$family[ltr], RVM$par[ltr], RVM$par2[ltr])
  copulas <- rev(apply(copDef,1, function(x) { 
                                   copulaFromFamilyIndex(x[1],x[2],x[3])
                                 }))
  
  new("vineCopula", copulas=copulas, dimension = as.integer(nrow(RVM$Matrix)),
      RVM=RVM, parameters = numeric(),
      param.names = character(), param.lowbnd = numeric(), 
      param.upbnd = numeric(), fullname = paste("RVine copula family."))
}

showVineCopula <- function(object) {
  dim <- object@dimension
  cat(object@fullname, "\n")
  cat("Dimension: ", dim, "\n")
  cat("Represented by the following",dim*(dim-1)/2, "copulas:\n")
  for (i in 1:length(object@copulas)) {
    cat("  ", class(object@copulas[[i]]), "with parameter(s)", 
        object@copulas[[i]]@parameters, "\n")
  }
}

setMethod("show", signature("vineCopula"), showVineCopula)

## density ##

dRVine <- function(u, copula, log=FALSE) {
  RVM <- copula@RVM
  vineLoglik <- RVineLogLik(u, RVM, separate=TRUE)$loglik
  if(log)
    return(vineLoglik)
  else
    return(exp(vineLoglik))
}

setMethod("dCopula", signature("numeric","vineCopula"), 
          function(u, copula, log, ...) {
            dRVine(matrix(u, ncol=copula@dimension), copula, log, ...)
          })
setMethod("dCopula", signature("matrix","vineCopula"), dRVine)
setMethod("dCopula", signature("data.frame","vineCopula"), 
          function(u, copula, log, ...) {
            dRVine(as.matrix(u), copula, log, ...)
          })

## simulation

rRVine <- function(n, copula) {
  RVM <- copula@RVM
  RVineSim(n, RVM)
}

setMethod("rCopula", signature("numeric","vineCopula"), rRVine)

# fitting using RVine
fitVineCop <- function(copula, data, method) {
  stopifnot(copula@dimension==ncol(data))
  if(!is.null(method[["familyset"]]))
    familyset <- method[["familyset"]]
  else
    familyset <- NA
  if("StructureSelect" %in% method)
    vineCop <- vineCopula(RVineStructureSelect(data, familyset, indeptest="indeptest" %in% method))
  else
    vineCop <- vineCopula(RVineCopSelect(data, familyset, copula@RVM$Matrix, 
                                         indeptest="indeptest" %in% method))
  
  return(new("fitCopula", estimate = vineCop@parameters, var.est = matrix(NA), 
             method = sapply(method,paste,collapse=", "), 
             loglik = RVineLogLik(data, vineCop@RVM)$loglik,
             fitting.stats=list(convergence = as.integer(NA)),
             nsample = nrow(data), copula=vineCop))
}

setMethod("fitCopula", signature=signature("vineCopula"), fitVineCop) 