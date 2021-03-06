\name{RVineSeqEst}       
\alias{RVineSeqEst}

\title{Sequential Estimation of an R-Vine Copula Model}

\description{
This function sequentially estimates the pair-copula parameters of a d-dimensional R-vine copula model as specified by the corresponding \code{\link{RVineMatrix}} object. 
}

\usage{
RVineSeqEst(data, RVM, method = "mle", se = FALSE, max.df = 30,
            max.BB = list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1)), 
            progress = FALSE, weights = NA)
}

\arguments{
  \item{data}{An N x d data matrix (with uniform margins).}
  \item{RVM}{An \code{\link{RVineMatrix}} object including the structure, the pair-copula families and the pair-copula parameters (if they are known).}  
  \item{method}{Character indicating the estimation method:
    either pairwise maximum likelihood estimation (\code{method = "mle"}; default) or inversion of Kendall's tau (\code{method = "itau"}; see \code{\link{BiCopEst}}.
		For \code{method = "itau"} only one parameter pair-copula families can be used
		(\code{family = 1}, \code{3}, \code{4}, \code{5}, \code{6}, \code{13}, \code{14}, \code{16}, \code{23}, \code{24}, \code{26}, \code{33}, \code{34} or \code{36}).}
  \item{se}{Logical; whether standard errors are estimated (default: \code{se = FALSE}).}
  \item{max.df}{Numeric; upper bound for the estimation of the degrees of freedom parameter of the t-copula
    (default: \code{max.df = 30}; for more details see \code{\link{BiCopEst}}).}
  \item{max.BB}{List; upper bounds for the estimation of the two parameters (in absolute values) of the BB1, BB6, BB7 and BB8 copulas \cr
    (default: \code{max.BB = list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1))}).}
  \item{progress}{Logical; whether the pairwise estimation progress is printed (default: \code{progress = FALSE}).}
  \item{weights}{Numerical; weights for each observation (opitional).}
}

\details{
The pair-copula parameter estimation is performed tree-wise,
i.e., for each R-vine tree the results from the previous tree(s) are used
to calculate the new copula parameters using \code{\link{BiCopEst}}.
}


\value{
  \item{RVM}{\code{\link{RVineMatrix}} object with the sequentially estimated parameters stored in \code{RVM$par} and \code{RVM$par2}.}
  \item{se}{Lower triangular d x d matrix with estimated standard errors of the (first) pair-copula parameters
    for each (conditional) pair defined in the \code{\link{RVineMatrix}} object (if \code{se = TRUE}).}
  \item{se2}{Lower triangular d x d matrix with estimated standard errors of the second parameters for pair-copula families with two parameters
    for each (conditional) pair defined in the \code{\link{RVineMatrix}} object (if \code{se = TRUE}).}
}    

\author{Ulf Schepsmeier, Jeffrey Dissmann}

\seealso{\code{\link{BiCopEst}}, \code{\link{BiCopHfunc}}, \code{\link{RVineLogLik}}, \code{\link{RVineMLE}}, \code{\link{RVineMatrix}}}

\examples{
# define 5-dimensional R-vine tree structure matrix
Matrix <- c(5, 2, 3, 1, 4,
            0, 2, 3, 4, 1,
            0, 0, 3, 4, 1,
            0, 0, 0, 4, 1,
            0, 0, 0, 0, 1)
Matrix <- matrix(Matrix, 5, 5)

# define R-vine pair-copula family matrix
family <- c(0, 1, 3, 4, 4,
            0, 0, 3, 4, 1,
            0, 0, 0, 4, 1,
            0, 0, 0, 0, 3,
            0, 0, 0, 0, 0)
family <- matrix(family, 5, 5)

# define R-vine pair-copula parameter matrix
par <- c(0, 0.2, 0.9, 1.5, 3.9,
         0, 0, 1.1, 1.6, 0.9,
         0, 0, 0, 1.9, 0.5,
         0, 0, 0, 0, 4.8,
         0, 0, 0, 0, 0)
par <- matrix(par, 5, 5)

# define second R-vine pair-copula parameter matrix
par2 <- matrix(0, 5, 5)

# define RVineMatrix object
RVM <- RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("V1", "V2", "V3", "V4", "V5"))

# simulate a sample of size 300 from the R-vine copula model
set.seed(123)
simdata <- RVineSim(300, RVM)

# sequential estimation
RVineSeqEst(simdata, RVM, method = "itau", se = TRUE)
RVineSeqEst(simdata, RVM, method = "mle", se = TRUE)
}
