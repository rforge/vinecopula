\name{RVineStdError}
\alias{RVineStdError}

\title{Standard Errors of an R-Vine Copula Model}

\description{
This function calculates the standard errors of a d-dimensional R-vine copula model given the Hessian matrix.
}

\usage{
RVineStdError(hessian, RVM)
}

\arguments{
  \item{hessian}{The Hessian matrix of the given R-vine.}
  \item{RVM}{An \code{\link{RVineMatrix}} object including the structure, the pair-copula families, and the parameters.}
}


\value{
  \item{se}{The calculated standard errors for the first parameter matrix. The entries are ordered with respect to the ordering of the \code{RVM$par} matrix.}
  \item{se2}{The calculated standard errors for the second parameter matrix.}
}

\note{
The negative Hessian matrix should be positive semidefinite. Otherwise NAs will be returned in some entries and the non-NA entries may be wrong.
If the negaive Hessian matrix is negative definite, then one could try a near positive matrix. The package \code{Matrix} provides a function called
\code{nearPD} to estimate a matrix which is positive definite and close to the given matrix.
}


\references{
Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
Selecting and estimating regular vine copulae and application to financial returns.
Computational Statistics & Data Analysis, 59 (1), 52-69.

Schepsmeier, U. and J. Stoeber (2012).
Derivatives and Fisher information of bivariate copulas.
Submitted for publication.
\url{http://mediatum.ub.tum.de/node?id=1106541}.

Stoeber, J. and U. Schepsmeier (2012).
Is there significant time-variation in multivariate copulas?
Submitted for publication.
\url{http://arxiv.org/abs/1205.4841}.
}

\author{Ulf Schepsmeier, Jakob Stoeber}

\seealso{\code{\link{BiCopDeriv}}, \code{\link{BiCopDeriv2}}, \code{\link{BiCopHfuncDeriv}}, \code{\link{BiCopHfuncDeriv2}}, \cr
\code{\link{RVineMatrix}}, \code{\link{RVineHessian}}, \code{\link{RVineGrad}}}

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

# compute the Hessian matrix of the first row of the data
out2 <- RVineHessian(simdata,RVM)

# get the standard errors
RVineStdError(out2$hessian, RVM)
}
