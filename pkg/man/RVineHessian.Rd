\name{RVineHessian}        
\alias{RVineHessian}                

\title{Hessian Matrix of the Log-Likelihood of an R-Vine Copula Model}

\description{
This function calculates the Hessian matrix of the log-likelihood of a d-dimensional R-vine copula model with respect to the copula parameter and evaluates it on a given copula data set.
}

\usage{
RVineHessian(data, RVM)
}

\arguments{
  \item{data}{An N x d data matrix (with uniform margins).}
  \item{RVM}{An \code{\link{RVineMatrix}} object including the structure, the pair-copula families, and the parameters. \cr
  Only the following copula families are allowed in \code{RVM$family} \cr
		\code{0} = independence copula \cr
		\code{1} = Gaussian copula \cr
		\code{2} = Student t copula (t-copula) (WARNING: see details)\cr
		\code{3} = Clayton copula \cr
		\code{4} = Gumbel copula \cr
		\code{5} = Frank copula \cr
		\code{6} = Joe copula \cr 
		\code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
		\code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
		\code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr 
		\code{23} = rotated Clayton copula (90 degrees) \cr
		\code{24} = rotated Gumbel copula (90 degrees) \cr
		\code{26} = rotated Joe copula (90 degrees) \cr
		\code{33} = rotated Clayton copula (270 degrees) \cr
		\code{34} = rotated Gumbel copula (270 degrees) \cr
		\code{36} = rotated Joe copula (270 degrees) \cr
  }
}

\note{
The Hessian matrix is not available for R-vine copula models with two parameter Archimedean copulas, i.e. BB1, BB6, BB7, BB8 and their rotated versions.
}

\value{
  \item{hessian}{The calculated Hessian matrix of the log-likelihood value of the R-vine copula model.}
  \item{der}{The product of the gradient vector with its transposed version.}
}


\references{
Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
Selecting and estimating regular vine copulae and application to financial returns.
Computational Statistics & Data Analysis, 59 (1), 52-69.

Schepsmeier, U. and J. Stoeber (2012).
Derivatives and Fisher information of bivariate copulas.
Statistical Papers.
\url{http://link.springer.com/article/10.1007/s00362-013-0498-x}.

Stoeber, J. and U. Schepsmeier (2013).
Estimating standard errors in regular vine copula models. Computational Statistics, 1-29
\url{http://link.springer.com/article/10.1007/s00180-013-0423-8#}.
}

\author{Ulf Schepsmeier, Jakob Stoeber}

\seealso{\code{\link{BiCopDeriv}}, \code{\link{BiCopDeriv2}}, \code{\link{BiCopHfuncDeriv}}, \code{\link{BiCopHfuncDeriv2}}, \cr
\code{\link{RVineMatrix}}, \code{\link{RVineMLE}}, \code{\link{RVineGrad}}}

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
out2 <- RVineHessian(simdata[1,], RVM)
out2$hessian
}
