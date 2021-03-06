\name{BiCopHfuncDeriv}           
\alias{BiCopHfuncDeriv}

\title{Derivatives of the h-Function of a Bivariate Copula}

\description{
This function evaluates the derivative of a given conditional parametric bivariate copula (h-function) with respect to its parameter(s) or one of its arguments.
}

\usage{
BiCopHfuncDeriv(u1, u2, family, par, par2 = 0, deriv = "par",
                obj = NULL, check.pars = TRUE)
}

\arguments{
  \item{u1,u2}{numeric vectors of equal length with values in [0,1].}
  \item{family}{integer; single number or vector of size \code{length(u1)}; defines the bivariate copula family: \cr
		\code{0} = independence copula \cr
	        \code{1} = Gaussian copula \cr
	        \code{2} = Student t copula (t-copula) \cr
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
  \item{par}{numeric; single number or vector of size \code{length(u1)}; copula parameter.}
  \item{par2}{integer; single number or vector of size \code{length(u1)}; second parameter for the t-Copula; default is \code{par2 = 0}, should be an positive integer for the Students's t copula \code{family = 2}.}
  \item{deriv}{Derivative argument \cr
  	\code{"par"} = derivative with respect to the first parameter (default)\cr
  	\code{"par2"} = derivative with respect to the second parameter (only available for the t-copula) \cr
  	\code{"u2"} = derivative with respect to the second argument \code{u2} \cr
  	}
  \item{obj}{\code{BiCop} object containing the family and parameter specification.}
  \item{check.pars}{logical; default is \code{TRUE}; if \code{FALSE}, checks for family/parameter-consistency are ommited (should only be used with care).}  
}

\value{
A numeric vector of the conditional bivariate copula derivative
\itemize{
    \item{of the copula \code{family}}
    \item{with  parameter(s) \code{par}, \code{par2}}
    \item{with respect to \code{deriv}}
    \item{evaluated at \code{u1} and \code{u2}.}
}
}

\details{
If the family and parameter specification is stored in a \code{\link{BiCop}} object \code{obj}, the alternative version \cr
\preformatted{BiCopHfuncDeriv(u1, u2, obj, deriv = "par")} 
can be used.
}

\references{
Schepsmeier, U. and J. Stoeber (2012).
Derivatives and Fisher information of bivariate copulas.
Statistical Papers.
\url{http://link.springer.com/article/10.1007/s00362-013-0498-x}.
}

\author{Ulf Schepsmeier}

\seealso{\code{\link{RVineGrad}}, \code{\link{RVineHessian}}, \code{\link{BiCopDeriv2}}, \code{\link{BiCopDeriv2}}, \code{\link{BiCopHfuncDeriv}}, \code{\link{BiCop}}}

\examples{
## simulate from a bivariate t-copula
simdata <- BiCopSim(300, 2, -0.7, par2 = 4)

## derivative of the conditional bivariate t-copula 
## with respect to the first parameter
u1 <- simdata[,1]
u2 <- simdata[,2]
BiCopHfuncDeriv(u1, u2, 2, -0.7, par2 = 4, deriv = "par")

## estimate a bivariate copula from the data and
## evaluate its derivative w.r.t. the parameter
cop <- BiCopEst(u1, u2, family = 2)
BiCopHfuncDeriv(u1, u2, cop, deriv = "par")
}

