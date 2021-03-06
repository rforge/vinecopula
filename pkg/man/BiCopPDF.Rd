\name{BiCopPDF}         
\alias{BiCopPDF}

\title{Density of a Bivariate Copula}

\description{
This function evaluates the probability density function (PDF) of a given parametric bivariate copula.
}

\usage{
BiCopPDF(u1, u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE)
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
		\code{7} = BB1 copula \cr/
		\code{8} = BB6 copula \cr
		\code{9} = BB7 copula \cr
		\code{10} = BB8 copula \cr
		\code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
		\code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
		\code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr 
		\code{17} = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
		\code{18} = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
		\code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
		\code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
		\code{23} = rotated Clayton copula (90 degrees) \cr
		\code{24} = rotated Gumbel copula (90 degrees) \cr
		\code{26} = rotated Joe copula (90 degrees) \cr
		\code{27} = rotated BB1 copula (90 degrees) \cr
		\code{28} = rotated BB6 copula (90 degrees) \cr
		\code{29} = rotated BB7 copula (90 degrees) \cr
		\code{30} = rotated BB8 copula (90 degrees) \cr
		\code{33} = rotated Clayton copula (270 degrees) \cr
		\code{34} = rotated Gumbel copula (270 degrees) \cr
		\code{36} = rotated Joe copula (270 degrees) \cr
		\code{37} = rotated BB1 copula (270 degrees) \cr
		\code{38} = rotated BB6 copula (270 degrees) \cr
		\code{39} = rotated BB7 copula (270 degrees) \cr
		\code{40} = rotated BB8 copula (270 degrees) \cr
    \code{104} = Tawn type 1 copula \cr
    \code{114} = rotated Tawn type 1 copula (180 degrees) \cr
    \code{124} = rotated Tawn type 1 copula (90 degrees)  \cr
    \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
    \code{204} = Tawn type 2 copula  \cr
    \code{214} = rotated Tawn type 2 copula (180 degrees) \cr
    \code{224} = rotated Tawn type 2 copula (90 degrees)  \cr
    \code{234} = rotated Tawn type 2 copula (270 degrees) \cr
		}
  \item{par}{numeric; single number or vector of size \code{length(u1)}; copula parameter.}
  \item{par2}{numeric; single number or vector of size \code{length(u1)}; second parameter for bivariate copulas with two parameters (t, BB1, BB6, BB7, BB8, Tawn type 1 and type 2; default: \code{par2 = 0}). \code{par2} should be an positive integer for the Students's t copula \code{family = 2}.}
  \item{obj}{\code{BiCop} object containing the family and parameter specification.}
\item{check.pars}{logical; default is \code{TRUE}; if \code{FALSE}, checks for family/parameter-consistency are ommited (should only be used with care).}
}

\value{
A numeric vector of the bivariate copula density 
\itemize{
    \item{of the copula \code{family}}
    \item{with  parameter(s) \code{par}, \code{par2}}
    \item{evaluated at \code{u1} and \code{u2}.}
}
}

\details{
If the family and parameter specification is stored in a \code{\link{BiCop}} object \code{obj}, the alternative version \cr
\preformatted{BiCopPDF(u1, u2, obj)} 
can be used.
}

\author{Eike Brechmann}

\seealso{\code{\link{BiCopCDF}}, \code{\link{BiCopHfunc}}, \code{\link{BiCopSim}}, \code{\link{BiCop}}}

\examples{
## simulate from a bivariate t-copula
simdata <- BiCopSim(300, 2, -0.7, par2 = 4)

## evaluate the density of the bivariate t-copula
u1 <- simdata[,1]
u2 <- simdata[,2]
BiCopPDF(u1, u2, 2, -0.7, par2 = 4)

## estimate a bivariate copula from the data and evaluate its PDF
cop <- BiCopSelect(u1, u2)
round(BiCopPDF(u1, u2, cop), 8)
}

