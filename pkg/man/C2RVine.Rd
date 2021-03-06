\name{C2RVine}            
\alias{C2RVine}

\title{Transform C-Vine to R-Vine Structure}

\description{
This function transforms a C-vine structure from the package CDVine to the corresponding R-vine structure.
}

\usage{
C2RVine(order, family, par, par2 = rep(0, length(family)))
}

\arguments{
  \item{order}{A d-dimensional vector specifying the order of the root nodes in the C-vine.}  
  \item{family}{A d*(d-1)/2 vector of pair-copula families with values\cr
		\code{0} = independence copula \cr
	        \code{1} = Gaussian copula \cr
	        \code{2} = Student t copula (t-copula) \cr
	        \code{3} = Clayton copula \cr
	        \code{4} = Gumbel copula \cr
	        \code{5} = Frank copula \cr
	        \code{6} = Joe copula \cr 
		\code{7} = BB1 copula \cr
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
  \item{par}{A d*(d-1)/2 vector of pair-copula parameters.}
  \item{par2}{A d*(d-1)/2 vector of second pair-copula parameters (optional; default:\cr \code{par2 = rep(0,length(family))}), necessary for the t-, BB1, BB6, BB7, BB8, Tawn type 1 and type 2 copulas.}
}

\value{
An \code{\link{RVineMatrix}} object.
}

\author{Ulf Schepsmeier, Eike Brechmann}

\seealso{\code{\link{RVineMatrix}}, \code{\link{D2RVine}}}

\examples{
# simulate a sample of size 500 from a 4-dimensional C-vine  
# copula model with mixed pair-copulas
# load package CDVine
library(CDVine)
d <- 4
dd <- d*(d-1)/2
order <- 1:d
family <- c(1, 2, 3, 4, 7, 3)
par <- c(0.5, 0.4, 2, 1.5, 1.2, 1.5)
par2 <- c(0, 5, 0, 0, 2, 0)
type <- 1
simdata <- CDVineSim(500, family, par, par2, type)

# determine log-likelihood
out <- CDVineLogLik(simdata, family, par, par2, type)
out$loglik

# transform to R-vine matrix notation
RVM <- C2RVine(order, family, par, par2)

# check that log-likelihood stays the same
out2 <- RVineLogLik(simdata,RVM)
out2$loglik
}

