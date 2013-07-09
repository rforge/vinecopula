\name{cor2pcor}
\Rdversion{1.1}
\alias{pcor2cor}
\alias{cor2pcor}

\title{correlations to partial correlations and vice versa for R-vines}
\description{
correlations to partial correlations and vice versa for R-vines
(C vine, D vine or general R vine) 
}

\usage{
cor2pcor(corMat, RVM)
pcor2cor(RVM)
}
\arguments{
  \item{corMat}{correlation matrix}
  \item{RVM}{\code{\link{RVineMatrix}} defining only the R-vine structure for \code{cor2pcor} and providing as well the partial correlations for \code{pcor2cor}.}
}
\value{
  \item{cor}{correlation matrix}
  \item{RVM}{RVineMatrix with transformed partial correlations}
}

\examples{
corMat <- matrix(c(1.00, 0.17, 0.15, 0.14, 0.13,
                   0.17, 1.00, 0.30, 0.28, 0.05,
                   0.15, 0.30, 1.00, 0.17, 0.05,
                   0.14, 0.28, 0.17, 1.00, 0.04,
                   0.13, 0.05, 0.05, 0.04, 1.00),5,5)

RVM <- list()
RVM$Matrix <- matrix(c(5,4,3,2,1,
                       0,4,3,2,1,
                       0,0,3,2,1,
                       0,0,0,2,1,
                       0,0,0,0,1),5,5)
class(RVM) <- "RVineMatrix"
cor2pcor(corMat, RVM)
 
newRVM <- cor2pcor(corMat, RVM)
newRVM$family <- matrix(1,5,5)
newRVM$par2 <- matrix(0,5,5)

# simulate a 5-dim normal from the vine and calculate its correlations
round(cor(qnorm(RVineSim(10000, newRVM))),2)

pcor2cor(cor2pcor(corMat, RVM))
pcor2cor(newRVM)
}

\keyword{vine}
\keyword{partial correlation}
