\name{dduCopula}
\alias{dduCopula}
\alias{ddvCopula}

\title{
Partial Derivatives of Copulas
}
\description{
Similar to \code{\link{dCopula}} and \code{\link{pCopula}} the function \code{dduCopula} evaluates the partial derivative \eqn{\frac{\partial}{\partial u} C(u,v)}{d C(u,v)/du} and the function \code{ddvCopula} evaluates the partial derivative \eqn{\frac{\partial}{\partial v} C(u,v)}{d C(u,v)/dv} of the provided copula. 
}
\usage{
dduCopula(u, copula, ...)
ddvCopula(u, copula, ...)
}

\arguments{
  \item{u}{
Pairs of values for which the partial derivative should be evaluated.
}
  \item{copula}{
The copula object representing the family member of interest.
}
  \item{\dots}{additional arguments can be passed on to the underlying functions.}
}
\value{
A vector of the evaluated partial derivatives of the same length as rows in \code{u}.
}
\author{
Benedikt Graeler
}

\examples{
library(copula)

BB1Cop <- BB1Copula()
BB1CopSmpl <- rCopula(100, BB1Cop)

# conditional probabilities of a Gaussian copula given u
BB1GivenU <- dduCopula(BB1CopSmpl, BB1Cop)

# vs. conditional probabilities of a Gaussian copula given v
BB1GivenV <- ddvCopula(BB1CopSmpl[,c(2,1)], BB1Cop)

plot(BB1GivenU, BB1GivenV)
abline(0,1)
}

\keyword{partial derivative}
\keyword{conditional probabilities}
