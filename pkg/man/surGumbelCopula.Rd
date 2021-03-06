\name{surGumbelCopula}
\alias{surGumbelCopula}
\alias{r90GumbelCopula}
\alias{r270GumbelCopula}

\title{
Survival and Rotated Gumbel Copulas
}
\description{
These are wrappers to functions from \code{\link{VineCopula-package}}
}
\usage{
surGumbelCopula(param)
r90GumbelCopula(param)
r270GumbelCopula(param)
}
\arguments{
  \item{param}{
A single parameter defining the Copula.
}
}
\value{
An object of class \code{\linkS4class{surGumbelCopula}}, \code{\linkS4class{r90GumbelCopula}} or \code{\linkS4class{r270GumbelCopula}} respectively.
}
\author{
Benedikt Graeler
}
\examples{
library(copula)

persp(surGumbelCopula(1.5), dCopula, zlim = c(0,10))
persp(r90GumbelCopula(-1.5), dCopula, zlim = c(0,10))
persp(r270GumbelCopula(-1.5), dCopula, zlim = c(0,10))
}
\keyword{copula}
