\name{plot.BiCop}
\alias{plot.BiCop}
\alias{contour.BiCop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting tools for BiCop objects
}
\description{
There are several options for plotting BiCop objects. The density of a bivariate copula density can be visualized as surface/perspective or contour plot. Optionally, the density can be coupled with standard normal margins (default for contour plots). Furthermore, a lambda-plot is available (c.f. \code{\link{BiCopLambda}}).
}
\usage{
\method{plot}{BiCop}(x, type = "surface", margins, size, ...)
\method{contour}{BiCop}(x, margins = "norm", size = 100L, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{BiCop object.}}
  \item{type}{
plot type; either \code{"surface"}, \code{"contour"}, or \code{"lambda"} (partial matching is activated); the latter is only implemented for a few families (c.f., \code{\link{BiCopLambda}}).
}
  \item{margins}{
only relevant for types \code{"contour"} and \code{"surface"}; either \code{"unif"} for the original copula density or \code{"norm"} for the transformed density with standard normal margins (partial matching is activated). Default is \code{"norm"} for \code{type = "contour"}, and \code{"unif"} for \code{type = "surface"}.
} 
  \item{size}{
integer; only relevant for types \code{"contour"} and \code{"surface"}; the plot is based on values on a \eqn{size x size} grid; default is 100 for \code{type = "contour"}, and 25 for \code{type = "surface"}.
}
  \item{\dots}{
optional arguments passed to \code{\link{contour}} or \code{\link{wireframe}}.}
}
\author{
Thomas Nagler
}

\seealso{
\code{\link{BiCop}}, \code{\link{contour}}, \code{\link{wireframe}}
}
\examples{
## construct BiCop object for a Tawn copula
obj <- BiCop(family = 104, par = 2.5, par2 = 0.4)

## plots
plot(obj)  # surface plot of copula density 
contour(obj)  # contour plot with standard normal margins
contour(obj, margins = "unif")  # contour plot of copula density
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }