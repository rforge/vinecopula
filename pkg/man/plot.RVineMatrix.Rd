\name{plot.RVineMatrix}
\alias{plot.RVineMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting \code{RVineMatrix} objects.
}
\description{
This function plots one or all trees of a given R-vine copula model. Edges can be labeld with information about the corresponding pair-copula.
}
\usage{
\method{plot}{RVineMatrix}(x, tree = "ALL", type = 0, edge.labels = NULL, interactive = FALSE,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{RVineMatrix} object.}
  \item{tree}{\code{"ALL"} or integer vector; specifies which trees are plotted.}
  \item{type}{integer; specifies how to make use of variable names: \cr
  \code{0} = variable names are ignored, \cr 
  \code{1} = variable names are used to annotate vertices, \cr
  \code{2} = uses numbers in plot and adds a legend for variable names.}
  \item{edge.labels}{character; either a vector of edge labels
  or one of the following: \cr
  \code{"family"} = pair-copula family abbreviation (see \code{\link[VineCopula:BiCopName]{BiCopName}}), \cr
  \code{"par"} = pair-copula parameters, \cr
  \code{"tau"} = pair-copula Kendall's tau (by conversion of parameters) \cr
  \code{"family-par"} = pair-copula family and parameters \cr
  \code{"family-tau"} = pair-copula family and Kendall's tau.
  }
  \item{interactive}{logical; if TRUE, the user is asked to adjust the positioning of 
  vertices with his mouse.}
  \item{\dots}{
Arguments passed to \code{\link[network:plot.network]{plot.network}}.
}
}
\author{
Thomas Nagler
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[VineCopula:RVineMatrix]{RVineMatrix}},
\code{\link[network:plot.network]{plot.network}}
\code{\link[VineCopula:BiCopName]{BiCopName}}
}
\examples{
## build vine model
strucmat <- matrix(c(3,1,2,0,2,1,0,0,1),3,3)
fammat <- matrix(c(0,1,6,0,0,3,0,0,0),3,3)
parmat <- matrix(c(0,0.3,3,0,0,1,0,0,0),3,3)
par2mat <- matrix(c(0,0,0,0,0,0,0,0,0),3,3)
RVM  <- RVineMatrix(Matrix=strucmat, family=fammat, par=parmat, par2=par2mat)

# plot trees
plot(RVM)

## build new model
# simulate from previous model
u <- RVineSim(500, RVM)
colnames(u) <- c("A", "B", "C")

# estimate new model
RVM2 <- RVineStructureSelect(u)

# plot new model with variable names ...
plot(RVM2, type = 1)

# annotate edge with pair-copula family and parameter
plot(RVM2, type = 1, edge.labels = "family-par")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}