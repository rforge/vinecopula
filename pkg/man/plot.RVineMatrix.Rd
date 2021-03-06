\name{plot.RVineMatrix}
\alias{plot.RVineMatrix}
\alias{contour.RVineMatrix}

\title{
Plotting \code{RVineMatrix} objects.
}


\description{
There are two plotting generics for \code{RVineMatrix} objects. \code{plot.RVineMatrix} plots one or all trees of a given R-vine copula model. Edges can be labeld with information about the corresponding pair-copula. \code{contour.RVineMatrix} produces a matrix of contour plots (using \code{\link[VineCopula:plot.BiCop]{plot.BiCop}}).
}


\usage{
\method{plot}{RVineMatrix}(x, tree = "ALL", type = 0, edge.labels = NULL, legend.pos = "bottomleft",
     interactive = FALSE,  ...)
\method{contour}{RVineMatrix}(x, tree = "ALL", xylim = NULL, cex.nums = 1, ...)
}


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
  \item{legend.pos}{the \code{x} argument for \code{\link[graphics:legend]{legend}}.}
  \item{interactive}{logical; if TRUE, the user is asked to adjust the positioning of 
  vertices with his mouse.}
  \item{xylim}{numeric vector of length 2; sets \code{xlim} and \code{ylim} for the contours}
  \item{cex.nums}{numeric; expansion factor for font of the numbers.}
  \item{\dots}{
Arguments passed to \code{\link[network:plot.network]{plot.network}} or \code{\link[VineCopula:plot.BiCop]{plot.BiCop}} respectively.}
}


\details{
If you want the contour boxes to be perfect sqaures, the plot height should be \code{1.25/length(tree)*(d - min(tree))} times the plot width.
}


\author{
Thomas Nagler, Nicole Barthel
}


\seealso{
\code{\link[VineCopula:RVineMatrix]{RVineMatrix}},
\code{\link[network:plot.network]{plot.network}},
\code{\link[VineCopula:plot.BiCop]{plot.BiCop}},
\code{\link[VineCopula:BiCopName]{BiCopName}},
\code{\link[graphics:legend]{legend}}
}


\examples{
## build vine model
strucmat <- matrix(c(3,   1, 2, 0, 2, 1, 0, 0, 1), 3, 3)
fammat   <- matrix(c(0,   1, 6, 0, 0, 3, 0, 0, 0), 3, 3)
parmat   <- matrix(c(0, 0.3, 3, 0, 0, 1, 0, 0, 0), 3, 3)
par2mat  <- matrix(c(0,   0, 0, 0, 0, 0, 0, 0, 0), 3, 3)
RVM  <- RVineMatrix(strucmat, fammat, parmat, par2mat)

# plot trees
plot(RVM)

## build new model
# simulate from previous model
u <- RVineSim(500, RVM)
colnames(u) <- c("X", "Y", "Z")

# estimate new model
RVM2 <- RVineStructureSelect(u)

\dontrun{
# plot new model with legend
plot(RVM2, type = 1)

# show contour plots
contour(RVM2)
}
}

\keyword{plot}