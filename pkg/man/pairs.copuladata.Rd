\name{pairs.copuladata}
\alias{pairs.copuladata}

\title{Pairs Plot of Copula Data}

\description{
  This function provides pair plots for copula data. It plots bivariate contour plots on the lower panel, scatter plots and correlations on the upper panel and histograms on the diagonal panel.
}

\usage{
\method{pairs}{copuladata}(x, labels = names(x), ...,
                  label.pos = 0.85, cex.labels = 1, gap = 0, axes = FALSE,
                  pch = ".", col = "grey", cex.points = 1,
                  method.cor = "kendall", col.cor = "red", digits.cor = 2, cex.cor = 1,
                  bw = 2, size = 100, levels = seq(0.01, 0.2, length.out = 30),
                  margins = "norm", margins.par = 0, xylim = NA,
                  col.contour = terrain.colors(length(levels)),
                  col.hist = "grey")
}

\arguments{
  \item{x}{\code{copuladata} object.}
  \item{labels}{variable names/labels.}
  \item{\dots}{other graphical parameters (see \code{\link[graphics]{par}}).}
  \item{label.pos}{y position of labels in the diagonal panel; default: \code{label.pos = 0.85}.}
  \item{cex.labels}{magnification to be used for the labels of the diagonal panel; default: \code{cex.labels = 1}.}
  \item{gap}{distance between subplots, in margin lines; default: \code{gap = 0}.}
  \item{axes}{a logical value indicating whether both axes should be drawn on the plot; default: \code{axes = FALSE}.}
  \item{pch}{plotting characters/symbols to be used for the points of the scatter plots; default: \code{pch = "."}.}
  \item{col}{colour to be used for the points of the scatter plots; default: \code{col = "grey"}.}
  \item{cex.points}{magnification to be used for the points of the scatter plots; default: \code{cex.points = 1}.}
  \item{method.cor}{a character string indicating which correlation coefficients are computed. One of \code{pearson}, \code{kendall} (default), or \code{spearman}}
  \item{col.cor}{colour to be used for the correlation coefficients in the scatter plots; default: \code{col.cor = "red"}.}
  \item{digits.cor}{digits to be used for the correlation coefficients in the scatter plots; default: \code{digits.cor = 2}.}
  \item{cex.cor}{magnification to be used for the correlation coefficients in the scatter plots; default: \code{cex.cor = 1}.}
  \item{bw}{bandwidth to be used for the contour plots (smoothing factor; default: \code{bw = 1}).}
  \item{size}{number of grid points to be used for the contour plots; default: \code{size = 100}.}
  \item{levels}{vector of contour levels to be used for the contour plots.
                For Gaussian, Student t or exponential margins the default value (\code{levels = seq(0.01, 0.2, length.out = 30)}) typically is a good choice.
                For uniform margins we recommend\cr
                \code{levels = seq(0.1, 1.5, length.out = 30)}\cr
                and for Gamma margins\cr
                \code{levels = seq(0.005, 0.09, length.out = 30)}.}
  \item{margins}{character; margins for the contour plots. Possible margins are:\cr
                 \code{"norm"} = standard normal margins (default)\cr
                 \code{"t"} = Student t margins with degrees of freedom as specified by \code{margins.par}\cr
                 \code{"gamma"} = Gamma margins with shape and scale as specified by \code{margins.par}\cr
                 \code{"exp"} = Exponential margins with rate as specified by \code{margins.par}\cr
                 \code{"unif"} = uniform margins}
  \item{margins.par}{parameter(s) of the distribution of the margins (of the contour plots) if necessary (default: \code{margins.par = 0}), i.e., 
                     \itemize{
                       \item a positive real number for the degrees of freedom of Student t margins (see \code{\link{dt}}),
                        \item a 2-dimensional vector of positive real numbers for the shape and scale parameters of Gamma margins (see \code{\link{dgamma}}),
                        \item a positive real number for the rate parameter of exponential margins (see \code{\link{dexp}}). 
                     }}
  \item{xylim}{2-dimensional vector of the x- and y-limits to be used for the contour plots.
               By default (\code{xylim = NA}) standard limits for the selected margins are used.} 
  \item{col.contour}{colour to be used for the contour plots; default: \code{col.contour = terrain.colors(length(levels))}.}
  \item{col.hist}{colour to be used for histograms of the diagonal panel; default: \code{col.hist = "grey"}.}
}

\author{Tobias Erhardt}

\seealso{\code{\link[graphics]{pairs}}, \code{\link{as.copuladata}}, \code{\link{BiCopMetaContour}}}

\examples{
data(daxreturns)

data <- as.copuladata(daxreturns)

## pairs plot with default settings
pairs(data[1:5])

## pairs plot with custom settings
nlevels <- 20
pairs(data[1:5], cex.labels = 2, gap = 1,
      pch = 20, col = "black", cex.points = 0.5,
      method.cor = "spearman", col.cor = "green",
      digits.cor = 3, cex.cor = 1.5,
      bw = 1.5, levels = seq(0.01, 0.2, length.out = nlevels),
      margins = "t", margins.par = 5, xylim = c(-1,2),
      col.contour = heat.colors(nlevels), col.hist = "white")
}
