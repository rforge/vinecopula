\name{pairs.copuladata}
\alias{pairs.copuladata}

\title{Pairs Plot of Copula Data}

\description{
  This function provides pair plots for copula data. Using default setting it plots bivariate contour plots on the lower panel, scatter plots and correlations on the upper panel and histograms on the diagonal panel.
}

\usage{
\method{pairs}{copuladata}(x, labels = names(x), ...,
      lower.panel = lp.copuladata,
      upper.panel = up.copuladata,
      diag.panel = dp.copuladata,
      label.pos = 0.85, cex.labels = 1, gap = 0)
}

\arguments{
  \item{x}{\code{copuladata} object.}
  \item{labels}{variable names/labels.}
  \item{\dots}{other graphical parameters (see \code{\link[graphics]{par}}).}
  \item{lower.panel}{panel function to be used on the lower diagonal panels}
  \item{upper.panel}{panel function to be used on the upper diagonal panels}
  \item{diag.panel}{panel function to be used on the diagonal panels}
  \item{label.pos}{y position of labels in the diagonal panel; default: \code{label.pos = 0.85}.}
  \item{cex.labels}{magnification to be used for the labels of the diagonal panel; default: \code{cex.labels = 1}.}
  \item{gap}{distance between subplots, in margin lines; default: \code{gap = 0}.}
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
pairs(data[1:5], cex = 2, pch = 1,
      diag.panel = NULL, label.pos = 0.5,
      cex.labels = 2.5, gap = 1)
      
## pairs plot with own panel functions
up.copuladata <- function(x, y) {
  # lower panel: empirical contour plot
  op <- par(usr = c(-3, 3, -3, 3), new = TRUE)
  BiCopMetaContour(x, y, bw = 2, levels = c(0.01, 0.05, 0.1, 0.15, 0.2),
                    # exponential margins
                    margins = "exp", margins.par = 1,
                    axes = FALSE)
  on.exit(par(op))
}

lp.copuladata <- function(x, y) {
  # upper panel: scatter plot (copula data) and correlation
  op <- par(usr = c(0, 1, 0, 1), new = TRUE)
  points(x, y, pch = 1, col = "black")
  r <- cor(x, y, method = "spearman") # Spearman's rho
  txt <- format(x = r, digits = 3, nsmall = 3)[1]
  text(x = 0.5, y = 0.5, labels = txt, cex = 1 + abs(r) * 2, col = "blue")
  on.exit(par(op))
}

dp.copuladata <- function(x) {
  # diagonal panel: histograms (copula data)
  op <- par(usr = c(0, 1, 0, 1.5), new = TRUE)
  hist(x, freq = FALSE, add = TRUE, col = "brown", border = "white", main = "")
  abline(h = 1, col = "black", lty = 2)
  on.exit(par(op))
}

nlevels <- 20
pairs(data[1:5],
      lower.panel = lp.copuladata, upper.panel = up.copuladata,
      diag.panel = dp.copuladata, gap = 0.5)
}
