pairs.copuladata <- function(x,
                             labels = names(x),
                             ...,
                             lower.panel = lp.copuladata,
                             upper.panel = up.copuladata,
                             diag.panel = dp.copuladata,
                             label.pos = 0.85, 
                             cex.labels = 1,
                             gap = 0) {
  ## pairs plot for 'copuladata'
  
  # provide input data and set default labels, panel functions, etc.
  default <- list(x = as.matrix(x),
                  labels = labels,
                  lower.panel = lower.panel,
                  upper.panel = upper.panel,
                  diag.panel = diag.panel,
                  label.pos = label.pos, 
                  cex.labels = cex.labels,
                  gap = gap
                  )
  
  # pairs plot (with panel functions as defined below or as provided by user)
  pars <- modifyList(list(xaxt = "n", yaxt = "n"), list(...))
  op <- do.call(par, pars)
  do.call(pairs, modifyList(default, list(...)))
  on.exit(par(op))
}


## lower panel: empirical contour plot
lp.copuladata <- function(x, y, ...) {
  # set default parameters
  pars <- list(u1 = x,
               u2 = y,
               bw = 2, 
               size = 100, 
               levels = seq(0.01, 0.2, length.out = 30),
               margins = "norm", 
               margins.par = 0, 
               xylim = NA, 
               col = terrain.colors(30),
               axes = FALSE,
               drawlabels = FALSE)
  # get non-default parameters
  pars <- modifyList(pars, list(...))
  op <- par(usr = c(-3, 3, -3, 3), new = TRUE)
  # call BiCopMetaContour
  do.call(BiCopMetaContour, pars)
  on.exit(par(op))
}


## upper panel: scatter plot (copula data) and correlation
up.copuladata <- function(x, y, ...) {
  # set default parameters
  pars <- list(x = x,
               y = y,
               pch = ".",
               cex = 1,
               col = "grey"
  )
  # get non-default parameters
  pars <- modifyList(pars, list(...))
  op <- par(usr = c(0, 1, 0, 1), new = TRUE)
  # call points (to produce scatter plot)
  do.call(points, pars)
  r <- cor(x = x, y = y, method = "kendall")
  txt <- format(x = r, digits = 2, nsmall = 2)[1]
  # call text
  do.call(text, modifyList(list(x = 0.5,
                                y = 0.5,
                                labels = txt,
                                cex = 1 + abs(r) * 3,
                                col = "red"),
                           list(...)
                           )
  )
  on.exit(par(op))
}


## diagonal panel: histograms (copula data)
dp.copuladata <- function(x, ...) {
  # set default parameters
  pars <- list(x = x,
               freq = FALSE,
               add = TRUE,
               col = "grey",
               border = "black", 
               main = "")
  # get non-default parameters
  pars <- modifyList(pars, list(...))
  op <- par(usr = c(0, 1, 0, 1.6), new = TRUE)
  # call hist
  do.call(hist, pars)
  if (pars$freq == FALSE)
    abline(h = 1, col = "black", lty = 3)
  on.exit(par(op))
}
