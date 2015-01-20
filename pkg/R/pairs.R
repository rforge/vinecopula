pairs.copuladata <- function(x, labels = names(x), ...,
                             label.pos = 0.85, cex.labels = 1, gap = 0, axes = FALSE,
                             pch = ".", col = "grey", cex.points = 1,
                             method.cor = "kendall", col.cor = "red", digits.cor = 2, cex.cor = 1,
                             bw = 2, size = 100, levels = seq(0.01, 0.2, length.out = 30),
                             margins = "norm", margins.par = 0, xylim = NA,
                             col.contour = terrain.colors(length(levels)),
                             col.hist = "grey"){
  ## pairs plot for 'copuladata'
  
  ## labeling of axes
  if(axes){
    xaxt <- "s"
    yaxt <- "s"
  } else {
    xaxt <- "n"
    yaxt <- "n"
  }
  
  ## lower panel: empirical contour plot
  lower.panel.copuladata <- function(x, y, lower.bw = bw, lower.size = size,
                                     lower.levels = levels, lower.margins = margins, 
                                     lower.margins.par = margins.par, lower.xylim = xylim, 
                                     col = col.contour, ...){
    op <- par(usr = c(-3, 3, -3, 3), new = TRUE)
    BiCopMetaContour(x, y, bw = lower.bw, size = lower.size,
                     levels = lower.levels, axes = FALSE,
                     margins = lower.margins, margins.par = lower.margins.par,
                     xylim = lower.xylim, col = col, drawlabels = FALSE, ...)
    on.exit(par(op))
  }
  
  ## upper panel: scatter plot (copula data) and correlation
  upper.panel.copuladata  <- function(x, y, method=method.cor, upper.pch = pch, upper.col = col,
                                      upper.col.text = col.cor, upper.cex = cex.points, 
                                      upper.digits = digits.cor, upper.cex.cor = cex.cor, ...){
    op <- par(usr = c(0, 1, 0, 1), new = TRUE)
    plot(x, y, pch = upper.pch, cex = upper.cex, col = upper.col, axes=FALSE, ...)
    r <- cor(x, y, method = method)
    txt <- format(r, digits = upper.digits, nsmall = upper.digits)[1]
    text(0.5, 0.5, txt, cex = upper.cex.cor + abs(r)*3, col = upper.col.text)
    on.exit(par(op))
  }
  
  ## diagonal panel: histograms (copula data)
  diag.panel.copuladata  <- function(x, diag.col=col.hist, ...){
    op <- par(usr = c(0, 1, 0, 1.6), new = TRUE)
    hist(x, freq = FALSE, add = TRUE, col = diag.col, border = "black", main = "")
    abline(h = 1, col = "black", lty=3)
    on.exit(par(op))
  }
  
  ## pairs plot (with panel functions as defined above)
  pairs.default(x, labels = labels, ...,
                lower.panel = lower.panel.copuladata,
                upper.panel = upper.panel.copuladata,
                diag.panel = diag.panel.copuladata,
                label.pos = label.pos, cex.labels = cex.labels,
                gap = gap, xaxt=xaxt, yaxt=yaxt)
}