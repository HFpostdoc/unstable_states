display.cell <- function(x, bgcol="#FFFFFF", ...) {
  opar <- par(bg=bgcol, mar=rep(0,4))
  plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
  text(.5, .5, x, ...)
  lines(c(-0.1,1.1), c(0,0),col='lightgrey')
  par(opar)
}

format.digits <- function(x) as.character(paste("$", as.character(x), sep="   "))

display.sparkline <- function(x, y, bgcol="#FFFFFF", ...) {
  opar <- par(bg=bgcol, mar=rep(0,4))
  plot(c(0,max(x)), c(0,max(y)), type="n", axes=FALSE, xlab="", ylab="")
  lines(x, y, col='black')
  par(opar)
}
