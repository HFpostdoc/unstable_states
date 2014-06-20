###Get the early warning signals for the sensitivity analysis
###MKLau 20 Jun 2014

sddSens <- function(x,cx,eval=100){
  ##1. standardize (xi - xbar) / sdx
  x <- (x - mean(x))/sd(x)
  cx <- (cx - mean(cx))/sd(cx)
  ##2. fit function using loess
  fit.x <- loess.smooth(x=(1:length(x)),y=x,evaluation=eval)
  fit.cx <- loess.smooth(x=(1:length(x)),y=cx,evaluation=eval)
  ##3. detrend and decycle by subtracting the loess of controls
  fit.x$y-fit.cx$y
}

