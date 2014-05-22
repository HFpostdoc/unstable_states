###Simulation to compare the dynamics of systems with and without stable states
###Idea = AMEllison
###Code = MKLau

##Using Bob May's logistic growth parameterization for ecosystems
##Using the quadratic recurrence equation (or quadratic map) which May develops by defining X = bN/a
##Here, we use the r notation that is more common now instead of a
##From dN/dt = rN(K-N) / K, setting x = N/K yields
##dx/dt = rx(1-x) and
##xt+1 = rxt (1-xt)

qre <- function(x0=0.1,r=2,tf=100,ts=1){
  t <- seq(0,tf,by=ts)
  out <- t*0
  out[1] <- x0
  for (i in 2:length(out)){out[i] <- r*out[i-1]*(1-out[i-1])}
  out <- cbind(Time=t,X=out)
  out
}

###plots the phase space xt+1 ~ xt
plot.phase <- function(x,col=1,pch=19,cex=0.75){
  plot(x[2:nrow(x),2]~x[1:(nrow(x)-1),2],
       col=col,pch=pch,cex=cex,
       xlab=expression('x'['t']),ylab=expression('x'['t+1']))
}

###plots lines for the phase space using spline
lines.phase <- function(x,col=1){
  lines(spline(x[2:nrow(x),2]~x[1:(nrow(x)-1),2]),col=col)
}

###plots points for the phase space
points.phase <- function(x,col=1){
  points(x[2:nrow(x),2]~x[1:(nrow(x)-1),2],col=col)
}

