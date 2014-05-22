###Simulation to compare the dynamics of systems with and without stable states
###Idea = AMEllison
###Code = MKLau

###Use Bob May's logistic growth parameterization for ecosystems

#analytical solution to the least difference equation
grow.lde <- function(x0=0.5,r=1,t0=0,tn=10,n=0.001){
  t <- seq(t0,tn,by=n)
  out <- 1 / (1 + ((1/x0) - 1) * exp(1)^(-r*t))
  out <- cbind(t,out)
  colnames(out) <- c('time','X')
  out
}

r <- 1;i <- 1;dx <- seq(0.05,1,by=0.001)
plot(grow.lde(x0=0,r=r),ylim=c(0,1),col=rainbow(length(dx))[i],type='l')
for (x0 in dx){lines(grow.lde(x0,r=r),col=rainbow(length(dx))[i]);i <- i + 1}

grow.log <- function(N=1,a=2,b=0,tf=5,step=0.1){
  time <- seq(step,tf,by=step)
  out <- time * 0
  out[1] <- N
  for (t in 2:length(time)){out[t] <- out[t-1]*(a - b*out[t-1])}
  out <- cbind(time,out)
  colnames(out) <- c('Time','N')
  out
}

N.a2.707 <- grow.log(a=2.707,b=0.1,step=0.01)
N.a3.414 <- grow.log(a=3.414,b=0.1,step=0.01)
N.a3.570 <- grow.log(a=3.570,b=0.1,step=0.01)
N.a3.571 <- grow.log(a=3.571,b=0.1,step=0.01)
N.a3.700 <- grow.log(a=3.700,b=0.1,step=0.01)

png('../results/alpha_chaos.png',width=1900,height=900,pointsize=25)
par(mfcol=c(2,5))
N <- N.a2.707;plot(N,cex=0.5,pch=19);plot(N[2:nrow(N),2]~N[1:(nrow(N)-1),2],xlab='Nt',ylab='Nt+1')
N <- N.a3.414;plot(N,cex=0.5,pch=19);plot(N[2:nrow(N),2]~N[1:(nrow(N)-1),2],xlab='Nt',ylab='Nt+1')
N <- N.a3.570;plot(N,cex=0.5,pch=19);plot(N[2:nrow(N),2]~N[1:(nrow(N)-1),2],xlab='Nt',ylab='Nt+1')
N <- N.a3.571;plot(N,cex=0.5,pch=19);plot(N[2:nrow(N),2]~N[1:(nrow(N)-1),2],xlab='Nt',ylab='Nt+1')
N <- N.a3.700;plot(N,cex=0.5,pch=19);plot(N[2:nrow(N),2]~N[1:(nrow(N)-1),2],xlab='Nt',ylab='Nt+1')
dev.off()

###Use the quadratic recurrence equation (or quadratic map) which May develops by defining X = bN/a
##Use the r notation that is more common now instead of a

##From dN/dt = rN(K-N) / K, setting x = N/K yields
##dx/dt = rx(1-x) and
##xt+1 = rxt (1-xt)

plot.phase <- function(x,col=1){
  plot(x[2:nrow(x),2]~x[1:(nrow(x)-1),2],col=col)
}

lines.phase <- function(x,col=1){
  lines(spline(x[2:nrow(x),2]~x[1:(nrow(x)-1),2]),col=col)
}

points.phase <- function(x,col=1){
  points(x[2:nrow(x),2]~x[1:(nrow(x)-1),2],col=col)
}


qre <- function(x0=0.1,r=2,tf=500,ts=1){
  t <- seq(0,tf,by=ts)
  out <- t*0
  out[1] <- x0
  for (i in 2:length(out)){out[i] <- r*out[i-1]*(1-out[i-1])}
  out <- cbind(Time=t,X=out)
  out
}

##Re-creating May 1976 Fig 1
png('../results/qre_phase.png',width=500,height=500,pointsize=15)
plot(c(0,1),c(0,1),type='l',xlab=expression('X'['t']),ylab=expression('X'['t+1']))
points.phase(qre(x0=0.0001,r=2.707,tf=100,ts=0.1))
lines.phase(qre(x0=0.0001,r=2.707,tf=100,ts=0.1))
points.phase(qre(x0=0.0001,r=3.414,tf=100,ts=0.1))
lines.phase(qre(x0=0.0001,r=3.414,tf=100,ts=0.1))
dev.off()
