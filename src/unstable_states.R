###Simulation to compare the dynamics of systems with and without stable states
###Idea = AMEllison
###Code = MKLau
###Initiated 21May2014

###Find the start and end points for days in the model output
find.days <- function(y){
  f <- 0
  s <- 1
  for (i in 1:(length(y)-1)){
    if (y[i]==0&y[i+1]!=0){s <- c(s,i)}else if (y[i]!=0&y[i+1]==0){f <- c(f,i)}
  }
  f <- f[-1]
  out <- apply(rbind(s,f),2,function(x) x[1]:x[2])
  return(out)
}

###This formulation was taken from Gotelli 4th ed.
disrupt.mc <- function(mu=0,sd=0.1,ri=1,rf=2,tf=1000,N=10,K=100,burn=100,dump=FALSE){
  if (burn > tf){burn <- round(tf*0.25,2)}
  r <- c(rep(ri,burn),seq(ri,rf,length=tf))
  r <- r + rnorm(length(r),mu,sd)
  for (t in 2:burn){
    N[t] <- N[t-1] + r[t]*N[t-1]*(1-N[t-1]/K)
  }
  for (t in burn:(burn+tf)){
    N[t] <- N[t-1] + r[t]*N[t-1]*(1-N[t-1]/K)
  }
  if (dump){N <- data.frame(t=1:(burn+tf),r=r,N=N)}
  return(N)
}


grow <- function(r=2.570,N=10,K=100,tf=100){
  for (t in 2:tf){
    N[t] <- N[t-1] + r*N[t-1]*(1-N[t-1]/K)
  }
  return(N)
}

disrupt <- function(ri=1,rf=2,tf=1000,N=10,K=100,burn=100,dump=FALSE){
  if (burn > tf){burn <- round(tf*0.25,2)}
  r <- c(rep(ri,burn),seq(ri,rf,length=tf))
  for (t in 2:burn){
    N[t] <- N[t-1] + r[t]*N[t-1]*(1-N[t-1]/K)
  }
  for (t in burn:(burn+tf)){
    N[t] <- N[t-1] + r[t]*N[t-1]*(1-N[t-1]/K)
  }
  if (dump){N <- data.frame(t=1:(burn+tf),r=r,N=N)}
  return(N)
}


##Using Bob May's logistic growth parameterization for ecosystems
##Using the quadratic recurrence equation (or quadratic map) which May develops by defining X = bN/a
##Here, we use the r notation that is more common now instead of a
##From dN/dt = rN(K-N) / K, setting x = N/K yields
##dx/dt = rx(1-x) and
##xt+1 = rxt (1-xt)

###lvq = lotka volterra equations
lvq <- function(N1=10,N2=10,r1=2,r2=2,K1=100,K2=90,a=0.5,b=0.5){
  dn1 <- r1*N1-(N1^2*r1/K1)-(N1*N2*a*r1/K1)
  dn2 <- r2*N2-(N2^2*r2/K2)-(N1*N2*b*r2/K2)
  return(c(N1+dn1,N2+dn2))
}

###run model
runmod <- function(Ni,t,FUN){
  N <- Ni
  for (i in 2:t){
    N[i] <- N[i-1] + FUN(N[i-1])
  }
  N
}

###allee affect
allee <- function(N,r=0.1,K=100,a=1){
  r*N*(1-(N/K)-(a*(N^2)/K))
}

###
qre <- function(x0=0.1,r=2,tf=100,ts=1){
  t <- seq(0,tf,by=ts)
  out <- t*0
  out[1] <- x0
  for (i in 2:length(out)){out[i] <- r*out[i-1]*(1-out[i-1])}
  out <- cbind(Time=t,X=out)
  out
}

###Basic simulation
###simSS = simulate system state changes given a logistic difference model
simSS <- function(x0,r1,r2,n,ramp=FALSE,step=0.1){
                                        #set the start points for state transitions
  nf1 <- round(n/2,0)
  nf2 <- nf1 + 1
  out <- list()
  if (ramp){}else{
    out[[1]] <- qre(x0=x0,r=r1,tf=(nf1 + 1),ts=1)
    out[[2]] <- qre(x0=out[[1]][nrow(out[[1]]),2],r=r2,tf=length(1:nf2),ts=1)[-1,]
    out[[1]][,1] <- out[[1]][,1] 
    out[[2]][,1] <- out[[2]][,1] + (out[[1]][nrow(out[[1]]),1])
    out <- do.call(rbind,out)
    colnames(out) <- c('t','x')
  }
  return(out)
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

