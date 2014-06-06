###Lab bench for the unstable states project
###Initiated 22May2014
###MKLau

source('../src/unstable_states.R')
library('fractal')
library('YaleToolkit')

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

###Chaos and not
n1 <- grow(r=2.50000001,N=10,K=300,tf=1000)
n2 <- grow(r=2.50000001,N=11,K=300,tf=1000)
c1 <- grow(r=2.59999999,N=10,K=300,tf=1000)
c2 <- grow(r=2.59999999,N=11,K=300,tf=1000)

par(mfrow=c(2,2))
plot(n1,type='l');lines(n2,col='red')
plot(c1,type='l');lines(c2,col='red')
plot(n1-n2,type='l')
plot(c1-c2,type='l')
par(mfrow=c(1,1))
hist(c1-c2,freq=FALSE)
lines(density(rnorm(length(c1),mean(c1-c2),sd(c1-c2))))
plot(density(c1-c2))
lines(density(rnorm(length(c1),mean(c1-c2),sd(c1-c2))),lty=2)


###Disrupt simulations
ri <- c(1.997,2.455,2.550,2.57)-(0.02/2)
rf <- c(1.997,2.455,2.550,2.57)+(0.02/2)

###With error in r
par(mfcol=c(2,4))
plot(disrupt.mc(sd=0.01,ri=ri[1],rf=rf[1]),xlab='time',ylab='N')
title(main=rf[1])
plot(disrupt.mc(sd=0.01,ri=ri[1],rf=rf[1],dump=TRUE)$r,xlab='time',ylab='r')
plot(disrupt.mc(sd=0.01,ri=ri[2],rf=rf[2]),xlab='time',ylab='N')
title(main=rf[2])
plot(disrupt.mc(sd=0.01,ri=ri[2],rf=rf[2],dump=TRUE)$r,xlab='time',ylab='r')
plot(disrupt.mc(sd=0.01,ri=ri[3],rf=rf[3]),xlab='time',ylab='N')
title(main=rf[3])
plot(disrupt.mc(sd=0.01,ri=ri[3],rf=rf[3],dump=TRUE)$r,xlab='time',ylab='r')
plot(disrupt.mc(sd=0.01,ri=ri[4],rf=rf[4]),xlab='time',ylab='N')
title(main=rf[4])
plot(disrupt.mc(sd=0.01,ri=ri[4],rf=rf[4],dump=TRUE)$r,xlab='time',ylab='r')

###phases
par(mfrow=c(2,2))
for (i in 1:4){
  n <- disrupt.mc(sd=0.01,ri=ri[i],rf=rf[i])
  plot(n[2:length(n)]~n[1:(length(n)-1)],type='l',col='grey')
}

###focus on transition to chaos
opar <- par()

par(mfrow=c(2,2),mai=opar$mai*0.01)
plot(n1 <- disrupt.mc(sd=0,ri=ri[4],rf=rf[4]),xlab='',ylab='',type='l',xaxt='null',yaxt='null',bty='n')
plot(n2 <- disrupt.mc(sd=0.01,ri=ri[4],rf=rf[4]),xlab='',ylab='',type='l',xaxt='null',yaxt='null',bty='n')
plot(n1[11:length(n1)]~n1[10:(length(n1)-1)],type='l',col='black',lwd=0.25,xlab='',ylab='',xaxt='null',yaxt='null',bty='n')
plot(n2[11:length(n2)]~n2[10:(length(n2)-1)],type='l',col='black',lwd=0.25,xlab='',ylab='',xaxt='null',yaxt='null',bty='n')

par(mfrow=c(10,10),mai=opar$mai*0.01)
for (i in 1:100){
  n1 <- disrupt.mc(sd=0.01,ri=ri[4],rf=rf[4])
  plot(n1[11:length(n1)]~n1[10:(length(n1)-1)],type='l',col='black',lwd=0.25,xlab='',ylab='',xaxt='null',yaxt='null',bty='n')
}

###Looking at early warning signals
n0 <- disrupt.mc(sd=0,ri=ri[4],rf=rf[4],dump=TRUE)
n1 <- disrupt.mc(sd=0.005,ri=ri[4],rf=rf[4],dump=TRUE)
n2 <- disrupt.mc(sd=0.01,ri=ri[4],rf=rf[4],dump=TRUE)

##Crossing the threshold of r
par(mfrow=c(1,3))
plot(n0$r~n0$t,type='l',ylab='r')
abline(h=2.571451,lty=2)
plot(n1$r~n1$t,xlab='t')
abline(h=2.571451,lty=2)
plot(n2$r~n2$t)
abline(h=2.571451,lty=2)

##EWS
library(earlywarnings)
ews0 <- generic_ews(n0$N)
ews1 <- generic_ews(n1$N)


###No error in r
par(mfcol=c(2,4))
plot(disrupt(ri=ri[1],rf=rf[1]),xlab='time',ylab='N')
title(main=rf[1])
plot(disrupt(ri=ri[1],rf=rf[1],dump=TRUE)$r,xlab='time',ylab='r')
plot(disrupt(ri=ri[2],rf=rf[2]),xlab='time',ylab='N')
title(main=rf[2])
plot(disrupt(ri=ri[2],rf=rf[2],dump=TRUE)$r,xlab='time',ylab='r')
plot(disrupt(ri=ri[3],rf=rf[3]),xlab='time',ylab='N')
title(main=rf[3])
plot(disrupt(ri=ri[3],rf=rf[3],dump=TRUE)$r,xlab='time',ylab='r')
plot(disrupt(ri=ri[4],rf=rf[4]),xlab='time',ylab='N')
title(main=rf[4])
plot(disrupt(ri=ri[4],rf=rf[4],dump=TRUE)$r,xlab='time',ylab='r')

###phases
par(mfrow=c(2,2))
for (i in 1:4){
  n <- disrupt(ri=ri[i],rf=rf[i])
  plot(n[2:length(n)]~n[1:(length(n)-1)])
}

###
dndt <- function(N,r,k){
  -r*N + (r/k)*N^2
}

###Luewis Kareiva model
lk <- function(N,r,k,A){
  r*N*(1-(N/k))*((N/k)-(A/k))
}

Vdndt <- function(N,r,k){
  (r/(3*k))*N^3 - (r/2)*N^2
}

Vlk <- function(N,r,k,A){
  (r/k)*((1/(4*k))*N^4 - (1/3)*(1+(A/k))*N^3 + (A/2)*N^2)
}

N <- 1:80
par(mfcol=c(2,2))
plot(-sapply(N,dndt,r=0.1,k=62)~N,ylim=c(-2,2))
abline(h=0,v=62,lty=2)
plot(sapply(N,Vdndt,r=0.1,k=62)~N,ylim=c(-70,10))
abline(h=0,v=62,lty=2)
plot(sapply(N,lk,r=0.1,k=62,A=20)~N,ylim=c(-0.6,0.6))
abline(h=0,v=62,lty=2)
plot(sapply(N,Vlk,r=0.1,k=62,A=20)~N,ylim=c(-12,4))
abline(h=0,v=62,lty=2)

fold <- function(x,a){
  x^3 + a*x
}

fold.sim <- sapply(seq(-10,10,by=0.1),function(a) fold(0:100,a))

par(mfrow=c(1,1))
for (i in 1:ncol(fold.sim)){
  plot(fold.sim[,i])
}

while(2<3){
t <- 100
out <- matrix(rep(0,t*2),nrow=t)
out[1,] <- c(sample(1:100,1),sample(1:100,1))
pert <- c(1)
r1 <- r2 <- 2
for (i in 2:t){
  if (sample(c(TRUE,FALSE),1,prob=c(0.01,0.99))){
    out[i,] <- lvq(out[i-1,1]-(out[i-1,1]*runif(1,0,1)),
                   out[i-1,2]-(out[i-1,2]*runif(1,0,1)),r1=r1,r2=r2)
    pert[i] <- 2
  }else{
    out[i,] <- lvq(out[i-1,1],out[i-1,2],r1=r1,r2=r2)
    pert[i] <- 1
  }
}
par(mfrow=c(1,1))
plot(out[,1],type='l',col=0,xlab='',ylab='',frame.plot=FALSE,xaxt='n',yaxt='n',ylim=c(0,max(out)))
for (i in 1:nrow(out)){
  points(i,out[i,1],col=pert[i])
  points(i,out[i,2],col=pert[i],pch=19)
}
plot(out,type='l',col=0,xlab='',ylab='',frame.plot=FALSE,xaxt='n',yaxt='n')
for (i in 1:nrow(out)){
  points(out[i,1],out[i,2],col=pert[i]-1,pch=19)
}
for (i in 1:nrow(out)){
  points(out[i,1],out[i,2],col=rainbow(nrow(out))[i],pch=19,cex=7)
}
for (i in 1:nrow(out)){
  points(out[i,1],out[i,2],col=pert[i]-1,pch=19,cex=2)
}
}


###
while(2<3){
  xy <- sample(1:3,2)
  plot(lorenz[,xy],xaxt='n',yaxt='n',col=c('black','grey','white')[sample(1:3,1)],xlab='',ylab='',frame.plot=FALSE)
  for (i in 2:nrow(lorenz)){
    points(lorenz[i,xy[1]],lorenz[i,xy[2]],col=rainbow(nrow(lorenz))[i],pch=19)
  }
}

###
runmod(1,100,FUN=function(x) allee(x,a=2))
runmod(1,100,FUN=function(x) qre(x))

###Conducting simulations
n.sim <- 150
n.x0 <- 25
x0 <- runif(n.x0)
sim.12 <- sim.24 <- sim.48 <- sim.8c <- sim.cc <- list()
for (i in 1:length(x0)){
  sim.12[[i]] <- simSS(x0[i],r1=(2.995 - 0.024),r2=(2.995 + 0.024),n=n.sim)
  sim.24[[i]] <- simSS(x0[i],r1=(3.447 - 0.024),r2=(3.447 + 0.024),n=n.sim)
  sim.48[[i]] <- simSS(x0[i],r1=(3.543 - 0.024),r2=(3.543 + 0.024),n=n.sim)
  sim.8c[[i]] <- simSS(x0[i],r1=3.545,r2=(3.545+0.048),n=n.sim)
  sim.cc[[i]] <- simSS(x0[i],r1=(3.545+0.048),r2=(3.545+0.048+0.048),n=n.sim)
}
sim.out <- list(sim.12,sim.24,sim.48,sim.8c,sim.cc)
dput(sim.out,file='../results/simulations_22May2014')
df.12 <- data.frame(do.call(cbind,lapply(sim.12,function(x) x[,2])))
df.24 <- data.frame(do.call(cbind,lapply(sim.24,function(x) x[,2])))
df.48 <- data.frame(do.call(cbind,lapply(sim.48,function(x) x[,2])))
df.8c <- data.frame(do.call(cbind,lapply(sim.8c,function(x) x[,2])))
df.cc <- data.frame(do.call(cbind,lapply(sim.cc,function(x) x[,2])))
###graphs
sparklines(df.12,overlap=TRUE,yscale=c(0,1),xaxis=FALSE,yaxis=FALSE)
sparklines(df.24,overlap=TRUE,yscale=c(0,1),xaxis=FALSE,yaxis=FALSE)
sparklines(df.48,overlap=TRUE,yscale=c(0,1),xaxis=FALSE,yaxis=FALSE)
sparklines(df.8c,overlap=TRUE,yscale=c(0,1),xaxis=FALSE,yaxis=FALSE)
sparklines(df.cc,overlap=TRUE,yscale=c(0,1),xaxis=FALSE,yaxis=FALSE)
###state-space
par(mfrow=c(1,5))
plot.phase(sim.12[[1]],col=c('blue',rep('blue',(n.sim/2)),rep('red',(n.sim/2)),'red','red'))
plot.phase(sim.24[[1]],col=c('blue',rep('blue',(n.sim/2)),rep('red',(n.sim/2)),'red','red'))
plot.phase(sim.48[[1]],col=c('blue',rep('blue',(n.sim/2)),rep('red',(n.sim/2)),'red','red'))
plot.phase(sim.8c[[1]],col=c('blue',rep('blue',(n.sim/2)),rep('red',(n.sim/2)),'red','red'))
plot.phase(sim.cc[[1]],col=c('blue',rep('blue',(n.sim/2)),rep('red',(n.sim/2)),'red','red'))

#0. Choose ri and rf given dr
r. <- r[basins<20];basins. <- basins[basins<20];plot(r.,basins.)
abline(v=3.57,lty=2)
axis(side=1,at=3.57,labels='3.57',font=2)
table(basins.)[table(basins.)>20]
range(r.[basins.==1]);range(r.[basins.==2]);range(r.[basins.==4]);range(r.[basins.==8])

##       rl    ru
## [1] 1.004 2.995
## [1] 3.003 3.447
## [1] 3.451 3.543
## [1] 3.545 3.69


