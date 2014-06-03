###Lab bench for the unstable states project
###Initiated 22May2014
###MKLau

source('../src/unstable_states.R')
library('fractal')
library('YaleToolkit')

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


