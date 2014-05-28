###Lab bench for the unstable states project
###Initiated 22May2014
###MKLau

source('../src/unstable_states.R')
library('YaleToolkit')

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


