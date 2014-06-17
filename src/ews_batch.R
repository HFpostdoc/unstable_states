###Run repeated simulations for ews time series
source('../src/unstable_states.R')
library(earlywarnings)
cb8.16 <- 2.57 #choatic boundary, 8-16 cycles
ri <- (cb8.16)-(0.02/2)
rf <- (cb8.16)+(0.02/2)
n <- 150
stats8.16 <- list()
stats16.8 <- list()
for (i in 1:n){
  print(i)
  drmc8.16 <- disrupt.mc(N=i,sd=0.01,ri=ri,rf=rf,dump=TRUE)
  drmc16.8 <- disrupt.mc(N=i,sd=0.01,ri=rf,rf=ri,dump=TRUE)
  drmc8.16$N[abs(drmc8.16$N)==Inf] <- 0;drmc8.16$N[(drmc8.16$N)<0] <- 0
  drmc16.8$N[abs(drmc16.8$N)==Inf] <- 0;drmc16.8$N[(drmc16.8$N)<0] <- 0
  if (sd(drmc8.16$N[-1])!=0){
    ews8.16 <- generic_ews(drmc8.16$N)  
    stats8.16[[i]] <- cor(ews8.16,method='ken')[1,-1]
    dev.off()
  }else{stats8.16[[i]] <- NA}
  if (sd(drmc16.8$N[-1])!=0){
    ews16.8 <- generic_ews(drmc16.8$N)
    stats16.8[[i]] <- cor(ews16.8,method='ken')[1,-1]
    dev.off()
  }else{stats16.8[[i]] <- NA}
}
###
dput(stats8.16,file='../results/stats816.rdata')
dput(stats16.8,file='../results/stats168.rdata')
