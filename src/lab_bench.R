###Lab bench for the unstable states project
###Initiated 22May2014
###MKLau

source('../src/unstable_states.R')

###Basic simulation

#0. Choose ri and rf given dr
library(YaleToolkit)
r <- seq(2.4,4,length=10)
test <- list()
for (i in 1:length(r)){test[[i]] <- qre(r=r[i],tf=100,ts=0.1)[,2]}
test <- data.frame(do.call(cbind,test))[1:500,]
colnames(test) <- paste('r = ',round(r,2),sep='')
png('../results/spark_r.png',width=550,height=950,pointsize=15)
sparklines(test,overlap=FALSE,ylab=colnames(test),xlab='Time')
dev.off()

###honing in on the chaos threshold
r <- seq(3.56,3.58,by=0.001)
test <- list()
for (i in 1:length(r)){test[[i]] <- qre(r=r[i],tf=100,ts=0.1)[,2]}
test <- data.frame(do.call(cbind,test))[1:500,]
colnames(test) <- paste('r = ',round(r,2),sep='')
png('../results/spark_rchoatic.png',width=550,height=950,pointsize=10)
sparklines(test,overlap=FALSE,ylab=colnames(test),xlab='Time')
dev.off()

##
r <- seq(1,4,by=0.001)
test <- list()
for (i in 1:length(r)){test[[i]] <- qre(r=r[i],tf=500,ts=0.1)[,2]}
test <- data.frame(do.call(cbind,test))
test <- test[,-1]
r <- r[-1]
colnames(test) <- paste('r = ',round(r,2),sep='')

###Look at effects of threshold on basins
par(mfrow=c(1,5))
sparklines(test[,c(1,62,120,180,250,301)],overlap=FALSE,
           ylab=colnames(test)[c(1,62,120,180,250,301)],
           xlab='Time')
end.trans <- 3000
basins <- unlist(lapply(apply(test[end.trans:nrow(test),],2,
                              function(x) as.numeric(names(table(round(x,7))))),length)
                 )
par(mfrow=c(1,1))
plot(r,basins)
abline(v=r[((1:length(basins))[basins>100])[1]],lty=2)
text(x=r[((1:length(basins))[basins>100])],y=-50,labels=r[((1:length(basins))[basins>100])],col='red',cex=1.75)


#1. pick a random starting value for x between 0.0000000000001 and 1
#2. Run for tn time steps
#3. Change ri to rf
#4. Run tn time steps
#5. Plot time series
#6. Plot phase space
#7. Plot ews

