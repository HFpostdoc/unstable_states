###Lab bench for the unstable states project
###Initiated 22May2014
###MKLau

source('../src/unstable_states.R')

###Basic simulation


#0. Choose ri and rf given dr

r. <- r[basins<20];basins. <- basins[basins<20];plot(r.,basins.)
abline(v=3.57,lty=2)
axis(side=1,at=3.57,labels='3.57',font=2)

#1. pick a random starting value for x between 0.0000000000001 and 1
#2. Run for tn time steps
#3. Change ri to rf
#4. Run tn time steps
#5. Plot time series
#6. Plot phase space
#7. Plot ews

