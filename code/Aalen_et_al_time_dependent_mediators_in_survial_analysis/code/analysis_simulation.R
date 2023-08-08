#######################################################
# DATA ANALYSIS
#######################################################
####################################################
## ----------------------------------------------------------------------------
## Filename    :	analysis_simulation.R
## ----------------------------------------------------------------------------
##   Project     :       BiomJ article "Time-dependent mediators in survival analysis: 
##                             Modelling direct and indirect
##                             effects with the additive hazards model                                                            
##  Authors     :       Odd O Aalen, Mats J Stensrud, Vanessa Didelez, Rhian Daniel,
##                       Kjetil Roysland and Susanne Strohmaier   
##  Contact (code): Susanne Strohmaier (su.strohmaier@gmail.com)
## ----------------------------------------------------------------------------
##  Date        :       19.12.2018
##  Purpose     :       to apply dynamic path analysis as decribed in
##                      the manuscript to the data simulated in 
##                      data.generation.R
##                      resulting plots are stored in /Results  
##  R Version   :       R - 3.5.1                                                           
##  
##  Input data files  :   mydata.fin
##  Used functions  :  my.lm.alt() ; my.additive.new()
##  Output :    results/Cumulative_hazards.pdf
##              results/Relative_survival.pdf
##   Required R packages :  ----
## ------------------------------------------------------------------
#####################################################################


# load functions for DPA 
source("functions/my.additive.new.R")
source("functions/my.lm.alt.R")
#load data
load("../data/mydata.fin.rda")

#################################################################3
## application of the linear model function
lmfit.N.2.alt<- my.lm.alt( ~ x,event="E",startt="start",
                           stopt="stopt", aggstop='alt_stop',mediator= "M", mydata.fin)

## application of the additive hazard model 
a.fit.N <- my.additive.new(~ x + M,event="E",startt="start",
                           stopt="stopt", mydata.fin)

##############################################################
# bootstrap loop for confidence intervals 
##############################################################
start.time.bootstrap <- Sys.time()

N=10 # number of boostrap samples
lm.effect=c() # object ot store the effects
#Bootstrap loop:
for(i in 1:N){
  print(i)
  #Create i-th bootstrap dataset   
  bootset<-mydata.fin[sample(nrow(mydata.fin),replace=T),]
  # unique event times in the bootstrap sample
  obstimes.boot = sort(unique(bootset$stopt[bootset$E==1]))
  #Regressions: functions calls as before
  # linear regression
  looplm.M<-my.lm.alt( ~ x,event="E",startt="start",
                       stopt="stopt", aggstop='alt_stop',mediator= "M",
                       bootset)
  
  # additive regression 
  loopadd.M<-my.additive.new(~ x + M,event="E",startt="start",
                             stopt="stopt", bootset)
  
  # effect of treatment on current mediator M
  INTonM<-looplm.M$coeff[,2]
  
  # direct effects of treamnent and M ( effects from additve regression)
  loopdirect.INT<-  loopadd.M$cumcoeff[,2]
  loopdirect.M<-  loopadd.M$cumcoeff[,3]
  
  # indirect treatment effects through M
  loopindirect.M<-cumsum( looplm.M$coeff[,2]*loopadd.M$coeff[,3])
  #total treatment effect
  looptotal<-cumsum(looplm.M$coeff[,2]*loopadd.M$coeff[,3]
                    + loopadd.M$coeff[,2])
  
  # estimates 
  obstimes.estimates<-cbind(obstimes.boot,  INTonM, 
                            loopdirect.INT,   loopdirect.M, 
                            loopindirect.M,
                            looptotal)
  # store effects:
  lm.effect<-rbind(lm.effect,obstimes.estimates)
  #Output:
  list(lm.effects = lm.effect)
}

end.time.bootstrap <- Sys.time()
time.taken.bootstrap <- end.time.bootstrap - start.time.bootstrap
time.taken.bootstrap

##############################################################

##############################################################
# event times
obstimes.imput = sort(unique(mydata.fin$stopt[mydata.fin$E==1] ))

# objects to store the estimates:

INTonM.l=c()
INTonM.u=c()

direct.INT.l=c()
direct.INT.u=c()

direct.M.l=c()
direct.M.u=c()

indirect.M.l=c()
indirect.M.u=c()

total.l=c()
total.u=c()

# pointwise CIs 

for (j in 1:length(obstimes.imput)){
  # select the corresponding columns form the 'lm.effect' object,
  #which is the output from the bootstrap loop
  subseti=subset(lm.effect,lm.effect[,1]==obstimes.imput[j])
  # pointwise CI for the effect of treatment on the mediators
  temp.INTonM.l=quantile(subseti[,2],0.025,names=FALSE)
  temp.INTonM.u=quantile(subseti[,2],0.975,names=FALSE)
  
  # pointwise CI for direct effect of treatment
  temp.direct.INT.l=quantile(subseti[,3],0.025,names=FALSE)
  temp.direct.INT.u=quantile(subseti[,3],0.975,names=FALSE)
  
  # pointwise CI for direct effect of the mediators 
  temp.direct.M.l=quantile(subseti[,4],0.025,names=FALSE)
  temp.direct.M.u=quantile(subseti[,4],0.975,names=FALSE)

  # pointwise CI for the indirect effectss
  temp.indirect.M.l=quantile(subseti[,5],0.025,names=FALSE)
  temp.indirect.M.u=quantile(subseti[,5],0.975,names=FALSE)
  
  # pointwise CI for the total effect
  temp.total.l=quantile(subseti[,6],0.025,names=FALSE)
  temp.total.u=quantile(subseti[,6],0.975,names=FALSE)
  
  # cobine pointswise values to one vector:
  INTonM.l=rbind(INTonM.l,temp.INTonM.l)
  INTonM.u=rbind(INTonM.u,temp.INTonM.u)
  
  direct.INT.l=rbind(direct.INT.l,temp.direct.INT.l)
  direct.INT.u=rbind(direct.INT.u,temp.direct.INT.u)
  
  direct.M.l=rbind(direct.M.l,temp.direct.M.l)
  direct.M.u=rbind(direct.M.u,temp.direct.M.u)
  
  indirect.M.l=rbind(indirect.M.l,temp.indirect.M.l)
  indirect.M.u=rbind(indirect.M.u,temp.indirect.M.u)
  
  total.l=rbind(total.l,temp.total.l)
  total.u=rbind(total.u,temp.total.u)
  
}

################################################################
#### cumulative hazards plots 
################################################################
par(mfrow=c(2,3))

# direct effect of treatment on outcome
plot(obstimes.imput,a.fit.N$cumcoeff[,2],type="l",
     main='Direct intervention effect on the outcome',
     xlab='Years since randomisation',
     ylab='Cumulative Effect', ylim=c(-0.05,0.30))
lines(obstimes.imput,rep(0,length(obstimes.imput)))
lines(obstimes.imput,direct.INT.l,type='l',col="grey")
lines(obstimes.imput,direct.INT.u,type='l',col="grey")

# indirect effects
ind.M<-cumsum(lmfit.N.2.alt$coeff[,2]*a.fit.N$coeff[,3])
plot(obstimes.imput,ind.M,type="l",
     main='Indirect effect through M',
     xlab='Years since randomisation',
     ylab='Cumulative Effect', ylim=c(-0.05,0.30))
lines(obstimes.imput,rep(0,length(obstimes.imput)))
lines(obstimes.imput, indirect.M.l,type='l',col="grey")
lines(obstimes.imput, indirect.M.u,type='l',col="grey")

# total effect 
total<-cumsum(lmfit.N.2.alt$coeff[,2]*a.fit.N$coeff[,3]+
                + a.fit.N$coeff[,2])

plot(obstimes.imput,total, type="l",
     main='Total effect of the intervention on the outcome',
     xlab='Years since randomisation',
     ylab='Cumulative Effect', ylim=c(-0.05,0.30))
lines(obstimes.imput,rep(0,length(obstimes.imput)))
lines(obstimes.imput, total.l,type='l',col="grey")
lines(obstimes.imput, total.u,type='l',col="grey")

# effect treatment on mediators; 
plot(obstimes.imput,lmfit.N.2.alt$coeff[,2],type="l",
     main='Intervention effect on M',
     xlab='Years since randomisation',
     ylab='Regression coefficient')
lines(obstimes.imput,INTonM.l,type='l',col="grey")
lines(obstimes.imput,INTonM.u,type='l',col="grey")

# direct effect of mediators on outcome
plot(obstimes.imput,a.fit.N$cumcoeff[,3],type="l",
     main='Direct effect of M on outcome',
     xlab='Years since randomisation',
     ylab='Cumulative Effect',ylim=c(-0.03,0.03))
lines(obstimes.imput,rep(0,length(obstimes.imput)))
lines(obstimes.imput,direct.M.l,type='l',col="grey")
lines(obstimes.imput,direct.M.u,type='l',col="grey")


####### proportion mediated
plot(obstimes.imput, ind.M/total, type="l",
     main='Proportion of effect mediated',
     xlab='Years since randomisation',
     ylab='Proportion mediated',
     ylim=c(-0.15,0.15))
lines(obstimes.imput,rep(0,length(obstimes.imput)))

##### save plots!!!! 

################################################################
#### relative survival plots 
################################################################

par(mfrow=c(2,3))
# direct effect of treatment on outcome
plot(obstimes.imput,exp(-a.fit.N$cumcoeff[,2]),type="l",
     main='Direct intervention effect on the outcome',
     xlab='Years since randomisation',
     ylab='Relative Survival', ylim=c(0.8,1.005))
lines(obstimes.imput,rep(1,length(obstimes.imput)))
lines(obstimes.imput,exp(-direct.INT.l),type='l',col="grey")
lines(obstimes.imput,exp(-direct.INT.u),type='l',col="grey")


# indirect effects
ind.M<-cumsum(lmfit.N.2.alt$coeff[,2]*a.fit.N$coeff[,3])
plot(obstimes.imput,exp(-ind.M),type="l",
     main='Indirect effect through M',
     xlab='Years since randomisation',
     ylab='Relative Survival', ylim=c(0.8,1.005))
lines(obstimes.imput,rep(1,length(obstimes.imput)))
lines(obstimes.imput, exp(-indirect.M.l),type='l',col="grey")
lines(obstimes.imput, exp(-indirect.M.u),type='l',col="grey")

# total effect 
total<-cumsum(lmfit.N.2.alt$coeff[,2]*a.fit.N$coeff[,3]+
                + a.fit.N$coeff[,2])

plot(obstimes.imput,exp(-total), type="l",
     main='Total effect of the intervention on the outcome',
     xlab='Years since randomisation',
     ylab='Relative Survival', ylim=c(0.8,1.005))
lines(obstimes.imput,rep(1,length(obstimes.imput)))
lines(obstimes.imput,exp(-total.l),type='l',col="grey")
lines(obstimes.imput, exp(-total.u),type='l',col="grey")

# effect treatment on mediators; 
plot(obstimes.imput,lmfit.N.2.alt$coeff[,2],type="l",
     main='Intervention effect on diastolic BP',
     xlab='Years since randomisation',
     ylab='Regression coefficient')
lines(obstimes.imput,INTonM.l,type='l',col="grey")
lines(obstimes.imput,INTonM.u,type='l',col="grey")

# direct effect of mediators on outcome
plot(obstimes.imput,exp(-a.fit.N$cumcoeff[,3]),type="l",
     main='Direct effect of diastolic BP on outcome',
     xlab='Years since randomisation',
     ylab='Relative Survival', ylim=c(0.95,1.03))
lines(obstimes.imput,rep(1,length(obstimes.imput)))
lines(obstimes.imput,exp(-direct.M.l),type='l',col="grey")
lines(obstimes.imput,exp(-direct.M.u),type='l',col="grey")

##### save plots!!!! 

# save.image(/results/workspace_simulations) 

