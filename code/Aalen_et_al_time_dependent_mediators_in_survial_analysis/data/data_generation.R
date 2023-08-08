#######################################################################################
#
#                                                                                     
#   Filename    :	data_generation.R    												  
#                                                                                     
#   Project     :       BiomJ article "Time-dependent mediators in survival analysis: 
#                             Modelling direct and indirect
#                             effects with the additive hazards model                                                            
#   Authors     :       Odd O Aalen, Mats J Stensrud, Vanessa Didelez, Rhian Daniel,
#                       Kjetil Roysland and Susanne Strohmaier                                                                
#   Date        :       19.12.2018
#   Purpose     :       to simulate data minicking important features of 
#                       of the data used in the real life example in the manuscript    
#                       Detailed explanations of the data generating system 
#                       can be found in Strohmaier et.al. (Stats in Med), Section 4.1
#   R Version   :       R - 3.5.1                                                           
#   
#
#   Input data files  :    ----                                                    
#   Output data files :    mydata.fin
#
#   Required R packages :  ----
#
#
########################################################################################

#####################################
##### data simulation #######
#####################################
set.seed(5)
# number of patients
n <- 10000; 
delta <- 1/52; 
time<-seq(delta, 5, by=delta)
# The weeks measurements are made after treatment
measure.week <- c(4,8, 12, 26 , 52, 78, 104, 156)
measure.time <- time[measure.week]
f.w <- measure.week[1]/52
l.o <- (length(time)-measure.week[1])
fw <- measure.week[1]
lw <- length(time) - 1
nM <- length(measure.week)
# Randomise treatment x and sample M0
x <- sample(x=c(0, 1), n, replace=TRUE)
M0 <- rnorm(n, 7.8, 1.5)
# Effect of treatment 0 on M
b0 <- spline(0:5, rep(0, 6), xout=seq(delta, 5, delta))
b0 <- b0$y
# effect of treatment 1 on M
b1 <- spline(0:5, c(-0.1, -3, -2.2, -3.3, -2.9,-2.9)/2, xout=seq(delta, 5, delta))
b1 <- b1$y
# Generate M matrix and noise
M <- x %*% t(b1) + abs(x-1) %*% t(b0) + M0
Sigma <- matrix(0, length(time), length(time))
diag(Sigma) <- 0.05
mu <- rep(0, length(time))
noise <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
M.mat <- M + noise

patient <- seq(1, n)
# combine into a data set
sim.dat <- as.data.frame(cbind(patient, x, M0, M.mat[, measure.week]))
names(sim.dat)[4:(4 + length(measure.week) - 1)] <-
  paste0("M", measure.week)

# Effect of x on hazard
beta1 <- spline(c(0, 0.2, 0.8, 1.1, 3.5, 5),
                c(0.6, 0.2, 0.13, 0.1, 0.1, 0.1)/2,
                xout=seq(delta, 5, delta))
beta1 <- beta1$y
# Matrix with effect of x on hazard
X.mat <- sim.dat$x %*% t(beta1)
# Effect of M on hazard
beta2 <- spline(c(0,1,3,5), c(-0.04, -0.02,-0.02, -0.02)/7, xout=seq(delta,5,delta))
beta2 <- beta2$y
beta2.mat <- matrix(rep(beta2, n), byrow=T, nrow=n)
# Matrix with effects of M on hazard - elementwise product
Meff.mat <- M.mat * beta2.mat


# remember to include a matrix of beta0 values if needed
hmat <- X.mat + Meff.mat
beta0 <- min(hmat)
hmat <- hmat + abs(beta0)
pmat <- hmat*delta
# adding id column
p.mat.org <- cbind(1:n, pmat)

sim.dat$Time <- rep(time[length(time)], n)
p.mat <- p.mat.org
for(t in (measure.week[1]):length(time)) {
  t.i <- time[t]
  p.vec <- p.mat[, (t+1)]
  E <- rep(0, length(p.vec))
  id <- rep(0, length(p.vec))
  for(j in 1:length(p.vec)) {
    id[j] <- p.mat[j, 1]
    E[j] <- sample(c(0,1), 1, prob=c(1-p.vec[j], p.vec[j]))
  }
  if(any(E==1)) {
    id.event <- id[E==1]
    sim.dat[sim.dat$patient %in% id.event, "Time"] <- time[t]
    p.mat <- p.mat[E==0, ]
  }
}

sim.dat$E <- 1*(sim.dat$Time<5)  

# Censoring
ff <- runif(n, 0, 25)
sim.dat$Time <- sim.dat$Time - f.w
sim.dat$C.Time <- ff
sim.dat$C <- 1*(sim.dat$C.Time < sim.dat$Time)
sim.dat[sim.dat$C==1, "E"] <- 0
sim.dat[sim.dat$C==1, "Time"] <- sim.dat[sim.dat$C==1, "C.Time"] 

# Reshaping data to suited format
mydata <- reshape(sim.dat, paste0("M", measure.week),
                  times = (measure.time - f.w),
                  direction="long", v.names="M")
mydata$M12 <- rep(mydata[1:n, "M"], nM)

# Ordering
mydata <- mydata[ order(mydata[,"patient"], mydata[,"time"]), ]
row.names(mydata) <- 1:nrow(mydata)
mydata <- mydata[, c("patient","x","M0","M12", "M","Time","C", "E","time")]
names(mydata)[9] <- "start"
mydata$stopt <- rep(c(measure.time[2:nM], max(time)) - f.w, n)

m.split <- split(mydata, mydata$patient)

hax <- function(df) {
  df$E <- 1*(df$Time <= df$stopt & sum(df$E) == nM )
  df$C <- 1*(df$Time <= df$stopt & sum(df$C) == nM )
  df$stopt <- pmin(c(df$stopt), df$Time)
  keep <- 1:min(c(which(df$E==1), nM, which(df$C==1)))
  df <- df[keep, ]
}

m.split.hax <- lapply(m.split, hax)
mydata.hax <- do.call("rbind", m.split.hax)
mydata.hax[mydata.hax$E==1, "stopt"] <- mydata.hax[mydata.hax$E==1, "stopt"] + 
  runif(length(mydata.hax[mydata.hax$E==1, "stopt"]), 0, 0.005)

sim.dat.long <- mydata.hax[, c("patient","x","M0","M12", "M","C","E","start","stopt")]
mydata.fin <- sim.dat.long


#####################################################################
# data management step mainly for `bookkeeping` reasons to 
# for discrete mediator measurments
mydata.fin$alt_stop[mydata.fin$stopt<=(1/52)*4]<-(1/52)*4
mydata.fin$alt_stop[mydata.fin$stopt>(1/52)*4 & mydata.fin$stopt<=(1/52)*8 ]<-(1/52)*8
mydata.fin$alt_stop[mydata.fin$stopt>(1/52)*8 & mydata.fin$stopt<=(1/52)*12 ]<-(1/52)*12
mydata.fin$alt_stop[mydata.fin$stopt>(1/52)*12 & mydata.fin$stopt<=(1/52)*26 ]<-(1/52)*26
mydata.fin$alt_stop[mydata.fin$stopt>(1/52)*26 & mydata.fin$stopt<=(1/52)*52 ]<-(1/52)*52
mydata.fin$alt_stop[mydata.fin$stopt>(1/52)*52 & mydata.fin$stopt<=(1/52)*78]<-(1/52)*78
mydata.fin$alt_stop[mydata.fin$stopt>(1/52)*78 & mydata.fin$stopt<=(1/52)*104 ]<-(1/52)*104
mydata.fin$alt_stop[mydata.fin$stopt>(1/52)*104 & mydata.fin$stopt<=(1/52)*156 ]<-(1/52)*156
mydata.fin$alt_stop[mydata.fin$stopt>(1/52)*156]<-(1/52)*300

# save data
save(mydata.fin, file="../data/mydata.fin.rda")


