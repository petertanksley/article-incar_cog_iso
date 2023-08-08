#### #################################################
# linear regression function
####################################################
  ## ----------------------------------------------------------------------------
  ## Title: my.lm.alt
  ## ----------------------------------------------------------------------------
  ##   Project     :       BiomJ article "Time-dependent mediators in survival analysis: 
  ##                             Modelling direct and indirect
  ##                             effects with the additive hazards model                                                            
  ##  Authors     :       Odd O Aalen, Mats J Stensrud, Vanessa Didelez, Rhian Daniel,
  ##                       Kjetil Roysland and Susanne Strohmaier   
  ##  Contact (code): Susanne Strohmaier (su.strohmaier@gmail.com)
  ## ----------------------------------------------------------------------------
  ## Description: Function running linear regressions at each 
  ##              mediator measurement time, but storing a linear 
  ##              regression estimate for each event time.
  ## ----------------------------------------------------------------------------
  ## Required Packages: ---
  ## ----------------------------------------------------------------------------
  ## Usage: my.lm.alt= function(regformula,event,startt, stopt, aggstop,
  ##                    mediator,dataset, w=1)
  ##                            
  ## requires data in counting process format!! 

  ##        regformula: independent variable input as for standard lm fuction in R 
  ##        event:   event indicator
  ##        startt:  start time of the observation intervals 
  ##        stopt:  end time of the observation intervals (actual event times)
  ##        aggstop: 'alternative end time' corresponding 
  ##              to the discret mediator measurement times. 
  ##        mediator: mediator of interest
  ##        dataset: dataset (in counting process format)
  ## ----------------------------------------------------------------------------
  ## Output: list of "coeff" for independet variables in "regformula" 
  ## ----------------------------------------------------------------------------
  
my.lm.alt= function(regformula,event,startt, stopt, aggstop, mediator,
                    dataset, w=1){
  if(length(w)==1){w = rep(1,length(dataset[,1]))}
  dataset = cbind(dataset,w)
  b.matrix = c() #Matrix to store b for each i
  # define and sort event times 
  obestime_x=dataset[[stopt]][dataset[[event]]==1]
  obstimes = sort(unique(obestime_x))
  #Loop doing a linear regression in each time point:
  for(i in obstimes){
    epsilon = 0.0001 #Ridge parameter epsilon
    # select the subset for subjects under risk at time i: 
    sta<-dataset[[startt]]
    aggst<-dataset[[aggstop]]
    obstimei = subset(dataset, sta<i & aggst>=i)
    
    weightsi = obstimei$w
    #Designmatrix Y
    Y = model.matrix(regformula,obstimei) 
    # define mediator
    E=obstimei[[mediator]]
    # linear regression coeff 'by hand'
    b = solve((t(Y) %*% (weightsi * Y) + 
                 epsilon * diag(length(Y[1,]))),t(Y) %*% (weightsi * E))
    b.matrix = rbind(b.matrix,t(b))
    
  }
  #Output:
  list(coeff = b.matrix)
}


#######################################################
# NEW addtive hazard function
#######################################################
my.additive.new = function(regformula,event,startt, stopt, dataset, w=1){
  if(length(w)==1){w = rep(1,length(dataset[,1]))}
  dataset = cbind(dataset,w)
  b.matrix = c() #Matrix to store b for each i
  #obstimes = sort(unique(dataset$stopt_itt))
  
  obestime_x=dataset[[stopt]][dataset[[event]]==1]
  obstimes = sort(unique(obestime_x))
  #Loop doing a linear regression in each time point:
  for(i in obstimes){
    epsilon = 0.0001 #Ridge parameter epsilon
    # select the subset for subjects under risk at time i: 
    sta<-dataset[[startt]]
    sto<-dataset[[stopt]]
    obstimei = subset(dataset, sta<i & sto>=i)
    weightsi = obstimei$w
    
    Y = model.matrix(regformula,obstimei)
    E = as.numeric(obstimei[[stopt]] == i & obstimei[[event]]==1) 
    #Dependent variable = indicator variable for event at time i
    b = solve((t(Y) %*% (weightsi * Y) + epsilon * diag(length(Y[1,]))),t(Y) %*% (weightsi * E))
    b.matrix = rbind(b.matrix,t(b))
  }
  #Find and return cummulative coefficients:
  cumcoeff = c()
  for(j in 1:(length(b.matrix[1,]))){
    cumcoeff = cbind(cumcoeff,cumsum(b.matrix[,j]))
  }
  colnames(cumcoeff) = colnames(b.matrix)
  #Output:
  list(coeff = b.matrix,cumcoeff = cumcoeff)
}

##############################################################
##############################################################