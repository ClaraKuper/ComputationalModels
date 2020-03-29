exp_growth <- function(t,a=1,b=2,ter=1){
  # This function returns the exponential after time t
  # parameters: 
  # t = time how long the function was growing (unit), no default
  # a = constant, start value of the the function (units observation at t0), default: 1
  # b = basis, growth factor of the function, default: 2 (doubling)
  # ter = time constant, time required to reach the growth factor, default: 1
  
  # Example: if this is modeling a rabbit population and 
  # we observe this population over 30 days, 
  # the initial population was 20 animals, 
  # we know that every 2 months, the population is 1.5 times the time it was before, 
  # we would give it the values:
  # t = 30 (days)
  # a = 20 (individuals)
  # b = 1.5 (factor)
  # ter = 60 (days)
  
  grown_population <- a*b^(t/ter)
  
  return(grown_population)
}

# generate data
x <- seq(1,100)
y <- exp_growth(x,a=1,b=2,ter=15)

d <- cbind(x,y)

# parameter estimates using the optim function
startParams <- c(1,2,10)
names(startParams) <- c("a","b","ter")

# current data and preditions are being plotted
getregpred <- function(params,data){
  # here, we calculate the data from our independent variable under the assumption of the parameters b0(intercept) and b1(slope). Note that the names are hard coded here and that this function will only work for a linear regression. This is similar to the function that generated the data, but without the error term.
  getregpred <- params["a"]*params["b"]^(data[,1]/params["ter"])
  
  # return the computed regression data so we can calculate the deviation to the measured data
  return(getregpred)
}

# Get the computed predictions and compute how much they deviate

rmsd <- function(params,data){
  # generate predictions
  preds <- getregpred(params,data)
  #evaluate the deviation, this computes the sum of squared deviations between data points and devides by the total number of data points to get the mean summed square deviation
  rmsd <- sqrt(sum((preds-data[,2])^2)/length(preds))
}

#the parameter estimation with 
# for model1
xout1 <- optim(startParams,rmsd,data=d)

# plot data plus fit
plot(d)
y_pred <- exp_growth(d[,1],a=xout1$par[1],b=xout1$par[2],ter=xout1$par[3])
lines(d[,1],y_pred,col="red")