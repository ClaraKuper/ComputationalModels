#--------------------------------------------------------------------------
# Computational modeling
# Parameter estimation using optim
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
library(minpack.lm)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Generate data
time  <- seq(0,25)
alpha <- 0.05
y <- exp(alpha*time)

y <- y + rnorm(length(y),0,2)
y <- ifelse(y <=0, 0, y)
df <- data.frame(time=time,y=y)

plot(time,y)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Fit estimation with 
m1 <- nlsLM(y ~ exp(alpha*time),
             data=df,
             start = list(alpha=0))
df$fitted <- fitted(m1)
points(df$time,df$fitted,col="red")
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Fit estimation with optim based on code in book

# parameter estimates using the optim function
startParams <- c(0.01)
names(startParams) <- c("alpha")

# current data and preditions are being plotted
getpred <- function(params,data){
  getpred <-  exp(params["alpha"]*data$time)
  
  # draw graph when key is pressed
#  par(ask=FALSE)
#  plot(data[,2],type="n", las = 1, ylim = c(25,35),xlim=c(-2,2),xlab="X",ylab="Y")
#  par(ask=FALSE)
#  points(data[,2],data[,1],pch=21,bg="gray")
#  lines(data[,2],getpred,lty="solid")
  
  # return the computed regression data so we can calculate the deviation to the measured data
  return(getpred)
}

# Get the computed predictions and compute how much they deviate

rmsd <- function(params,data1){
  # generate predictions
  preds <- getpred(params,data1)
  #evaluate the deviation, this computes the sum of squared deviations between data points and devides by the total number of data points to get the mean summed square deviation
  rmsd <- sqrt(sum((preds-data1[,2])^2)/length(preds))
}

#the parameter estimation with 
# for model1
xout1 <- optim(startParams,rmsd,data1=df,method="BFGS")

df$fitted_optim <- exp(xout1$par*df$time)
points(df$time,df$fitted_optim,col="green")

# same model with simplex instead of BFGS
xout2 <- optim(startParams,rmsd, data1 = df)
df$fitted_smplx <- exp(xout2$par*df$time)
points(df$time,df$fitted_smplx, col = "orange")

#--------------------------------------------------------------------------