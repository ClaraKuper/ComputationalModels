#--------------------------------------------------------------------------
# Computational modeling
# MLE of bimodal distribution
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
library(maxLik)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Generate the data (the data could be saccades and express saccades, for example)
lat_sac <- 150  # mean latency for a saccade
lat_exp <- 80  # mean latency for an express saccade

width   <- 5   # width of the gaussian

# generate data with 1500 data points for normal saccades and 500 data points for express saccades
lat<-c(rnorm(1500,lat_sac,width),rnorm(500,lat_exp,width))
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# 1: Define likelihood
# 5 parameter model
# omega: weight of first gaussian
# mu1: mean of first gaussian
# sigma1: sd of first gaussian
# mu2: mean of 2nd gaussian
# sigma2: sd of 2nd gaussian
  LLfnc <- function(param) {
    omega  <- param[1]
    if(omega<0 || omega>1){return(NA) }
    mu1    <- param[2]
    sigma1 <- param[3]
    if(sigma1<=0){return(NA)}
    mu2    <- param[4]
    sigma2 <- param[5]
    if(sigma2<=0){return(NA)}    
    
    #LogLikelihood LL
    log((omega*mapply(dnorm,lat,mu1,sigma1,log=FALSE)) + ((1-omega)*mapply(dnorm,lat,mu2,sigma2,log=FALSE)))
   
  }
  
  # Run estimation
  mle3 <- maxLik(logLik = LLfnc, start = c(omega= runif(1,0,1), 
                                           mu1    = rnorm(1,100,10),
                                           sigma1 = runif(1,1,20),
                                           mu2    = rnorm(1,100,10),
                                           sigma2 = runif(1,1,20)),
                                           method="NM")
#--------------------------------------------------------------------------
  
#--------------------------------------------------------------------------
# plot results
hist(lat,freq=FALSE,xlim=c(0,400),ylim=c(0,0.1),breaks=50)
curve((1-coef(mle3)[1])*dnorm(x,coef(mle3)[2],coef(mle3)[3])+
      (coef(mle3)[1])*dnorm(x,coef(mle3)[4],coef(mle3)[5]),col="red",add=TRUE)
#--------------------------------------------------------------------------