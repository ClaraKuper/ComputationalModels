#--------------------------------------------------------------------------
# Computational modeling
# MLE of bimodal distribution
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
library(maxLik)
library(ggplot2)
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
    
    ##### CLARA QUESTION #####
    # dnorm(lat[1],mu1,sigma1, log = FALSE) returns the likelihood for value lat[1] given parameters mu1 
    # and sigma1 under a normal distribution? 
    ##########################
  }
  
  # Run estimation
  mle3 <- maxLik(logLik = LLfnc, start = c(omega= runif(1,0,1), 
                                           mu1    = rnorm(1,100,10),
                                           sigma1 = runif(1,1,20),
                                           mu2    = rnorm(1,100,10),
                                           sigma2 = runif(1,1,20)),
                                           method="NM")  # method is Nelder-Mead
  
  
#--------------------------------------------------------------------------
  
#--------------------------------------------------------------------------
# plot results
hist(lat,freq=FALSE,xlim=c(0,400),ylim=c(0,0.1),breaks=50)
curve((1-coef(mle3)[1])*dnorm(x,coef(mle3)[2],coef(mle3)[3])+
      (coef(mle3)[1])*dnorm(x,coef(mle3)[4],coef(mle3)[5]),col="red",add=TRUE)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# MLE by hand
# 0. simulate data
# 1. define 2 parameter model
# 2. evaluate model for random parameter estimation by estimation log likelihood
#    save values for each iteration
#    find value closest do zero
# 3. plot likelihood function
# 4. plot results from model

# 0 -----------------------------------------------------------------------
# Generate the data (the data could be saccades and express saccades, for example)
mean1 <- 130  # mean latency for a saccade
mean2 <- 40  # mean latency for an express saccade

sd1    <- 10   # width of the gaussian
sd2    <- 10   
weight <- 0.5  # weight of the first mean

n     <- 3000 # number of data points

# generate data
simdat <-c(rnorm(weight*n,mean1,sd1),rnorm((1-weight)*n,mean2,sd2))

# 1 -----------------------------------------------------------------------
get_LL_bm <- function(data,mu1,mu2,sigma1=10,sigma2=10,omega = 0.5){
  # This function computes and return log likelihood for a bimodal distribution.
  # sd and weight values are assigned in the function:
  # sigma1: 10
  # sigma2: 10
  # omega: 0.7
  # distribution means mu1 and mu2 are expected as input
  data <- sort(data)
  
  return(sum(log((omega*(dnorm(data,mu1,sigma1,log=FALSE))) 
                  + ((1-omega)*(dnorm(data,mu2,sigma2,log=FALSE))))))
}

# 2 ----------------------------------------------------------------------

fit_means <- function(data,mu1_range,mu2_range,n_refine = 100, n_gen = 200){
  # this function fits two mean values of a bimodal distribution to the given data.
  # mu1_range and mu2_range are given as sequence start:end of a reasonable estimate for
  # the mean values. The start value is drawn randomly.
  # The model works in two stages:
  # Matrix stage: random value for mu1 and mu2 are drawn from a normal distribution and log likelihood is estimated
  # Loop stage: the best fit is chosen and defined as new mean for the normal distribution to generate new 
  # random estimates of mu1 and mu2. At the same time, the width of the curve is narrowed.
  # if the last stage of parameter estimation contains a better estimate than the current stage,
  # the values are not updated.
  
  params   <- c(mu1_mean = sample(mu1_range,1),mu2_mean = sample(mu2_range,1),mu_width = n_refine/2)
  parameter_space <- matrix(nrow = 0, ncol = 3)
  colnames(parameter_space) <- c("mu1", "mu2", "LL")
  print(params)

  
  for(r in 1:n_refine){
    randmat <- matrix(nrow = n_gen, ncol = 3)
    randmat[,1] <- rnorm(n_gen,params['mu1_mean'],params['mu_width'])
    randmat[,2] <- rnorm(n_gen,params['mu2_mean'],params['mu_width'])
    randmat[,3] <- mapply(get_LL_bm,randmat[,1],randmat[,2],MoreArgs = list(simdat))

    best_fit    <- which(randmat[,3]==max(randmat[,3]))
    print(randmat[best_fit,3])
    
    params['mu1_mean']   <- randmat[best_fit,1]
    params['mu2_mean']   <- randmat[best_fit,2]
    
    if(params['mu_width']>5){
      params['mu_width']   <- params['mu_width']-1
    }  
    parameter_space <- rbind(parameter_space,randmat)
    
    
  }
  return(parameter_space)
}


fit <- fit_means(simdat,1:50,100:150)

params <- c(fit[which(fit[,3]==max(fit[,3])),])

# 3 -------------------------------------------------------------------
vismat <- matrix(nrow=0,ncol=3)
for (m1 in 1:200){
  for (m2 in 1:200){
    writeMat     <- matrix(nrow=1,ncol=3)
    writeMat[1,1]<- m1
    writeMat[1,2]<- m2
    writeMat[1,3]<- get_LL_bm(simdat,m1,m2)
    vismat<-rbind(vismat,writeMat)
  }
}
colnames(vismat) <- c('mu1','mu2','LL')

fitLL<-as.data.frame(vismat)

fitLL <- fitLL[which(fitLL$LL>-20000),]

LL_surface <- ggplot(data = fitLL, aes(mu1, mu2)) + 
  geom_tile(aes(fill = LL), colour = "white") + 
  scale_fill_gradient(high = "white",low = "blue")

LL_surface

## -------------------------------
## this plot doesn't fill in colors... why?
##--------------------------------

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

hist(simdat,freq=FALSE,xlim=c(0,400),ylim=c(0,0.1),breaks=50)
curve(dnorm(x,params['mu1'],sd1),0,400,add=TRUE)
curve(dnorm(x,params['mu2'],sd2),0,400,add=TRUE)
