#written by Roger Mundry
#version from Oct. 29, 2011
#function returning the correction factor needed to be added to an AIC-value in
#case of a small sample
#arguments are:
#N: the sample size
#k: number of estimated parameters (needs to include the intercept and potentially
#dispersion parameter
#returns the correction factor to be added to the uncorrected AIC
aic.c.fac<-function(N, k){
  return(2*k*(k+1)/(N-k-1))
}
