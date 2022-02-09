##-------------------------------##
## Martingale type residual (MT) ##
## left censoring
##-------------------------------##

library(PerformanceAnalytics)

resMT <- function(theta,y,X,cc,family="ST")
{
  n <- length(y)
  p <- ncol(X)
  resm <- resmt <- S <- c()
  
  
  if(family=="SN")
  {
    betas <- theta[1:p]
    sigma2 <- theta[(p+1)]
    lambda <- theta[(p+2)]
    
    Delta <- sqrt(sigma2)*(lambda/sqrt(lambda^2 + 1))
    mu <- X%*%betas # - sqrt(2/pi)*Delta
    S <- 1 - cdfSNI(y,mu,sigma2,lambda,nu=NULL,type=family)
  }
  if(family=="ST")
  {
    betas <- theta[1:p]
    sigma2 <- theta[(p+1)]
    lambda <- theta[(p+2)]
    nu <- theta[(p+3)]
    
    Delta <- sqrt(sigma2)*(lambda/sqrt(lambda^2 + 1))
    k1 <- sqrt(nu/2)*(gamma((nu-1)/2)/gamma(nu/2))
    mu <- X%*%betas# - sqrt(2/pi)*k1*Delta
    S <- 1 - cdfSNI(y,mu,sigma2,lambda,nu,type=family)
  }
    
  for (i in 1:n)
  {
    resm[i] <-  1-cc[i]+log(S[i])
    resmt[i] <-  sqrt(-2*((1-cc[i])*log(1-cc[i]-resm[i])+resm[i]))*sign(resm[i])
  }
  
  return(list(resm=resm , resmt=resmt)) 
}

## ---------------------------------------------------- ##
## Envelopes of the MT residuals for the SMSN-CR models ##
## ---------------------------------------------------- ##

EnvelopeRMT <- function(theta,y,X,cc,family="ST")
{
  res <- resMT(theta,y,X,cc,family=family)
  #windows()
  if(family=="SN")
  {
    chart.QQPlot(res$resmt,distribution='norm',envelope=0.95,ylab="Martingale-type residuals",xlab="Standard normal quantile",main=" ",line="quartile", col=c(20,1),pch=19)
  }
  if(family=="ST")
  {
    chart.QQPlot(res$resmt,distribution='norm',envelope=0.95,ylab="Martingale-type residuals",xlab="Standard normal quantile",main=" ",line="quartile", col=c(20,1),pch=19)
  }
  }


#EnvelopeRMT(est_T$theta,y,x,cc,cens="left",dist="ST")
