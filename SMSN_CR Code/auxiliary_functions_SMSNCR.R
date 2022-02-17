### =================================================================== ###
### Likelihood Based Inference for Censored Linear Regression Models    ###
###         with Scale Mixtures of Skew-Normal Distributions            ###
###  Authors: Thalita do Bem Mattos, Aldo M. Garay and Victor H. Lachos ###
### =================================================================== ###

### --------------------------------###
### AUXILIARY FUNCTIONS FOR SMSN-CR ###
### --------------------------------###

if(!require(mnormt)) install.packages("mnormt");  library(mnormt)
if(!require(moments)) install.packages("moments");library(moments)
if(!require(tcltk2)) install.packages("tcltk2");library(tcltk2)
if(!require(sn)) install.packages("sn");library(sn)
if(!require(Matrix)) install.packages("Matrix");library(Matrix)







## -------------------------------- ##
## CUMULATIVE DISTRIBUTION FUNCTION ##
## -------------------------------- ##

cdfSNI<-function(x,mu,sigma2,lambda,nu,type="SN")
{
  n <- length(x)
  resp<-matrix(0,n,1)
  if(type=="Normal")
  {
    resp <- pnorm((x-mu)/sqrt(sigma2))
    return(resp)
  }
  
  if(type=="T")
  {
    resp <- pt((x-mu)/sqrt(sigma2),df=nu)
    return(resp)
  }
  
  if(type=="SN")
  {
    delta <- lambda/sqrt(1+lambda^2)
    SIGMA <- matrix(c(sigma2,-delta*sqrt(sigma2),-delta*sqrt(sigma2),1),byrow=TRUE,ncol=2,nrow=2)
    if(length(mu)==1)
    {
      MU <- cbind(rep(mu,n),0)
    }
    if(length(mu)==n)
    {
      MU <- cbind(mu,0)
    }
    Y <- cbind(x,0)
    for(i in 1:n)
    {
      resp[i] <- 2*pmnorm(x=Y[i,],mean=MU[i,],varcov=SIGMA)
    }
    return(resp)
  }
  
  if(type=="ST")
  {
    delta <- lambda/sqrt(1+lambda^2)
    SIGMA <- matrix(c(sigma2,-delta*sqrt(sigma2),-delta*sqrt(sigma2),1),byrow=TRUE,ncol=2,nrow=2)
    if(length(mu)==1)
    {
      MU <- cbind(rep(mu,n),0)
    }
    if(length(mu)==n)
    {
      MU <- cbind(mu,0)
    }
    Y <- cbind(x,0)
    nu <- round(nu)
    for(i in 1:n)
    {
      resp[i] <- 2*pmt(x=Y[i,], mean = MU[i,], S=SIGMA, df=nu)
    }  
    return(resp)
    
  }
  
  if(type=="SSL")
  {
    cdf<- function(y)
    {
      f <- function(u) 2*nu*u^(nu - 1)*dnorm(y,mu,sqrt(u^(-1)*sigma2))*pnorm(u^(1/2)*lambda*(y-mu)/sqrt(sigma2))
      cdf <- integrate(Vectorize(f),0,1)$value
    }
    densidade <- as.numeric(cdf(x))
    resp<-as.numeric(integrate(Vectorize(cdf),-Inf,x)$value)
    return(list(pdf=densidade,cdf=resp))
  }
  
  if(type=="SCN")
  {
    if(length(mu)==1)
    {
      MU <- cbind(rep(mu,n))
    }
    if(length(mu)==n)
    {
      MU <- cbind(mu)
    }
    
    dSNC <- function(y, mu, sigma2, lambda, nu)
    {
      dens <- 2*(nu[1]*dnorm(y, mu, sqrt(sigma2/nu[2]))*pnorm(sqrt(nu[2])*lambda*sigma2^(-1/2)*(y-mu)) + (1 - nu[1])*dnorm(y, mu, sqrt(sigma2))*pnorm(lambda*sigma2^(-1/2)*(y-mu)))
      return(dens)
    }
    
    for(i in 1:n)
    {
      resp[i] <- as.numeric(integrate(dSNC,-Inf,x[i],MU[i],sigma2,lambda,nu)$value)
    }
    return(resp)
  }
  
}

## ---------------------------- ##
## PROBABILITY DENSITY FUNCTION ##
## ---------------------------- ##

pdfSNI <- function(y,mu,sigma2,lambda,nu,type="SN")
{
  # para SCN : nu = (nu,gamma)
  resp<-matrix(0,length(y),1) 
  if(type=="Normal")
  {
    resp <- dnorm((y-mu)/sqrt(sigma2))/sqrt(sigma2)
  }
  
  if(type=="T")
  {
    resp <- dt((y-mu)/sqrt(sigma2),df=nu)/sqrt(sigma2)
  }
  
  if(type=="SN")
  {
    resp<-2*dnorm((y-mu)/sqrt(sigma2))*pnorm(lambda*(y-mu)/sqrt(sigma2))/sqrt(sigma2)
  }
  if(type=="ST")
  {
    z=(y-mu)/sqrt(sigma2)
    resp=2*dt(z,df=nu)*pt(sqrt(nu+1)*lambda*z/sqrt(nu+z^2),df=nu+1)/sqrt(sigma2)
  }
  if(type=="SSL")
  {
    f <- function(u){ 2*nu*u^(nu - 1)*dnorm(y,mu,sqrt(u^(-1)*sigma2))*pnorm(u^(1/2)*lambda*(y-mu)/sqrt(sigma2))}
    resp <- integrate(Vectorize(f),0,1)$value
  }
  if(type=="SCN")
  {
    resp <- 2*(nu[1]*dnorm(y, mu, sqrt(sigma2/nu[2]))*pnorm(sqrt(nu[2])*lambda*sigma2^(-1/2)*(y-mu)) + (1 - nu[1])*dnorm(y, mu, sqrt(sigma2))*pnorm(lambda*sigma2^(-1/2)*(y-mu)))
  }
  return(resp)
}

## ----------------------- ##
## GENERATE VARIABLES SMSN ##
## ----------------------- ##

rSMSN <- function(n,mu,sigma2,lambda,nu,dist){
  if(dist=="SN"|dist=="SSL"|dist=="ST"|dist=="SCN"){
    if(length(lambda) == 0) stop("lambda must be provided.")
  }
  if(dist=="ST"|dist=="SSL"){
    if(length(nu) == 0) stop("nu must be provided.") 
  }
  if(dist=="SCN"){
    if(length(nu) != 2) stop("nu must be a vector of size 2.") 
  }
  y <- rep(0,n)
  if(dist=="SN")
  {
    u <- rep(1,n)
  }
  if(dist=="ST")
  {
    u <- rgamma(n=n,shape=nu/2,rate=nu/2)
  }  
  if(dist=="SSL")
  {
    u <- rbeta(n=n,shape1=nu,shape2=1)	
  }
  if(dist=="SCN")
  {
    p <- runif(n)
    u <- rep(1,n)
    u[p<nu[1]] <- nu[2]
  }
  deltinha <- lambda/sqrt(1+lambda^2)
  Delta <-  sqrt(sigma2)*deltinha
  tau <- sigma2*(1-deltinha^2)
  
  T0 <- rnorm(n)
  T1 <- rnorm(n)
  T2 <- abs(T0)*u^(-1/2)
  y <-  mu + Delta*T2 + u^(-1/2)*sqrt(tau)*T1
  
  return(y)
}

## ---------------------------------- ##
## GENERATE VARIABLE TRUNCATED NORMAL ##
## ---------------------------------- ##

rTN <- function(n,mu,sigma2,a,b)
{
  #a: lower
  #b: upper
  u <- runif(n)
  sigma <- sqrt(sigma2)
  aux <- u*(pnorm(b,mean=mu,sd=sigma)-pnorm(a,mean=mu,sd=sigma))+pnorm(a,mean=mu,sd=sigma) 
  amostra.x <- qnorm(aux,mean=mu,sd=sigma)
  return(amostra.x)
}

## ----------------------------- ##
## GENERATE VARIABLE TRUNCATED t ##
## ----------------------------- ##

rTt <- function(n,mu,sigma2,nu,a,b)
{
  u <- runif(n)
  aux <- u*(cdfSNI(b,mu,sigma2,lambda=NULL,nu,type="T") - cdfSNI(a,mu,sigma2,lambda=NULL,nu,type="T")) + cdfSNI(a,mu,sigma2,lambda=NULL,nu,type="T")     
  #aux <- u*(pt((b-mu)/sqrt(sigma2),df=nu) - pt((a-mu)/sqrt(sigma2),df=nu)) +  pt((a-mu)/sqrt(sigma2),df=nu) 
  amostra.x <- mu + sqrt(sigma2)*qt(aux, df=nu)
  return(amostra.x)
}

## --------------------------------------- ##
## GENERATE VARIABLE TRUNCATED SKEW-NORMAL ##
## --------------------------------------- ##

rTSN <- function(n,mu,sigma2,lambda,a,b,cens) 
{
  m <- n*20
  if(cens=="left"){cens = 1}
  if(cens=="right"){cens = 2}
  if(cens=="interval"){cens = 3}
  y <- rTN(m,mu,sigma2,a,b)
  if(cens=="1") 
  {
    densTSN <-  pdfSNI(y,mu,sigma2,lambda,nu=NULL,type="SN")/c(cdfSNI(b,mu,sigma2,lambda,nu=NULL,type="SN"))
  } 
  if(cens=="2")
  {
    densTSN <-  pdfSNI(y,mu,sigma2,lambda,nu=NULL,type="SN")/c(1-cdfSNI(a,mu,sigma2,lambda,nu=NULL,type="SN"))
  }
  if(cens=="3")
  {
    densTSN <-  pdfSNI(y,mu,sigma2,lambda,nu=NULL,type="SN")/c(cdfSNI(b,mu,sigma2,lambda,nu=NULL,type="SN")-cdfSNI(a,mu,sigma2,lambda,nu=NULL,type="SN"))
  }
  densTN <- dnorm(y,mu,sqrt(sigma2))/(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
  W <- densTSN/densTN
  pesos <- W/sum(W)
  x <- sample(y,m/20,replace = FALSE, prob = pesos)
  return(x)
}

## ---------------------------------- ##
## GENERATE VARIABLE TRUNCATED SKEW-t ##
## ---------------------------------- ##

rTSt <- function(n,mu,sigma2,lambda,nu,a,b,cens)
{
  m <- n*20 
  if(cens=="left"){cens = 1}
  if(cens=="right"){cens = 2}
  if(cens=="interval"){cens = 3}
  y <- rTt(m,mu,sigma2,nu,a,b)
  if(cens=="1") 
  {
    densTSt <- pdfSNI(y,mu,sigma2,lambda,nu,type="ST")/c(cdfSNI(b,mu,sigma2,lambda,nu,type="ST")) 
    densTt <-  pdfSNI(y,mu,sigma2,lambda=NULL,nu,type="T")/c(cdfSNI(b,mu,sigma2,lambda=NULL,nu,type="T")) 
  } 
  if(cens=="2")
  {
    densTSt <- pdfSNI(y,mu,sigma2,lambda,nu,type="ST")/c(1 - cdfSNI(a,mu,sigma2,lambda,nu,type="ST"))
    densTt <-  pdfSNI(y,mu,sigma2,lambda=NULL,nu,type="T")/c(1-cdfSNI(a,mu,sigma2,lambda=NULL,nu,type="T"))
  }
  if(cens=="3")
  {
    densTSt <- pdfSNI(y,mu,sigma2,lambda,nu,type="ST")/c(cdfSNI(b,mu,sigma2,lambda,nu,type="ST") - cdfSNI(a,mu,sigma2,lambda,nu,type="ST"))
    densTt <- pdfSNI(y,mu,sigma2,lambda=NULL,nu,type="T")/c(cdfSNI(b,mu,sigma2,lambda=NULL,nu,type="T") - cdfSNI(a,mu,sigma2,lambda=NULL,nu,type="T"))  
  }
  W <- densTSt/densTt
  pesos <- W/sum(W)
  x <- sample(y,m/20,replace = FALSE, prob = pesos)
  return(x)
}

### ---------------------------------------------------- ###
### GENERATE VARIABLE TRUNCATED SKEW-NORMAL CONTAMINATED ###
### ---------------------------------------------------- ###

rTSCN <- function(n,mu,sigma2,lambda,nu,a,b,cens)
{
  z <- c()
  p <- runif(n)
  u <- rep(1,n)
  u[p<nu[1]] <- nu[2]
  a1 <- (a-mu)*sqrt(u)
  b1 <- (b-mu)*sqrt(u)
  for(i in 1:n)
  {
    z[i] <- rTSN(n=1,mu=0,sigma2,lambda,a1[i],b1[i],cens) 
  }
  y <- mu + u^(-1/2)*z
  return(y)
}

### -------------------------------------- ### 
### GENERATE VARIABLE TRUNCATED SKEW-SLASH ###
### -------------------------------------- ###

rTSSl <- function(n,mu,sigma2,lambda,nu,a,b,cens)
{
  z <- c()
  u <- rbeta(n=n,shape1=nu,shape2=1)
  a1 <- (a-mu)*sqrt(u)
  b1 <- (b-mu)*sqrt(u)
  for(i in 1:n)
  {
    z[i] <- rTSN(n=1,mu=0,sigma2,lambda,a1[i],b1[i],cens) 
  }
  y <- mu + u^(-1/2)*z
  return(y)
}

## ------------------------ ##
## EXPECTATION NOT CENSORED ##
## ------------------------ ##

NCensEsp <- function(y,X,betas,sigma2,lambda,nu,dist)
{            
  deltinha <- lambda/sqrt(1+(lambda^2))
  Delta <-  sqrt(sigma2)*deltinha
  tau <- sigma2 - Delta^2 
  n <- length(y)
  
  if(dist=="SN")
  {
    eta <- -sqrt(2/pi)
    mu <- X%*%betas + Delta*eta
    
    muT <- (Delta/(tau+Delta^2))*(y - mu) 
    M2T <- tau/(tau + Delta^2)
    MT <- sqrt(M2T)
    Wphi <- dnorm(muT/MT)/pnorm(muT/MT)
    
    EU <- cbind(rep(1,n))
    EUY <- y
    EUY2 <- y^2
    EUT <- (muT + eta) + MT*Wphi
    EUT2 <- (muT + eta)^2 + M2T + MT*(muT + 2*eta)*Wphi
    EUYT <- y*((muT + eta) + MT*Wphi)
    ElogU <- log(EU) 
  }
  
  if(dist=="ST")
  {
    k1 <- sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
    eta <- -sqrt(2/pi)*k1
    mu <- X%*%betas + Delta*eta
    
    d <- ((y - mu)/sqrt(sigma2))^2
    M2T <- 1/(1 + (Delta^2)*(tau^(-1)))
    MT <- sqrt(M2T)
    muT <- M2T*Delta*(tau^(-1))*(y - mu)
    muT_eta <- muT + eta
    A <- muT/ MT
    
    E <-(2*(nu)^(nu/2)*gamma((2+nu)/2)*((d + nu + A^2))^(-(2+nu)/2)) / (gamma(nu/2)*pi*sqrt(sigma2)*pdfSNI(y,mu,sigma2,lambda,nu,type="ST"))
    U <- ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(d + nu)^(-(nu+3)/2))/(gamma(nu/2)*sqrt(pi)*sqrt(sigma2)*pdfSNI(y,mu,sigma2,lambda,nu,type="ST")))*pt(sqrt((3+nu)/(d+nu))*A,3+nu)
    
    EU <- U
    EUY <- y*U
    EUY2 <- (y^2)*U
    EUT <- U*muT_eta + MT*E
    EUT2 <- U*(muT_eta^2) + M2T + MT*(muT_eta + eta)*E
    EUYT <- y*(U*muT_eta + MT*E)
    ElogU <- log(U)
  }
  
  if(dist=="SSL")
  {
    k1 <- 2*nu/(2*nu-1)
    eta <- -sqrt(2/pi)*k1
    mu <- X%*%betas + Delta*eta
    
    d <- ((y - mu)/sqrt(sigma2))^2
    M2T <- 1/(1 + (Delta^2)*(tau^(-1)))
    MT <- sqrt(M2T) 
    muT <- M2T*Delta*(tau^(-1))*(y - mu)
    muT_eta <- muT + eta
    A <- muT/ MT
    
    u <- vector(mode="numeric",length=n)
    E <- vector(mode="numeric",length=n)
    for(i in 1:n)
    {
      E[i] <- (((2^(nu + 1))*nu*gamma(nu + 1))/(pdfSNI(y[i],mu[i],sigma2,lambda,nu,type="SSL")*pi*sqrt(sigma2)))* ((d[i]+A[i]^2)^(-nu-1))*pgamma(1,nu+1,(d[i]+A[i]^2)/2)
      faux <- function(u) u^(nu+0.5)*exp(-u*d[i]/2)*pnorm(u^(1/2)*A[i])
      aux22 <- integrate(faux,0,1)$value
      u[i] <- ((sqrt(2)*nu) / (pdfSNI(y[i],mu[i],sigma2,lambda,nu,type="SSL")*sqrt(pi)*sqrt(sigma2)))*aux22
    } 
    
    EU <- u
    EUY <- y*u
    EUY2 <- (y^2)*u
    EUT <- u*muT_eta + MT*E
    EUT2 <- u*(muT_eta^2) + M2T + MT*(muT_eta + eta)*E
    EUYT <- y*(u*muT_eta + MT*E)
    ElogU <- log(u)
  }
  
  if(dist=="SCN")
  {
    k1 <- nu[1]/nu[2]^(1/2)+1-nu[1]
    eta <- -sqrt(2/pi)*k1
    mu <- X%*%betas + Delta*eta
    
    d <- ((y - mu)/sqrt(sigma2))^2
    M2T <- 1/(1 + (Delta^2)*(tau^(-1)))
    MT <- sqrt(M2T) 
    muT <- M2T*Delta*(tau^(-1))*(y - mu)
    muT_eta <- muT + eta
    A <- muT/ MT
    
    u <- (2/pdfSNI(y,mu,sigma2,lambda,nu,type="SCN"))*(nu[1]*nu[2]*dnorm(y,mu,sqrt(sigma2/nu[2]))*pnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dnorm(y,mu,sqrt(sigma2))*pnorm(A,0,1))
    E <- (2/pdfSNI(y,mu,sigma2,lambda,nu,type="SCN"))*(nu[1]*sqrt(nu[2])*dnorm(y,mu,sqrt(sigma2/nu[2]))*dnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dnorm(y,mu,sqrt(sigma2))*dnorm(A,0,1))
    
    EU <- u
    EUY <- y*u
    EUY2 <- (y^2)*u
    EUT <- u*(muT_eta) + MT*E
    EUT2 <- u*(muT_eta^2) + M2T + MT*(muT_eta + eta)*E
    EUYT <- y*(u*muT_eta + MT*E)
    ElogU <- log(u)
  }
  
  return(list(EU=EU,EUY=EUY,EUY2=EUY2,EUT=EUT,EUT2=EUT2,EUYT=EUYT,ElogU = ElogU))
}

## -------------------- ##
## EXPECTATION CENSORED ##
## -------------------- ##

CensEsp <- function(y,X,betas,sigma2,lambda,nu,cc,cens,LS,N,dist)
{
  # cc: vetor 0 ou 1
  
  deltinha <- lambda/sqrt(1+(lambda^2))
  Delta <-  sqrt(sigma2)*deltinha
  tau <- sigma2 - Delta^2 
  
  if(cens=="left"){cens = 1}
  if(cens=="right"){cens = 2}
  if(cens=="interval"){cens = 3}
  
  if(dist=="SN")
  {
    eta <- -sqrt(2/pi)
    mu <- X%*%betas + Delta*eta
    mu1 <- mu[cc==1]
    ys <- matrix(ncol=length(mu1),nrow=N)
    
    if(cens=="1")
    {
      a <- rep(-Inf,length(mu1))
      b <- y[cc==1]
    }
    if(cens=="2")
    {
      a <- y[cc==1]
      b <- rep(Inf,length(mu1))
    }
    if(cens=="3")
    {
      a <- y[cc==1]
      b <- LS[cc==1]
    }
    
    for(i in 1:length(mu1))
    {             
      ys[,i] <- rTSN(N,mu1[i],sigma2,lambda,a[i],b[i],cens)
    }  
    
    muT <- (Delta/(tau+Delta^2))*(ys -kronecker(matrix(mu1,1,ncol(ys)),matrix(1,nrow(ys),1)))
    M2T <- tau/(tau + Delta^2)
    MT <- sqrt(M2T)
    Wphi <- dnorm(muT/MT)/pnorm(muT/MT)
    
    EUY <- apply(ys,2,mean)
    EUY2 <- apply(ys^2,2,mean)
    
    auxT1<-(muT + eta) + MT*Wphi
    EUT <- apply(auxT1,2,mean)
    auxT2 <- (muT + eta)^2 + M2T + MT*(muT + 2*eta)*Wphi
    EUT2 <- apply(auxT2,2,mean)
    auxTY <- ys*((muT + eta) + MT*Wphi)
    EUYT <- apply(auxTY,2,mean)
    
    EU <- cbind(rep(1,length(mu1)))
    ElogU <- log(EU)
  }
  
  if(dist=="ST")
  {
    k1 <- sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
    eta <- -sqrt(2/pi)*k1
    mu <- X%*%betas + Delta*eta
    mu1 <- mu[cc==1]
    
    ys <- matrix(ncol=length(mu1),nrow=N)
    
    if(cens=="1")
    {
      a <- rep(-Inf,length(mu1))
      b <- y[cc==1]
    }
    if(cens=="2")
    {
      a <- y[cc==1]
      b <- rep(Inf,length(mu1))
    }
    if(cens=="3")
    {
      a <- y[cc==1]
      b <- LS[cc==1]
    }
    
    for(i in 1:length(mu1))
    {             
      ys[,i] <- rTSt(N,mu1[i],sigma2,lambda,nu,a[i],b[i],cens)
    }
    
    M2T <- 1/(1 + (Delta^2)*(tau^(-1)))
    MT <- sqrt(M2T)
    muT <- M2T*Delta*(tau^(-1))*(ys - kronecker(matrix(mu1,1,ncol(ys)),matrix(1,nrow(ys),1)))
    muT_eta <- muT + eta
    A <- muT/MT
    auxd <- (ys - kronecker(matrix(mu1,1,ncol(ys)),matrix(1,nrow(ys),1)))
    d <- (auxd^2)/sigma2  
    
    E <- U <- matrix(0,nrow=dim(ys)[1], ncol=dim(ys)[2])
    for(i in 1:length(mu1))
    {
      E[,i] <- (2*(nu)^(nu/2)*gamma((2+nu)/2)*((d[,i] + nu + (A[,i])^2))^(-(2+nu)/2))/(gamma(nu/2)*pi*sqrt(sigma2)*pdfSNI(ys[,i],mu1[i],sigma2,lambda,nu,type="ST"))
      U[,i] <- ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(d[,i] + nu)^(-(nu+3)/2))/(gamma(nu/2)*sqrt(pi)*sqrt(sigma2)*pdfSNI(ys[,i],mu1[i],sigma2,lambda,nu, type="ST")))*pt(sqrt((3+nu)/(d[,i]+nu))*A[,i],3+nu)
    }
    
    auxUY <- ys*U
    auxUY2 <-  (ys^2)*U
    auxUT <-  U*(muT_eta) +  MT*E
    auxUT2 <-  U*(muT_eta^2) + M2T + MT*(muT_eta + eta)*E
    auxUYT <- ys*(U*muT_eta + MT*E)
    auxlogU <- log(U)
    
    EU <- apply(U,2,mean)
    EUY <- apply(auxUY,2,mean)
    EUY2 <- apply(auxUY2,2,mean)
    EUT <- apply(auxUT,2,mean)
    EUT2 <- apply(auxUT2,2,mean)
    EUYT <- apply(auxUYT,2,mean)
    ElogU <- apply(auxlogU,2,mean)
  }
  
  if(dist=="SSL")
  {
    k1 <- 2*nu/(2*nu-1)
    eta <- -sqrt(2/pi)*k1
    mu <- X%*%betas + Delta*eta
    mu1 <- mu[cc==1]
    
    ys <- matrix(ncol=length(mu1),nrow=N)
    
    if(cens=="1")
    {
      a <- rep(-Inf,length(mu1))
      b <- y[cc==1]
    }
    if(cens=="2")
    {
      a <- y[cc==1]
      b <- rep(Inf,length(mu1))
    }
    if(cens=="3")
    {
      a <- y[cc==1]
      b <- LS[cc==1]
    }
    
    for(i in 1:length(mu1))
    {
      ys[,i] <- rTSSl(N,mu1[i],sigma2,lambda,nu,a[i],b[i],cens)
    } 
    
    auxd <- (ys - kronecker(matrix(mu1,1,ncol(ys)),matrix(1,nrow(ys),1)))
    d <- (auxd^2)/sigma2 
    M2T <- 1/(1 + (Delta^2)*(tau^(-1)))
    MT <- sqrt(M2T) 
    muT <- M2T*Delta*(tau^(-1))*(ys - kronecker(matrix(mu1,1,ncol(ys)),matrix(1,nrow(ys),1)))
    muT_eta <- muT + eta
    A <- muT/ MT
    
    E <- u <- matrix(0,nrow=nrow(ys), ncol=ncol(ys))
    for(i in 1:ncol(ys))
    {
      for(j in 1:nrow(ys))
      {
        E[j,i] <- (((2^(nu + 1))*nu*gamma(nu + 1))/(pdfSNI(ys[j,i],mu1[i],sigma2,lambda,nu,type="SSL")*pi*sqrt(sigma2)))* ((d[j,i]+A[j,i]^2)^(-nu-1))*pgamma(1,nu+1,(d[j,i]+A[j,i]^2)/2) 
        faux <- function(u) u^(nu+0.5)*exp(-u*d[j,i]/2)*pnorm(u^(1/2)*A[j,i])
        aux22 <- integrate(faux,0,1)$value
        u[j,i] <- ((sqrt(2)*nu) / (pdfSNI(ys[j,i],mu1[i],sigma2,lambda,nu,type="SSL")*sqrt(pi)*sqrt(sigma2)))*aux22
      }
    }
    
    auxUY <- ys*u
    auxUY2 <-  (ys^2)*u
    auxUT <-  u*(muT_eta) +  MT*E
    auxUT2 <-  u*(muT_eta^2) + M2T + MT*(muT_eta + eta)*E
    auxUYT <- ys*(u*muT_eta + MT*E)
    auxlogU <- log(u) 
    
    EU <- apply(u,2,mean)
    EUY <- apply(auxUY,2,mean)
    EUY2 <- apply(auxUY2,2,mean)
    EUT <- apply(auxUT,2,mean)
    EUT2 <- apply(auxUT2,2,mean)
    EUYT <- apply(auxUYT,2,mean)
    ElogU <- apply(auxlogU,2,mean)
    
  }
  
  if(dist=="SCN")
  {
    k1 <- nu[1]/nu[2]^(1/2)+1-nu[1]
    eta <- -sqrt(2/pi)*k1
    mu <- X%*%betas + Delta*eta
    mu1 <- mu[cc==1]
    
    ys <- matrix(ncol=length(mu1),nrow=N)
    
    if(cens=="1")
    {
      a <- rep(-Inf,length(mu1))
      b <- y[cc==1]
    }
    if(cens=="2")
    {
      a <- y[cc==1]
      b <- rep(Inf,length(mu1))
    }
    if(cens=="3")
    {
      a <- y[cc==1]
      b <- LS[cc==1]
    }
    
    for(i in 1:length(mu1))
    {
      ys[,i] <- rTSCN(N,mu1[i],sigma2,lambda,nu,a[i],b[i],cens)
    }  
    
    auxd <- (ys - kronecker(matrix(mu1,1,ncol(ys)),matrix(1,nrow(ys),1)))
    d <- (auxd^2)/sigma2 
    M2T <- 1/(1 + (Delta^2)*(tau^(-1)))
    MT <- sqrt(M2T) 
    muT <- M2T*Delta*(tau^(-1))*(ys - kronecker(matrix(mu1,1,ncol(ys)),matrix(1,nrow(ys),1)))
    muT_eta <- muT + eta
    A <- muT/ MT
    
    u <- E <- matrix(0,nrow=nrow(ys), ncol=ncol(ys))
    for(i in 1:ncol(ys))
    {
      u[,i] <- (2/pdfSNI(ys[,i],mu1[i],sigma2,lambda,nu,type="SCN"))*(nu[1]*nu[2]*dnorm(ys[,i],mu1[i],sqrt(sigma2/nu[2]))*pnorm(sqrt(nu[2])*A[,i],0,1) + (1-nu[1])*dnorm(ys[,i],mu1[i],sqrt(sigma2))*pnorm(A[,i],0,1))
      E[,i] <- (2/pdfSNI(ys[,i],mu1[i],sigma2,lambda,nu,type="SCN"))*(nu[1]*sqrt(nu[2])*dnorm(ys[,i],mu1[i],sqrt(sigma2/nu[2]))*dnorm(sqrt(nu[2])*A[,i],0,1) + (1-nu[1])*dnorm(ys[,i],mu1[i],sqrt(sigma2))*dnorm(A[,i],0,1))
    } 
    
    auxUY <- ys*u
    auxUY2 <-  (ys^2)*u
    auxUT <-  u*(muT_eta) +  MT*E
    auxUT2 <-  u*(muT_eta^2) + M2T + MT*(muT_eta + eta)*E
    auxUYT <- ys*(u*muT_eta + MT*E)
    
    EU <- apply(u,2,mean)
    EUY <- apply(auxUY,2,mean)
    EUY2 <- apply(auxUY2,2,mean)
    EUT <- apply(auxUT,2,mean)
    EUT2 <- apply(auxUT2,2,mean)
    EUYT <- apply(auxUYT,2,mean)
    ElogU <- log(EU)
  }
  
  return(list(EU=EU,EUY=EUY,EUY2=EUY2,EUT=EUT,EUT2=EUT2,EUYT=EUYT, ElogU = ElogU)) 
}

## --------------------- ##
## LOG-LIKELIHOOD - SMSN ##
## --------------------- ##

logVerosCens <- function(theta,y,X,cc,cens,dist)
{
  
  if(cens=="left"){cens = 1}
  if(cens=="right"){cens = 2}
  
  p <- ncol(X)
  n <- nrow(X)
  betas <- theta[1:p]
  sigma2 <- theta[(p+1)]
  lambda <- theta[(p+2)]
  
  deltinha <- lambda/sqrt(1+(lambda^2))
  Delta <-  sqrt(sigma2)*deltinha
  
  ver <- matrix(0,n,1)
  auxy <- matrix(0,n,1)
  
  if(dist=="SN")
  {
    eta <- -sqrt(2/pi)
    mu <- X%*%betas + eta*Delta
    
    ver[cc==0] <- pdfSNI(y[cc==0],mu[cc==0],sigma2,lambda,nu=NULL,type="SN")
    if(cens=="1")
    {
      if(sum(cc)>0)
      {
        ver[cc==1] <- cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu=NULL,type="SN")
      }
    }
    if(cens=="2")
    {
      if(sum(cc)>0)
      {
        ver[cc==1] <- (1-cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu=NULL,type="SN"))
        
      }
    }
    if(length(which(ver == 0)) > 0) ver[which(ver == 0)] <- .Machine$double.xmin
  }
  
  if(dist=="ST")
  {
    nu <- theta[(p+3)]
    k1 <- sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
    eta <- -sqrt(2/pi)*k1
    
    mu <- X%*%betas + eta*Delta
    
    ver[cc==0] <- pdfSNI(y[cc==0],mu[cc==0],sigma2,lambda,nu,type="ST")
    if(cens==1)
    {
      if(sum(cc)>0)
      {
        ver[cc==1] <- cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu,type="ST")
      }
    }
    if(cens==2)
    {
      if(sum(cc)>0)
      {
        ver[cc==1] <- (1-cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu,type="ST"))
        
      }
    }
    if(length(which(ver == 0)) > 0) ver[which(ver == 0)] <- .Machine$double.xmin
  }
  
  if(dist=="SSL")
  {
    nu <- theta[(p+3)]
    k1 <- 2*nu/(2*nu-1)
    eta <- -sqrt(2/pi)*k1
    
    mu <- X%*%betas + Delta*eta
    
    for(i in 1:length(y[cc==0]))
    {
      ver[cc==0][i] <- pdfSNI(y[cc==0][i],mu[cc==0][i],sigma2,lambda,nu,type="SSL")
    }
    
    if(cens=="1")
    {
      if(sum(cc)>0)
      {
        for(i in 1:length(y[cc==1]))
        {
          ver[cc==1][i] <- cdfSNI(y[cc==1][i],mu[cc==1][i],sigma2,lambda,nu,type="SSL")$cdf
        }
      }
    }
    if(cens=="2")
    {
      if(sum(cc)>0)
      {
        for(i in 1:length(y[cc==1]))
        {
          ver[cc==1][i] <- (1-cdfSNI(y[cc==1][i],mu[cc==1][i],sigma2,lambda,nu,type="SSL")$cdf)
        }
      }
    }
    if(length(which(ver == 0)) > 0) ver[which(ver == 0)] <- .Machine$double.xmin
    
  }
  
  if(dist=="SCN")
  {
    nu1 <- theta[(p+3)]
    nu2 <- theta[(p+4)]
    
    nu <- rbind(nu1,nu2)
    
    k1 <- (nu1/(nu2^(1/2))) +1-nu1
    eta <- -sqrt(2/pi)*k1
    
    mu <- X%*%betas + Delta*eta
    
    ver[cc==0] <- pdfSNI(y[cc==0], mu[cc==0], sigma2, lambda, nu, type="SCN") 
    if(cens==1)
    {
      if(sum(cc)>0)
      {
        ver[cc==1] <- cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu,type="SCN") 
      }
    }
    if(cens==2)
    {
      if(sum(cc)>0)
      {
        ver[cc==1] <- (1-cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu,type="SCN")) 
      }
    }
    if(length(which(ver == 0)) > 0) ver[which(ver == 0)] <- .Machine$double.xmin
  }  
  
  logver <- sum(log(ver))
  return(logver)
}

## ------------------------------------- ##
## LOG-LIKELIHOOD - SMSN - ESTIMATION NU ##
## ------------------------------------- ##

logVerosCens.nu <- function(nu,theta,y,X,cc,cens,dist)
{
  
  if(cens=="left"){cens = 1}
  if(cens=="right"){cens = 2}
  
  p <- ncol(X)
  n <- nrow(X)
  betas <- theta[1:p]
  sigma2 <- theta[(p+1)]
  lambda <- theta[(p+2)]
  
  deltinha <- lambda/sqrt(1+(lambda^2))
  Delta <-  sqrt(sigma2)*deltinha
  
  ver <- matrix(0,n,1)
  auxy <- matrix(0,n,1)
  
  #if(dist=="SN")
  #{
  #  eta <- -sqrt(2/pi)
  #  mu <- X%*%betas + eta*Delta
  #  
  #  ver[cc==0] <- pdfSNI(y[cc==0],mu[cc==0],sigma2,lambda,nu=NULL,type="SN")
  #  if(cens=="1")
  #  {
  #    if(sum(cc)>0)
  #    {
  #      ver[cc==1] <- cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu=NULL,type="SN")
  #    }
  #  }
  #  if(cens=="2")
  #  {
  #    if(sum(cc)>0)
  #    {
  #     ver[cc==1] <- (1-cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu=NULL,type="SN"))
  #      
  #    }
  #  }
  #  if(length(which(ver == 0)) > 0) ver[which(ver == 0)] <- .Machine$double.xmin
  #}
  
  if(dist=="ST")
  {
    nu <- nu
    k1 <- sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
    eta <- -sqrt(2/pi)*k1
    
    mu <- X%*%betas + eta*Delta
    
    ver[cc==0] <- pdfSNI(y[cc==0],mu[cc==0],sigma2,lambda,nu,type="ST")
    if(cens==1)
    {
      if(sum(cc)>0)
      {
        ver[cc==1] <- cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu,type="ST")
      }
    }
    if(cens==2)
    {
      if(sum(cc)>0)
      {
        ver[cc==1] <- (1-cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu,type="ST"))
        
      }
    }
    if(length(which(ver == 0)) > 0) ver[which(ver == 0)] <- .Machine$double.xmin
  }
  
  if(dist=="SSL")
  {
    nu <- nu
    k1 <- 2*nu/(2*nu-1)
    eta <- -sqrt(2/pi)*k1
    
    mu <- X%*%betas + Delta*eta
    
    for(i in 1:length(y[cc==0]))
    {
      ver[cc==0][i] <- pdfSNI(y[cc==0][i],mu[cc==0][i],sigma2,lambda,nu,type="SSL")
    }
    
    if(cens=="1")
    {
      if(sum(cc)>0)
      {
        for(i in 1:length(y[cc==1]))
        {
          ver[cc==1][i] <- cdfSNI(y[cc==1][i],mu[cc==1][i],sigma2,lambda,nu,type="SSL")$cdf
        }
      }
    }
    if(cens=="2")
    {
      if(sum(cc)>0)
      {
        for(i in 1:length(y[cc==1]))
        {
          ver[cc==1][i] <- (1-cdfSNI(y[cc==1][i],mu[cc==1][i],sigma2,lambda,nu,type="SSL")$cdf)
        }
      }
    }
    if(length(which(ver == 0)) > 0) ver[which(ver == 0)] <- .Machine$double.xmin
    
  }
  
  if(dist=="SCN")
  {
    nu1 <- nu[1]
    nu2 <- nu[2]
    
    k1 <- nu1/nu2^(1/2)+1-nu1
    eta <- -sqrt(2/pi)*k1
    
    mu <- X%*%betas + Delta*eta
    
    ver[cc==0] <- pdfSNI(y[cc==0], mu[cc==0], sigma2, lambda, nu, type="SCN") 
    if(cens==1)
    {
      if(sum(cc)>0)
      {
        ver[cc==1] <- cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu,type="SCN") 
      }
    }
    if(cens==2)
    {
      if(sum(cc)>0)
      {
        ver[cc==1] <- (1-cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu,type="SCN")) 
      }
    }
    if(length(which(ver == 0)) > 0) ver[which(ver == 0)] <- .Machine$double.xmin
  }  
  
  logver <- sum(log(ver))
  return(-logver)
}

## ----------------------------------- ##
## GENERATION OF CENSORED SMSN SAMPLES ##
## ----------------------------------- ##

generate_SMSNCR <- function(X,betas,sigma2,lambda,n,cens,perc,dist,nu)
{
  deltinha <- lambda/(sqrt(1 + lambda^2))
  Delta <- sqrt(sigma2)*deltinha
  
  if(dist=="SN")
  {
    eta <- -sqrt(2/pi)
  }
  if(dist=="ST")
  { 
    k1 <- sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
    eta <- -sqrt(2/pi)*k1
  }
  if(dist=="SSL")
  { 
    k1 <- 2*nu/(2*nu-1)
    eta <- -sqrt(2/pi)*k1
  }
  if(dist=="SCN")
  {
    k1 <- (nu[1]/(nu[2]^(1/2))) + 1-nu[1]
    eta <- -sqrt(2/pi)*k1
  }
  
  mu <-  eta*Delta 
  error <- rSMSN(n=n,mu=mu,sigma2=sigma2,lambda=lambda,nu=nu,dist=dist)
  y <- X%*%betas + error
  
  yc <- y
  
  if(perc==0) cc <- rep(0,n)
  
  if(perc > 0)
  {
    if(cens=="left")
    {
      aa <- sort(yc, decreasing=FALSE)
      cutof <- aa[ceiling(perc*n)]
      cc <- matrix(1,n,1)*(yc <= cutof)
      yc[cc==1] <- cutof
    }
    if(cens=="right")
    {
      aa <- sort(yc, decreasing=TRUE)
      cutof <- aa[ceiling(perc*n)]
      cc <- matrix(1,n,1)*(yc >= cutof)
      yc[cc==1] <- cutof
    }
  }    
  return(list(y=y,yc=yc,cc=cc))  
}

##-------------------------------##
## Martingale type residual (MT) ##
##-------------------------------##

library(PerformanceAnalytics)

resMT <- function(theta,y,X,cc,cens,dist)
{
  n <- length(y)
  p <- ncol(X)
  resm <- resmt <- S <- c()
  
  loglikelihood <- logVerosCens(theta,y,X,cc,cens,dist)
  npar <- (p + 2)
  AIC <- -(2)*loglikelihood + 2*npar
  BIC <- -(2)*loglikelihood + (npar)*log(n)
  EDC <- -(2)*loglikelihood + (npar)*0.2*sqrt(n) 
  
  cat('-----------\n')
  cat('log-likelihood =',round(loglikelihood,4), '\n')
  cat('-----------\n')
  cat('AIC =', round(AIC,4), '\n')
  cat('-----------\n')
  cat('BIC =', round(BIC,4), '\n')
  cat('-----------\n')
  cat('EDC =', round(EDC,4), '\n')
  cat('-----------\n')
  
  if(dist=="SN")
  {
    betas <- theta[1:p]
    sigma2 <- theta[(p+1)]
    lambda <- theta[(p+2)]
    
    Delta <- sqrt(sigma2)*(lambda/sqrt(lambda^2 + 1))
    mu <- X%*%betas  - sqrt(2/pi)*Delta
    S <- 1 - cdfSNI(y,mu,sigma2,lambda,nu=NULL,type=dist)
  }
  if(dist=="ST")
  {
    betas <- theta[1:p]
    sigma2 <- theta[(p+1)]
    lambda <- theta[(p+2)]
    nu <- theta[(p+3)]
    
    Delta <- sqrt(sigma2)*(lambda/sqrt(lambda^2 + 1))
    k1 <- sqrt(nu/2)*(gamma((nu-1)/2)/gamma(nu/2))
    mu <- X%*%betas - sqrt(2/pi)*k1*Delta
    S <- 1 - cdfSNI(y,mu,sigma2,lambda,nu,type=dist)
  }
  if(dist=="SSL")
  {
    betas <- theta[1:p]
    sigma2 <- theta[(p+1)]
    lambda <- theta[(p+2)]
    nu <- theta[(p+3)]
    
    Delta <- sqrt(sigma2)*(lambda/sqrt(lambda^2 + 1))
    k1 <- (2*nu)/(2*nu - 1)
    mu <- X%*%betas - sqrt(2/pi)*k1*Delta
    for(i in 1:n)
    {
      S[i] <- 1 - cdfSNI(y[i],mu[i],sigma2,lambda,nu,type="SSL")$cdf
    }
  }
  if(dist=="SCN")
  {
    betas <- theta[1:p]
    sigma2 <- theta[(p+1)]
    lambda <- theta[(p+2)]
    nu1 <- theta[(p+3)]
    nu2 <- theta[(p+4)]
    
    nu <- rbind(nu1,nu2)
    
    Delta <- sqrt(sigma2)*(lambda/sqrt(lambda^2 + 1))
    k1 <- (nu1/nu2^(1/2))+1-nu1
    mu <- X%*%betas - sqrt(2/pi)*k1*Delta
    S <- 1 - cdfSNI(y,mu,sigma2,lambda,nu,type="SCN")
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

EnvelopeRMT <- function(theta,y,X,cc,cens,dist)
{
  res <- resMT(theta,y,X,cc,cens,dist)
  windows()
  if(dist=="SN")
  {
    chart.QQPlot(res$resmt,distribution='norm',envelope=0.95,ylab="Deviance residuals",xlab="Standard normal quantile",main="SN-CR",line="quartile", col=c(20,1),pch=19)
  }
  if(dist=="ST")
  {
    chart.QQPlot(res$resmt,distribution='norm',envelope=0.95,ylab="Deviance residuals",xlab="Standard normal quantile",main="ST-CR",line="quartile", col=c(20,1),pch=19)
  }
  if(dist=="SSL")
  {
    chart.QQPlot(res$resmt,distribution='norm',envelope=0.95,ylab="Deviance residuals",xlab="Standard normal quantile",main="SSL-CR",line="quartile", col=c(20,1),pch=19)
  }
  if(dist=="SCN")
  {
    chart.QQPlot(res$resmt,distribution='norm',envelope=0.95,ylab="Deviance residuals",xlab="Standard normal quantile",main="SCN-CR",line="quartile", col=c(20,1),pch=19)
  }
}