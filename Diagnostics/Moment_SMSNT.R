### Remove all
#rm(list=ls())

##################################################
### Aldo M. Garay - Victor Hugo Lachos Davila
### 09/10    
#install.packages("mvtnorm")
#install.packages("mvnorm")
#install.packages("hbmem")
#install.packages("truncdist")
#install.packages("mnormt")
#install.packages("sn")
#library(mvtnorm)
#library(mnormt)
#library(hbmem)
#library(truncdist)             
#library(mnormt)
#library(sn)

if(!require(mvtnorm))   install.packages("mvtnorm");   library(mvtnorm)
if(!require(mnormt))    install.packages("mnormt");    library(mnormt)
if(!require(hbmem))     install.packages("hbmem");     library(hbmem)
if(!require(truncdist)) install.packages("truncdist"); library(truncdist)
if(!require(sn))        install.packages("sn");        library(sn)
if(!require(moments))   install.packages("moments");   library(moments)

##################################################
##################################################
######### E_phi e E_Phi - Simetricas #############
##################################################
##################################################

E_Phi <- function(r,a,nu,delta,type=type)
{
  n <- length(a)
  if(type=="N")
  {
    resp <- pnorm(a)
  }        
  if(type=="T")
  {
    Aux0 <- gamma(0.5*(nu+(2*r)))
    Aux1 <- gamma(nu/2) 
    Aux2 <- Aux0/Aux1
    Aux3 <- (0.5*nu)^(-r)
    Aux4 <- AcumPearsonVII(a,0,1,nu+(2*r),nu)
    resp <- Aux2*Aux3*Aux4
  }
  if(type=="PVII")
  {
    Aux0 <- gamma(0.5*(nu+(2*r)))
    Aux1 <- gamma(nu/2) 
    Aux2 <- Aux0/Aux1
    Aux3 <- (0.5*delta)^(-r)
    Aux4 <- AcumPearsonVII(a,0,1,nu+(2*r),delta)
    resp <- Aux2*Aux3*Aux4  
  }
  if(type=="SL")
  {
    Aux0 <- nu/(nu+r)
    Aux1 <- AcumSlash(a,0,1,nu+r)
    resp <- Aux0*Aux1
  }        
  if(type=="CN")
  {
    Aux0 <- (nu[2]^r)*AcumNormalC(a,0,1,nu)
    Aux1 <- (1-nu[2]^r)*(1-nu[1])*pnorm(a)
    resp <- Aux0 + Aux1
  } 
  return(resp)
}

################################################################################
################################################################################

E_phi <- function(r,a,nu,delta,type=type)
{
  n <- length(a)
  b <- rep(Inf,n)
  b1<- rep(-Inf,n)
  
  if(setequal(a,b)== TRUE | setequal(a,b1)== TRUE)
  {
    resp <- rep(0,n)
  }
  else
  {
    if(type=="N")
    {
      resp <- dnorm(a)
    }
    if(type=="T")
    {
      Aux0 <- gamma(0.5*(nu+2*r))
      Aux1 <- gamma(nu/2)*sqrt(2*pi) 
      Aux2 <- Aux0/Aux1
      Aux3 <- (0.5*nu)^(nu/2)
      Aux4 <- (0.5*(a^2+nu))^(-0.5*(nu+2*r))
      resp <- Aux2*Aux3*Aux4 
    }
    if(type=="PVII")
    {
      Aux0 <- gamma(0.5*(nu+2*r))
      Aux1 <- gamma(nu/2)*sqrt(2*pi) 
      Aux2 <- Aux0/Aux1
      Aux3 <- (0.5*delta)^(nu/2)
      Aux4 <- (0.5*(a^2+delta))^(-0.5*(nu+2*r))
      resp <- Aux2*Aux3*Aux4 
    }
    if(type=="SL")
    {
      Aux0 <- nu/sqrt(2*pi) 
      Aux1 <- (0.5*a^2)^(-(nu+r))
      Aux2 <- GamaInc(nu+r,0.5*a^2)
      resp <- Aux0*Aux1*Aux2 
    } 
    if(type=="CN")
    {
      Aux0 <- nu[1]*(nu[2]^r)*dnorm(a*sqrt(nu[2]))
      Aux1 <- (1-nu[1])*dnorm(a)
      resp <- Aux0 + Aux1
    }
  }      
  return(resp)
}

##################################################
##################################################
######### E_phi e E_Phi - Assimetricas ###########
##################################################
##################################################

E_phiSNI <- function(r,a,nu,delta,lambda,type=type)
{
  n <- length(a)
  b <- rep(Inf,n)
  b1<- rep(-Inf,n)
  if(setequal(a,b)== TRUE | setequal(a,b1)== TRUE)
  {
    resp <- rep(0,n)
  }
  else
  {
    if(type=="SN")
    {
      resp <- 2*dnorm(a)*pnorm(lambda*a)
    }
    if(type=="ST")
    {
      Aux0 <- (2^(r+1))*(nu^(nu/2))*gamma(0.5*nu+r)
      Aux1 <- gamma(nu/2)*sqrt(2*pi) 
      Aux2 <- (a^2+nu)^(0.5*nu+r)
      Aux3 <- Aux0/(Aux1*Aux2)
      Aux4 <- sqrt((2*r+nu)/(a^2+nu))*lambda*a
      #Aux5 <- cdfSNI(Aux4,mu=0,sigma2=1,lambda=0,nu=(2*r+nu),type="ST")
      Aux5 <- pt(Aux4,df=(2*r+nu))
      resp <- Aux3*Aux5
    }
    if(type=="SPVII")
    {
      Aux0 <- gamma(0.5*(nu+2*r))
      Aux1 <- gamma(nu/2)*sqrt(2*pi) 
      Aux2 <- Aux0/Aux1
      Aux3 <- (0.5*delta)^(nu/2)
      Aux4 <- (0.5*(a^2+delta))^(-0.5*(nu+2*r))
      resp <- Aux2*Aux3*Aux4 
    }
    if(type=="SSL")     ##Observ "mu" nao pode ser igual a "Lim1"
    {
      Aux0 <- (2^(nu+r+1))*nu
      Aux1 <- gamma(nu+r)
      Aux2 <- (a^(2*(nu+r)))*(sqrt(2*pi))
      resp <- (Aux0*Aux1/Aux2)*Mont_SSL(m=20000,r,a,lambda,nu, case="1")*pgamma(1, shape=(r+nu), scale = 1/(0.5*a^2)) 
    } 
    if(type=="SCN")
    {
      Aux0 <- nu[1]*(nu[2]^r)*(2*dnorm(a*sqrt(nu[2]))*pnorm(lambda*a*sqrt(nu[2])))        
      Aux1 <- (1-nu[1])*(2*dnorm(a)*pnorm(lambda*a))        
      resp <- Aux0 + Aux1
    }
  }      
  return(resp)
}

E_PhiSNI <- function(r,a,nu,delta,lambda,type=type)
{
  deltinha <- lambda/sqrt(1+lambda^2)
  n <- length(a)
  #         b <- rep(Inf,n)
  #         b1<- rep(-Inf,n)
  #  if(setequal(a,b)== TRUE | setequal(a,b1)== TRUE)
  #  {
  #    resp <- rep(0,n)
  #  }
  #  else
  #  {
  if(type=="SN")
  {
    resp <- cdfSNI(a,mu=0,sigma2=1,lambda,nu=NUL,type=type)
  }        
  if(type=="ST")
  {
    Aux0 <- 2^(r+1)
    Aux1 <- gamma(0.5*nu+r)
    Aux2 <- gamma(nu/2) 
    Aux3 <- nu^(r)
    Aux4 <- (Aux0*Aux1)/(Aux2*Aux3)
    Aux5 <- sqrt((2*r+nu)/nu)
    valor <- Aux5*c(a,0)
    mean1 <- c(0,0)
    Sigma1 <- matrix(c(1,-deltinha,-deltinha,1), 2, 2)
    GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
    Aux6 <- pmvt(lower = -Inf, upper = valor-mean1, sigma = Sigma1,  df=round(2*r+nu), algorithm = GB)[1]
    #GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
    #Aux6 <- pmvt(lower = -Inf, upper = valor-mean1, sigma = Sigma1,  df=2*r+nu, algorithm = GB)[1]
    #Aux6 <- pmt(valor, mean1, Sigma1, df=2*r+nu)
    #Aux6 <- acumt2(NULL,valor,mean1,Sigma1,2*r+nu)
    resp <- Aux4*Aux6
  }
  if(type=="SPVII")
  {
    Aux0 <- gamma(0.5*(nu+(2*r)))
    Aux1 <- gamma(nu/2) 
    Aux2 <- Aux0/Aux1
    Aux3 <- (0.5*delta)^(-r)
    Aux4 <- AcumPearsonVII(a,0,1,nu+(2*r),delta)
    resp <- Aux2*Aux3*Aux4  
  }
  if(type=="SSL")
  {
    Aux0 <- 2*nu
    Aux1 <- (r+nu)
    resp <- (Aux0/Aux1)*Mont_SSL(m=20000,r,a,lambda,nu, case="2") 
  }        
  if(type=="SCN")
  {
    Aux0 <- (nu[2]^r)*cdfSNI(a,mu=0,sigma2=1,lambda,nu,type="SCN")
    Aux1 <- 2*(1-nu[2]^r)*(1-nu[1])
    valor <- c(a,0)
    mean1 <- c(0,0)
    Sigma1 <- matrix(c(1,-deltinha,-deltinha,1), 2, 2)
    if(a==Inf)
    {
      valor=c(500,0)
    }  
    Aux3 <- pmnorm(valor, mu=c(0,0), varcov=Sigma1)
    resp <- Aux0 + (Aux1*Aux3)
  } 
  #  }
  return(resp)
}

# E_PhiSNI(2,2,nu = 3,1,2,type="ST")
# E_PhiSNI(2,2,nu = 3.8,1,2,type="ST")
# E_PhiSNI(2,2,nu = 4,1,2,type="ST")

################################################################
##########          Densidades das SNI              ############

## Densidade/CDF da SN com locação escala #######
dSN <- function(y, mu = 0, sigma2 = 1, shape=1){
  dens <- 2*dnorm(y, mu, sqrt(sigma2))*pnorm(shape*((y - mu)/sqrt(sigma2)))
  return(dens)
}



## Densidade/CDF da ST com locação escala #######
dt.ls <- function(x, loc = 0, sigma2 = 1,shape=1, nu = 4){
  d <- (x - loc)/sqrt(sigma2)
  dens <- 2*dt(d, df = nu)*pt(sqrt((1+nu)/(d^2+nu))*d*shape,1+nu)/sqrt(sigma2)
  return(dens)
}



## Densidade/CDF da Skew Normal Contaminada #######
  dSNC <- function(y, mu, sigma2, shape, nu){
    dens <- 2*(nu[1]*dnorm(y, mu, sqrt(sigma2/nu[2]))*pnorm(sqrt(nu[2])*shape*sigma2^(-1/2)*(y-mu)) + (1 - nu[1])*dnorm(y, mu, sqrt(sigma2))*pnorm(shape*sigma2^(-1/2)*(y-mu)))
    return(dens)
  }


### Densidade da Skew Slash  ######
dSS <- function(y, mu, sigma2, shape,nu){
  resp <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y)) {
    f <- function(u) 2*nu*u^(nu - 1)*dnorm(y[i],mu,sqrt(sigma2/u))*pnorm(u^(1/2)*shape*(sigma2^(-1/2))*(y[i]-mu))
    resp[i] <- integrate(f,0,1)$value
  }
  return(resp)
}


dT <- function(cc, y, mu, sigma2 = 1, nu=4){
  densN<- vector(mode = "numeric", length = length(y))
  aux<-(y-mu)/sqrt(sigma2)
  densN[cc==0] <- dt(aux[cc==0],nu)/sqrt(sigma2)
  densN[cc==1]<- pt(aux[cc==1],nu)
  if(length(which(densN == 0)) > 0) densN[which(densN == 0)] <- .Machine$double.xmin
  return(densN)
}


################################################################
##########          Acumuladas das SNI              ############
################################################################
 
cdfSNI<- function(x, mu, sigma2, lambda, nu, type = "SN"){
  n <- length(x)
  resp <- matrix(0, n, 1)
  if (type == "N") {
    resp <- pnorm((x - mu)/sqrt(sigma2))
    return(resp)
  }
  
  if (type == "T") {
    resp <- pt((x - mu)/sqrt(sigma2), df = nu)
    return(resp)
  }
  
  if (type == "SN") {
    delta <- lambda/sqrt(1 + lambda^2)
    SIGMA <- matrix(c(sigma2, -delta * sqrt(sigma2), -delta * sqrt(sigma2), 
                      1), byrow = TRUE, ncol = 2, nrow = 2)
    if (length(mu) == 1) {
      MU <- cbind(rep(mu, n), 0)
    }
    if (length(mu) == n) {
      MU <- cbind(mu, 0)
    }
    Y <- cbind(x, 0)
    for (i in 1:n) {
      resp[i] <- 2 * pmnorm(x = Y[i, ], mean = MU[i, ], varcov = SIGMA)
    }
    return(resp)
  }
  
  if (type == "ST") {
    delta <- lambda/sqrt(1 + lambda^2)
    SIGMA <- matrix(c(sigma2, -delta * sqrt(sigma2), -delta * sqrt(sigma2), 
                      1), byrow = TRUE, ncol = 2, nrow = 2)
    if (length(mu) == 1) {
      MU <- cbind(rep(mu, n), 0)
    }
    if (length(mu) == n) {
      MU <- cbind(mu, 0)
    }
    Y <- cbind(x, 0)
    #nu <- round(nu)
    if(nu%%1 == 0){
      #nu integer
      for (i in 1:n){
        resp[i] <- 2 * pmt(x = Y[i, ], mean = MU[i, ], S = SIGMA, df = nu)
      }
    }else{
      #nu decimal
      for (i in 1:n) {
        resp[i] <- 2 * acumt2(a = NULL,b = Y[i, ],mu = MU[i, ],Sigma = SIGMA,nu = nu)
      }
    }
    return(resp)
    
  }
  if (type == "SSL") {
    cdf <- function(y) {
      f <- function(u) 2 * nu * u^(nu - 1) * dnorm(y, mu, sqrt(u^(-1) * 
                                                                 sigma2)) * pnorm(u^(1/2) * lambda * (y - mu)/sqrt(sigma2))
      cdf <- integrate(Vectorize(f), 0, 1)$value
    }
    densidade <- as.numeric(cdf(x))
    resp <- as.numeric(integrate(Vectorize(cdf), -Inf, x)$value)
    return(list(pdf = densidade, cdf = resp))
  }
}

# cdfSNI(c(1,2,3), 2, 3, 2.5, nu = 4, type = "ST")
# cdfSNI(c(1,2,3), 2, 3, 2.5, nu = 4.56, type = "ST")
# cdfSNI(c(1,2,3), 2, 3, 2.5, nu = 5, type = "ST")


################################################################################
####   Acumuladas das NI
################################################################################
AcumPearsonVII <- function(y,mu,sigma2,nu,delta)
{
  Acum <- z <- vector(mode = "numeric", length = length(y))
  sigma2a <- sigma2*(delta/nu)
  for (i in 1:length(y))
  { 
    z[i] <- (y[i]-mu)/sqrt(sigma2a)
    Acum[i] <- pt(z[i],df=nu)
  }
  return(Acum)
}

P <- function(y,mu,sigma2,nu,delta)
{
  A <- z <- vector(mode = "numeric", length = length(y))
  sigma2a <- sigma2*(delta/nu)
  n <- length(y)
  i <- 0
  while (i<n)
  { 
    i <- i +1      
    z[i] <- (y[i]-mu)/sqrt(sigma2a)
    A[i] <- pt(z[i],df=nu)
  }
  return(A)
}

AcumSlash <- function(y,mu,sigma2,nu)
{
  Acum <- z <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y))
  { 
    z[i] <- (y[i]-mu)/sqrt(sigma2)
    f1 <- function(u) nu*u^(nu-1)*pnorm(z[i]*sqrt(u)) 
    Acum[i]<- integrate(f1,0,1)$value	 	
  }
  return(Acum)
}

AcumNormalC <- function(y,mu,sigma2,nu)
{
  Acum <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y))
  { 
    eta  <- nu[1]    
    gama <- nu[2]
    Acum[i] <- eta*pnorm(y[i],mu,sqrt(sigma2/gama)) + (1-eta)*pnorm(y[i],mu,sqrt(sigma2))
  }
  return(Acum)
}

################################################################################
####   Funcoes Auxiliares
################################################################################

GamaInc <- function(a,x)
{
  res <- vector(mode = "numeric", length = length(x))
  f <-function(t) exp(-t)*t^(a-1)
  for  (i in 1:length(x))
  {
    res[i] <- integrate(f,0,x[i])$value
  }
  return(res)
}

Mont_SSL <- function(m,r,a,lambda,nu, case="1")
{
  deltinha <- lambda/sqrt(1+lambda^2)
  Sigma1 <- matrix(c(1,-deltinha,-deltinha,1), 2, 2)  
       y <- d <- vector(mode = "numeric", length = m)
  if(case=="1")
  {
   #y <- rtgamma(1,shape=(r+nu),scale=(0.5*a^2),a=0,b=1)
    y <- rtrunc(m, "gamma", a=0, b=1, shape=(r+nu), scale=1/(0.5*a^2))
    d <- pnorm(sqrt(y)*lambda*a)
  }
  for (i in 1:m)
  {
    if(case=="2")
    {
      y[i] <- rbeta(1, shape1=(r+nu), shape2=1)
      if(a==Inf)
      { 
      a=500
      }
     valor <- c(a,0)*sqrt(y[i])
      d[i] <- pmnorm(valor, mu=c(0,0), varcov=Sigma1)
    }
  }
  emp <- (1/m)*sum(d)
  return(emp=emp)
}

##################################################
##################################################
########  MOMENTOS DAS SNI - TRUNCADAS  ##########
##################################################
##################################################
#install.packages("mnormt")

MomenSNI <- function(mu,sigma2,lambda,nu,delta,Lim1,Lim2,type)    ### Por em quanto valido so para STT
  {
    Lim11 <- (Lim1-mu)/sqrt(sigma2)        
    Lim21 <- (Lim2-mu)/sqrt(sigma2)
        b <- sqrt(2/pi)
        n <- length(Lim11) 
  lambda1 <- sqrt(1+lambda^2)
 Slambda1 <- b*lambda/(lambda1)^2
    if(type=="SN")
    {
     type1 <- "N"
      EU <-  1
    #  FNIb   <-  cdfSNI(Lim21,mu=0,sigma2=1,lambda,nu=NULL,type="SN")
    #  FNIa   <-  cdfSNI(Lim11,mu=0,sigma2=1,lambda,nu=NULL,type="SN")
    }
    if(type=="ST")
    {
   type1 <- "T"
    #  FNIb   <-  cdfSNI(Lim21,mu=0,sigma2=1,lambda,nu,type="ST")
    #  FNIa   <-  cdfSNI(Lim11,mu=0,sigma2=1,lambda,nu,type="ST")
    }
    if(type=="SSL")
    {
       type1 <- "SL"      
    # FNIb   <- cdfSNI(Lim21,mu=0,sigma2=1,lambda,nu,type="SSL")
    #  FNIa   <- cdfSNI(Lim11,mu=0,sigma2=1,lambda,nu,type="SSL")
    }
    if(type=="SCN")
    { 
       type1 <- "CN"      
    #      EU <- (nu[1]*nu[2]) + (1-nu[1])
    #  FNIb   <- cdfSNI(Lim21,mu=0,sigma2=1,lambda,nu,type="SCN")
    #  FNIa   <- cdfSNI(Lim11,mu=0,sigma2=1,lambda,nu,type="SCN")
    }
    if(Lim11==-Inf)
     {
      FNIa <- 0
        S4 <- (Lim21)*E_phiSNI(-0.5,Lim21,nu,delta=NULL,lambda,type=type)
   #     S7 <- (Lim21^2)*E_phiSNI(-0.5,Lim21,nu,delta=NULL,lambda,type=type)
   #     S9 <- ((b*lambda/lambda1^2))*(Lim21*E_phi(-1,Lim21*lambda1,nu,delta,type=type1))  
   #    S15 <- (Lim21)*E_phiSNI(-1.5,Lim21,nu,delta=NULL,lambda,type=type)
   #    S16 <- (Lim21)^3*E_phiSNI(-0.5,Lim21,nu,delta=NULL,lambda,type=type)
   #    S18 <- ((b*lambda/lambda1^2))*((Lim21)^2*E_phi(-1,Lim21*lambda1,nu,delta,type=type1))  
     }  
    else
     {
      FNIa <- cdfSNI(Lim11,mu=0,sigma2=1,lambda,nu=nu,type=type)
     }
    if(Lim21==Inf)
     {
       FNIb <- 1
         S4 <- - (Lim11)*E_phiSNI(-0.5,Lim11,nu,delta=NULL,lambda,type=type)
 #        S7 <- - (Lim11^2)*E_phiSNI(-0.5,Lim11,nu,delta=NULL,lambda,type=type)
  #       S9 <- - -((b*lambda/lambda1^2))*(Lim11*E_phi(-1,Lim11*lambda1,nu,delta,type=type1))         
  #      S15 <- - (Lim11)*E_phiSNI(-1.5,Lim11,nu,delta=NULL,lambda,type=type)
  #      S16 <- - (Lim11)^3*E_phiSNI(-0.5,Lim11,nu,delta=NULL,lambda,type=type)
  #      S18 <- - -((b*lambda/lambda1^2))*((Lim11)^2*E_phi(-1,Lim11*lambda1,nu,delta,type=type1))  
     }  
    else
     {
        FNIb <- cdfSNI(Lim21,mu=0,sigma2=1,lambda,nu,type=type)
     }  
     prob <- FNIb-FNIa
          if(length(which(prob == 0)) > 0) prob[which(prob == 0)] <- .Machine$double.xmin
           
       K <- 1/prob
      S0 <- (b*lambda/lambda1)*(E_Phi(-0.5,Lim21*lambda1, nu,delta,type=type1)- E_Phi(-0.5,Lim11*lambda1, nu,delta,type=type1))
      S1 <- E_phiSNI(-0.5,Lim21,nu,delta=NULL,lambda,type=type)- E_phiSNI(-0.5,Lim11,nu,delta=NULL,lambda,type=type)
      S2 <- E_PhiSNI(-1,Lim21,nu,delta,lambda,type=type)- E_PhiSNI(-1,Lim11,nu,delta,lambda,type=type)
      S3 <- Slambda1*(E_phi(-1,Lim21*lambda1,nu,delta,type=type1)-E_phi(-1,Lim11*lambda1,nu,delta,type=type1))
###########
 #     S5 <- E_phiSNI(-1.5,Lim21,nu,delta=NULL,lambda,type=type)- E_phiSNI(-1.5,Lim11,nu,delta=NULL,lambda,type=type)
 #     S6 <- (b*lambda/lambda1)*(E_Phi(-1.5,Lim21*lambda1, nu,delta,type=type1)- E_Phi(-1.5,Lim11*lambda1, nu,delta,type=type1))  
#      S8 <- ((b*lambda/lambda1^3))*(E_Phi(-1.5,Lim21*lambda1, nu,delta,type=type1)- E_Phi(-1.5,Lim11*lambda1, nu,delta,type=type1))      
###########
  #   S13<- E_PhiSNI(-2,Lim21,nu,delta,lambda,type=type)- E_PhiSNI(-2,Lim11,nu,delta,lambda,type=type)
  #   S14<- (b*lambda/lambda1^2)*(E_phi(-2,Lim21*lambda1,nu,delta,type=type1)-E_phi(-2,Lim11*lambda1,nu,delta,type=type1))
  #   S17<-  ((b*lambda/lambda1^4))*(E_phi(-2,Lim21*lambda1,nu,delta,type=type1)-E_phi(-2,Lim11*lambda1,nu,delta,type=type1))
###########
if(setequal(Lim11,-Inf)== FALSE & setequal(Lim21,Inf)== FALSE)
{
  S4 <- (Lim21)*E_phiSNI(-0.5,Lim21,nu,delta=NULL,lambda,type=type)- (Lim11)*E_phiSNI(-0.5,Lim11,nu,delta=NULL,lambda,type=type)  
 # S7 <- (Lim21^2)*E_phiSNI(-0.5,Lim21,nu,delta=NULL,lambda,type=type)- (Lim11^2)*E_phiSNI(-0.5,Lim11,nu,delta=NULL,lambda,type=type)
 # S9 <- ((b*lambda/lambda1^2))*(Lim21*E_phi(-1,Lim21*lambda1,nu,delta,type=type1)-Lim11*E_phi(-1,Lim11*lambda1,nu,delta,type=type1))  
# S15 <- (Lim21)*E_phiSNI(-1.5,Lim21,nu,delta=NULL,lambda,type=type)-(Lim11)*E_phiSNI(-1.5,Lim11,nu,delta=NULL,lambda,type=type)
 #S16 <- (Lim21)^3*E_phiSNI(-0.5,Lim21,nu,delta=NULL,lambda,type=type)-(Lim11)^3*E_phiSNI(-0.5,Lim11,nu,delta=NULL,lambda,type=type)
# S18 <- ((b*lambda/lambda1^2))*((Lim21)^2*E_phi(-1,Lim21*lambda1,nu,delta,type=type1)-(Lim11)^2*E_phi(-1,Lim11*lambda1,nu,delta,type=type1))  
}  
# Saux1<- 3*(S13-S14-S15)
# Saux2<- 2*S17+S18
    #### 
    EUX1 <- K*(S0-S1) 
    EUX2 <- K*(S2-S3-S4)
   # EUX3 <- K*(2*(-S5+S6)-S7+S8-S9)
   # EUX4 <- K*(Saux1-S16-Saux2)
   sigma <- sqrt(sigma2)
    EUY1 <- mu + sigma*EUX1
    EUY2 <- mu^2 + 2*mu*sigma*EUX1 + (sigma^2)*EUX2 
    #EUY3 <- mu^3 + 3*(mu^2)*sigma*EUX1 + 3*mu*(sigma^2)*EUX2 + (sigma^3)*EUX3
    #EUY4 <- (mu^4) +  4*(mu^3)*sigma*EUX1 + 6*(mu^2)*(sigma^2)*EUX2 + 4*mu*(sigma^3)*EUX3 + (sigma^4)*EUX4
   #Skewn <- (EUY3-3*EUY1*EUY2+2*EUY1^3)/(EUY2-EUY1^2)^(1.5)
   # Kurt <- (EUY4-4*EUY1*EUY3+6*EUY2*EUY1^2-3*EUY1^4)/(EUY2-EUY1^2)^2
    return(list(EUY1=EUY1,EUY2=EUY2))
  }

## ----------------------------------- ##
## GENERATION OF CENSORED SMSN SAMPLES ##
## ----------------------------------- ##

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
  
  #mu <-  eta*Delta
  mu<-0 
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