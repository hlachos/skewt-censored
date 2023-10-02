#library(mvtnorm)

#location-scale student-T distribution pdf
dent<-function(x,mu,sigma2,nu){
  z<-(x-mu)/sqrt(sigma2)
  return(1/sqrt(sigma2)*dt(x = z,df = nu))
}

#########################################################################

#location-scale student-T distribution pcf
pent<-function(x,mu,sigma2,nu){
  return(pt((x-mu)/sqrt(sigma2),nu))
}

##########################################################################

acumt2 = function(a = NULL,b,mu,Sigma,nu){
  Sigma21 = c(Sigma[-1,-1] - Sigma[-1,1]%*%solve(Sigma[1,1])%*%Sigma[1,-1])
  expab = function(x,y){
    delta1 = c((x - mu[1])^2/Sigma[1,1])
    mu21 = mu[-1] + Sigma[-1,1]%*%solve(Sigma[1,1])%*%(x - mu[1])
    return(pent(y,mu21,Sigma21*(nu+delta1)/(nu+1),nu+1))
  }
  func = function(x,val){
    val*dent(x,mu[1],Sigma[1,1],nu)
  }
  if(all(is.infinite(a)) | is.null(a)){
    return(integrate(f = function(x) func(x,expab(x,b[2])),lower = -Inf,upper = b[1])$value)
    #return(integrate(f = function(x) func(x,expab(x,b[2])),lower = mu[1]-37*sqrt(Sigma[1,1]),upper = b[1])$value)
  }else{
    return(integrate(f = function(x) func(x,expab(x,b[2])),lower = a[1],upper = b[1])$value - 
             integrate(f = function(x) func(x,expab(x,a[2])),lower = a[1],upper = b[1])$value)
  }
}

##########################################################################
#TEST
##########################################################################

# a = c(-1,1)
# b = c(2,3)
# mu = as.matrix(c(1,2))
# Sigma = matrix(c(2,-0.5,-0.5,1),2,2)
 
# nu = 3.2
# acumt2(a,b,mu,Sigma,nu)
# tlrmvnmvt::pmvt(lower = c(a-mu),upper = c(b-mu),df = nu,sigma = Sigma)[1]

# 
# nu = 3.27
# acumt2(a,b,mu,Sigma,nu)
# pmvt(lower = c(a-mu),upper = c(b-mu),df = nu,sigma = Sigma)[1]
# 
# #CDF
# a = c(-Inf,-Inf)
# acumt2(a,b,mu,Sigma,nu)
# #or a = NULL
# acumt2(b = b,mu = mu,Sigma = Sigma,nu = nu)
# tlrmvnmvt::pmvt(lower = c(a-mu),upper = c(b-mu),df = nu,sigma = Sigma)[1]
