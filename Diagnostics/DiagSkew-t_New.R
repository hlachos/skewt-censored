### ========================================================================
###    Influence diagnostics for skew-t censored linear regression models   
### Authors: Marcos S. Oliveira, Daniela C. R. Oliveira and Victor H. Lachos
### ========================================================================

### --------------------------------
### HESSIAN Q                       
### --------------------------------

HessianaQ <- function(cc, x, y, theta, E00, E01, E02, E10, E20, E11, family = "ST"){
  
  p <- ncol(x)
  n <- nrow(x)
  
  if(family=="ST"){
    beta    <- theta[1:p]
    sigma2  <- theta[p+1]
    lambda  <- theta[p+2]
    nu      <- theta[p+3] 
    mu      <- x%*%beta
    sigma   <- sqrt(sigma2)
      
    suma1   <- matrix(0, p, p)  # betabeta
    suma2   <- matrix(0, p, 1)  # betasigma2
    suma3   <- matrix(0, p, 1)  # betalambda
    suma4   <- 0                # sigma2sigma2
    suma5   <- 0                # sigma2lambda
    suma6   <- 0                # lambdalambda
    
    sbeta   <- matrix(0, n, p)
    ssigma2 <- matrix(0, n, 1)
    slambda <- matrix(0, n, 1)
      
    for(i in 1:n){
      suma1          <- suma1 - ((1 + lambda^2)/sigma2)*E00[i]*(x[i,])%*%t(x[i,])
      aux2           <- x[i,]*E01[i] - E00[i]*(x[i,])%*%t(x[i,])%*%beta
      aux3           <- x[i,]*E10[i]
      aux4           <- E02[i] - 2*E01[i]*mu[i] + E00[i]*mu[i]^2
      aux5           <- E11[i] - E10[i]*mu[i]
      suma2          <- suma2 - ((1 + lambda^2)/sigma2)*(1/sigma2*aux2 - 0.5*lambda/(sigma*sqrt(1 + lambda^2))*aux3)
      suma3          <- suma3 + 2*lambda/sigma2*aux2 - (2*lambda^2 + 1)/(sigma*sqrt(1 + lambda^2))*aux3
      suma4          <- suma4 + 1/(2*sigma2^2) - (1 + lambda^2)/sigma2^3*aux4 + 3/4*lambda*sqrt(lambda^2 + 1)/sigma^5*aux5
      suma5          <- suma5 + lambda^2/sigma2^2*aux4 - (2*lambda^2 + 1)/(2*sigma^3*sqrt(1 + lambda^2))*aux5
      suma6          <- suma6 - (lambda^2 - 1)/(lambda^2 + 1)^2 - 1/sigma2*aux4 + 1/sigma*(lambda*(2*lambda^2 -3)/(lambda^2 + 1)^1.5)*aux5 - E20[i]
        
      sbeta[i,]      <- ((1 + lambda^2)/sigma2)*(E01[i]*t(x[i,]) - E00[i]*t(x[i,])*mu[i] - sigma*lambda/sqrt(1 + lambda^2)*E10[i]*t(x[i,]))
      ssigma2[i]     <- -1/(2*sigma2) + ((1 + lambda^2)/(2*sigma2^2))*(E02[i] - 2*E01[i]*mu[i] + (t(mu[i])%*%mu[i])*E00[i]) - ((lambda*sqrt(1 + lambda^2))/(2*sigma^3))*(E11[i] - E10[i]*mu[i])
      slambda[i]     <- lambda/(1 + lambda^2) - (lambda/sigma2)*(E02[i] - 2*E01[i]*mu[i] + E00[i]*(t(mu[i])%*%mu[i])) +  ((1 + 2*lambda^2)/(sigma*sqrt(1 + lambda^2)))*(E11[i] - E10[i]*mu[i]) - lambda*E20[i]
    }
      
    derbetabeta      <- suma1
    derbetasigma2    <- suma2
    derbetalambda    <- suma3
    dersigma2sigma2  <- suma4
    dersigma2lambda  <- suma5
    derlambdalambda  <- suma6
      
    MatrizQ          <- matrix(0, nrow = (p+2), ncol = (p+2))
    MatrizQ[1:p,1:p] <- derbetabeta
    MatrizQ[p+1,1:p] <- t(derbetasigma2)
    MatrizQ[1:p,p+1] <- derbetasigma2
    MatrizQ[1:p,p+2] <- derbetalambda
    MatrizQ[p+2,1:p] <- t(derbetalambda)
    MatrizQ[p+1,p+2] <- dersigma2lambda
    MatrizQ[p+1,p+1] <- dersigma2sigma2
    MatrizQ[p+2,p+1] <- dersigma2lambda
    MatrizQ[p+2,p+2] <- derlambdalambda
      
    dQtheta          <- matrix(0, p+2, n)
    thetai           <- matrix(0, p+2, n)
      
    for(i in 1:n){
      dQtheta[1:p,i] <- apply(sbeta[-i,], 2, sum)
      dQtheta[p+1,i] <- sum(ssigma2[-i])
      dQtheta[p+2,i] <- sum(slambda[-i])
      thetai[,i]     <- theta[1:(p+2)] + solve(-MatrizQ)%*%dQtheta[,i]
    }
  }
    
  if(family=="SN"){
    beta    <- theta[1:p]
    sigma2  <- theta[p+1]
    lambda  <- theta[p+2] 
    mu      <- x%*%beta
    sigma   <- sqrt(sigma2)
    
    suma1   <- matrix(0, p, p)  # betabeta
    suma2   <- matrix(0, p, 1)  # betasigma2
    suma3   <- matrix(0, p, 1)  # betalambda
    suma4   <- 0                # sigma2sigma2
    suma5   <- 0                # sigma2lambda
    suma6   <- 0                # lambdalambda
    
    sbeta   <- matrix(0, n, p)
    ssigma2 <- matrix(0, n, 1)
    slambda <- matrix(0, n, 1)

    for(i in 1:n){
      suma1          <- suma1 - ((1 + lambda^2)/sigma2)*E00[i]*(x[i,])%*%t(x[i,])
      aux2           <- x[i,]*E01[i] - E00[i]*(x[i,])%*%t(x[i,])%*%beta
      aux3           <- x[i,]*E10[i]
      aux4           <- E02[i] - 2*E01[i]*mu[i] + E00[i]*mu[i]^2
      aux5           <- E11[i] - E10[i]*mu[i]
      suma2          <- suma2 - ((1 + lambda^2)/sigma2)*(1/sigma2*aux2 - 0.5*lambda/(sigma*sqrt(1 + lambda^2))*aux3)
      suma3          <- suma3 + 2*lambda/sigma2*aux2 - (2*lambda^2 + 1)/(sigma*sqrt(1 + lambda^2))*aux3
      suma4          <- suma4 + 1/(2*sigma2^2) - (1 + lambda^2)/sigma2^3*aux4 + 3/4*lambda*sqrt(lambda^2 + 1)/sigma^5*aux5
      suma5          <- suma5 + lambda^2/sigma2^2*aux4 - (2*lambda^2 + 1)/(2*sigma^3*sqrt(1 + lambda^2))*aux5
      suma6          <- suma6 - (lambda^2 - 1)/(lambda^2 + 1)^2 - 1/sigma2*aux4 + 1/sigma*(lambda*(2*lambda^2 -3)/(lambda^2 + 1)^1.5)*aux5 - E20[i]
      
      sbeta[i,]      <- ((1 + lambda^2)/sigma2)*(E01[i]*t(x[i,]) - E00[i]*t(x[i,])*mu[i] - sigma*lambda/sqrt(1 + lambda^2)*E10[i]*t(x[i,]))
      ssigma2[i]     <- -1/(2*sigma2) + ((1 + lambda^2)/(2*sigma2^2))*(E02[i] - 2*E01[i]*mu[i] + (t(mu[i])%*%mu[i])*E00[i]) - ((lambda*sqrt(1 + lambda^2))/(2*sigma^3))*(E11[i] - E10[i]*mu[i])
      slambda[i]     <- lambda/(1 + lambda^2) - (lambda/sigma2)*(E02[i] - 2*E01[i]*mu[i] + E00[i]*(t(mu[i])%*%mu[i])) +  ((1 + 2*lambda^2)/(sigma*sqrt(1 + lambda^2)))*(E11[i] - E10[i]*mu[i]) - lambda*E20[i]
    }
    
    derbetabeta      <- suma1
    derbetasigma2    <- suma2
    derbetalambda    <- suma3
    dersigma2sigma2  <- suma4
    dersigma2lambda  <- suma5
    derlambdalambda  <- suma6
    
    MatrizQ          <- matrix(0, nrow = (p+2), ncol = (p+2))
    MatrizQ[1:p,1:p] <- derbetabeta
    MatrizQ[p+1,1:p] <- t(derbetasigma2)
    MatrizQ[1:p,p+1] <- derbetasigma2
    MatrizQ[1:p,p+2] <- derbetalambda
    MatrizQ[p+2,1:p] <- t(derbetalambda)
    MatrizQ[p+1,p+2] <- dersigma2lambda
    MatrizQ[p+1,p+1] <- dersigma2sigma2
    MatrizQ[p+2,p+1] <- dersigma2lambda
    MatrizQ[p+2,p+2] <- derlambdalambda
    
    dQtheta          <- matrix(0, p+2, n)
    thetai           <- matrix(0, p+2, n)
    
    for(i in 1:n){
      dQtheta[1:p,i] <- apply(sbeta[-i,], 2, sum)
      dQtheta[p+1,i] <- sum(ssigma2[-i])
      dQtheta[p+2,i] <- sum(slambda[-i])
      thetai[,i]     <- theta[1:(p+2)] + solve(-MatrizQ)%*%dQtheta[,i] 
    }
  }
  
  if (family=="N"){
    beta    <- theta[1:p]
    sigma2  <- theta[p+1]
    lambda  <- 0                # theta[p+2] 
    mu      <- x%*%beta
    sigma   <- sqrt(sigma2)
    
    suma1   <- matrix(0, p, p)  # betabeta
    suma2   <- matrix(0, p, 1)  # betasigma2
    #suma3   <- matrix(0, p, 1)  # betalambda
    suma4   <- 0                # sigma2sigma2
    #suma5   <- 0                # sigma2lambda
    #suma6   <- 0                # lambdalambda
    
    sbeta   <- matrix(0, n, p)
    ssigma2 <- matrix(0, n, 1)
    #slambda <- matrix(0, n, 1)
    
    for(i in 1:n){
      suma1          <- suma1 - ((1 + lambda^2)/sigma2)*E00[i]*(x[i,])%*%t(x[i,])
      aux2           <- x[i,]*E01[i] - E00[i]*(x[i,])%*%t(x[i,])%*%beta
      aux3           <- x[i,]*E10[i]
      aux4           <- E02[i] - 2*E01[i]*mu[i] + E00[i]*mu[i]^2
      aux5           <- E11[i] - E10[i]*mu[i]
      suma2          <- suma2 - ((1 + lambda^2)/sigma2)*(1/sigma2*aux2 - 0.5*lambda/(sigma*sqrt(1 + lambda^2))*aux3)
      #suma3          <- suma3 + 2*lambda/sigma2*aux2 - (2*lambda^2 + 1)/(sigma*sqrt(1 + lambda^2))*aux3
      suma4          <- suma4 + 1/(2*sigma2^2) - (1 + lambda^2)/sigma2^3*aux4 + 3/4*lambda*sqrt(lambda^2 + 1)/sigma^5*aux5
      #suma5          <- suma5 + lambda^2/sigma2^2*aux4 - (2*lambda^2 + 1)/(2*sigma^3*sqrt(1 + lambda^2))*aux5
      #suma6          <- suma6 - (lambda^2 - 1)/(lambda^2 + 1)^2 - 1/sigma2*aux4 + 1/sigma*(lambda*(2*lambda^2 -3)/(lambda^2 + 1)^1.5)*aux5 - E20[i]
      
      sbeta[i,]      <- ((1 + lambda^2)/sigma2)*(E01[i]*t(x[i,]) - E00[i]*t(x[i,])*mu[i] - sigma*lambda/sqrt(1 + lambda^2)*E10[i]*t(x[i,]))
      ssigma2[i]     <- -1/(2*sigma2) + ((1 + lambda^2)/(2*sigma2^2))*(E02[i] - 2*E01[i]*mu[i] + (t(mu[i])%*%mu[i])*E00[i]) - ((lambda*sqrt(1 + lambda^2))/(2*sigma^3))*(E11[i] - E10[i]*mu[i])
      #slambda[i]     <- lambda/(1 + lambda^2) - (lambda/sigma2)*(E02[i] - 2*E01[i]*mu[i] + E00[i]*(t(mu[i])%*%mu[i])) +  ((1 + 2*lambda^2)/(sigma*sqrt(1 + lambda^2)))*(E11[i] - E10[i]*mu[i]) - lambda*E20[i]
    }
    
    derbetabeta      <- suma1
    derbetasigma2    <- suma2
    #derbetalambda    <- suma3
    dersigma2sigma2  <- suma4
    #dersigma2lambda  <- suma5
    #derlambdalambda  <- suma6
    
    MatrizQ          <- matrix(0, nrow = (p+1), ncol = (p+1))
    MatrizQ[1:p,1:p] <- derbetabeta
    MatrizQ[p+1,1:p] <- t(derbetasigma2)
    MatrizQ[1:p,p+1] <- derbetasigma2
    #MatrizQ[1:p,p+2] <- derbetalambda
    #MatrizQ[p+2,1:p] <- t(derbetalambda)
    #MatrizQ[p+1,p+2] <- dersigma2lambda
    MatrizQ[p+1,p+1] <- dersigma2sigma2
    #MatrizQ[p+2,p+1] <- dersigma2lambda
    #MatrizQ[p+2,p+2] <- derlambdalambda
    
    dQtheta          <- matrix(0, p+1, n)
    thetai           <- matrix(0, p+1, n)
    
    for(i in 1:n){
      dQtheta[1:p,i] <- apply(sbeta[-i,], 2, sum)
      dQtheta[p+1,i] <- sum(ssigma2[-i])
      #dQtheta[p+2,i] <- sum(slambda[-i])
      thetai[,i]     <- theta[1:(p+1)] + solve(-MatrizQ)%*%dQtheta[,i] 
    }
  }
  
  if (family=="T"){
    beta    <- theta[1:p]
    sigma2  <- theta[p+1]
    nu      <- theta[p+2] 
    lambda  <- 0
    mu      <- x%*%beta
    sigma   <- sqrt(sigma2)
    
    suma1   <- matrix(0, p, p)  # betabeta
    suma2   <- matrix(0, p, 1)  # betasigma2
    #suma3   <- matrix(0, p, 1)  # betalambda
    suma4   <- 0                # sigma2sigma2
    #suma5   <- 0                # sigma2lambda
    #suma6   <- 0                # lambdalambda
    
    sbeta   <- matrix(0, n, p)
    ssigma2 <- matrix(0, n, 1)
    #slambda <- matrix(0, n, 1)
    
    for(i in 1:n){
      suma1          <- suma1 - ((1 + lambda^2)/sigma2)*E00[i]*(x[i,])%*%t(x[i,])
      aux2           <- x[i,]*E01[i] - E00[i]*(x[i,])%*%t(x[i,])%*%beta
      aux3           <- x[i,]*E10[i]
      aux4           <- E02[i] - 2*E01[i]*mu[i] + E00[i]*mu[i]^2
      aux5           <- E11[i] - E10[i]*mu[i]
      suma2          <- suma2 - ((1 + lambda^2)/sigma2)*(1/sigma2*aux2 - 0.5*lambda/(sigma*sqrt(1 + lambda^2))*aux3)
      #suma3          <- suma3 + 2*lambda/sigma2*aux2 - (2*lambda^2 + 1)/(sigma*sqrt(1 + lambda^2))*aux3
      suma4          <- suma4 + 1/(2*sigma2^2) - (1 + lambda^2)/sigma2^3*aux4 + 3/4*lambda*sqrt(lambda^2 + 1)/sigma^5*aux5
      #suma5          <- suma5 + lambda^2/sigma2^2*aux4 - (2*lambda^2 + 1)/(2*sigma^3*sqrt(1 + lambda^2))*aux5
      #suma6          <- suma6 - (lambda^2 - 1)/(lambda^2 + 1)^2 - 1/sigma2*aux4 + 1/sigma*(lambda*(2*lambda^2 -3)/(lambda^2 + 1)^1.5)*aux5 - E20[i]
      
      sbeta[i,]      <- ((1 + lambda^2)/sigma2)*(E01[i]*t(x[i,]) - E00[i]*t(x[i,])*mu[i] - sigma*lambda/sqrt(1 + lambda^2)*E10[i]*t(x[i,]))
      ssigma2[i]     <- -1/(2*sigma2) + ((1 + lambda^2)/(2*sigma2^2))*(E02[i] - 2*E01[i]*mu[i] + (t(mu[i])%*%mu[i])*E00[i]) - ((lambda*sqrt(1 + lambda^2))/(2*sigma^3))*(E11[i] - E10[i]*mu[i])
      #slambda[i]     <- lambda/(1 + lambda^2) - (lambda/sigma2)*(E02[i] - 2*E01[i]*mu[i] + E00[i]*(t(mu[i])%*%mu[i])) +  ((1 + 2*lambda^2)/(sigma*sqrt(1 + lambda^2)))*(E11[i] - E10[i]*mu[i]) - lambda*E20[i]
    }
    
    derbetabeta      <- suma1
    derbetasigma2    <- suma2
    #derbetalambda    <- suma3
    dersigma2sigma2  <- suma4
    #dersigma2lambda  <- suma5
    #derlambdalambda  <- suma6
    
    MatrizQ          <- matrix(0, nrow = (p+1), ncol = (p+1))
    MatrizQ[1:p,1:p] <- derbetabeta
    MatrizQ[p+1,1:p] <- t(derbetasigma2)
    MatrizQ[1:p,p+1] <- derbetasigma2
    #MatrizQ[1:p,p+2] <- derbetalambda
    #MatrizQ[p+2,1:p] <- t(derbetalambda)
    #MatrizQ[p+1,p+2] <- dersigma2lambda
    MatrizQ[p+1,p+1] <- dersigma2sigma2
    #MatrizQ[p+2,p+1] <- dersigma2lambda
    #MatrizQ[p+2,p+2] <- derlambdalambda
    
    dQtheta          <- matrix(0, p+1, n)
    thetai           <- matrix(0, p+1, n)
    
    for(i in 1:n){
      dQtheta[1:p,i] <- apply(sbeta[-i,], 2, sum)
      dQtheta[p+1,i] <- sum(ssigma2[-i])
      #dQtheta[p+2,i] <- sum(slambda[-i])
      thetai[,i]     <- theta[1:(p+1)] + solve(-MatrizQ)%*%dQtheta[,i]
    }
  }
  
  obj.out <- list(MatrizQ = MatrizQ, dQtheta = dQtheta, thetai=thetai)

  return(obj.out)
}

### --------------------------------
### Q FUNCTION                      
### --------------------------------

Q.function <- function(cc, x, y, theta, E00, E01, E02,  E10,  E20,  E11,  family = "ST"){
  
  p <- ncol(x)
	n <- nrow(x)
	
	if(family=="ST"){
	  beta    <- theta[1:p]
    sigma2  <- theta[p+1]
    lambda  <- theta[p+2]
    nu      <- theta[p+3]
    mu      <- x%*%beta
    sigma   <- sqrt(sigma2)
    delta   <- lambda / (sqrt(1 + lambda^2))
    Delta   <- sqrt(sigma2)*delta
    Gama    <- sigma2 - Delta^2   
    suma    <- 0
    
    for (i in 1:n){
      suma  <- suma - 0.5*log(Gama) - 1/(2*Gama)*(E02[i] - 2*E01[i]*mu[i] + mu[i]^2*E00[i] + Delta^2*E20[i] - (2*Delta)*(E11[i] - E10[i]*mu[i]))
    } 
  }
	
	if(family=="SN"){
	  beta    <- theta[1:p]
	  sigma2  <- theta[p+1]
	  lambda  <- theta[p+2]
	  #nu      <- theta[p+3]
	  mu      <- x%*%beta
	  sigma   <- sqrt(sigma2)
	  delta   <- lambda / (sqrt(1 + lambda^2))
	  Delta   <- sqrt(sigma2)*delta
	  Gama    <- sigma2 - Delta^2   
	  suma    <- 0
	  
	  for (i in 1:n){
	    suma  <- suma - 0.5*log(Gama) - 1/(2*Gama)*(E02[i] - 2*E01[i]*mu[i] + mu[i]^2*E00[i] + Delta^2*E20[i] - (2*Delta)*(E11[i] - E10[i]*mu[i]))
	  }
	}
	
	if(family=="N"){
	  beta    <- theta[1:p]
	  sigma2  <- theta[p+1]
	  lambda  <- 0
	  #nu      <- theta[p+3]
	  mu      <- x%*%beta
	  sigma   <- sqrt(sigma2)
	  delta   <- lambda / (sqrt(1 + lambda^2))
	  Delta   <- sqrt(sigma2)*delta
	  Gama    <- sigma2 - Delta^2   
	  suma    <- 0
	  
	  for (i in 1:n){
	    suma  <- suma - 0.5*log(Gama) - 1/(2*Gama)*(E02[i] - 2*E01[i]*mu[i] + mu[i]^2*E00[i] + Delta^2*E20[i] - (2*Delta)*(E11[i] - E10[i]*mu[i]))
	  } 
	}
	
	if(family=="T"){
	  beta    <- theta[1:p]
	  sigma2  <- theta[p+1]
	  nu      <- theta[p+2]
	  lambda  <- 0
	  mu      <- x%*%beta
	  sigma   <- sqrt(sigma2)
	  delta   <- lambda / (sqrt(1 + lambda^2))
	  Delta   <- sqrt(sigma2)*delta
	  Gama    <- sigma2 - Delta^2   
	  suma    <- 0
	  
	  for (i in 1:n){
	    suma  <- suma - 0.5*log(Gama) - 1/(2*Gama)*(E02[i] - 2*E01[i]*mu[i] + mu[i]^2*E00[i] + Delta^2*E20[i] - (2*Delta)*(E11[i] - E10[i]*mu[i]))
	  } 
	}
	
	return(suma)
}

### --------------------------------
### CASE-DELETION                   
### --------------------------------

CaseDele <- function(cc, x, y, theta, E00,  E01,  E02,  E10,  E20,  E11,  family = "ST"){
  
  HQ <- HessianaQ(cc, x, y, theta, E00, E01,  E02,  E10,  E20,  E11, family = family)
  
  p  <- ncol(x)
	n  <- nrow(x)
	
	if(family=="ST"){
	  beta     <- theta[1:p]
    sigma2   <- theta[p+1]
    lambda   <- theta[p+2]
    nu       <- theta[p+3]
    mu       <- x%*%beta
    sigma    <- sqrt(sigma2)
    
    GD       <- matrix(0, n, 1)
    QD       <- matrix(0, n, 1)
    
    QComp    <- Q.function(cc, x, y, theta, E00, E01,  E02,  E10,  E20,  E11,  family = family)
    
    for (i in 1:n){
      GD[i]  <- abs(t(HQ$dQtheta[,i])%*%solve(-HQ$MatrizQ)%*%HQ$dQtheta[,i])
      Qexc   <- Q.function(cc, x, y, c(HQ$thetai[,i], nu), E00, E01, E02, E10, E20, E11,  family = family)
      QD[i]  <- abs(2*(QComp-Qexc))
    }
  }
	
	if(family=="SN"){
	  beta     <- theta[1:p]
	  sigma2   <- theta[p+1]
	  lambda   <- theta[p+2]
	  #nu       <- theta[p+3]
	  mu       <- x%*%beta
	  sigma    <- sqrt(sigma2)
	  
	  GD       <- matrix(0, n, 1)
	  QD       <- matrix(0, n, 1)
	  
	  QComp    <- Q.function(cc, x, y, theta, E00, E01,  E02,  E10,  E20,  E11,  family = family)
	  
	  for (i in 1:n){
	    GD[i]  <- abs(t(HQ$dQtheta[,i])%*%solve(-HQ$MatrizQ)%*%HQ$dQtheta[,i])
	    Qexc   <- Q.function(cc, x, y, HQ$thetai[,i], E00, E01, E02, E10, E20, E11,  family = family)
	    QD[i]  <- abs(2*(QComp-Qexc))
	  }
	}
	  
	if(family=="N"){
	  beta     <- theta[1:p]
	  sigma2   <- theta[p+1]
	  #lambda   <- theta[p+2]
	  #nu       <- theta[p+3]
	  mu       <- x%*%beta
	  sigma    <- sqrt(sigma2)
	  
	  GD       <- matrix(0, n, 1)
	  QD       <- matrix(0, n, 1)
	  
	  QComp    <- Q.function(cc, x, y, theta, E00, E01,  E02,  E10,  E20,  E11,  family = family)
	  
	  for (i in 1:n){
	    GD[i]  <- abs(t(HQ$dQtheta[,i])%*%solve(-HQ$MatrizQ)%*%HQ$dQtheta[,i])
	    Qexc   <- Q.function(cc, x, y, HQ$thetai[,i], E00, E01, E02, E10, E20, E11,  family = family)
	    QD[i]  <- abs(2*(QComp-Qexc))
	  }
	}
	 
	if(family=="T"){ 
	  beta     <- theta[1:p]
	  sigma2   <- theta[p+1]
	  nu       <- theta[p+2]
	  lambda   <- 0
	  mu       <- x%*%beta
	  sigma    <- sqrt(sigma2)
	  
	  GD       <- matrix(0, n, 1)
	  QD       <- matrix(0, n, 1)
	  
	  QComp    <- Q.function(cc, x, y, theta, E00, E01,  E02,  E10,  E20,  E11,  family = family)
	  
	  for (i in 1:n){
	    GD[i]  <- abs(t(HQ$dQtheta[,i])%*%solve(-HQ$MatrizQ)%*%HQ$dQtheta[,i])
	    Qexc   <- Q.function(cc, x, y, c(HQ$thetai[,i], nu), E00, E01, E02, E10, E20, E11,  family = family)
	    QD[i]  <- abs(2*(QComp-Qexc))
	  }
	}
	
	obj.out <- list(GD = GD, QD = QD)

  return(obj.out)
}
 
### --------------------------------
### INFLUENCE DIAGNOSTICS           
### --------------------------------

If <- function(cc, x, y, theta, E00,  E01,  E02,  E10,  E20,  E11,  Perturb, family = "ST"){
  
  HQ <- HessianaQ(cc, x, y, theta, E00, E01,  E02,  E10,  E20,  E11, family = family)
  
  p  <- ncol(x)
  n  <- nrow(x)

  if(family=="ST"){
    beta     <- theta[1:p]
    sigma2   <- theta[p+1]
    lambda   <- theta[p+2]
    nu       <- theta[p+3]
    mu       <- x%*%beta
    sigma    <- sqrt(sigma2)
    
    A        <- E00*mu
    B        <- E02 - 2*E01*mu + E00*mu^2
    D        <- E11 - E10*mu
    delta1   <- c()
    delta2   <- c()
    delta3   <- c()
    mdelta   <- matrix(0, (p+2), n)
    
    if(Perturb==1){
      delta1 <- ((1 + lambda^2)/sigma2)*(t(x)%*%diag(E01[1:n]) - t(x)%*%diag(A[1:n]) - sigma*lambda/sqrt(1 + lambda^2)*t(x)%*%diag(E10[1:n])) # der_betaw
      delta2 <- (-1/(2*sigma2))*(rep(1,n) - ((1 + lambda^2)/sigma2)*t(B) + ((lambda*sqrt(1+lambda^2))/sigma)*t(D))                            # der_sigma2w
      delta3 <- (lambda/(1 + lambda^2))*rep(1,n) - (lambda/sigma2)*t(B) + ((1 + 2*lambda^2)/(sigma*sqrt(1+lambda^2)))*t(D) - lambda*t(E20)    # der_lambdaw
        
      mdelta[1:p,] <- delta1
      mdelta[p+1,] <- delta2
      mdelta[p+2,] <- delta3
    }
    if(Perturb==2){
      delta1 <- ((1 + lambda^2)/sigma2)*(t(x)%*%diag(E01[1:n]) - t(x)%*%diag(A[1:n]) - sigma*lambda/(2*sqrt(1 + lambda^2))*t(x)%*%diag(E10[1:n])) # der_betaw
      delta2 <- (1/(2*sigma2))*( ((1 + lambda^2)/sigma2)*t(B) - ((lambda*sqrt(1+lambda^2))/(2*sigma))*t(D))                                       # der_sigma2w
      delta3 <-  - (lambda/sigma2)*t(B) + ((1 + 2*lambda^2)/(2*sigma*sqrt(1+lambda^2)))*t(D)                                                      # der_lambdaw
      
      mdelta[1:p,] <- delta1
      mdelta[p+1,] <- delta2
      mdelta[p+2,] <- delta3
    }
  }  
  
  if(family=="SN"){
    beta     <- theta[1:p]
    sigma2   <- theta[p+1]
    lambda   <- theta[p+2]
    ##nu       <- theta[p+3]
    mu       <- x%*%beta
    sigma    <- sqrt(sigma2)
    
    A        <- E00*mu
    B        <- E02 - 2*E01*mu + E00*mu^2
    D        <- E11 - E10*mu
    delta1   <- c()
    delta2   <- c()
    delta3   <- c()
    mdelta   <- matrix(0, (p+2), n)
    
    if(Perturb==1){
      delta1 <- ((1 + lambda^2)/sigma2)*(t(x)%*%diag(E01[1:n]) - t(x)%*%diag(A[1:n]) - sigma*lambda/sqrt(1 + lambda^2)*t(x)%*%diag(E10[1:n])) # der_betaw
      delta2 <- (-1/(2*sigma2))*(rep(1,n) - ((1 + lambda^2)/sigma2)*t(B) + ((lambda*sqrt(1+lambda^2))/sigma)*t(D))                            # der_sigma2w
      delta3 <- (lambda/(1 + lambda^2))*rep(1,n) - (lambda/sigma2)*t(B) + ((1 + 2*lambda^2)/(sigma*sqrt(1+lambda^2)))*t(D) - lambda*t(E20)    # der_lambdaw
      
      mdelta[1:p,] <- delta1
      mdelta[p+1,] <- delta2
      mdelta[p+2,] <- delta3
    }
    if(Perturb==2){
      delta1 <- ((1 + lambda^2)/sigma2)*(t(x)%*%diag(E01[1:n]) - t(x)%*%diag(A[1:n]) - sigma*lambda/(2*sqrt(1 + lambda^2))*t(x)%*%diag(E10[1:n])) # der_betaw
      delta2 <- (1/(2*sigma2))*( ((1 + lambda^2)/sigma2)*t(B) - ((lambda*sqrt(1+lambda^2))/(2*sigma))*t(D))                                       # der_sigma2w
      delta3 <-  - (lambda/sigma2)*t(B) + ((1 + 2*lambda^2)/(2*sigma*sqrt(1+lambda^2)))*t(D)                                                      # der_lambdaw
      
      mdelta[1:p,] <- delta1
      mdelta[p+1,] <- delta2
      mdelta[p+2,] <- delta3
    }
  }  
  
  if(family=="T"){
    beta     <- theta[1:p]
    sigma2   <- theta[p+1]
    nu       <- theta[p+2]
    lambda   <- 0 ##theta[p+2]
    mu       <- x%*%beta
    sigma    <- sqrt(sigma2)
    
    A        <- E00*mu
    B        <- E02 - 2*E01*mu + E00*mu^2
    D        <- E11 - E10*mu
    delta1   <- c()
    delta2   <- c()
    ##delta3   <- c()
    mdelta   <- matrix(0, (p+1), n)
    
    if(Perturb==1){
      delta1 <- ((1 + lambda^2)/sigma2)*(t(x)%*%diag(E01[1:n]) - t(x)%*%diag(A[1:n]) - sigma*lambda/sqrt(1 + lambda^2)*t(x)%*%diag(E10[1:n])) # der_betaw
      delta2 <- (-1/(2*sigma2))*(rep(1,n) - ((1 + lambda^2)/sigma2)*t(B) + ((lambda*sqrt(1+lambda^2))/sigma)*t(D))                            # der_sigma2w
      ##delta3 <- (lambda/(1 + lambda^2))*rep(1,n) - (lambda/sigma2)*t(B) + ((1 + 2*lambda^2)/(sigma*sqrt(1+lambda^2)))*t(D) - lambda*t(E20)  # der_lambdaw
      
      mdelta[1:p,] <- delta1
      mdelta[p+1,] <- delta2
      ##mdelta[p+2,] <- delta3
    }
    if(Perturb==2){
      delta1 <- ((1 + lambda^2)/sigma2)*(t(x)%*%diag(E01[1:n]) - t(x)%*%diag(A[1:n]) - sigma*lambda/(2*sqrt(1 + lambda^2))*t(x)%*%diag(E10[1:n])) # der_betaw
      delta2 <- (1/(2*sigma2))*( ((1 + lambda^2)/sigma2)*t(B) - ((lambda*sqrt(1+lambda^2))/(2*sigma))*t(D))                                       # der_sigma2w
      ##delta3 <-  - (lambda/sigma2)*t(B) + ((1 + 2*lambda^2)/(2*sigma*sqrt(1+lambda^2)))*t(D)                                                    # der_lambdaw
      
      mdelta[1:p,] <- delta1
      mdelta[p+1,] <- delta2
      ##mdelta[p+2,] <- delta3
    }
  }
  
  if(family=="N"){
    beta     <- theta[1:p]
    sigma2   <- theta[p+1]
    lambda   <- 0 ##theta[p+2]
    ##nu       <- theta[p+3]
    mu       <- x%*%beta
    sigma    <- sqrt(sigma2)
    
    A        <- E00*mu
    B        <- E02 - 2*E01*mu + E00*mu^2
    D        <- E11 - E10*mu
    delta1   <- c()
    delta2   <- c()
    ##delta3   <- c()
    mdelta   <- matrix(0, (p+1), n)
    
    if(Perturb==1){
      delta1 <- ((1 + lambda^2)/sigma2)*(t(x)%*%diag(E01[1:n]) - t(x)%*%diag(A[1:n]) - sigma*lambda/sqrt(1 + lambda^2)*t(x)%*%diag(E10[1:n])) # der_betaw
      delta2 <- (-1/(2*sigma2))*(rep(1,n) - ((1 + lambda^2)/sigma2)*t(B) + ((lambda*sqrt(1+lambda^2))/sigma)*t(D))                            # der_sigma2w
      ##delta3 <- (lambda/(1 + lambda^2))*rep(1,n) - (lambda/sigma2)*t(B) + ((1 + 2*lambda^2)/(sigma*sqrt(1+lambda^2)))*t(D) - lambda*t(E20)  # der_lambdaw
      
      mdelta[1:p,] <- delta1
      mdelta[p+1,] <- delta2
      ##mdelta[p+2,] <- delta3
    }
    if(Perturb==2){
      delta1 <- ((1 + lambda^2)/sigma2)*(t(x)%*%diag(E01[1:n]) - t(x)%*%diag(A[1:n]) - sigma*lambda/(2*sqrt(1 + lambda^2))*t(x)%*%diag(E10[1:n])) # der_betaw
      delta2 <- (1/(2*sigma2))*( ((1 + lambda^2)/sigma2)*t(B) - ((lambda*sqrt(1+lambda^2))/(2*sigma))*t(D))                                       # der_sigma2w
      ##delta3 <-  - (lambda/sigma2)*t(B) + ((1 + 2*lambda^2)/(2*sigma*sqrt(1+lambda^2)))*t(D)                                                    # der_lambdaw
      
      mdelta[1:p,] <- delta1
      mdelta[p+1,] <- delta2
      ##mdelta[p+2,] <- delta3
    }
  }
  
  If        <- t(mdelta)%*%solve(HQ$MatrizQ)%*%mdelta
  medida    <- diag(If)/sum(diag(If))
  benchmark <- mean(medida) + 3.5*sd(medida)
  
  obs <- c(rep(0,n))
  
  for(i in 1:n){
    if(medida[i] > benchmark){ obs[i] <- i }
  }
  
  influent <- obs[obs!=0]
  
  if(length(influent) == 0){influent <- NULL}
  
  obj.out <- list(If = If, Medida = medida, Benchmark = benchmark, Influent = influent)
  
  return(obj.out)
}