#################################################################################################
## EM algorithm for estimating the parameters of ST-CR Model  (left censored)
## For right censored, use -y
## Last version date: 11/08/2022
#################################################################################################

EM.skewCens <- function(cc, x, y, beta = NULL, sigma2 = NULL, shape = NULL, nu = NULL, get.init = TRUE, show.envelope = "FALSE", error = 0.000001, iter.max = 500, family = "ST"){

  ## cc is a vector nx1 of left-censoring (0 = uncensoring or 1 = censoring)
  ## x is the design matrix of dimension nxp
  ## y vector of responses nx1
  ## Options of family: "ST", "T" , "SN" , "N"

  if ((family != "ST") && (family != "T") && (family != "SN") && (family != "N")) stop(paste("Family", family, "not recognized.\n", sep=" "))
    
  p <- ncol(x)
  n <- nrow(x)
  
  if (get.init == FALSE){
    if((length(beta) != p) || (length(sigma2) == 0) || (length(shape) == 0) || (length(nu) == 0))  stop("The model is not specified correctly.\n")
    
    delta  <- shape / (sqrt(1 + shape^2))
    Delta  <- sqrt(sigma2)*delta
    Gama   <- sigma2 - Delta^2
    mu     <- x%*%beta
  }
  
  if (get.init == TRUE){             # Least square
    reg        <- lm(y ~ x[ , 2:p])
    beta       <- as.vector(coefficients(reg), mode = "numeric")
    sigma2     <- sum((y - x%*%beta)^2)/(n - p)
    shape      <- as.numeric(sign(skewness(y - x%*%beta))*3)
    delta      <- shape / (sqrt(1 + shape^2))
    Delta      <- sqrt(sigma2)*delta
    Gama       <- sigma2 - Delta^2
    mu         <- x%*%beta
    nu         <- 3
  }

  #########################
  ##   Skew-t
  #########################
  
  if(family=="ST"){
    
    cont       <- 0
    criterio   <- 1
    lkante     <- 1
    ychap      <- y
    n.par      <- p + 2
    
    while(criterio > error){
      
      cont     <- (cont + 1)
  
      dj       <- ((y - mu)/sqrt(sigma2))^2
      Mtij2    <- 1/(1 + (Delta^2)*(Gama^(-1)))
      Mtij     <- sqrt(Mtij2)
      mutij    <- Mtij2*Delta*(Gama^(-1))*(y - mu)
      A        <- mutij / Mtij
      cnu      <- 2*gamma((nu+1)/2)/(gamma(nu/2)*sqrt(nu*pi*(1 + shape^2)))
      E        <- (2*(nu)^(nu/2)*gamma((2+nu)/2)*((dj + nu + A^2))^(-(2+nu)/2)) / (gamma(nu/2)*pi*sqrt(sigma2)*dt.ls(y, mu, sigma2, shape, nu))
      u        <- ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(dj + nu)^(-(nu+3)/2)) / (gamma(nu/2)*sqrt(pi)*sqrt(sigma2)*dt.ls(y, mu, sigma2, shape, nu)) )*pt(sqrt((3+nu)/(dj+nu))*A,3+nu)
      
      S1       <- u
      S2       <- (mutij*u + Mtij*E)
      S3       <- (mutij^2*u + Mtij2 + Mtij*mutij*E)
      
      E00      <- S1
      E01      <- y*S1
      E02      <- y^2*S1
      E10      <- S2
      E20      <- S3
      E11      <- y*S2
      sigma2s  <- nu/(nu + 2)*sigma2
      sigma2ss <- nu/((nu + 1)*(1 + shape^2))*sigma2
      
      Aux1     <- cdfSNI(y, mu,  sigma2s, shape, nu+2, type = "ST") 
      Aux11    <- cdfSNI(y, mu, sigma2ss,     0, nu+1, type = "ST")
      Aux2     <- cdfSNI(y, mu,   sigma2, shape,   nu, type = "ST")
      
      mu1      <- mu[cc==1]
      y1       <- y[cc==1]
      np       <- length(mu1)
      aux1MomW <- aux2MomS <- aux3MomS <- matrix(0, np, 2)
      
      for(j in 1:np){      
        A1a          <- MomenSNI(mu1[j],  sigma2s, shape, nu+2, delta=NULL, Lim1=-Inf, Lim2=y1[j], type = "ST") # CalMom(mu,sigma2,nu,y,type)
        A2a          <- MomenSNI(mu1[j], sigma2ss,     0, nu+1, delta=NULL, Lim1=-Inf, Lim2=y1[j], type = "ST")
        A3a          <- MomenSNI(mu1[j],   sigma2, shape,   nu, delta=NULL, Lim1=-Inf, Lim2=y1[j], type = "ST")
        aux1MomW[j,] <- c(A1a$EUY1, A1a$EUY2)
        aux2MomS[j,] <- c(A2a$EUY1, A2a$EUY2)
        aux3MomS[j,] <- c(A3a$EUY1, A3a$EUY2)
      }
      
      P1aux <- P2aux <- WW <- u
      P1aux[cc==1]   <- aux1MomW[,1]
      P2aux[cc==1]   <- aux1MomW[,2]
      WW[cc==1]      <- aux2MomS[,1]
      
      Wphi           <- Aux11/Aux2
      E00Aux         <- Aux1/Aux2
      E01Aux         <- Aux1/Aux2*P1aux
      E02Aux         <- Aux1/Aux2*P2aux
      E10Aux         <- Delta/(Gama + Delta^2)*(E01Aux - E00Aux*mu) + Mtij*cnu*Wphi
      E20Aux         <- (Delta/(Gama + Delta^2))^2*(E02Aux - 2*E01Aux*mu + mu^2*E00Aux) + Mtij*(Delta/(Gama + Delta^2))*cnu*Wphi*(WW - mu) + Mtij2
      E11Aux         <- Delta/(Gama + Delta^2)*(E02Aux - E01Aux*mu) + Mtij*Wphi*cnu*WW
      
      ychap[cc==1]   <- aux3MomS[,1]
      E00[cc==1]     <- E00Aux[cc==1]        # E[U]
      E01[cc==1]     <- E01Aux[cc==1]        # E[UY]  
      E02[cc==1]     <- E02Aux[cc==1]        # E[UY2]   
      E10[cc==1]     <- E10Aux[cc==1]        # E[UT]
      E20[cc==1]     <- E20Aux[cc==1]        # E[UT2]
      E11[cc==1]     <- E11Aux[cc==1]        # E[UTY]  
      
      beta           <- solve(t(x)%*%diag(c(E00))%*%x)%*%(t(x)%*%matrix(E01 - E10*Delta, n, 1))
      mu             <- x%*%beta
      Delta          <- sum(E11 - E10*mu)/sum(E20)
      Gama           <- sum(E02 - 2*E01*mu + E00*mu^2 + Delta^2*E20 - 2*Delta*E11 + 2*Delta*E10*mu)/n		
      sigma2         <- Gama + Delta^2
      shape          <- ((sigma2^(-1/2))*Delta)/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
      
      f <- function(nu){
        sum(log(dt.ls(y[cc==0], mu[cc==0], sigma2, shape, nu))) + sum(log(cdfSNI(y[cc==1], mu[cc==1], sigma2, shape, nu, type = "ST")))
      }
      
      nu             <- optimize(f, c(2.01,150), tol = 0.0000000001, maximum = TRUE)$maximum
      lk             <- f(nu)
      logver         <- lk
      criterio       <- abs((lk/lkante - 1))
      lkante         <- lk
     
      if (cont==iter.max){
        criterio     <- error/10
      }
    }
    
    teta_novo        <- matrix(c(beta, sigma2, shape, nu), ncol = 1) # to compute the number of parameters
    #teta_novo1       <- matrix(c(beta, sigma2, shape, nu), ncol = 1)
    
    #########################
    ##  Information Matrix
    #########################
    
    sbeta            <- c()
    ssigma2          <- c()
    slambda          <- c()
    MIE              <- matrix(0, p+2, p+2)
    S                <- matrix(0,   1, p+2)
    sigma            <- sqrt(sigma2)
    lambda           <- shape
    
    for(i in 1:n){
      sbeta          <- ((1 + lambda^2)/sigma2)*(E01[i]*t(x[i,]) - E00[i]*t(x[i,])*mu[i] - Delta*E10[i]*t(x[i,]))
      ssigma2        <- -1/(2*sigma2) + ((1 + lambda^2)/(2*sigma2^2))*(E02[i] - 2*E01[i]*mu[i] + (t(mu[i])%*%mu[i])*E00[i]) - ((lambda*sqrt(1+lambda^2))/(2*sigma^3))*(E11[i] - E10[i]*mu[i])
      slambda        <- lambda/(1 + lambda^2) - (lambda/sigma2)*(E02[i] - 2*E01[i]*mu[i] + E00[i]*(t(mu[i])%*%mu[i])) + ((1 + 2*lambda^2)/(sigma*sqrt(1+lambda^2)))*(E11[i] - E10[i]*mu[i]) - lambda*E20[i]      
      S              <- c(sbeta, ssigma2, slambda)
      MIE1           <- S%*%t(S)
      ind            <- lower.tri(MIE1)
      MIE1[ind]      <- t(MIE1)[ind]
      MIE            <- MIE1 + MIE
    }
    
    se               <- sqrt(diag(solve(MIE)))    # standard errors
    
    if(show.envelope=="TRUE"){
      envelop        <- EnvelopeRMT(teta_novo, y, x, cc, family = "ST")
    }
  }
  
  #########################
  ## Student-t
  #########################
  
  if(family=="T"){
    
    cont       <- 0
    criterio   <- 1
    lkante     <- 1
    shape      <- 0
    ychap      <- y
    n.par      <- p + 1
    
    while(criterio > error){
      
      cont     <- (cont+1)
      
      dj       <- ((y - mu)/sqrt(sigma2))^2
      Mtij2    <- 1/(1 + (Delta^2)*(Gama^(-1)))
      Mtij     <- sqrt(Mtij2)
      mutij    <- Mtij2*Delta*(Gama^(-1))*(y - mu)
      A        <- mutij / Mtij
      cnu      <- 2*gamma((nu + 1)/2)/(gamma(nu/2)*sqrt(nu*pi*(1 + shape^2)))
      E        <- (2*(nu)^(nu/2)*gamma((2+nu)/2)*((dj + nu + A^2))^(-(2+nu)/2)) / (gamma(nu/2)*pi*sqrt(sigma2)*dt.ls(y, mu, sigma2, shape, nu))
      u        <- ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(dj + nu)^(-(nu+3)/2)) / (gamma(nu/2)*sqrt(pi)*sqrt(sigma2)*dt.ls(y, mu, sigma2, shape, nu)) )*pt(sqrt((3+nu)/(dj+nu))*A,3+nu)
      
      S1       <- u
      S2       <- (mutij*u + Mtij*E)
      S3       <- (mutij^2*u + Mtij2 + Mtij*mutij*E)
      
      E00      <- S1
      E01      <- y*S1
      E02      <- y^2*S1
      E10      <- S2
      E20      <- S3
      E11      <- y*S2
      sigma2s  <- nu/(nu + 2)*sigma2
      sigma2ss <- nu/((nu + 1)*(1 + shape^2))*sigma2
      
      Aux1     <- cdfSNI(y, mu,  sigma2s, shape, nu+2, type = "T") 
      Aux11    <- cdfSNI(y, mu, sigma2ss,     0, nu+1, type = "T")
      Aux2     <- cdfSNI(y, mu,   sigma2, shape,   nu, type = "T")
      
      mu1      <- mu[cc==1]
      y1       <- y[cc==1]
      np       <- length(mu1)
      aux1MomW <- aux2MomS <- aux3MomS <- matrix(0, np, 2)
      
      for(j in 1:np){      
        A1a          <- MomenSNI(mu1[j],  sigma2s, 0, nu+2, delta=NULL, Lim1=-Inf, Lim2=y1[j], type = "ST") # CalMom(mu,sigma2,nu,y,type)
        A2a          <- MomenSNI(mu1[j], sigma2ss, 0, nu+1, delta=NULL, Lim1=-Inf, Lim2=y1[j], type = "ST")
        A3a          <- MomenSNI(mu1[j],   sigma2, 0,   nu, delta=NULL, Lim1=-Inf, Lim2=y1[j], type = "ST")
        aux1MomW[j,] <- c(A1a$EUY1, A1a$EUY2)
        aux2MomS[j,] <- c(A2a$EUY1, A2a$EUY2)
        aux3MomS[j,] <- c(A3a$EUY1, A3a$EUY2)
      }
      
      P1aux <- P2aux <- WW <- u
      P1aux[cc==1]   <- aux1MomW[,1]
      P2aux[cc==1]   <- aux1MomW[,2]
      WW[cc==1]      <- aux2MomS[,1]
      
      Wphi           <- Aux11/Aux2
      E00Aux         <- Aux1/Aux2
      E01Aux         <- Aux1/Aux2*P1aux
      E02Aux         <- Aux1/Aux2*P2aux
      E10Aux         <- Delta/(Gama + Delta^2)*(E01Aux - E00Aux*mu) + Mtij*cnu*Wphi
      E20Aux         <- (Delta/(Gama + Delta^2))^2*(E02Aux - 2*E01Aux*mu + mu^2*E00Aux) + Mtij*(Delta/(Gama + Delta^2))*cnu*Wphi*(WW - mu) + Mtij2
      E11Aux         <- Delta/(Gama + Delta^2)*(E02Aux - E01Aux*mu) + Mtij*Wphi*cnu*WW
      
      ychap[cc==1]   <- aux3MomS[,1] 
      E00[cc==1]     <- E00Aux[cc==1]        # E[U]
      E01[cc==1]     <- E01Aux[cc==1]        # E[UY]  
      E02[cc==1]     <- E02Aux[cc==1]        # E[UY2]   
      E10[cc==1]     <- E10Aux[cc==1]        # E[UT]
      E20[cc==1]     <- E20Aux[cc==1]        # E[UT2]
      E11[cc==1]     <- E11Aux[cc==1]        # E[UTY]  
      
      beta           <- solve(t(x)%*%diag(c(E00))%*%x)%*%(t(x)%*%matrix(E01 - E10*Delta, n, 1))
      mu             <- x%*%beta
      Delta          <- 0
      Gama           <- sum(E02 - 2*E01*mu + E00*mu^2 + Delta^2*E20 - 2*Delta*E11 + 2*Delta*E10*mu)/n		
      sigma2         <- Gama + Delta^2
      shape          <- ((sigma2^(-1/2))*Delta)/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
  
      f <- function(nu){
        sum(log(dT(cc, y, mu, sigma2, nu)))
      }
        
      nu             <- optimize(f, c(2.01, 150), tol = 0.0000000001, maximum = TRUE)$maximum 
      lk             <- f(nu)
      logver         <- lk
      criterio       <- abs((lk/lkante - 1))
      lkante         <- lk
      
      if (cont==iter.max){
        criterio     <- error/10
      }
    }
   
    teta_novo        <- matrix(c(beta, sigma2,    nu), ncol = 1)
    teta_novo1       <- matrix(c(beta, sigma2, 0, nu), ncol = 1)  
    
    #########################
    ##  Information Matrix
    #########################
      
    sbeta            <- c()
    ssigma2          <- c()
    MIE              <- matrix(0, p+1, p+1)
    S                <- matrix(0,   1, p+1)
    sigma            <- sqrt(sigma2)
    lambda           <- 0
    
    for(i in 1:n){
      sbeta          <- ((1 + lambda^2)/sigma2)*(E01[i]*t(x[i,]) - E00[i]*t(x[i,])*mu[i] - Delta*E10[i]*t(x[i,]))
      ssigma2        <- -1/(2*sigma2) + ((1 + lambda^2)/(2*sigma2^2))*(E02[i] - 2*E01[i]*mu[i] + (t(mu[i])%*%mu[i])*E00[i]) - ((lambda*sqrt(1+lambda^2))/(2*sigma^3))*(E11[i] - E10[i]*mu[i])
      S              <- c(sbeta, ssigma2)
      MIE1           <- S%*%t(S)
      ind            <- lower.tri(MIE1)
      MIE1[ind]      <- t(MIE1)[ind]
      MIE            <- MIE1 + MIE
    }
    
    se               <- sqrt(diag(solve(MIE)))    # standard errors
    
    if(show.envelope=="TRUE"){
      envelop        <- EnvelopeRMT(teta_novo1, y, x, cc, family = "ST")
    }
  }
  
  #########################
  ## Skew-Normal
  #########################
  
  if(family=="SN"){
    
    cont       <- 0
    criterio   <- 1
    lkante     <- 1
    ychap      <- y
    n.par      <- p + 2
    
    while(criterio > error){
      
      cont     <- (cont+1)
      
      dj       <- ((y - mu)/sqrt(sigma2))^2
      Mtij2    <- 1/(1 + (Delta^2)*(Gama^(-1)))
      Mtij     <- sqrt(Mtij2)
      mutij    <- Mtij2*Delta*(Gama^(-1))*(y - mu)
      A        <- mutij / Mtij
      cnu      <- 2/sqrt(2*pi*(1 + shape^2))
      prob     <- pnorm(mutij/Mtij)
      
      if(length(which(prob == 0)) > 0) prob[which(prob == 0)] <- .Machine$double.xmin
      
      E        <- dnorm(mutij/Mtij)/prob
      u        <- rep(1, n)
      
      S1       <- u
      S2       <- (mutij*u + Mtij*E)
      S3       <- (mutij^2*u + Mtij2 + Mtij*mutij*E)
      
      E00      <- S1
      E01      <- y*S1
      E02      <- y^2*S1
      E10      <- S2
      E20      <- S3
      E11      <- y*S2
      #ychap    <-y
      sigma2s  <- (1/(1 + shape^2))*sigma2
      
      Aux11    <- cdfSNI(y, mu, sigma2s,     0, NULL, type = "SN")
      Aux2     <- cdfSNI(y, mu,  sigma2, shape, NULL, type = "SN")
      
      mu1      <- mu[cc==1]
      y1       <- y[cc==1]
      np       <- length(mu1)
      aux1MomW <- aux2MomS <- aux3MomS <- matrix(0, np, 2)
       
      for(j in 1:np){
        A1a          <- MomenSNI(mu1[j],  sigma2, shape, NULL, delta=NULL, Lim1=-Inf, Lim2=y1[j], type = "SN") # CalMom(mu,sigma2,nu,y,type)
        A2a          <- MomenSNI(mu1[j], sigma2s,     0, NULL, delta=NULL, Lim1=-Inf, Lim2=y1[j], type = "SN")
        A3a          <- MomenSNI(mu1[j],  sigma2,     0, NULL, delta=NULL, Lim1=-Inf, Lim2=y1[j], type = "SN")
        aux1MomW[j,] <- c(A1a$EUY1, A1a$EUY2)
        aux2MomS[j,] <- c(A2a$EUY1, A2a$EUY2)
        aux3MomS[j,] <- c(A3a$EUY1, A3a$EUY2)
      }
      
      Wphi           <- Aux11/Aux2
      E00Aux         <- E01Aux <- E02Aux <- WW1 <- u
      E01Aux[cc==1]  <- aux1MomW[,1]
      E02Aux[cc==1]  <- aux1MomW[,2]
      WW1[cc==1]     <- aux2MomS[,1]
      E10Aux         <- Delta/(Gama + Delta^2)*(E01Aux - E00Aux*mu) + Mtij*cnu*Wphi
      E20Aux         <- (Delta/(Gama + Delta^2))^2*(E02Aux - 2*E01Aux*mu + mu^2*E00Aux) + Mtij*(Delta/(Gama + Delta^2))*cnu*Wphi*(WW1 - mu) + Mtij2
      E11Aux         <- Delta/(Gama + Delta^2)*(E02Aux - E01Aux*mu) + Mtij*Wphi*cnu*WW1
      
      ychap[cc==1]   <- aux3MomS[,1]
      E00[cc==1]     <- E00Aux[cc==1]        # E[U]
      E01[cc==1]     <- E01Aux[cc==1]        # E[UY]  
      E02[cc==1]     <- E02Aux[cc==1]        # E[UY2]   
      E10[cc==1]     <- E10Aux[cc==1]        # E[UT]
      E20[cc==1]     <- E20Aux[cc==1]        # E[UT2]
      E11[cc==1]     <- E11Aux[cc==1]        # E[UTY]  
      
      beta           <- solve(t(x)%*%diag(E00)%*%x)%*%(t(x)%*%(E01 - E10*Delta))
      mu             <- x%*%beta
      Delta          <- sum(E11 - E10*mu)/sum(E20)
      Gama           <- sum(E02 - 2*E01*mu + E00*mu^2 + Delta^2*E20 - 2*Delta*E11 + 2*Delta*E10*mu)/n		
      sigma2         <- Gama + Delta^2
      shape          <- ((sigma2^(-1/2))*Delta)/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
      
      auxpdf         <- dSN(   y[cc==0], mu[cc==0], sigma2, shape)
      auxcdf         <- cdfSNI(y[cc==1], mu[cc==1], sigma2, shape, NULL, type = "SN")
      lk             <- sum(log(auxpdf)) + sum(log(auxcdf))
      logver         <- lk
      criterio       <- abs((lk/lkante - 1))
      lkante         <- lk
      
      if (cont==iter.max){
        criterio     <- error/10
      }
    }
    
    ychap            <- E01
    
    teta_novo        <- matrix(c(beta, sigma2, shape), ncol = 1)
    #teta_novo1       <-matrix(c(beta,sigma2,shape),ncol=1)
    
    #########################
    ##  Information Matrix
    #########################
    
    sbeta            <- c()
    ssigma2          <- c()
    slambda          <- c()
    MIE              <- matrix(0, p+2, p+2)
    S                <- matrix(0,   1, p+2)
    sigma            <- sqrt(sigma2)
    lambda           <- shape
    
    for(i in 1:n){
      sbeta          <- ((1 + lambda^2)/sigma2)*(E01[i]*t(x[i,]) - E00[i]*t(x[i,])*mu[i] - Delta*E10[i]*t(x[i,]))
      ssigma2        <- -1/(2*sigma2) + ((1 + lambda^2)/(2*sigma2^2))*(E02[i] - 2*E01[i]*mu[i] + (t(mu[i])%*%mu[i])*E00[i]) - ((lambda*sqrt(1+lambda^2))/(2*sigma^3))*(E11[i] - E10[i]*mu[i])
      slambda        <- lambda/(1 + lambda^2) - (lambda/sigma2)*(E02[i] - 2*E01[i]*mu[i] + E00[i]*(t(mu[i])%*%mu[i])) + ((1 + 2*lambda^2)/(sigma*sqrt(1+lambda^2)))*(E11[i] - E10[i]*mu[i]) - lambda*E20[i]      
      S              <- c(sbeta, ssigma2, slambda)
      MIE1           <- S%*%t(S)
      ind            <- lower.tri(MIE1)
      MIE1[ind]      <- t(MIE1)[ind]
      MIE            <- MIE1 + MIE
    }
    
    se               <- sqrt(diag(solve(MIE)))    # standard errors
    
    if(show.envelope=="TRUE"){
      envelop        <- EnvelopeRMT(teta_novo, y, x, cc, family = "SN")
    }   
  }

  #########################
  ## Normal
  #########################
  
  if(family=="N"){
    
    cont       <- 0
    criterio   <- 1
    lkante     <- 1
    shape      <- 0
    ychap      <- y
    n.par      <- p + 1
    
    while(criterio > error){
      
      cont     <- (cont+1)
      
      dj       <- ((y - mu)/sqrt(sigma2))^2
      Mtij2    <- 1/(1 + (Delta^2)*(Gama^(-1)))
      Mtij     <- sqrt(Mtij2)
      mutij    <- Mtij2*Delta*(Gama^(-1))*(y - mu)
      A        <- mutij / Mtij
      cnu      <- 2/sqrt(2*pi*(1 + shape^2))
      prob     <- pnorm(mutij/Mtij)
      
      if(length(which(prob == 0)) > 0) prob[which(prob == 0)] <- .Machine$double.xmin
      
      E        <- dnorm(mutij/Mtij)/prob
      u        <- rep(1, n)
      
      S1       <- u
      S2       <- (mutij*u + Mtij*E)
      S3       <- (mutij^2*u + Mtij2 + Mtij*mutij*E)
      
      E00      <- S1
      E01      <- y*S1
      E02      <- y^2*S1
      E10      <- S2
      E20      <- S3
      E11      <- y*S2
      
      sigma2s  <- (1/(1 + shape^2))*sigma2
      
      Aux11    <- cdfSNI(y, mu, sigma2s,     0, NULL, type = "N")
      Aux2     <- cdfSNI(y, mu,  sigma2, shape, NULL, type = "N")
       
      mu1      <- mu[cc==1]
      y1       <- y[cc==1]
      np       <- length(mu1)
      aux1MomW <- aux2MomS <- matrix(0, np, 2)
       
      for(j in 1:np){
        A1a          <- MomenSNI(mu1[j],  sigma2, shape, NULL, delta=NULL, Lim1=-Inf, Lim2=y1[j], type = "SN")
        A2a          <- MomenSNI(mu1[j], sigma2s,     0, NULL, delta=NULL, Lim1=-Inf, Lim2=y1[j], type = "SN")
        aux1MomW[j,] <- c(A1a$EUY1, A1a$EUY2)
        aux2MomS[j,] <- c(A2a$EUY1, A2a$EUY2)
      }
      
      Wphi           <- Aux11/Aux2
      E00Aux         <- E01Aux <- E02Aux <- WW1 <- u
      E01Aux[cc==1]  <- aux1MomW[,1]
      E02Aux[cc==1]  <- aux1MomW[,2]
      WW1[cc==1]     <- aux2MomS[,1]
      E10Aux         <- Delta/(Gama + Delta^2)*(E01Aux - E00Aux*mu) + Mtij*cnu*Wphi
      E20Aux         <- (Delta/(Gama + Delta^2))^2*(E02Aux - 2*E01Aux*mu + mu^2*E00Aux) + Mtij*(Delta/(Gama + Delta^2))*cnu*Wphi*(WW1 - mu) + Mtij2
      E11Aux         <- Delta/(Gama + Delta^2)*(E02Aux - E01Aux*mu) + Mtij*Wphi*cnu*WW1
      
      E00[cc==1]     <- E00Aux[cc==1]        # E[U]
      E01[cc==1]     <- E01Aux[cc==1]        # E[UY]  
      E02[cc==1]     <- E02Aux[cc==1]        # E[UY2]   
      E10[cc==1]     <- E10Aux[cc==1]        # E[UT]
      E20[cc==1]     <- E20Aux[cc==1]        # E[UT2]
      E11[cc==1]     <- E11Aux[cc==1]        # E[UTY]  
      
      beta           <- solve(t(x)%*%diag(E00)%*%x)%*%(t(x)%*%(E01 - E10*Delta))
      mu             <- x%*%beta
      Delta          <- 0 # sum(E11 - E10*mu)/sum(E20)
      Gama           <- sum(E02 - 2*E01*mu + E00*mu^2 + Delta^2*E20 - 2*Delta*E11 + 2*Delta*E10*mu)/n		
      sigma2         <- Gama + Delta^2
      shape          <- ((sigma2^(-1/2))*Delta)/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
      
      auxpdf         <- dSN(   y[cc==0], mu[cc==0], sigma2, shape)
      auxcdf         <- cdfSNI(y[cc==1], mu[cc==1], sigma2, shape, NULL, type = "SN")
      lk             <- sum(log(auxpdf)) + sum(log(auxcdf))
      logver         <- lk
      criterio       <- abs((lk/lkante - 1))
      lkante         <- lk
      
      if (cont==iter.max){     
        criterio     <- error/10
      }
    }
    
    ychap            <- E01
    #ychap[cc==0]     <- mu[cc==0] 
    teta_novo        <- matrix(c(beta, sigma2   ), ncol = 1)
    teta_novo1       <- matrix(c(beta, sigma2, 0), ncol = 1)
    
    #########################
    ##  Information Matrix
    #########################
    
    sbeta            <- c()
    ssigma2          <- c()
    #slambda          <- c()
    MIE              <- matrix(0, p+1, p+1)
    S                <- matrix(0,   1, p+1)
    sigma            <- sqrt(sigma2)
    lambda           <- 0
    
    for(i in 1:n){
      sbeta          <- ((1 + lambda^2)/sigma2)*(E01[i]*t(x[i,]) - E00[i]*t(x[i,])*mu[i] - Delta*E10[i]*t(x[i,]))
      ssigma2        <- -1/(2*sigma2) + ((1 + lambda^2)/(2*sigma2^2))*(E02[i] - 2*E01[i]*mu[i] + (t(mu[i])%*%mu[i])*E00[i]) - ((lambda*sqrt(1 + lambda^2))/(2*sigma^3))*(E11[i] - E10[i]*mu[i])
      S              <- c(sbeta, ssigma2)
      MIE1           <- S%*%t(S)
      ind            <- lower.tri(MIE1)
      MIE1[ind]      <- t(MIE1)[ind]
      MIE            <- MIE1 + MIE
    }
    
    se               <- sqrt(diag(solve(MIE)))    # standard errors
    
    if(show.envelope=="TRUE"){
      envelop        <- EnvelopeRMT(teta_novo1, y, x, cc, family = "SN")
    }
  }
 
  aic                <- -2*lk +             2*length(teta_novo)
  bic                <- -2*lk +        log(n)*length(teta_novo)
  caic               <- -2*lk +  (log(n) + 1)*length(teta_novo)
  hq                 <- -2*lk + 2*log(log(n))*length(teta_novo)
 
  #####################################################
  ### Summary of the model results.
  #####################################################

  namesx <- ('(Intercept)     ')
  
  if(ncol(as.matrix(x)) > 1){
    for(i in 2:ncol(as.matrix(x))){
      namesx <- cbind(namesx, paste("x", i-1,"     ", sep = ""))
    }
  }
  
  t.est            <- teta_novo[1:n.par]/se[1:n.par]
  p.est            <- 2*(1 - pt(abs(t.est), n-1))
  
  param            <- round(cbind(teta_novo[1:n.par], se[1:n.par], t.est, p.est), digits = 3)
  namespar         <- colnames(x)
  colx             <- ncol(as.matrix(x))
  
  if(length(namespar) == 0){
    namespar       <- namesx[1:colx]
  }
  
  if(family=="N" || family=="T" ){
    dimnames(param)  <- list(c(namespar, expression(sigma^2)), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  }
  
  if(family=="SN" || family=="ST" ){
    dimnames(param)  <- list(c(namespar, expression(sigma^2), expression(lambda)), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  }
  
  cat('\n')
  cat('-------------------------------------------------------\n')
  cat('                  EM estimates and SE                  \n')
  cat('                                                       \n')
  cat('Coefficients:                                          \n')
  print(param)
  cat('-------------------------------------------------------\n')
  
  if(family=="T" || family=="ST" ){
    cat("Degree of freedom: ", round(nu, 3), "\n")
  }
  
  cat('-------------------------------------------------------\n')
  cat('\r \n')
  critFin           <- c(logver, aic, bic, caic, hq)
  critFin           <- round(t(as.matrix(critFin)), digits = 3)
  dimnames(critFin) <- list(c("Value"), c("Loglik", "AIC", "BIC","CAIC", "HQ"))
  cat('\n')
  cat('Model selection criteria                               \n')
  cat('-------------------------------------------------------\n')
  print(critFin)
  cat('-------------------------------------------------------\n')
  cat('\r \n')
  
  return(list(theta=teta_novo, Se=se, nu=nu, iter=cont, logver=logver, AIC=aic,
              BIC=bic, CAIC=caic, HQ=hq, yhat=ychap, E00=E00, E01=E01, E02=E02,
              E10=E10, E20=E20, E11=E11))	
}