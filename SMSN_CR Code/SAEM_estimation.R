### =================================================================== ###
### Likelihood Based Inference for Censored Linear Regression Models    ###
###         with Scale Mixtures of Skew-Normal Distributions            ###
###  Authors: Thalita do Bem Mattos, Aldo M. Garay and Victor H. Lachos ###
### =================================================================== ###

## -------------- ##
## SAEM - SMSN-CR ##
## -------------- ##

SAEM_EST <- function(y,X,cc,cens,LS=NULL,precisao=0.0001,MaxIter=200, M=20, pc=0.30, dist="SN",nu.fixed=FALSE,nu=NULL,show.convergence="TRUE")
{ 
  
  # y : variable response
  # X: matrix design
  # cc: vector of 0 or 1's of censors ; 1 = censored / 0=non-censored
  # cens: censorship: right -> "right" ; left -> "left" ; interval -> "interval"
  # LS: If the censorship is interval, it contains the upper limits of the observed intervals
  # precisao: precision
  # MaxIter : maximum iteration
  # M: Monte Carlo sample size (10 to 20)
  # pc: cutoff point (0 < c < 1)
  # dist: distribution: Skew-Normal -> "SN" ; Skew-T -> "ST" ; Skew-Slash -> "SSL"; Skew-Contaminated Normal -> "SCN"
  # nu.fixed: if FALSE, the parameter nu is estimated; TRUE, the parameter nu is fixed.
  # nu : if ST and SSL : numeric ; if SCN: a vector c(x,y). 
  
  start.time <- Sys.time()

  p <- ncol(X)
  n <- nrow(X)
  
  ## n = n1 + n2
  n1 <- length(y[cc==0])
  n2 <- length(y[cc==1])
  
  if(cens=="left"){cens = 1}
  if(cens=="right"){cens = 2}
  if(cens=="interval"){cens = 3}
  
  # INITIAL VALUES
  
  betas <- solve(t(X)%*%X)%*%t(X)%*%y
  sigma2 <- sum((y-X%*%betas)^2)/(n-p)
  lambda <- as.numeric(sign(skewness(y-X%*%betas))*3) 
  
  deltinha <- lambda/sqrt(1+(lambda^2))
  Delta <-  sqrt(sigma2)*deltinha
  tau <- sigma2 - Delta^2 
  
  if(dist=="SN")
  {
    nu <- nu
    teta0 <- rbind(betas,sigma2,lambda)
  }
  if(dist=="ST"|dist=="SSL")
  {
    nu <- nu
    teta0 <- rbind(betas,sigma2,lambda,nu)
  }
  if(dist=="SCN")
  {
    nu1 <- nu[1]
    nu2 <- nu[2] 
    teta0 <- rbind(betas,sigma2,lambda,nu1,nu2) 
  }
  
  #### 
  
  critval <- critval2 <- 1
  
  delta1 <- 0.001
  delta2 <- precisao
  
  count <- 0
  
  teta.vec <- matrix(nrow=length(teta0),ncol=MaxIter+1)
  teta.vec[,1] <- teta0
  
  log0 <- logVerosCens(teta0,y,X,cc,cens,dist)
  log.vec <- matrix(nrow=MaxIter+1, ncol=1)
  log.vec[1,] <- log0
  
  if(pc==1)
  {
    seqq=rep(1,pc*MaxIter)
  }else
  {
    seqq = c(rep(1,pc*MaxIter),(1/((((pc*MaxIter)+1):MaxIter)-(pc*MaxIter))))
    seqq = c(rep(1,MaxIter-length(seqq)),seqq)
  }
  
  # MEMORY SAEM
  S1 <- S2 <- S3 <- S4 <- S5 <- S6 <- S7 <- matrix(0,nrow=MaxIter+1,ncol=n2)
  
  pb <- tkProgressBar(title = "progress bar", min = 0, max = MaxIter, width = 300)
  setTkProgressBar(pb, 0, label=paste("Iter ",0,"/",MaxIter,"     -     ",0,"% done",sep = ""))
  
  ####
  while (critval < 3 && critval2 < 3)
  {
    count <- count + 1
    
    #STEP E: 
    
    EU <- EUY <- EUY2 <- EUT <- EUT2 <- EUYT <- ElogU <- matrix(0,n,1)
    sum_ui <- sum_uiyi <- sum_uiy2i <- sum_uiti <- sum_uit2i <- sum_uiyiti <- sum_logui <- matrix(0,n2,1)
    
    # IF ci = 0
    NCen <- NCensEsp(y=y,X=X,betas=betas,sigma2=sigma2,lambda=lambda,nu=nu,dist=dist)
    
    EU <- NCen$EU
    EUY <- NCen$EUY
    EUY2 <- NCen$EUY2
    EUT <- NCen$EUT
    EUT2 <- NCen$EUT2
    EUYT <- NCen$EUYT 
    ElogU <- NCen$ElogU
    
    # IF ci = 1
    if(sum(cc)>0)
    { 
      Cens <- CensEsp(y=y,X=X,betas=betas,sigma2=sigma2,lambda=lambda,nu=nu,cc=cc,cens=cens,LS=LS,N=M,dist=dist)
      
      E_ui <- Cens$EU
      E_uiyi <- Cens$EUY
      E_uiy2i <- Cens$EUY2
      E_uiti <- Cens$EUT
      E_uit2i <- Cens$EUT2
      E_uiyiti <- Cens$EUYT
      E_logui <- Cens$ElogU
      
      S1[count+1,] <- S1[count,] + seqq[count]*(E_ui - S1[count,])
      S2[count+1,] <- S2[count,] + seqq[count]*(E_uiyi - S2[count,])
      S3[count+1,] <- S3[count,] + seqq[count]*(E_uiy2i - S3[count,])
      S4[count+1,] <- S4[count,] + seqq[count]*(E_uiti - S4[count,])
      S5[count+1,] <- S5[count,] + seqq[count]*(E_uit2i - S5[count,])
      S6[count+1,] <- S6[count,] + seqq[count]*(E_uiyiti - S6[count,])
      S7[count+1,] <- S7[count,] + seqq[count]*(E_logui - S7[count,])
      
      EU[cc==1] <- S1[count+1,]
      EUY[cc==1] <- S2[count+1,]
      EUY2[cc==1] <- S3[count+1,]
      EUT[cc==1] <- S4[count+1,]
      EUT2[cc==1] <- S5[count+1,] 
      EUYT[cc==1] <- S6[count+1,]
      ElogU[cc==1] <- S7[count+1,] 
    }
    
    #STEP M: 
    
    aux1 <- Diagonal(n,EU)
    aux2 <- t(X)%*%aux1%*%X
    aux3 <- t(X)%*%EUY
    aux4 <- t(X)%*%EUT 
    
    betas <- as.matrix(solve(aux2)%*%(aux3 - Delta*aux4))
    
    aux5 <- sum(EUT2)
    auxbeta <- X%*%betas
        
    tau <- as.numeric(sum(EUY2 - 2*auxbeta*EUY + EU*auxbeta*auxbeta -2*Delta*EUYT + 2*Delta*auxbeta*EUT + (Delta^2)*EUT2)/n)
    
    Delta <- as.numeric(sum(EUYT - auxbeta*EUT)/aux5)
    
    sigma2 <- tau + Delta^2
    lambda <- Delta/sqrt(tau)

    ## Estimation NU
    
    if(nu.fixed==FALSE)
    {
      if(dist=="ST")
      {
        nu <- optimize(f=logVerosCens.nu, interval=c(1.01,150),tol=0.00001,theta=c(betas,sigma2,lambda),y=y,X=X,cc=cc,cens=cens,dist=dist)$minimum
      }
      if(dist=="SSL")
      {
        nu <- optimize(f=logVerosCens.nu, interval=c(1.1,150),tol=0.00001,theta=c(betas,sigma2,lambda),y=y,X=X,cc=cc,cens=cens,dist=dist)$minimum
      }
      if(dist=="SCN")
      {
        nu <- optim(c(nu1,nu2), method = "L-BFGS-B", logVerosCens.nu, lower = rep(0.01, 2), upper = rep(0.99,2), theta=c(betas,sigma2,lambda), y=y,X=X,cc=cc,cens=cens,dist=dist,hessian=TRUE)$par
        nu1 <- nu[1]
        nu2 <- nu[2]
      }
    }
    

    

    ####### MATRIZ DE INFORMAÇÃO EMPÍRICA
    sbeta <- c()
    ssigma2 <- c()
    slambda <- c()
    MIE <- matrix(0,p+2,p+2)
    S <- matrix(0,1,p+2)
    sigma <- sqrt(sigma2)
    for(i in 1:n)
    {
      sbeta <- ((1+lambda^2)/sigma2)*(EUY[i]*t(X[i,]) - EU[i]*t(X[i,])*auxbeta[i] - Delta*EUT[i]*t(X[i,]))
      ssigma2 <- -1/(2*sigma2) + ((1+lambda^2)/(2*sigma2^2))*(EUY2[i] - 2*EUY[i]*auxbeta[i] + (t(auxbeta[i])%*%auxbeta[i])*EU[i]) - ((lambda*sqrt(1+lambda^2))/(2*sigma^3))*(EUYT[i] - EUT[i]*auxbeta[i])
      slambda <- lambda/(1+lambda^2) - (lambda/sigma2)*(EUY2[i] - 2*EUY[i]*auxbeta[i] + EU[i]*(t(auxbeta[i])%*%auxbeta[i])) +  ((1+ 2*lambda^2)/(sigma*sqrt(1+lambda^2)))*(EUYT[i] - EUT[i]*auxbeta[i]) - lambda*EUT2[i]      
      S <- c(sbeta,ssigma2,slambda)
      MIE1 <- S%*%t(S)
      ind <- lower.tri(MIE1)
      MIE1[ind] <- t(MIE1)[ind]
      MIE <- MIE1 + MIE
    }
    se <- sqrt(diag(solve(MIE)))
    
    ########### 
    
    if(dist=="SN")
    {
      teta <- rbind(betas,sigma2,lambda)
    }
    if(dist=="ST"|dist=="SSL")
    {
      teta <- rbind(betas,sigma2,lambda,nu)
    }
    if(dist=="SCN")
    {
      teta <- rbind(betas,sigma2,lambda,nu1,nu2)
    }
 
    teta.vec[,count+1] <- teta
    
    criterio  <- abs(teta-teta0)/(abs(teta0)+ delta1)
    #criterio2 <- abs(teta-teta0)/(se+0.0001)
    criterio2 <- sqrt(sum(teta-teta0)^2)
    if(max(criterio) < delta2){critval <- critval+1}else{critval <- 0}
    if(max(criterio2) < 0.0002){critval2 <- critval2+1}else{critval2 <- 0}
    
    if(count == MaxIter){critval <- 10}
    
    teta0 <- teta
    
    if(dist=="SN")
    {
      teta.out <- rbind(betas,sigma2,lambda)
    }
    if(dist=="ST"|dist=="SSL")
    {
      teta.out <- rbind(betas,sigma2,lambda,nu)
    }
    if(dist=="SCN")
    {
      teta.out <- rbind(betas,sigma2,lambda,nu1,nu2)
    }
    
    Sys.sleep(0.1)
    setTkProgressBar(pb, count, label=paste("Iter ",count,"/",MaxIter,"     -     ",round(count/MaxIter*100,0),"% done",sep = ""))
  }
  
  if(show.convergence=="TRUE")
  {
    pcl <- pc*MaxIter
    p <- ncol(X)
    
    labels = list()
    for(i in 1:p){labels[[i]] <- bquote(beta[.(i)])}
    labels[[p+1]] <- bquote(sigma^2)
    labels[[p+2]] <- bquote(lambda)
    
    if(dist=="SN"){npar <- p+2}
    if(dist=="ST"|dist=="SSL")
    {
      npar <- p+3
      labels[[p+3]] <- bquote(nu)
    }
    if(dist=="SCN")
    {
      npar <- p+4
      labels[[p+3]] <- bquote(nu)
      labels[[p+4]] <- bquote(gamma)
    }
 
    par(mar=c(4, 4.5, 1, 0.5))
    op <- par(mfrow=c(ifelse(npar%%3==0,npar%/%3,(npar%/%3)+1),3))
    
    for(i in 1:npar)
    {
      plot.ts(teta.vec[i,],xlab="Iteration",ylab=labels[[i]])
      abline(v=pcl,lty=2)
    }
    
  }
  
  close(pb)
  
  cat('-----------\n')
  cat('Estimates\n')
  cat('-----------\n')
  cat('\n')
  print(round(teta.out[,1],4))
  cat('\n')
  cat('-----------\n')
  cat('Standard error\n')
  cat('-----------\n')
  cat('\n')
  print(round(se,4))
  
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  out <- list(teta=teta.out, time=time.taken, iter=count, EP=se, criterio=max(criterio,criterio2),EU=EU,EUY=EUY,EUY2=EUY2,EUT=EUT,EUT2=EUT2,EUYT=EUYT)
  return(out)   
}
  
  
  
  
  