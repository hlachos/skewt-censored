rm(list=ls(all=TRUE))

source("EMSkew-t.R")
source("Moment_SMSNT.R")
source("Integral_nu_float.R")
source("envelopes.R")

################################################################################
### Simulated Data
################################################################################
n<-200
x<-cbind(1, runif(n,1,5), runif(n,0,1))
y<-matrix(0,n,1)
betas<-c(1,2,3)
sigma2 <- 1
shape<- -2
nu<- 4
perc<- 0.2 ## censoring level


sdata<-generate_SMSNCR(x,betas,sigma2,shape,n,"left",perc,"ST",nu)
y<-sdata$yc
cc<-sdata$cc


################################################################################
## EM algorithm with OLS starting values get.init =TRUE (see Subsection 3.1)
################################################################################

est_ST.EM <- EM.skewCens(cc, x, y, get.init = TRUE, show.envelope="FALSE", error = 0.00001, iter.max = 300, family="ST")
est_SN.EM <- EM.skewCens(cc, x, y, get.init = TRUE,  show.envelope="FALSE", error = 0.00001, iter.max = 300, family="SN")
est_T.EM <- EM.skewCens(cc, x, y, get.init = TRUE,  show.envelope="FALSE", error = 0.00001, iter.max = 300, family="T")
est_N.EM <- EM.skewCens(cc, x, y, get.init = TRUE,  show.envelope="FALSE", error = 0.00001, iter.max = 300, family="N")



################################################################################
###  Random Starting values in the EM algorithm for the ST case
################################################################################


est_ST.EM.Random <- EM.skewCens(cc, x, y, beta = rnorm(3,0,5), sigma2 = rgamma(1,0.1,0.1), shape = rnorm(0,10),  
                                nu=3, show.envelope="FALSE", error = 0.00001, iter.max = 300, family="ST")

## Comparision between different starting values (OLS and Random)

      df.comp <- data.frame(
        ThetaTrue    = c(betas,sigma2,shape,nu),
        est_ST.EM.OLS  = est_ST.EM$theta,
        est_ST.EM.Random  = est_ST.EM.Random$theta
         )
    
    df.comp
    

################################################################################
## Comparing SAEM and EM for simulated data (time)
## SAEM proposed by Mattos et al. (2018)
################################################################################
                                                        
source("./SMSN_CR Code/auxiliary_functions_SMSNCR.R")
source("./SMSN_CR Code/SAEM_estimation.R")

start_time <- Sys.time()
est_ST.EM <-  EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 400, family="ST")
end_time <- Sys.time()
T.EM<-difftime(end_time, start_time, units = "min")

start_time <- Sys.time()
est_ST.SAEM <- SAEM_EST(y,x,cc,cens="left",LS=NULL,precisao=0.00001,MaxIter=400, M=20, pc=0.3, dist="ST",nu.fixed=FALSE,nu=3,show.convergence="FALSE")
end_time <- Sys.time()
T.SAEM<-difftime(end_time, start_time, units = "min")


################################################################################
###                     REAL DATA EXAMPLES
################################################################################


################################################################################
## CASE STUDY I:      Wage rate data
#  See paper Mattos et al. (2018) and Massuia et al (2017) for details 
################################################################################

if(!require(SMNCensReg)) install.packages("SMNCensReg");  library(SMNCensReg)
## library   SMNCensReg firs a CR model with symmetric distributions N, T

data(wage.rates)
y <- wage.rates$wage
x <- cbind(1,wage.rates$age,wage.rates$educ,wage.rates$kidslt6,wage.rates$kidsge6)
cc<- c(rep(0,428),rep(1,325))


est_N <-EM.skewCens(cc, x, y, get.init = TRUE, show.envelope="FALSE", error = 0.00001, iter.max = 300, family="N")
est_T <-EM.skewCens(cc, x, y, get.init = TRUE, show.envelope="FALSE", error = 0.00001, iter.max = 300, family="T")
est_SN <-EM.skewCens(cc, x, y, get.init = TRUE, show.envelope="FALSE", error = 0.00001, iter.max = 300, family="SN")
est_ST <-EM.skewCens(cc, x, y, get.init = TRUE, show.envelope="FALSE", error = 0.00001, iter.max = 300, family="ST")



## library SMNCensReg  : symmetric distribution only

fitT <- CensReg.SMN(cc,x,y,nu=3,cens="left",dist="T",show.envelope="TRUE", error=0.000001 ,iter.max=300)
fitN <- CensReg.SMN(cc,x,y,cens="left",dist="Normal",show.envelope="TRUE")


################################################################################
## CASE STUDY II:    Stellar abundances data
## See paper Mattos et al. (2018) for details
################################################################################

if(!require(astrodatR)) install.packages("astrodatR");  library(astrodatR)


data(censor_Be)
dados <- censor_Be
y <- dados[,5]
cc <- dados[,4]
x1 <- dados[,3]/1000
x <- cbind(1,x1)


est_N <-EM.skewCens(cc, x, y, get.init = TRUE, show.envelope="FALSE", error = 0.00001, iter.max = 300, family="N")
est_T <-EM.skewCens(cc, x, y, get.init = TRUE,show.envelope="FALSE", error = 0.00001, iter.max = 300, family="T")
est_ST <-EM.skewCens(cc, x, y, get.init = TRUE,show.envelope="FALSE", error = 0.0001, iter.max = 300, family="ST")
est_SN <-EM.skewCens(cc, x, y, get.init = TRUE, show.envelope="FALSE", error = 0.00001, iter.max = 300, family="SN")

### library   SMNCensReg: symmetric distributions only

fitT <- CensReg.SMN(cc,x,y,nu=3,cens="left",dist="T",show.envelope="TRUE", error=0.000001 ,iter.max=300)
fitN <- CensReg.SMN(cc,x,y,cens="left",dist="Normal",show.envelope="TRUE")

################################################################################
## CASE STUDY III: 
################################################################################

data(wage.rates)
y <- wage.rates$hours/1000
lattice:: densityplot(y)
x <- cbind(1,wage.rates$educ,wage.rates$age,wage.rates$exper,wage.rates$expersq)
#censure
cc <- c(rep(0,428),rep(1,325))
n <- length(y)


est_N <-EM.skewCens(cc, x, y, get.init = TRUE,show.envelope="FALSE", error = 0.00001, iter.max = 300, family="N")
est_T <-EM.skewCens(cc, x, y, get.init = TRUE,show.envelope="FALSE", error = 0.00001, iter.max = 300, family="T")
est_ST <-EM.skewCens(cc, x, y, get.init = TRUE,show.envelope="FALSE", error = 0.00001, iter.max = 300, family="ST")

### library   SMNCensReg: symmetric distributions only

fitT <- CensReg.SMN(cc,x,y,nu=3,cens="left",dist="T",show.envelope="TRUE", error=0.000001 ,iter.max=300)
fitN <- CensReg.SMN(cc,x,y,cens="left",dist="Normal",show.envelope="TRUE")

