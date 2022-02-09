rm(list=ls(all=TRUE))
folderSNI<-"C:/Users/vid09002/Dropbox (Uconn)/PaperMixCensurasST/ProgramsSTnu/FinalCodes-EM"
setwd(folderSNI)
source("EMSkew-t.R")
source("Moment_SMSNT.R")
source("Integral_nu_float.R")
source("Envelopes.R")

################################################################################
### Simulated Data
################################################################################
library(mixsmsn)
n<-1000
x<-cbind(1, runif(n,1,5), runif(n,0,1))
y<-matrix(0,n,1)
betas1<-c(1,2,3)
Sigma1 <- 1
mu1 <- x%*%betas1
shape1<- -2
nu<- 4
perc<- 0.2 ## censoring level

for (i in 1:n){
arg1 = list(mu=mu1[i], sigma2=Sigma1, shape=shape1, nu=nu)
y[i]= rmix(1, 1, "Skew.t", list(arg1))         # change to ("t", "Skew.t", "Skew.cn", "Skew.slash",
#"Skew.normal", "Normal" from the mixsmsns package
}

aa=sort(y,decreasing=FALSE)
cutof<-aa[ceiling(perc*n)]
cc=matrix(1,n,1)*(y<=cutof)
y[cc==1]=cutof


est_ST.EM <- EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 300, family="ST")
est_SN.EM <- EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 300, family="SN")
est_T.EM <- EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 300, family="T")
est_N.EM <- EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 300, family="N")

par(mfrow = c(2, 2))
EnvelopeRMT(est_ST.EM$theta,y,x,cc,family="ST")
EnvelopeRMT(c(est_SN.EM$theta),y,x,cc,family="SN")
EnvelopeRMT(c(est_T.EM$theta),y,x,cc,family="ST")
EnvelopeRMT(c(est_N.EM$theta,0),y,x,cc,family="SN")


################################################################################
## Comparing SAEM and EM 
## SAEM proposed by Mattos et al. (2018)
################################################################################
                                                        

folder<-"C:/Users/vid09002/Dropbox (Uconn)/PaperMixCensurasST/ProgramsSTnu/SMSN_CR Code"
setwd(folder)
source("auxiliary_functions_SMSNCR.R")
source("SAEM_estimation.R")

start_time <- Sys.time()
est_ST.EM <-  EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 300, family="ST"
end_time <- Sys.time()
T.EM<-difftime(end_time, start_time, units = "min")

start_time <- Sys.time()
est_ST.SAEM <- SAEM_EST(y,x,cc,cens="left",LS=NULL,precisao=0.000001,MaxIter=400, M=20, pc=0.3, dist="ST",nu.fixed=FALSE,nu=3,show.convergence="FALSE")
end_time <- Sys.time()
T.SAEM<-difftime(end_time, start_time, units = "min")

################################################################################
## CASE STUDY I 
################################################################################

library(SMNCensReg)
data(wage.rates)
y <- wage.rates$wage
x <- cbind(wage.rates$age,wage.rates$educ,wage.rates$kidslt6,wage.rates$kidsge6)
cc<- c(rep(0,428),rep(1,325))
##Fits a left censored Student-t model to the data

est_N <-EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 300, family="N")
EnvelopeRMT(c(est_N$theta,0),y,x,cc,family="SN")
est_T <-EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 300, family="T")
EnvelopeRMT(est_T$theta,y,x,cc,family="ST")
est_ST <-EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 300, family="ST")
EnvelopeRMT(est_ST$theta,y,x,cc,family="ST")


## library SMNCensReg  : symmetric distribution only

fitT <- CensReg.SMN(cc,x,y,nu=3,cens="left",dist="T",show.envelope="TRUE", error=0.000001 ,iter.max=300)
fitN <- CensReg.SMN(cc,x,y,cens="left",dist="Normal",show.envelope="TRUE")


################################################################################
## CASE STUDY II 
################################################################################

library(astrodatR)

data(censor_Be)
dados <- censor_Be
y <- dados[,5]
cc <- dados[,4]
x1 <- dados[,3]/1000
x <- cbind(1,x1)


est_N <-EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 300, family="N")
EnvelopeRMT(c(est_N$theta,0),y,x,cc,family="SN")
est_T <-EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 300, family="T")
EnvelopeRMT(est_T$theta,y,x,cc,family="ST")
est_ST <-EM.skewCens(cc, x, y, get.init = TRUE, error = 0.0001, iter.max = 300, family="ST")
EnvelopeRMT(est_ST$theta,y,x,cc,family="ST")
est_SN <-EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 300, family="SN")
EnvelopeRMT(est_SN$theta,y,x,cc,family="SN")

### symmetric distribution only

fitT <- CensReg.SMN(cc,x,y,nu=3,cens="left",dist="T",show.envelope="TRUE", error=0.000001 ,iter.max=300)
fitN <- CensReg.SMN(cc,x,y,cens="left",dist="Normal",show.envelope="TRUE")

################################################################################
## CASE STUDY III
################################################################################

library(SMNCensReg)
data(wage.rates)
y <- wage.rates$hours/1000
lattice:: densityplot(y)
x <- cbind(1,wage.rates$educ,wage.rates$age,wage.rates$exper,wage.rates$expersq)
#censure
cc <- c(rep(0,428),rep(1,325))
n <- length(y)


est_N <-EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 300, family="N")
EnvelopeRMT(c(est_N$theta,0),y,x,cc,family="SN")
est_T <-EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 300, family="T")
EnvelopeRMT(est_T$theta,y,x,cc,family="ST")
est_ST <-EM.skewCens(cc, x, y, get.init = TRUE, error = 0.00001, iter.max = 300, family="ST")
EnvelopeRMT(est_ST$theta,y,x,cc,family="ST")


fitT <- CensReg.SMN(cc,x,y,nu=3,cens="left",dist="T",show.envelope="TRUE", error=0.000001 ,iter.max=300)
fitN <- CensReg.SMN(cc,x,y,cens="left",dist="Normal",show.envelope="TRUE")

