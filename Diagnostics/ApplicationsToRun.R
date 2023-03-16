#######################################################################
###                 DIAGNOSTICS - REAL DATA EXAMPLES
#######################################################################

rm(list=ls(all=TRUE))

source("EMSkew-t_New.R")
source("DiagSkew-t_New.R")
source("Envelopes_New.R")
source("Moment_SMSNT.R")
source("Integral_nu_float.R")

#######################################################################
## CASE STUDY I:      Wage rate data
##  See paper Mattos et al. (2018) and Massuia et al (2017) for details 
#######################################################################

library(AER)
data("PSID1976")

y  <- PSID1976$wage
x  <- cbind(1,PSID1976$age,PSID1976$education,PSID1976$youngkids,PSID1976$experience)
cc <- ifelse(PSID1976$participation=="yes", 0, 1)

est_N  <- EM.skewCens(cc, x, y, show.envelope = "FALSE", family = "N")
est_SN <- EM.skewCens(cc, x, y, show.envelope = "FALSE", family = "SN")
est_T  <- EM.skewCens(cc, x, y, show.envelope = "FALSE", family = "T")
est_ST <- EM.skewCens(cc, x, y, show.envelope = "FALSE", family = "ST")

# Case-deletion
CD_N   <- CaseDele(cc, x, y,  est_N$theta,  est_N$E00,  est_N$E01,  est_N$E02,  est_N$E10,  est_N$E20,  est_N$E11, family =  "N")
CD_SN  <- CaseDele(cc, x, y, est_SN$theta, est_SN$E00, est_SN$E01, est_SN$E02, est_SN$E10, est_SN$E20, est_SN$E11, family = "SN")
CD_T   <- CaseDele(cc, x, y,  est_T$theta,  est_T$E00,  est_T$E01,  est_T$E02,  est_T$E10,  est_T$E20,  est_T$E11, family =  "T")
CD_ST  <- CaseDele(cc, x, y, est_ST$theta, est_ST$E00, est_ST$E01, est_ST$E02, est_ST$E10, est_ST$E20, est_ST$E11, family = "ST")

### Benchmarks
p <- ncol(x)
n <- nrow(x)

benchmark_GD_N  <- 2*(p+1)/n 
benchmark_GD_SN <- 2*(p+2)/n 
benchmark_GD_T  <- 2*(p+1)/n 
benchmark_GD_ST <- 2*(p+2)/n 

# GD Plots
par(mfrow=c(2,2))
plot(CD_N$GD, ylim = c(0, 1.4), main="N", xlab="Index", ylab="GD", type="h")
abline(h = benchmark_GD_N, lty = 3)
plot(CD_SN$GD, ylim = c(0, 0.1), main="SN", xlab="Index", ylab="GD", type="h")
abline(h = benchmark_GD_SN, lty = 3)
plot(CD_T$GD, ylim = c(0, 0.1), main="T", xlab="Index", ylab="GD", type="h")
abline(h = benchmark_GD_T, lty = 3)
plot(CD_ST$GD, ylim = c(0, 0.1), main="ST", xlab="Index", ylab="GD", type="h")
abline(h = benchmark_GD_ST, lty = 3)

##  Mahalanobis versus Weights
maha_T  <- (est_T$yhat -x%*%est_T$theta[1:p])^2/est_T$theta[p+1]
maha_ST <- (est_ST$yhat-x%*%est_ST$theta[1:p])^2/est_ST$theta[p+1]

### Figure. Mahalanobis versus Weights
par(mfrow=c(1,2))
plot(maha_T[cc==0], est_T$E00[cc==0], xlim = c(min(maha_T),max(maha_T)), main="T", xlab = "Mahalanobis distance", ylab ="Weights", pch=3)
points(x = maha_T[cc==1], y = est_T$E00[cc==1], pch = 1, col = "red")
abline(h = 1)
plot(maha_ST[cc==0], est_ST$E00[cc==0], xlim = c(min(maha_ST),max(maha_ST)), main="ST", xlab = "Mahalanobis distance", ylab ="Weights", pch=3)
points(x = maha_ST[cc==1], y = est_ST$E00[cc==1], pch = 1, col = "red")
abline(h = 1)

## Pertubartion 1 - Case-weight
If_N  <- If(cc,x,y,est_N$theta ,est_N$E00 ,est_N$E01 ,est_N$E02 ,est_N$E10 ,est_N$E20 ,est_N$E11 , Perturb = 1, family = "N")
If_SN <- If(cc,x,y,est_SN$theta,est_SN$E00,est_SN$E01,est_SN$E02,est_SN$E10,est_SN$E20,est_SN$E11, Perturb = 1, family = "SN")
If_T  <- If(cc,x,y,est_T$theta ,est_T$E00 ,est_T$E01 ,est_T$E02 ,est_T$E10 ,est_T$E20 ,est_T$E11 , Perturb = 1, family = "T")
If_ST <- If(cc,x,y,est_ST$theta,est_ST$E00,est_ST$E01,est_ST$E02,est_ST$E10,est_ST$E20,est_ST$E11, Perturb = 1, family = "ST")

## Pertubartion 2 - Scale perturbation
If_N2  <- If(cc,x,y,est_N$theta ,est_N$E00 ,est_N$E01 ,est_N$E02 ,est_N$E10 ,est_N$E20 ,est_N$E11 , Perturb = 2, family = "N")
If_SN2 <- If(cc,x,y,est_SN$theta,est_SN$E00,est_SN$E01,est_SN$E02,est_SN$E10,est_SN$E20,est_SN$E11, Perturb = 2, family = "SN")
If_T2  <- If(cc,x,y,est_T$theta ,est_T$E00 ,est_T$E01 ,est_T$E02 ,est_T$E10 ,est_T$E20 ,est_T$E11 , Perturb = 2, family = "T")
If_ST2 <- If(cc,x,y,est_ST$theta,est_ST$E00,est_ST$E01,est_ST$E02,est_ST$E10,est_ST$E20,est_ST$E11, Perturb = 2, family = "ST")

### Figure. Case-weight and Scale perturbation
par(mfrow=c(4,2))
plot(If_N$Medida, ylim = c(0, 0.20), main ="N", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_N$Benchmark, lty = 3)
plot(If_N2$Medida, ylim = c(0, 0.20), main ="N", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_N2$Benchmark, lty = 3)
plot(If_SN$Medida, ylim = c(0, 0.10), main ="SN", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_SN$Benchmark, lty = 3)
plot(If_SN2$Medida, ylim = c(0, 0.10), main ="SN", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_SN2$Benchmark, lty = 3)
plot(If_T$Medida, ylim = c(0, 0.04), main ="T", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_T$Benchmark, lty = 3)
plot(If_T2$Medida, ylim = c(0, 0.04), main ="T", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_T2$Benchmark, lty = 3)
plot(If_ST$Medida, ylim = c(0, 0.04), main ="ST", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_ST$Benchmark, lty = 3)
plot(If_ST2$Medida, ylim = c(0, 0.04), main ="ST", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_ST2$Benchmark, lty = 3)

#######################################################################
## CASE STUDY II:    Stellar abundances data
## See paper Mattos et al. (2018) for details
#######################################################################

if(!require(astrodatR)) install.packages("astrodatR");  library(astrodatR)

data(censor_Be)
dados <- censor_Be
y     <- dados$logN_Be
cc    <- ifelse(dados$Ind_Be==1, 0, 1) # cc == 1 represents the censored!
x1    <- dados$Teff/1000
x     <- cbind(1, x1)

est_N  <- EM.skewCens(cc, x, y, show.envelope = "FALSE", family = "N")
est_SN <- EM.skewCens(cc, x, y, show.envelope = "FALSE", family = "SN")
est_T  <- EM.skewCens(cc, x, y, show.envelope = "FALSE", family = "T")
est_ST <- EM.skewCens(cc, x, y, show.envelope = "FALSE", family = "ST")

# Case-deletion
CD_N   <- CaseDele(cc, x, y,  est_N$theta,  est_N$E00,  est_N$E01,  est_N$E02,  est_N$E10,  est_N$E20,  est_N$E11, family =  "N")
CD_SN  <- CaseDele(cc, x, y, est_SN$theta, est_SN$E00, est_SN$E01, est_SN$E02, est_SN$E10, est_SN$E20, est_SN$E11, family = "SN")
CD_T   <- CaseDele(cc, x, y,  est_T$theta,  est_T$E00,  est_T$E01,  est_T$E02,  est_T$E10,  est_T$E20,  est_T$E11, family =  "T")
CD_ST  <- CaseDele(cc, x, y, est_ST$theta, est_ST$E00, est_ST$E01, est_ST$E02, est_ST$E10, est_ST$E20, est_ST$E11, family = "ST")

### Benchmarks
p <- ncol(x)
n <- nrow(x)

benchmark_GD_N  <- 2*(p+1)/n 
benchmark_GD_SN <- 2*(p+2)/n 
benchmark_GD_T  <- 2*(p+1)/n 
benchmark_GD_ST <- 2*(p+2)/n 

round(c(min(CD_N$GD),max(CD_N$GD)),4)
round(c(min(CD_SN$GD),max(CD_SN$GD)),4)
round(c(min(CD_T$GD),max(CD_T$GD)),4)
round(c(min(CD_ST$GD),max(CD_ST$GD)),4)

# GD Plots
par(mfrow=c(2,2))
plot(CD_N$GD, ylim = c(0,0.12), main="N", xlab="Index", ylab="GD", type="h")
abline(h = benchmark_GD_N, lty = 3)
plot(CD_SN$GD, ylim = c(0, 0.12), main="SN", xlab="Index", ylab="GD", type="h")
abline(h = benchmark_GD_SN, lty = 3)
plot(CD_T$GD, ylim = c(0, 0.12), main="T", xlab="Index", ylab="GD", type="h")
abline(h = benchmark_GD_T, lty = 3)
plot(CD_ST$GD, ylim = c(0, 0.12), main="ST", xlab="Index", ylab="GD", type="h")
abline(h = benchmark_GD_ST, lty = 3)

##  Mahalanobis versus Weights
maha_T  <- (est_T$yhat -x%*%est_T$theta[1:p])^2/est_T$theta[p+1]
maha_ST <- (est_ST$yhat-x%*%est_ST$theta[1:p])^2/est_ST$theta[p+1]

### Figure. Mahalanobis versus Weights
par(mfrow=c(1,2))
plot(maha_T[cc==0], est_T$E00[cc==0], xlim = c(min(maha_T),max(maha_T)), ylim=c(min(est_T$E00), max(est_T$E00)), main="T", xlab = "Mahalanobis distance", ylab ="Weights", pch = 3)
points(x = maha_T[cc==1], y = est_T$E00[cc==1], pch = 1, col = "red")
abline(h = 1)
plot(maha_ST[cc==0], est_ST$E00[cc==0], xlim = c(min(maha_ST),max(maha_ST)), ylim=c(min(est_ST$E00), max(est_T$E00)), main="ST", xlab = "Mahalanobis distance", ylab ="Weights", pch = 3)
points(x = maha_ST[cc==1], y = est_ST$E00[cc==1], pch = 1, col = "red")
abline(h = 1)

## Pertubartion 1 - Case-weight
If_N  <- If(cc,x,y,est_N$theta ,est_N$E00 ,est_N$E01 ,est_N$E02 ,est_N$E10 ,est_N$E20 ,est_N$E11 , Perturb = 1, family = "N")
If_SN <- If(cc,x,y,est_SN$theta,est_SN$E00,est_SN$E01,est_SN$E02,est_SN$E10,est_SN$E20,est_SN$E11, Perturb = 1, family = "SN")
If_T  <- If(cc,x,y,est_T$theta ,est_T$E00 ,est_T$E01 ,est_T$E02 ,est_T$E10 ,est_T$E20 ,est_T$E11 , Perturb = 1, family = "T")
If_ST <- If(cc,x,y,est_ST$theta,est_ST$E00,est_ST$E01,est_ST$E02,est_ST$E10,est_ST$E20,est_ST$E11, Perturb = 1, family = "ST")

## Pertubartion 2 - Scale perturbation
If_N2  <- If(cc,x,y,est_N$theta ,est_N$E00 ,est_N$E01 ,est_N$E02 ,est_N$E10 ,est_N$E20 ,est_N$E11 , Perturb = 2, family = "N")
If_SN2 <- If(cc,x,y,est_SN$theta,est_SN$E00,est_SN$E01,est_SN$E02,est_SN$E10,est_SN$E20,est_SN$E11, Perturb = 2, family = "SN")
If_T2  <- If(cc,x,y,est_T$theta ,est_T$E00 ,est_T$E01 ,est_T$E02 ,est_T$E10 ,est_T$E20 ,est_T$E11 , Perturb = 2, family = "T")
If_ST2 <- If(cc,x,y,est_ST$theta,est_ST$E00,est_ST$E01,est_ST$E02,est_ST$E10,est_ST$E20,est_ST$E11, Perturb = 2, family = "ST")

### Figure. Case-weight and Scale perturbation
par(mfrow=c(4,2))
plot(If_N$Medida, ylim = c(0, 0.60), main ="N", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_N$Benchmark, lty = 3)
plot(If_N2$Medida, ylim = c(0, 0.60), main ="N", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_N2$Benchmark, lty = 3)
plot(If_SN$Medida, ylim = c(0, 0.70), main ="SN", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_SN$Benchmark, lty = 3)
plot(If_SN2$Medida, ylim = c(0, 0.30), main ="SN", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_SN2$Benchmark, lty = 3)
plot(If_T$Medida, ylim = c(0, 0.20), main ="T", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_T$Benchmark, lty = 3)
plot(If_T2$Medida, ylim = c(0, 0.20), main ="T", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_T2$Benchmark, lty = 3)
plot(If_ST$Medida, ylim = c(0, 0.30), main ="ST", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_ST$Benchmark, lty = 3)
plot(If_ST2$Medida, ylim = c(0, 0.20), main ="ST", xlab = "Index", ylab="M(0)", type="h")
abline(h = If_ST2$Benchmark, lty = 3)
