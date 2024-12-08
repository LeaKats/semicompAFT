rm(list = ls())
a<-Sys.time()
print(a)

#Load packages
suppressPackageStartupMessages(library(SemiCompRisks))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(readxl))

###### Load data ######
data <- as.data.frame(read_excel("breast.xlsx"))
data<-data[data$nodes>0,]

# data$V<-pmin(data$rtime,data$dtime)
# data$delta1<-data$recur
# data$W<-data$delta1*data$dtime
# data$delta2<-ifelse(data$recur==0&data$death==1,1,0)
# data$delta3<-ifelse(data$recur==1&data$death==1,1,0)
data$size20_50<-ifelse(data$size=="20-50",1,0)
data$size50<-ifelse(data$size==">50",1,0)
data$grade3<-ifelse(data$grade==3,1,0)
data$log_pgr<-log(1+data$pgr)
data$log_er<-log(1+data$er)
data$age_10<-data$age/10
data$log_nodes<-log(data$nodes)

n<-dim(data)[1]


data$rtime_years<-data$rtime/365.25
data$rtime_10<-data$rtime/10
data$rtime_years_10<-data$rtime/365.25/10
data$dtime_10<-data$dtime/10
data$dtime_years<-data$dtime/365.25
#dealing with cases when recur time=death time (added 0.5 day):
data$dtime_years[(data$recur==1)&(data$rtime==data$dtime)]<-data$dtime_years[(data$recur==1)&(data$rtime==data$dtime)]+0.5/365.25
data$age_at_relapse_10<-data$age_10+data$rtime_years_10

set.seed(1)

form <- Formula(rtime + recur | dtime + death ~ age_10+log_nodes+log_er+log_pgr+meno+size20_50+size50+hormon+chemo+grade3 | 
                  age_10+log_nodes+log_er+log_pgr+meno+size20_50+size50+hormon+chemo+grade3 | 
                  age_at_relapse_10+log_nodes+log_er+log_pgr+meno+size20_50+size50+hormon+chemo+grade3)

#####################
## Hyperparameters ##
#####################

## Subject-specific frailty variance component
##  - prior parameters for 1/theta
##
## Subject-specific frailty variance component
##  - prior parameters for 1/theta
##
theta.ab <- c(0.7, 0.7)

## Weibull baseline hazard function: alphas, kappas
##
WB.ab1 <- c(0.5, 0.01) # prior parameters for alpha1
WB.ab2 <- c(0.5, 0.01) # prior parameters for alpha2
WB.ab3 <- c(0.5, 0.01) # prior parameters for alpha3
##
WB.cd1 <- c(0.5, 0.05) # prior parameters for kappa1
WB.cd2 <- c(0.5, 0.05) # prior parameters for kappa2
WB.cd3 <- c(0.5, 0.05) # prior parameters for kappa3

## PEM baseline hazard function
##
PEM.ab1 <- c(0.7, 0.7) # prior parameters for 1/sigma_1^2
PEM.ab2 <- c(0.7, 0.7) # prior parameters for 1/sigma_2^2
PEM.ab3 <- c(0.7, 0.7) # prior parameters for 1/sigma_3^2
##
PEM.alpha1 <- 10 # prior parameters for K1
PEM.alpha2 <- 10 # prior parameters for K2
PEM.alpha3 <- 10 # prior parameters for K3

## MVN cluster-specific random effects
##
Psi_v <- diag(1, 3)
rho_v <- 100

## DPM cluster-specific random effects
##
Psi0  <- diag(1, 3)
rho0  <- 10
aTau  <- 1.5
bTau  <- 0.0125

##
hyperParams <- list(theta=theta.ab,
                    WB=list(WB.ab1=WB.ab1, WB.ab2=WB.ab2, WB.ab3=WB.ab3,
                            WB.cd1=WB.cd1, WB.cd2=WB.cd2, WB.cd3=WB.cd3),
                    PEM=list(PEM.ab1=PEM.ab1, PEM.ab2=PEM.ab2, PEM.ab3=PEM.ab3,
                             PEM.alpha1=PEM.alpha1, PEM.alpha2=PEM.alpha2, PEM.alpha3=PEM.alpha3),
                    MVN=list(Psi_v=Psi_v, rho_v=rho_v),
                    DPM=list(Psi0=Psi0, rho0=rho0, aTau=aTau, bTau=bTau))

###################
## MCMC SETTINGS ##
###################

## Setting for the overall run
##
numReps    <- 20000000
thin       <- 10000
burninPerc <- 0.5

## Settings for storage
##
nGam_save <- 0
storeV    <- rep(TRUE, 3)

## Tuning parameters for specific updates
##
##  - those common to all models
mhProp_theta_var  <- 0.05
mhProp_Vg_var     <- c(0.05, 0.05, 0.05)
##
## - those specific to the Weibull specification of the baseline hazard functions
mhProp_alphag_var <- c(0.01, 0.01, 0.01)
##
## - those specific to the PEM specification of the baseline hazard functions
Cg        <- c(0.2, 0.2, 0.2)
delPertg  <- c(0.5, 0.5, 0.5)
rj.scheme <- 1
Kg_max    <- c(50, 50, 50)
sg_max    <- c(max(data$rtime[data$recur == 1]),
               max(data$dtime[data$recur == 0 & data$death == 1]),
               max(data$dtime[data$recur == 1 & data$death == 1]))

time_lambda1 <- seq(1, sg_max[1], 1)
time_lambda2 <- seq(1, sg_max[2], 1)
time_lambda3 <- seq(1, sg_max[3], 1)               

##
mcmc.WB  <- list(run=list(numReps=numReps, thin=thin, burninPerc=burninPerc),
                 storage=list(nGam_save=nGam_save, storeV=storeV),
                 tuning=list(mhProp_theta_var=mhProp_theta_var,
                             mhProp_Vg_var=mhProp_Vg_var, mhProp_alphag_var=mhProp_alphag_var))

##
mcmc.PEM <- list(run=list(numReps=numReps, thin=thin, burninPerc=burninPerc),
                 storage=list(nGam_save=nGam_save, storeV=storeV),
                 tuning=list(mhProp_theta_var=mhProp_theta_var,
                             mhProp_Vg_var=mhProp_Vg_var, Cg=Cg, delPertg=delPertg,
                             rj.scheme=rj.scheme, Kg_max=Kg_max,
                             time_lambda1=time_lambda1, time_lambda2=time_lambda2,
                             time_lambda3=time_lambda3))

#####################
## Starting Values ##
#####################

##
Sigma_V <- diag(0.1, 3)
Sigma_V[1,2] <- Sigma_V[2,1] <- -0.05
Sigma_V[1,3] <- Sigma_V[3,1] <- -0.06
Sigma_V[2,3] <- Sigma_V[3,2] <- 0.07


#############
## Data analysis: PEM-DPM ##
#PEM: non-parametric mixture of piecewise exponential models 
#DPM: non-parametric Dirichlet process mixture of multivariate normals 
#############

##
myModel <- c("Markov", "PEM")
myPath  <- "Output/02-Results-PEM/"

startValues      <- initiate.startValues_HReg(form, data, model=myModel, nChain=2)

##
fit_PEM <- BayesID_HReg(form, data, id=NULL, model=myModel,
                        hyperParams, startValues, mcmc.PEM, path=myPath)

fit_PEM
summ.fit_PEM <- summary(fit_PEM); names(summ.fit_PEM)
summ.fit_PEM

write.csv(coefficients(fit_PEM),file="lee2015_10M.csv")
