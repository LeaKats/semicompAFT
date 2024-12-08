rm(list = ls())
a<-Sys.time()
print(a)

#load packages
library(SemiCompRisks)
library(Formula)
library(readxl)

#read data
data <- as.data.frame(read_excel("breast.xlsx"))
data<-data[data$nodes>0,]

#dealing with cases when recurring time is equal to death time.
data$dtime[data$rtime==data$dtime&data$recur==1&data$death==1]<-data$dtime[data$rtime==data$dtime&data$recur==1&data$death==1]+0.5
data$dtime[data$rtime==data$dtime&data$recur==1&data$death==0]<-data$dtime[data$rtime==data$dtime&data$recur==1&data$death==0]+0.5

#define variables
data$age_10<-data$age/10
data$size20_50<-ifelse(data$size=="20-50",1,0)
data$size50<-ifelse(data$size==">50",1,0)
data$log_nodes<-log(data$nodes)
data$log_pgr<-log(1+data$pgr)
data$log_er<-log(1+data$er)
data$grade3<-ifelse(data$grade==3,1,0)
data$age_since_relapse_10<-data$age_10+data$rtime/10

#define the times used in the format of package SemiCompRisks (in years)
data$y1L<-data$y1U<-data$rtime/365.25
data$y1U[which(data$recur == 0)] <- Inf
data$y2L<-data$y2U<-data$dtime/365.25
data$y2U[which(data$death == 0)] <- Inf
data$LT <- rep(0, dim(data)[1])


#define the formula that includes all transitions
form <- Formula(LT | y1L + y1U | y2L + y2U  ~ 
                  age_10+meno+size20_50+size50+log_nodes+log_er+log_pgr+hormon+chemo+grade3 | 
                  age_10+meno+size20_50+size50+log_nodes+log_er+log_pgr+hormon+chemo+grade3 | 
                  age_since_relapse_10+meno+size20_50+size50+log_nodes+log_er+log_pgr+hormon+chemo+grade3)

#####################
## Define the Hyperparameters ##
#####################

## Subject-specific random effects variance component
##
theta.ab <- c(0.5, 0.05)

## log-Normal model
##
LN.ab1 <- c(0.3, 0.3)
LN.ab2 <- c(0.3, 0.3)
LN.ab3 <- c(0.3, 0.3)

## DPM model
##
DPM.mu1 <- log(12)
DPM.mu2 <- log(12)
DPM.mu3 <- log(12)

DPM.sigSq1 <- 100
DPM.sigSq2 <- 100
DPM.sigSq3 <- 100

DPM.ab1 <-  c(2, 1)
DPM.ab2 <-  c(2, 1)
DPM.ab3 <-  c(2, 1)

Tau.ab1 <- c(1.5, 0.0125)
Tau.ab2 <- c(1.5, 0.0125)
Tau.ab3 <- c(1.5, 0.0125)

##
hyperParams <- list(theta=theta.ab,
                    LN=list(LN.ab1=LN.ab1, LN.ab2=LN.ab2, LN.ab3=LN.ab3),
                    DPM=list(DPM.mu1=DPM.mu1, DPM.mu2=DPM.mu2, DPM.mu3=DPM.mu3, DPM.sigSq1=DPM.sigSq1,
                             DPM.sigSq2=DPM.sigSq2, DPM.sigSq3=DPM.sigSq3, DPM.ab1=DPM.ab1, DPM.ab2=DPM.ab2,
                             DPM.ab3=DPM.ab3, Tau.ab1=Tau.ab1, Tau.ab2=Tau.ab2, Tau.ab3=Tau.ab3))

###################
## Define the MCMC SETTINGS ##
###################

## Setting for the overall run
##
numReps    <- 1000000
thin       <- 1000
burninPerc <- 0.5

## Setting for storage
##
nGam_save <- 10
nY1_save <- 10
nY2_save <- 10
nY1.NA_save <- 10

## Tuning parameters for specific updates
##
##  - those common to all models
betag.prop.var	<- c(0.01,0.01,0.01)
mug.prop.var	<- c(0.1,0.1,0.1)
zetag.prop.var	<- c(0.1,0.1,0.1)
gamma.prop.var	<- 0.01

##
mcmcParams	<- list(run=list(numReps=numReps, thin=thin, burninPerc=burninPerc),
                   storage=list(nGam_save=nGam_save, nY1_save=nY1_save, nY2_save=nY2_save, nY1.NA_save=nY1.NA_save),
                   tuning=list(betag.prop.var=betag.prop.var, mug.prop.var=mug.prop.var,
                               zetag.prop.var=zetag.prop.var, gamma.prop.var=gamma.prop.var))

#################################################################
## Analysis ############
#################################################################

#########
## DPM (Semiparametric Dirichlet process mixture model)##
#########

##
myModel <- "DPM"
myPath  <- paste0("Output/02-Results-DPM/",numReps,"/")

#initiate start values
startValues      <- initiate.startValues_AFT(form, data, model=myModel, nChain=3)

#fit
fit_DPM <- BayesID_AFT(form, data, model=myModel, hyperParams,
                       startValues, mcmcParams, path=myPath)

#summary
summ.fit_DPM <- summary(fit_DPM); names(summ.fit_DPM)
summ.fit_DPM

write.csv(summ.fit_DPM$theta,file=paste0(myPath,r,"_","theta_DPM.csv"))
write.csv(summ.fit_DPM$coef,file=paste0(myPath,r,"_","coef_DPM.csv"))
write.csv(summ.fit_DPM$psrf,file=paste0(myPath,r,"_","psrf_DPM.csv"))
write.csv(summ.fit_DPM$h0,file=paste0(myPath,r,"_","h0_DPM.csv"))
write.csv(summ.fit_DPM$setup[6:11],file=paste0(myPath,r,"_","setup_DPM.csv"))
