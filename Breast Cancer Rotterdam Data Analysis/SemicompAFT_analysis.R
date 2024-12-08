rm(list = ls())

##### set working directory
setwd("H:/My Drive/Education/Phd/biom_paper/R/scripts_biometrics")
##### Load the SemicompAFT package; it takes time to load
source("./SemicompAFT/SemicompAFT.R")


###### Load data ######
data <- as.data.frame(read_excel("breast.xlsx"))
data<-data[data$nodes>0,]

###### Prepare the data
#define variables
data$age_10<-data$age/10
data$log_nodes<-log(data$nodes)
data$log_er<-log(1+data$er)
data$log_pgr<-log(1+data$pgr)
data$size20_50<-ifelse(data$size=="20-50",1,0)
data$size50<-ifelse(data$size==">50",1,0)
data$grade3<-ifelse(data$grade==3,1,0)

#define the observed times and indicators
V<-data$V<-pmin(data$rtime,data$dtime)/365.25
delta1<-data$delta1<-data$recur
W<-data$W<-data$delta1*data$dtime/365.25
delta2<-data$delta2<-ifelse(data$recur==0&data$death==1,1,0)
delta3<-data$delta3<-ifelse(data$recur==1&data$death==1,1,0)
#dealing with cases when W=V (added 0.5 day):
W<-data$W<-ifelse(data$W==data$V,data$W+0.5/365.25,data$W)
data$delta1<-data$delta1==1
data$delta2<-data$delta2==1
data$delta3<-data$delta3==1

#define variable age at relapse
data$age_since_relapse_10<-data$age_10 + data$V/10

###### Set Parameters ######
## initial values
initial_sigma<-2 #initial sigma

## bandwidths settings
zeta_beta<-65
zeta_h<-0.01

## Convergence bounds
conv_betas_bound<-0.00001
conv_Hs_bound<-0.0001
conv_sigma_bound<-0.0001

## Number of replications and stopping iteration 
stop_iter_num<-1000

#number of bootstrap iterations
B<-100


#### Analysis ####
#choose the covariates at each transition
X_names<-c("age_10","log_nodes","log_er","log_pgr","meno","size20_50","size50","hormon","chemo",
           "grade3")
X_names12<-c("age_since_relapse_10","log_nodes","log_er","log_pgr","meno","size20_50","size50",
             "hormon","chemo","grade3")

X01<-as.matrix((data[,X_names]))
X02<-as.matrix((data[,X_names]))
X12<-as.matrix((data[,X_names12]))


results_with_frailty<-estimation_with_frailty(X01=X01,X02=X02,X12=X12,V=V,W=W,delta1=delta1,
                                              delta2=delta2,delta3=delta3,
                                              B=B,print=T,
                                              zeta_beta=zeta_beta,zeta_h=zeta_h,initial_sigma=initial_sigma,
                                              stop_iter_num=stop_iter_num,conv_betas_bound=conv_betas_bound,
                                              conv_Hs_bound=conv_Hs_bound,conv_sigma_bound=conv_sigma_bound)
#summary of the results
results_with_frailty$results
# estimated values
results_with_frailty$est
# estimated se
results_with_frailty$ESE
# show all bootstrap results
results_with_frailty$all_boot_res



# without frailty
results_without_frailty<-estimation_without_frailty(X01=X01,X02=X02,X12=X12,V=V,W=W,delta1=delta1,
                                                    delta2=delta2,delta3=delta3,
                                                    zeta_beta=zeta_beta,zeta_h=zeta_h,B=B,print=T)

#summary of the results
results_without_frailty