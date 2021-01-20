# examples for estimation with or without frailty
data<-read.csv(file="simulated_data.csv")
var_index01<-c(1,2)
var_index02<-c(2,3)
var_index12<-c(1,2,4)
X01<-as.matrix(data[,var_index01])
X02<-as.matrix(data[,var_index02])
X12<-as.matrix(data[,var_index12])
delta1<-data[,"delta1"]
delta2<-data[,"delta2"]
delta3<-data[,"delta3"]
V<-data[,"V"]
W<-delta1*data[,"W"]

results_with_frailty<-estimation_with_frailty_f(X01=X01,X02=X02,X12=X12,V=V,W=W,delta1=delta1,delta2=delta2,delta3=delta3,B=100,print=T)

results_without_frailty<-estimation_without_frailty_f(X01=X01,X02=X02,X12=X12,V=V,W=W,delta1=delta1,delta2=delta2,delta3=delta3,B=100,print=T)