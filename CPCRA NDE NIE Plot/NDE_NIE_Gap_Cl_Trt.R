
library(MASS)
library(survival)
library(dplyr)
###load in self-defined functions###
source("Lambdainv.R")
source("mysimRec.R")
source("mysimRec_gap.R")
source("datagenQ13_M_gap.R")
source("myprod.R")
source("S.R")
source("mymed_Q13.R")


####check real data
####read in data
covdata=read.csv("covdata.csv")

covdata$trt=covdata$RANDGRP-1
covdata$gender=covdata$GENDER-1
covdata$hemobl=covdata$HEMOBL-12
covdata$stratum=covdata$STRATUM
covdata$cd4bl=covdata$CD4BL/100
covdata$prevoi=covdata$PREVOI
Z=covdata[,"trt"]
Q1=0
Q3=1

tseq=seq(1,15,by=1)

X=as.matrix(covdata[,c("prevoi","gender","stratum","hemobl","cd4bl")])
tau=15
n=nrow(X)

############################Import the Estimation######################
 
 est=read.table(paste("estI_realgap",".txt",sep=""),header=FALSE,row.names=1)
 est = as.data.frame(t(est))
 
 betaz=list(est$betaz,est$betaz)
 betax=list(c(est$betax1,est$betax2,est$betax3,est$betax4,est$betax5),
            c(est$betax1,est$betax2,est$betax3,est$betax4,est$betax5))
 lambda.m=list(exp(c(est$log0_r1,est$log0_r2,est$log0_r3)),exp(c(est$log1_r1,est$log1_r2,est$log1_r3)))
 t.jump.m=list(c(2.5,5.983),c(6.6,10.033))
 
 lambda.y=exp(c(est$log_h1,est$log_h2,est$log_h3,est$log_h4,
                est$log_h5,est$log_h6,est$log_h7,est$log_h8,
                est$log_h9,est$log_h10))
 t.jump.y=c(2.4,3.633,5.433,7.3,8.6,10.233,11.133,12.267,14.867)
 etaz=est$etaz
 etax=c(est$etax1,est$etax2,est$etax3,est$etax4,est$etax5)
 etam=est$etam
 delta1=est$delta1
 delta2=0
 sigmav=sqrt(exp(est$log_varc))
 
 
 cov1=read.table("covI_realgap.txt",row.names=1)
#######################est1 from the real Estimates###############
 
 est1=c(est$log0_r1,est$log0_r2,est$log0_r3,
        est$log1_r1,est$log1_r2,est$log1_r3,
        est$log_h1,est$log_h2,est$log_h3,est$log_h4,est$log_h5,
        est$log_h6,est$log_h7,est$log_h8,est$log_h9,est$log_h10,
        est$betaz,est$betax1,est$betax2,est$betax3,est$betax4,est$betax5,
      
        est$etaz,est$etax1,est$etax2,est$etax3,est$etax4,est$etax5,
        est$etam,est$delta1,est$log_varc)
 
 
 ###parametric Bootstrap
 
 B=100000/20
 
##Mediation Analysis##
for (batch in 1:20){

  for (simseed in 0:999){
 set.seed(1111+simseed)
 Bres1_new=NULL
 Best=mvrnorm(n = 1, mu=est1, Sigma=cov1)
 Bsigmav=sqrt(exp(Best[31]))
 Blambda.m=list(exp(Best[1:3]),exp(Best[4:6]))
 Blambda.y=exp(Best[7:16])
 Bbetaz=list(Best[17],Best[17])
 Bbetax=list(Best[18:22],Best[18:22])
 Betaz=Best[23]
 Betax=Best[c(24:28)]
 Bdelta1=Best[30]
 Bdelta2=0
 Betam=Best[29]
 Bres1_new=mymed_Q13(X,n,Z,Q1,Q3,Bsigmav,tau,Bbetaz,Bbetax,Blambda.m,t.jump.m,Betaz,Betax,Betam,Bdelta1,Bdelta2,Blambda.y,t.jump.y,tseq,B,seed=1111+batch*B)
 
 
 write.csv(Bres1_new,paste("Bres_Gap_Trt",as.character(1111+simseed),"batch=",as.character(batch),".csv",sep=""),row.names=FALSE)
}
}
 