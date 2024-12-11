##Load the package

library(MASS)
library(sas7bdat)
library(reda)
library(survival)
library(dplyr)
###load in all self defined functions###
source("Lambdainv.R")
source("myprod.R")
source("mysimRec.R")
source("mysimRec_gap.R")
source("datagen_X.R")


################################################################################################
####################################Setting###################################################
n=400
tau=4
betaz=list(-0.5,-0.4,-0.3)
betax=list(0.5,0.5,0.5)
lambda.m=list(c(3,3.2),c(3.6,4.5),c(4.2,5))
t.jump.m=list(0.1,0.1,0.1)
lambda.y=c(0.9,1,1.1)
t.jump.y=c(1,2)
etaz=-1
etax=1
delta1=1
delta2=0
etam=0.25
sigmav=1

#################################Generate 200 datasets from Setting I#########################
for (iter in 1:200){
  set.seed(8888+iter)
  vi=rnorm(n)*sigmav
  data=datagen_X_gap(n,vi,tau,betaz,betax,lambda.m,t.jump.m,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,mod=1,cen=TRUE)
  
  write.csv(data$alldata,paste("sim1_gapdata",as.character(iter),".csv",sep=""),row.names=F)

  }


#################################Generate 200 datasets from Setting II#########################
delta2=0.15
for (iter in 1:200){
  set.seed(8888+iter)
  vi=rnorm(n)*sigmav
  data=datagen_X_gap(n,vi,tau,betaz,betax,lambda.m,t.jump.m,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,mod=2,cen=TRUE)
  write.csv(data$alldata,paste("sim2_gapdata",as.character(iter),".csv",sep=""),row.names=F)
  
}


#################################Generate 200 datasets from Setting III#########################
for (iter in 1:200){
  set.seed(8888+iter)
  v_gamma=rgamma(n,shape = 1, scale = 1)
  vi=log(v_gamma)
  data=datagen_X_gap(n,vi,tau,betaz,betax,lambda.m,t.jump.m,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,mod=2,cen=TRUE)
  write.csv(data$alldata,paste("sim3_gapdata",as.character(iter),".csv",sep=""),row.names=F)
  
}

#################################Generate 200 datasets from Setting V#########################
####New Settings#####
n=400
tau=4
betaz=list(-0.5,-0.4,-0.3)
betax=list(0.5,0.5,0.5)
lambda.m=list(c(1,3),c(5,2),c(6,1))
t.jump.m=list(2,2,2)
lambda.y=c(0.9,1,1.1)
t.jump.y=c(1,2)
etaz=-1
etax=1
delta1=1
delta2=0
etam=0.25
sigmav=1


for (iter in 1:200){
  set.seed(8888+iter)
  vi=rnorm(n)*sigmav
  data=datagen_X_gap(n,vi,tau,betaz,betax,lambda.m,t.jump.m,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,mod=1,cen=TRUE)
  
  write.csv(data$alldata,paste("sim5_gapdata",as.character(iter),".csv",sep=""),row.names=F)

  }