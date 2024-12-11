library(MASS)
library(survival)
library(dplyr)
###load in all self defined functions ##############
source("Lambdainv.R")
source("mysimRec.R")
source("mysimRec_gap.R")
source("myprod.R")
source("datagen_M_gap.R")
source("S00.R")
source("S01.R")
source("S10.R")
source("S11.R")
source("mymed.R")

###setting running parameters###

for (simseed in 0:999){

bootbatch=simseed%%5+1
simseed=simseed%/%5

####read in simulated datasets by simseed####

data=read.csv(paste("sim5_gapdata",as.character(simseed+1),".csv",sep=""),header=TRUE)
X2 <-data %>% group_by(ID) %>% filter(row_number(X2)==1) %>% select(ID,X2)
X=X2$X2
n=length(X)

####read in estimated  parameters by simseed####

est=read.table(paste("estVgap",as.character(simseed+1),".txt",sep=""),header=FALSE,row.names=1)

est = cbind(as.data.frame(t(est)), delta2=0)
est_m=c(est$log0_r1,est$log0_r2,est$log1_r1,est$log1_r2,est$log2_r1,est$log2_r2,
        est$log_h1,est$log_h2,est$log_h3,
        est$betax1,est$betax2,est$betax3,est$betaz1,est$betaz2,est$betaz3,
        est$etaz,est$etax,est$etam,est$delta1,est$log_varc)

####read in estimated variance by simseed####

estvar=read.table(paste("covVgap",as.character(simseed+1),".txt",sep=""),header=FALSE,row.names=1,na.strings = ".")


###define parameters
tseq=seq(1,3,by=1)
###number of sampling 
B=10000
###number of bootstrap 
M=20

###get point estimate (save in a vector form)
myest=mymed(X,n,est,tseq,B)

myBest=matrix(data=NA,nrow=M,ncol=length(myest$NDE)*2)


###compute variance via Bootstrap
if (!is.na(sum(estvar))){
for (m in 1:M){
  
	set.seed(8888+(bootbatch-1)*20+m)
	estb=mvrnorm(1,mu=est_m,Sigma=estvar)
	est=cbind(as.data.frame(t(estb)),delta2=0)
	mymed_m=mymed(X,n,est,tseq,B)
	myBest[m,]=c(mymed_m$NDE,mymed_m$NIE)
	
}
}
###save results

write.csv(myBest,paste("Boot5_gaplres",as.character(simseed+1),"batch=", as.character(bootbatch),".csv",sep=""),row.names=FALSE)

}


