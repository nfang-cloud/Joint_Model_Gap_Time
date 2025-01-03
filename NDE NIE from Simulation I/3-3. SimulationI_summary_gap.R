library(MASS)
library(survival)
library(dplyr)
library(reda)
library(clusterGeneration)
###load in self-defined functions#####
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
source("datagen_X_gap.R")


#######################compute true effects using large data#####################################
#################################################################################################
##### TRUE parameter ######

n=10000
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

#####Generate a large dataset#######

set.seed(8888)
vi=rnorm(n)*sigmav
data=datagen_X_gap(n,vi,tau,betaz,betax,lambda.m,t.jump.m,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,mod=1,cen=TRUE)$alldata

X2 <-data %>% group_by(ID) %>% filter(row_number(X2)==1) %>% select(ID,X2)
X=X2$X2
n=length(X)
tseq=seq(1,3,by=1)

#####Get the true value####

est=cbind(log0_r1=log(lambda.m[[1]][1]),log0_r2=log(lambda.m[[1]][2]),log1_r1=log(lambda.m[[2]][1]),
          log1_r2=log(lambda.m[[2]][2]),log2_r1=log(lambda.m[[3]][1]),log2_r2=log(lambda.m[[3]][2]),
          log_h1=log(lambda.y[1]),log_h2=log(lambda.y[2]),log_h3=log(lambda.y[3]),
          betax1=betax[[1]],betax2=betax[[2]],betax3=betax[[3]],
          betaz1=betaz[[1]],betaz2=betaz[[2]],betaz3=betaz[[3]],
          etaz=etaz,etax=etax,etam=etam,delta1=delta1,delta2=0,log_varc=log(sigmav^2))

est=as.data.frame(est)

truevalue0=mymed(X,n,est,tseq,B=1000)
truevalue=c(truevalue0$NDE,truevalue0$NIE)
truevalue1=truevalue$x[c(2,4,6,9,11,13)]

###handle each simulated dataset

allest=NULL
allsd=NULL
myest=NULL

for (i in 1:200){
  myBest=NULL
  myest1=read.csv(paste("est_gaplres",as.character(i),".csv",sep=""))[c(2,4,6),]
	myest=c(myest1$NDE,myest1$NIE)
	for (batch in 1:5){
	  for (m in 1:4){
	  skip_to_next <- FALSE
	  tryCatch({myBest_temp=read.csv(paste("Boot_gaplres1",as.character(i),"batch=",as.character(batch),"M=",as.character(m*5),".csv",sep=""))[c(1,2,3,4,5),]}, error = function(e) { skip_to_next <<- TRUE})
	  
	  if(skip_to_next) { next }
	  if (!is.na(sum(myBest_temp))){
	  myBest=rbind(myBest,myBest_temp)
	  }
	  }
	}
	
	
	if (!is.null(myBest)){
	myBest[abs(myBest)>1]<-NA
	myest.sd=apply(myBest,2,sd,na.rm=TRUE)
	allest=rbind(allest,myest)
	allsd=rbind(allsd,myest.sd)
	}
}

bias=apply(allest,2,mean,na.rm=TRUE)-truevalue1
ese=apply(allsd*is.finite(allsd),2,mean,na.rm=TRUE)
meSE=apply(allsd,2,median, na.rm=TRUE)
sd=apply(allest,2,sd,na.rm=TRUE)
cr=apply((abs(allest-rep(1,nrow(allest))%o%truevalue1)/allsd)<1.96,2,mean,na.rm=TRUE)

###Save the result###
write.csv(rbind(truevalue1,bias,ese,sd,meSE,cr),"sim1gapres_log.csv")



