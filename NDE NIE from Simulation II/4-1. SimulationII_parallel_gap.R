
library(MASS)
library(survival)
library(dplyr)
###load in all self defined functions (please add)
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


####read in simulated datasets by simseed####
for (simseed in 0:199){
data=read.csv(paste("sim2_gapdata",as.character(simseed+1),".csv",sep=""),header=TRUE)
X2 <-data %>% group_by(ID) %>% filter(row_number(X2)==1) %>% select(ID,X2)
X=X2$X2


####read in estimated parameter by simseed####


est=read.table(paste("estIIgap",as.character(simseed+1),".txt",sep=""),header=FALSE,row.names=1)
est=as.data.frame(t(est))


###define parameters
tseq=seq(1,3,by=1)
###number of sampling
B=10000
###number of bootstrap 
M=100
n=length(X)

###get point estimate (save in a vector form)

myest=mymed(X,n,est,tseq,B)

###save results
write.csv(myest,paste("est2_gaplres",as.character(simseed+1),".csv",sep=""),row.names=FALSE)
}
