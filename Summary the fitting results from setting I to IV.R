
####################################################################
##########################True parameters###########################
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

################Computing log() for lambda.m, lambda.y and sigmav###############

log_r0_m1=log(lambda.m[[1]][1])
log_r0_m2=log(lambda.m[[1]][2])
log_r1_m1=log(lambda.m[[2]][1])
log_r1_m2=log(lambda.m[[2]][2])
log_r2_m1=log(lambda.m[[3]][1])
log_r2_m2=log(lambda.m[[3]][2])

log_lam_y1=log(lambda.y[1])
log_lam_y2=log(lambda.y[2])
log_lam_y3=log(lambda.y[3])

log_varc=log(sigmav^2)

##############################################################################
###################################Setting I###################################
True_est=cbind(log_lam_m1,log_lam_m2,log_lam_m3,log_lam_m4,
            log_lam_y1,log_lam_y2,log_lam_y3,log_lam_y4,
            betax,betaz,etaz,etax,etam,delta1,log_varc)

#####Import the estimation and covariance######


est=NULL
var=NULL
for (iter in 1:200){
  est.next=NULL
  var.next=NULL
  cov=NULL
  skip_to_next <- FALSE
  tryCatch({myBest_temp=read.table(paste("covIgap",as.character(iter),".txt",sep=""))}, error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  est.next=read.table(paste("estIgap",as.character(iter),".txt",sep=""))
  est.next=as.data.frame(t(est.next))
  cov=read.table(paste("covIgap",as.character(iter),".txt",sep=""),row.names=1,na.strings=c("."))
  var.next=t(diag(as.matrix(cov)))
  est=rbind(est,est.next[2,])
  var=rbind(var,var.next)
}

names(est)<-c("log_r0_m1","log_r0_m2",
              "log_r1_m1","log_r1_m2",
              "log_r2_m1","log_r2_m2",
              "log_lam_y1","log_lam_y2","log_lam_y3",
              "betax1","betax2","betax3","betaz1","betaz2","betaz3",
              "etaz","etax","etam","delta1",
              
              "log_varc"
              )

names(var)<-c("log_r0_m1","log_r0_m2",
              "log_r1_m1","log_r1_m2",
              "log_r2_m1","log_r2_m2",
              "log_lam_y1","log_lam_y2","log_lam_y3",
              "betax1","betax2","betax3","betaz1","betaz2","betaz3",
              "etaz","etax","etam","delta1",
              
              "log_varc"
              ) 

est=as.data.frame(sapply(est, as.numeric))


################Computing bias, empirical SD, median SE, CR###############

est_mean <- colMeans(est)
est_mean=apply(est,2,mean)
est_bias=apply(est-rep(1,nrow(est))%o%c(True_est),2,mean,na.rm=TRUE)
est_SD=apply(est,2,sd)
est_CR=apply((abs(est-rep(1,nrow(est))%o%c(True_est))/sqrt(var))<1.96,2,mean,na.rm=TRUE)
est_meSE=apply(sqrt(var),2,median, na.rm=TRUE)
est_ESE=apply(sqrt(var),2,mean,na.rm=TRUE)
my_eva<-data.frame(Mean=est_mean,Bias=est_bias,SD=est_SD,meSE=est_meSE, ESE=est_ESE,CR=est_CR)
write.csv(my_eva,paste("sim1_estbiagap",".csv",sep=""),row.names = TRUE)


##############################################################################
###################################Setting II###################################
delta2=0.15
True_est=cbind(log_lam_m1,log_lam_m2,log_lam_m3,log_lam_m4,
          log_lam_y1,log_lam_y2,log_lam_y3,log_lam_y4,
          betax,betaz,etaz,etax,etam,delta1,delta2,log_varc)
#####Import the estimation and covariance######

est=NULL
var=NULL
for (iter in 1:200){
  est.next=NULL
  var.next=NULL
  cov=NULL
  skip_to_next <- FALSE
  tryCatch({myBest_temp=read.table(paste("covIIgap",as.character(iter),".txt",sep=""))}, error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  est.next=read.table(paste("estIIgap",as.character(iter),".txt",sep=""))
  est.next=as.data.frame(t(est.next))
  cov=read.table(paste("covIIgap",as.character(iter),".txt",sep=""),row.names=1,na.strings=c("."))
  var.next=t(diag(as.matrix(cov)))
  est=rbind(est,est.next[2,])
  var=rbind(var,var.next)
}

names(est)<-c("log_r0_m1","log_r0_m2",
              "log_r1_m1","log_r1_m2",
              "log_r2_m1","log_r2_m2",
              "log_lam_y1","log_lam_y2","log_lam_y3",
              "betax1","betax2","betax3","betaz1","betaz2","betaz3",
              "etaz","etax","etam","delta1","delta2",
              
              "log_varc"
)

names(var)<-c("log_r0_m1","log_r0_m2",
              "log_r1_m1","log_r1_m2",
              "log_r2_m1","log_r2_m2",
              "log_lam_y1","log_lam_y2","log_lam_y3",
              "betax1","betax2","betax3","betaz1","betaz2","betaz3",
              "etaz","etax","etam","delta1","delta2",
              
              "log_varc"
) 

est=as.data.frame(sapply(est, as.numeric))

################Computing bias, empirical SD, median SE, CR###############
est_mean <- colMeans(est)
est_mean=apply(est,2,mean)
est_bias=apply(est-rep(1,nrow(est))%o%c(True_est),2,mean)
est_SD=apply(est,2,sd)
est_CR=apply((abs(est-rep(1,nrow(est))%o%c(True_est))/sqrt(var))<1.96,2,mean,na.rm=TRUE)
est_meSE=apply(sqrt(var),2,median, na.rm=TRUE)
est_ESE=apply(sqrt(var),2,mean,na.rm=TRUE)
my_eva<-data.frame(Mean=est_mean,Bias=est_bias,SD=est_SD,meSE=est_meSE, ESE=est_ESE,CR=est_CR)
write.csv(my_eva,paste("sim2_estbiagap",".csv",sep=""),row.names=TRUE)



##############################################################################
###################################Setting III###################################

True_est=cbind(log_lam_m1,log_lam_m2,log_lam_m3,log_lam_m4,
          log_lam_y1,log_lam_y2,log_lam_y3,log_lam_y4,
          betax,betaz,etaz,etax,etam,delta1,delta2,log_varc)

#####Import the estimation and covariance######

est=NULL
var=NULL
for (iter in 1:200){
  est.next=NULL
  var.next=NULL
  cov=NULL
  skip_to_next <- FALSE
  tryCatch({myBest_temp=read.table(paste("covIIIgap",as.character(iter),".txt",sep=""))}, error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  est.next=read.table(paste("estIIIgap",as.character(iter),".txt",sep=""))
  est.next=as.data.frame(t(est.next))
  cov=read.table(paste("covIIIgap",as.character(iter),".txt",sep=""),row.names=1,na.strings=c("."))
  var.next=t(diag(as.matrix(cov)))
  est=rbind(est,est.next[2,])
  var=rbind(var,var.next)
}

names(est)<-c("log_lam_m1","log_lam_m2","log_lam_m3","log_lam_m4","log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4","betax","betaz","etaz","etax","etam","delta1","delta2","log_varc") 
names(var)<-c("log_lam_m1","log_lam_m2","log_lam_m3","log_lam_m4","log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4","betax","betaz","etaz","etax","etam","delta1","delta2","log_varc") 


est=as.data.frame(sapply(est, as.numeric))

################Computing bias, empirical SD, median SE, CR###############
est_mean <- colMeans(est)
est_mean=apply(est,2,mean)
est_bias=apply(est-rep(1,nrow(est))%o%c(True_est),2,mean)
est_SD=apply(est,2,sd)
est_CR=apply((abs(est-rep(1,nrow(est))%o%c(True_est))/sqrt(var))<1.96,2,mean,na.rm=TRUE)
est_meSE=apply(sqrt(var),2,median, na.rm=TRUE)
est_ESE=apply(sqrt(var),2,mean,na.rm=TRUE)
my_eva<-data.frame(Mean=est_mean,Bias=est_bias,SD=est_SD,meSE=est_meSE, ESE=est_ESE,CR=est_CR)
write.csv(my_eva,paste("sim3_estbialog",".csv",sep=""),row.names=TRUE)



##############################################################################
###################################Setting IV#################################

True_est=cbind(log_lam_m1,log_lam_m2,log_lam_m3,log_lam_m4,
          log_lam_y1,log_lam_y2,log_lam_y3,log_lam_y4,
          betax,betaz,etaz,etax,etam,delta1,log_varc)
#####Import the estimation and covariance######

est=NULL
var=NULL
for (iter in 1:200){
  est.next=NULL
  var.next=NULL
  cov=NULL
  skip_to_next <- FALSE
  tryCatch({myBest_temp=read.table(paste("covIVgap",as.character(iter),".txt",sep=""))}, error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  est.next=read.table(paste("estIVgap",as.character(iter),".txt",sep=""))
  est.next=as.data.frame(t(est.next))
  cov=read.table(paste("covIVgap",as.character(iter),".txt",sep=""),row.names=1,na.strings=c("."))
  var.next=t(diag(as.matrix(cov)))
  est=rbind(est,est.next[2,])
  var=rbind(var,var.next)
}

names(est)<-c("log_r0_m1","log_r0_m2",
              "log_r1_m1","log_r1_m2",
              "log_r2_m1","log_r2_m2",
              "log_lam_y1","log_lam_y2","log_lam_y3",
              "betax1","betax2","betax3","betaz1","betaz2","betaz3",
              "etaz","etax","etam","delta1",
              
              "log_varc"
)

names(var)<-c("log_r0_m1","log_r0_m2",
              "log_r1_m1","log_r1_m2",
              "log_r2_m1","log_r2_m2",
              "log_lam_y1","log_lam_y2","log_lam_y3",
              "betax1","betax2","betax3","betaz1","betaz2","betaz3",
              "etaz","etax","etam","delta1",
              
              "log_varc"
) 

est=as.data.frame(sapply(est, as.numeric))

################Computing bias, empirical SD, median SE, CR###############
est_mean <- colMeans(est)
est_mean=apply(est,2,mean)
est_bias=apply(est-rep(1,nrow(est))%o%c(True_est),2,mean)
est_SD=apply(est,2,sd)
est_CR=apply((abs(est-rep(1,nrow(est))%o%c(True_est))/sqrt(var))<1.96,2,mean,na.rm=TRUE)
est_meSE=apply(sqrt(var),2,median, na.rm=TRUE)
est_ESE=apply(sqrt(var),2,mean,na.rm=TRUE)
my_eva<-data.frame(Mean=est_mean,Bias=est_bias,SD=est_SD,meSE=est_meSE,ESE=est_ESE,CR=est_CR)
write.csv(my_eva,paste("sim4_estbiagap",".csv",sep=""),row.names=TRUE)

