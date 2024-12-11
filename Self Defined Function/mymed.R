

mymed<-function(X,n,est,tseq,B){
  delta2=est$delta2
  tau=4
  betaz=list(est$betaz1,est$betaz2,est$betaz3)
  betax=list(est$betax1,est$betax2,est$betax3)
  
  lambda.m=list(c(exp(est$log0_r1),exp(est$log0_r2)),
                c(exp(est$log1_r1),exp(est$log1_r2)),c(exp(est$log2_r1),exp(est$log2_r2)))
    

  t.jump.m=list(0.1,0.1,0.1)
  t.jump.y=c(1,2)
  lambda.y=c(exp(est$log_h1),exp(est$log_h2),exp(est$log_h3))
  #Generate the M0/M1
  
  avgres=matrix(data=0,nrow=n,ncol=4*length(tseq))
  avgressq= matrix(data=0,nrow=n,ncol=4*length(tseq))
  for (b in 1:B){
    
    set.seed(1111+b)
    vi=rnorm(n,0,sqrt(exp(est$log_varc)))
    M_T=datagen_M_gap(X,vi,betaz,betax,tau,lambda.m,t.jump.m)
    R0=M_T[which((M_T$M==0)&(M_T$event==1)),][,-c(3,4)]
    R1=M_T[which((M_T$M==1)&(M_T$event==1)),][,-c(3,4)]
    tmpres=avgres
    for (i in 1:n){
      Xi=X[i]
      R0i=R0[which(R0$ID==i),]
      R1i=R1[which(R1$ID==i),]
      S00b=S00(Xi,vi[i],R0i,est$etaz,est$etax,est$etam,est$delta1,delta2,lambda.y,t.jump.y,tseq)
      S01b=S01(Xi,vi[i],R0i,est$etaz,est$etax,est$etam,est$delta1,delta2,lambda.y,t.jump.y,tseq)
      S10b=S10(Xi,vi[i],R1i,est$etaz,est$etax,est$etam,est$delta1,delta2,lambda.y,t.jump.y,tseq)
      S11b=S11(Xi,vi[i],R1i,est$etaz,est$etax,est$etam,est$delta1,delta2,lambda.y,t.jump.y,tseq)
      tmpres[i,]= c(S00b,S01b,S10b,S11b)
      
    }
    
    tmpressq=tmpres^2
    ###please check if the following update code is correct
    avgres=avgres+(tmpres-avgres)/b
    avgressq=avgressq+(tmpressq-avgressq)/b
  }
  avgsd=sqrt(avgressq-avgres^2)
  
  avg_S=apply(avgres,2,mean,na.rm=TRUE)
  sd_S=apply(avgres,2,sd,na.rm=TRUE)
  
  S00_n=avg_S[1:length(tseq)]
  S01_n=avg_S[(length(tseq)+1):(length(tseq)*2)]
  S10_n=avg_S[(length(tseq)*2+1):(length(tseq)*3)]
  S11_n=avg_S[(length(tseq)*3+1):(length(tseq)*4)]
  
  NDE01=S01_n-S00_n
  NDE10=S10_n-S11_n
  NIE01=S11_n-S01_n
  NIE10=S00_n-S10_n
  
  NDE=(NDE01-NDE10)/2
  NIE=(NIE01-NIE10)/2
  return(list(NDE=NDE,NIE=NIE))
  
}
