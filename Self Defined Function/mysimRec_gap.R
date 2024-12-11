mysimRec_gap<-function(t.jump_list,lambda_list,tau,risk_list,recurrent=TRUE){
  n=lengths(risk_list)[1]
  ID=1:n
  s=rep(0,n)
  ###use to record the data
  dat=NULL
  ###use this to record the total time
  T=rep(0,n)
  sss=1:n
  count=0
  
  while(length(sss)>0){
    count=count+1
    if (count<=length(lambda_list)){
      lambda=lambda_list[[count]]
      risk=risk_list[[count]]
      t.jump=t.jump_list[[count]]
    }
    if (count>length(lambda_list)){
      lambda=lambda_list[[length(lambda_list)]]
      risk=risk_list[[length(lambda_list)]]
      t.jump=t.jump_list[[length(lambda_list)]]
    }
    ###generate gap times
    u=-log(runif(n))
    ID=sss
    D=Lambdainv(u[sss],t.jump,lambda,tau,risk[sss])
    ###update time
    T[sss]=T[sss]+D
    sss=which(T<tau)
    if (!recurrent){sss=NULL}
    event=as.numeric(T<tau)
    tmp=data.frame(ID=ID,stoptime=T[ID],event=event[ID])
    dat=rbind(dat,tmp)
  }
  dat[order(dat$ID,dat$stoptime),]
  return(dat[order(dat$ID,dat$stoptime),])
}