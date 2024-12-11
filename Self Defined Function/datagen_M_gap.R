#Function to generate M0 and M1
datagen_M_gap<-function(X,vi,betaz,betax,tau,lambda.m,t.jump.m){
  
  ####Simulate M^z(t)
  ID=1:n
  beta.all=risk0=risk1=list()
  for (j in 1:length(betaz)){
    beta.all[[j]]=c(betaz[[j]],betax[[j]],1)
    risk0[[j]]=c(exp(cbind(0,X,vi)%*%beta.all[[j]]))
    risk1[[j]]=c(exp(cbind(1,X,vi)%*%beta.all[[j]]))
  }
  ####generate survival time until all individual has follow-up at lest time tau
  M0=mysimRec_gap(t.jump.m,lambda.m,tau,risk0,recurrent=TRUE)
  M1=mysimRec_gap(t.jump.m,lambda.m,tau,risk1,recurrent=TRUE)
  M0$M=0
  M1$M=1
  rbind(M0,M1)
}  