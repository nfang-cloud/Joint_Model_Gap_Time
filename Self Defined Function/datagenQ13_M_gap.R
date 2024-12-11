####Function used to generate recurrent events M_Q1 and M_Q3
datagenQ13_M_gap<-function(X,Q1,Q3,vi,betaz,betax,tau,lambda.m,t.jump.m){
  
  ####Simulate M^z(t)
  ID=1:n
  beta.all=risk_Q1=risk_Q3=list()
  for (j in 1:length(betaz)){
    beta.all[[j]]=c(betaz[[j]],betax[[j]],1)
    risk_Q1[[j]]=c(exp(cbind(Q1,X,vi)%*%beta.all[[j]]))
    risk_Q3[[j]]=c(exp(cbind(Q3,X,vi)%*%beta.all[[j]]))
  }
  
  ####generate survival time until all individual has follow-up at lest time tau
  M_Q1=mysimRec_gap(t.jump.m,lambda.m,tau,risk_Q1,recurrent=TRUE)
  M_Q3=mysimRec_gap(t.jump.m,lambda.m,tau,risk_Q3,recurrent=TRUE)
  M_Q1$M=Q1
  M_Q3$M=Q3
  rbind(M_Q1,M_Q3)
}