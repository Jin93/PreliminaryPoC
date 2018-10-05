########### This code is for the implement of the propose Bayesian approach after adding an interim analysis for futility stop. 
########### The approach was considered but not adopted in the end. The results were not shown in the manuscript.
load("simulationdata.RData")
K=3
q0=0.2 # historical benchmark response rate
q1=0.3 # target response rate
num.sim=2000 # number of simulations per setting
alpha=0.1 #level of false positive rate we wish to control.
p.sce=t(sapply(0:K,FUN=function(x){c(rep(q1,x),rep(q0,K-x))}))
######## obtained from the minimax Simon's two-stage design from 
######## http://cancer.unc.edu/biostatistics/program/ivanova/SimonsTwoStageDesign.aspx
bayesinterimtable=matrix(NA,12,5+2*K+2)
bayesinterimtable[1,1:5]=c(17,35,4,11,21.4)
bayesinterimtable[2,1:5]=c(13,27,2,9,20)
bayesinterimtable[3,1:5]=c(12,27,2,9,18.6)
bayesinterimtable[4,1:5]=c(14,24,3,8,17)
bayesinterimtable[5,1:5]=c(8,24,1,8,15.9)
bayesinterimtable[6,1:5]=c(10,20,2,7,13.2)
bayesinterimtable[7,1:5]=c(11,17,3,6,12)
bayesinterimtable[8,1:5]=c(6,17,1,6,9.8)
bayesinterimtable[9,1:5]=c(3,13,0,5,7.9)
bayesinterimtable[10,1:5]=c(5,10,1,4,6.3)
bayesinterimtable[11,1:5]=c(3,7,0,3,5)
bayesinterimtable[12,1:5]=c(2,7,0,3,3.8)

Decision.bayesinterim=list()
for (x in 1:nrow(bayesinterimtable)
{
  Ni1=bayesinterimtable[x,1]
  Ni=bayesinterimtable[x,2]
  r1=bayesinterimtable[x,3]
  r2=bayesinterimtable[x,4]
  samplesize.bayesinterim=ceiling(bayesinterimtable[x,5])
  nik=matrix(NA,2,K) # number of patients in indication k at state i
  rik=matrix(NA,2,K) # number of responders in indication k at state i
  nik[1,]=rep(Ni1,K) # number of patients enrolled at stage 1
  Decision.bayesinterim[[samplesize.bayesinterim]]=rep(0,K+1)
  
  ###### tuning:
  decision.bayesinterim=Tstat=numeric()
  p0=rep(q0,K)
  
  tp=which(p0>=q1)
  tn=which(p0<q1)
  
  #samplesize=0
  for (sim in 1:num.sim)
  {
    ##### Stage 1:
    rik[1,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[1,x],prob=p0[x])})
    
    ## Futility stop:
    stage2.stop=which(rik[1,]<=r1)
    stage2.cont=which(rik[1,]>r1)
    nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
    if (length(stage2.stop)>0) 
    {
      Tstat[sim]=0
    }
    
    ## Stage 2:
    if (length(stage2.cont)>0)
    {
      rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
      ri=colSums(as.matrix(rik[,stage2.cont]))
      ni=colSums(as.matrix(nik[,stage2.cont]))
      K1=length(stage2.cont)
      
      q=rep(log(((q0+q1)/2)/(1-(q0+q1)/2)),K) ## can consider different settings based on histological data
      ############ Jags model for BHM:
      jags.data <- list("n"=ni, "Y"=ri, "K"=K1, "q"=q)
      jags.fit <- jags.model(file = "C:/Users/e0359820/Desktop/Jin/Signal Detection Project/a-ina_p.txt",data = jags.data,
                             n.adapt=1000,n.chains=1,quiet=T)
      update(jags.fit, 4000)
      bayes.out <- coda.samples(jags.fit,variable.names = c("p","d","pi","delta","tausq","mu1","mu2"),n.iter=10000)
      
      ### Interim analysis:
      if (K1 == 1)
      {
        Tstat[sim]=sum(bayes.out[[1]][,"d"]>0)/nrow(bayes.out[[1]])
      }
      if (K1 > 1)
      {
        Tstat[sim]=sum(apply(bayes.out[[1]][,sapply(1:K1,FUN=function(x){paste("d[",x,"]",sep="")})],1,sum)>0)/nrow(bayes.out[[1]])
      }
    }
    print(sim)
  }
  c.bayesinterim=quantile(Tstat,1-alpha)
  
  ###################### Simulations:
  for (scenario in 2:(nrow(p.sce)))
  {
    p0=p.sce[scenario,]
    decision.bayesinterim=numeric()
    
    tp=which(p0>=q1)
    tn=which(p0<q1)
    
    decision.bayesinterim=matrix(NA,num.sim,K)
    samplesize=0
    for (sim in 1:num.sim)
    {
      ##### Stage 1:
      rik[1,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[1,x],prob=p0[x])})
      
      ## Futility stop:
      stage2.stop=which(rik[1,]<=r1)
      stage2.cont=which(rik[1,]>r1)
      nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
      if (length(stage2.stop)>0) 
      {
        Tstat[sim]=0
        decision.bayesinterim[sim]=0
      }
      
      ## Stage 2:
      if (length(stage2.cont)>0)
      {
        rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
        ri=colSums(as.matrix(rik[,stage2.cont]))
        ni=colSums(as.matrix(nik[,stage2.cont]))
        K1=length(stage2.cont)
        
        q=rep(log(((q0+q1)/2)/(1-(q0+q1)/2)),K) ## can consider different settings based on histological data
        ############ Jags model for BHM:
        jags.data <- list("n"=ni, "Y"=ri, "K"=K1, "q"=q)
        jags.fit <- jags.model(file = "C:/Users/e0359820/Desktop/Jin/Signal Detection Project/a-ina_p.txt",data = jags.data,
                               n.adapt=1000,n.chains=1,quiet=T)
        update(jags.fit, 4000)
        bayes.out <- coda.samples(jags.fit,variable.names = c("p","d","pi","delta","tausq","mu1","mu2"),n.iter=10000)
        
        ### Interim analysis:
        if (K1 == 1)
        {
          Tstat[sim]=sum(bayes.out[[1]][,"d"]>0)/nrow(bayes.out[[1]])
        }
        if (K1 > 1)
        {
          Tstat[sim]=sum(apply(bayes.out[[1]][,sapply(1:K1,FUN=function(x){paste("d[",x,"]",sep="")})],1,sum)>0)/nrow(bayes.out[[1]])
        }
        decision.bayesinterim[sim]=ifelse(Tstat[sim]>c.bayesinterim,1,0)
      }
      Decision.bayesinterim[[samplesize.bayesinterim]][scenario]=Decision.bayesinterim[[samplesize.bayesinterim]][scenario] + ifelse((length(tp)>0)&(decision.bayesinterim[sim]==0),1,0)
      
      print(sim)
      print(Decision.bayesinterim[[samplesize.bayesinterim]][scenario])
      samplesize=samplesize+sum(nik)
    }
    bayesinterimtable[x,5+scenario]=Decision.bayesinterim[[samplesize.bayesinterim]][scenario]/num.sim
    bayesinterimtable[x,9+scenario]=samplesize/num.sim
  }
  save(bayesinterimtable,Decision.bayesinterimtable,file="bayesinterim_simu_fn_q0=0.3_k3.RData")
}



K=6
q0=0.2
q1=0.3
p.sce=t(sapply(0:K,FUN=function(x){c(rep(q1,x),rep(q0,K-x))}))
######## obtained from the minimax Simon's two-stage design from 
######## http://cancer.unc.edu/biostatistics/program/ivanova/SimonsTwoStageDesign.aspx
bayesinterimtable=matrix(NA,13,5+2*K+2)
bayesinterimtable[1,1:5]=c(17,36,4,12,21.6)
bayesinterimtable[2,1:5]=c(15,32,3,11,21)
bayesinterimtable[3,1:5]=c(11,32,2,11,19)
bayesinterimtable[4,1:5]=c(11,28,2,10,17.5)
bayesinterimtable[5,1:5]=c(13,25,3,9,16)
bayesinterimtable[6,1:5]=c(11,21,2,8,14.8)
bayesinterimtable[7,1:5]=c(10,21,2,8,13.5)
bayesinterimtable[8,1:5]=c(7,18,1,7,11.7)
bayesinterimtable[9,1:5]=c(3,15,0,6,8.9)
bayesinterimtable[10,1:5]=c(5,15,1,6,7.6)
bayesinterimtable[11,1:5]=c(2,12,0,5,5.6)
bayesinterimtable[12,1:5]=c(2,9,0,4,4.5)
bayesinterimtable[13,1:5]=c(2,6,0,3,3.4)
colnames(bayesinterimtable)=c("n1","n2","r1","r2","ss1","type1",
                    sapply(2:(K+1),FUN=function(x){paste("type2-",x,sep="")}),
                    sapply(1:(K+1),FUN=function(x){paste("ss",x,sep="")}))

Decision.bayesinterim=list()
for (x in 1:nrow(bayesinterimtable))
{
  Ni1=bayesinterimtable[x,1]
  Ni=bayesinterimtable[x,2]
  r1=bayesinterimtable[x,3]
  r2=bayesinterimtable[x,4]
  samplesize.bayesinterim=ceiling(bayesinterimtable[x,5])
  nik=matrix(NA,2,K) # number of patients in indication k at state i
  rik=matrix(NA,2,K) # number of responders in indication k at state i
  nik[1,]=rep(Ni1,K) # number of patients enrolled at stage 1
  Decision.bayesinterim[[samplesize.bayesinterim]]=rep(0,K+1)
  
  ###### tuning:
  decision.bayesinterim=Tstat=numeric()
  p0=rep(q0,K)
  
  tp=which(p0>=q1)
  tn=which(p0<q1)
  
  #samplesize=0
  for (sim in 1:num.sim)
  {
    ##### Stage 1:
    rik[1,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[1,x],prob=p0[x])})
    
    ## Futility stop:
    stage2.stop=which(rik[1,]<=r1)
    stage2.cont=which(rik[1,]>r1)
    nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
    if (length(stage2.stop)>0) 
    {
      Tstat[sim]=0
    }
    
    ## Stage 2:
    if (length(stage2.cont)>0)
    {
      rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
      ri=colSums(as.matrix(rik[,stage2.cont]))
      ni=colSums(as.matrix(nik[,stage2.cont]))
      K1=length(stage2.cont)
      
      q=rep(log(((q0+q1)/2)/(1-(q0+q1)/2)),K) ## can consider different settings based on histological data
      ############ Jags model for BHM:
      jags.data <- list("n"=ni, "Y"=ri, "K"=K1, "q"=q)
      jags.fit <- jags.model(file = "C:/Users/e0359820/Desktop/Jin/Signal Detection Project/a-ina_p.txt",data = jags.data,
                             n.adapt=1000,n.chains=1,quiet=T)
      update(jags.fit, 4000)
      bayes.out <- coda.samples(jags.fit,variable.names = c("p","d","pi","delta","tausq","mu1","mu2"),n.iter=10000)
      
      ### Interim analysis:
      if (K1 == 1)
      {
        Tstat[sim]=sum(bayes.out[[1]][,"d"]>0)/nrow(bayes.out[[1]])
      }
      if (K1 > 1)
      {
        Tstat[sim]=sum(apply(bayes.out[[1]][,sapply(1:K1,FUN=function(x){paste("d[",x,"]",sep="")})],1,sum)>0)/nrow(bayes.out[[1]])
      }
    }
    print(sim)
  }
  c.bayesinterim=quantile(Tstat,1-alpha)
  ###################### Simulations:
  for (scenario in 2:(nrow(p.sce)))
  {
    p0=p.sce[scenario,]
    decision.bayesinterim=numeric()
    
    tp=which(p0>=q1)
    tn=which(p0<q1)
    
    decision.bayesinterim=matrix(NA,num.sim,K)
    samplesize=0
    for (sim in 1:num.sim)
    {
      ##### Stage 1:
      rik[1,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[1,x],prob=p0[x])})
      
      ## Futility stop:
      stage2.stop=which(rik[1,]<=r1)
      stage2.cont=which(rik[1,]>r1)
      nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
      if (length(stage2.stop)>0) 
      {
        Tstat[sim]=0
        decision.bayesinterim[sim]=0
      }
      
      ## Stage 2:
      if (length(stage2.cont)>0)
      {
        rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
        ri=colSums(as.matrix(rik[,stage2.cont]))
        ni=colSums(as.matrix(nik[,stage2.cont]))
        K1=length(stage2.cont)
        
        q=rep(log(((q0+q1)/2)/(1-(q0+q1)/2)),K) ## can consider different settings based on histological data
        ############ Jags model for BHM:
        jags.data <- list("n"=ni, "Y"=ri, "K"=K1, "q"=q)
        jags.fit <- jags.model(file = "C:/Users/e0359820/Desktop/Jin/Signal Detection Project/a-ina_p.txt",data = jags.data,
                               n.adapt=1000,n.chains=1,quiet=T)
        update(jags.fit, 4000)
        bayes.out <- coda.samples(jags.fit,variable.names = c("p","d","pi","delta","tausq","mu1","mu2"),n.iter=10000)
        
        ### Interim analysis:
        if (K1 == 1)
        {
          Tstat[sim]=sum(bayes.out[[1]][,"d"]>0)/nrow(bayes.out[[1]])
        }
        if (K1 > 1)
        {
          Tstat[sim]=sum(apply(bayes.out[[1]][,sapply(1:K1,FUN=function(x){paste("d[",x,"]",sep="")})],1,sum)>0)/nrow(bayes.out[[1]])
        }
        decision.bayesinterim[sim]=ifelse(Tstat[sim]>c.bayesinterim,1,0)
      }
      Decision.bayesinterim[[samplesize.bayesinterim]][scenario]=Decision.bayesinterim[[samplesize.bayesinterim]][scenario] + ifelse((length(tp)>0)&(decision.bayesinterim[sim]==0),1,0)
      
      print(sim)
      print(Decision.bayesinterim[[samplesize.bayesinterim]][scenario])
      samplesize=samplesize+sum(nik)
    }
    bayesinterimtable[x,5+scenario]=Decision.bayesinterim[[samplesize.bayesinterim]][scenario]/num.sim
    bayesinterimtable[x,12+scenario]=samplesize/num.sim
  }
  save(bayesinterimtable,Decision.bayesinterim,file="bayesinterim_simu_fn_q0=0.3_k6.RData")
}
