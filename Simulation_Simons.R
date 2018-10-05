K=3 # total number of indications
q0=0.2 # historical benchmark response rate
q1=0.3 # target response rate
#### obtained from minimax Simon's two-stage design
#### http://cancer.unc.edu/biostatistics/program/ivanova/SimonsTwoStageDesign.aspx
Simonstable=matrix(NA,18,5+K+1)
Simonstable[1,1:5]=c(21,28,6,9,21.8)
Simonstable[2,1:5]=c(17,35,4,11,21.4)
Simonstable[3,1:5]=c(13,27,2,9,20)
Simonstable[4,1:5]=c(12,27,2,9,18.6)
Simonstable[5,1:5]=c(15,23,3,8,17.8)
Simonstable[6,1:5]=c(14,24,3,8,17)
Simonstable[7,1:5]=c(8,24,1,8,15.9)
Simonstable[8,1:5]=c(11,20,2,7,14.4)
Simonstable[9,1:5]=c(10,20,2,7,13.2)
Simonstable[10,1:5]=c(11,17,3,6,12)
Simonstable[11,1:5]=c(7,16,1,6,10.8)
Simonstable[12,1:5]=c(6,17,1,6,9.8)
Simonstable[13,1:5]=c(6,13,1,5,8.4)
Simonstable[14,1:5]=c(3,13,0,5,7.9)
Simonstable[15,1:5]=c(5,10,1,4,6.3)
Simonstable[16,1:5]=c(2,10,0,4,4.9)
Simonstable[17,1:5]=c(2,7,0,3,3.8)
Simonstable[18,1:5]=c(2,4,0,2,2.7)

p.sce=t(sapply(0:K,FUN=function(x){c(rep(q1,x),rep(q0,K-x))}))

Decision.ind=list()
for (x in 1:nrow(Simonstable))
{
  Ni1=Simonstable[x,1]
  Ni=Simonstable[x,2]
  r1=Simonstable[x,3]
  r2=Simonstable[x,4]
  samplesize.ind=ceiling(Simonstable[x,5])
  nik=matrix(NA,2,K) # number of patients in indication k at state i
  rik=matrix(NA,2,K) # number of responders in indication k at state i
  nik[1,]=rep(Ni1,K) # number of patients enrolled at stage 1
  Decision.ind[[samplesize.ind]]=rep(0,K+1)
  
  ###################### Simulations:
  for (scenario in 1:(nrow(p.sce)))
  {
    p0=p.sce[scenario,]
    decision.independent=numeric()
    
    tp=which(p0>=q1)
    tn=which(p0<q1)
    
    decision.ind=matrix(NA,num.sim,K)
    for (sim in 1:num.sim)
    {
      ##### Stage 1:
      rik[1,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[1,x],prob=p0[x])})

      ## Futility stop:
      stage2.stop=which(rik[1,]<=r1)
      stage2.cont=which(rik[1,]>r1)
      nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
      if (length(stage2.stop)>0) 
        {decision.ind[sim,stage2.stop]=0}
      ## Stage 2:
      if (length(stage2.cont)>0)
      {
        rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
        ri=colSums(as.matrix(rik[,stage2.cont]))
        ni=colSums(as.matrix(nik[,stage2.cont]))
        K1=length(stage2.cont)
        
        ## Futility stop:
        decision.ind[sim,stage2.cont]=ifelse(ri>r2,1,0)
      }
      decision.independent[sim]=ifelse(sum(decision.ind[sim,])>0,1,0)
      if (length(tp)>0)
      {
        Decision.ind[[samplesize.ind]][scenario]=Decision.ind[[samplesize.ind]][scenario] + ifelse((length(tp)>0)&(decision.independent[sim]==0),1,0)
      }
      if (length(tp)==0)
      {
        Decision.ind[[samplesize.ind]][scenario]=Decision.ind[[samplesize.ind]][scenario] + ifelse((length(tp)==0)&(decision.independent[sim]==1),1,0)
      }
    }
    print(Decision.ind[[samplesize.ind]][scenario])
    Simonstable[x,5+scenario]=Decision.ind[[samplesize.ind]][scenario]/num.sim
  }
}



K=6
q0=0.2
q1=0.3
#### obtained from minimax Simon's two-stage design
#### http://cancer.unc.edu/biostatistics/program/ivanova/SimonsTwoStageDesign.aspx
Simonstable=matrix(NA,17,5+(nrow(p.sce)))
Simonstable[1,1:5]=c(17,36,4,12,21.6)
Simonstable[2,1:5]=c(15,32,3,11,21)
Simonstable[3,1:5]=c(14,32,3,11,19.4)
Simonstable[4,1:5]=c(11,32,2,11,19)
Simonstable[5,1:5]=c(11,28,2,10,17.5)
Simonstable[6,1:5]=c(13,28,3,10,16.8)
Simonstable[7,1:5]=c(13,25,3,9,16)
Simonstable[8,1:5]=c(11,21,2,8,14.8)
Simonstable[9,1:5]=c(10,21,2,8,13.5)
Simonstable[10,1:5]=c(7,21,1,8,12.9)
Simonstable[11,1:5]=c(7,18,1,7,11.7)
Simonstable[12,1:5]=c(8,18,2,7,10)
Simonstable[13,1:5]=c(3,15,0,6,8.9)
Simonstable[14,1:5]=c(5,15,1,6,7.6)
Simonstable[15,1:5]=c(2,12,0,5,5.6)
Simonstable[16,1:5]=c(2,9,0,4,4.5)
Simonstable[17,1:5]=c(2,6,0,3,3.4)

p.sce=t(sapply(0:K,FUN=function(x){c(rep(q1,x),rep(q0,K-x))}))

Decision.ind=list()
for (x in 1:nrow(Simonstable))
{
  Ni1=Simonstable[x,1]
  Ni=Simonstable[x,2]
  r1=Simonstable[x,3]
  r2=Simonstable[x,4]
  samplesize.ind=ceiling(Simonstable[x,5])
  nik=matrix(NA,2,K) # number of patients in indication k at state i
  rik=matrix(NA,2,K) # number of responders in indication k at state i
  nik[1,]=rep(Ni1,K) # number of patients enrolled at stage 1
  Decision.ind[[samplesize.ind]]=rep(0,K+1)
  
  ###################### Simulations:
  for (scenario in 1:(nrow(p.sce)))
  {
    p0=p.sce[scenario,]
    decision.independent=numeric()
    
    tp=which(p0>=q1)
    tn=which(p0<q1)
    
    decision.ind=matrix(NA,num.sim,K)
    for (sim in 1:num.sim)
    {
      ##### Stage 1:
      rik[1,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[1,x],prob=p0[x])})
      
      ## Futility stop:
      stage2.stop=which(rik[1,]<=r1)
      stage2.cont=which(rik[1,]>r1)
      nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
      if (length(stage2.stop)>0) 
      {decision.ind[sim,stage2.stop]=0}
      ## Stage 2:
      if (length(stage2.cont)>0)
      {
        rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
        ri=colSums(as.matrix(rik[,stage2.cont]))
        ni=colSums(as.matrix(nik[,stage2.cont]))
        K1=length(stage2.cont)
        
        ## Futility stop:
        decision.ind[sim,stage2.cont]=ifelse(ri>r2,1,0)
      }
      decision.independent[sim]=ifelse(sum(decision.ind[sim,])>0,1,0)
      if (length(tp)>0)
      {
        Decision.ind[[samplesize.ind]][scenario]=Decision.ind[[samplesize.ind]][scenario] + ifelse((length(tp)>0)&(decision.independent[sim]==0),1,0)
      }
      if (length(tp)==0)
      {
        Decision.ind[[samplesize.ind]][scenario]=Decision.ind[[samplesize.ind]][scenario] + ifelse((length(tp)==0)&(decision.independent[sim]==1),1,0)
      }
    }
    print(c(samplesize.ind,scenario))
    print(Decision.ind[[samplesize.ind]][scenario])
    Simonstable[x,5+scenario]=Decision.ind[[samplesize.ind]][scenario]/num.sim
  }
}


