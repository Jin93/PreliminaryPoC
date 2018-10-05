Decision.bayes=list()
for (samplesize.bayes in 4:22) ### assuming that the sample size is between 4 and 22 for each indication group.
{
  Decision.bayes[[samplesize.bayes]]=rep(0,K+1)
  #################### Calibrate for Bayesian method to control overall type I error rate:
  nik.bayes=rep(samplesize.bayes,K) # number of patients in indication k
  decision.bayes=Tstat=numeric()
  p0=rep(q0,K) # set all true response rates to be equal to standard of care response rate
  for (sim in 1:2000)
  {
    rik.bayes=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik.bayes[x],prob=p0[x])})
    q=rep(log(((q0+q1)/2)/(1-(q0+q1)/2)),K) ## can consider different settings based on histological data
    ############ Jags model for BHM:
    jags.data <- list("n"=nik.bayes, "Y"=rik.bayes, "K"=K, "q"=q,"g1"=log(q1/(1-q1))-log((q1+q0)/2/(1-(q1+q0)/2)),"g0"=-log((q1+q0)/2/(1-(q1+q0)/2))+log(q0/(1-q0)))
    jags.fit <- jags.model(file = "C:/Users/e0359820/Desktop/Jin/Signal Detection Project/a-ina.txt",data = jags.data,
                           n.adapt=1000,n.chains=1,quiet=T)
    update(jags.fit, 4000)
    bayes.out <- coda.samples(jags.fit,variable.names = c("p","d","pi","delta","tausq","mu1","mu2"),n.iter=10000)
    #mcmc1=mcmc(data=bayes.out[[1]],start=1,end=nrow(bayes.out[[1]]))
    #mcmc2=mcmc(data=bayes.out[[2]],start=1,end=nrow(bayes.out[[2]]))
    #mcmc=mcmc.list(mcmc1,mcmc2)
    #gelman.diag(mcmc) # 1.1
    
    ### Interim analysis:
    Tstat[sim]=sum(apply(bayes.out[[1]][,sapply(1:K,FUN=function(x){paste("d[",x,"]",sep="")})],1,sum)>0)/nrow(bayes.out[[1]])
    print(sim)
  }
  c.bayes=quantile(Tstat,1-alpha) # probability cutoff for the final decision that is calibrated to control the type I error rate to be at alpha
  
  ###################### Simulations:
  for (scenario in 2:(nrow(p.sce))) # different scenario: different number of truely active indications
  {
    p0=p.sce[scenario,] # set true response rates
    decision.bayes=numeric()
    tp=which(p0>=q1) # indicators of the true positives
    tn=which(p0<q1) # indicators of the true negatives
    Tstat=numeric() # posterior probability for the final decision
    
    for (sim in 1:num.sim)
    {
      #################### Bayesian Method: ####################
      rik.bayes=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik.bayes[x],prob=p0[x])})
      q=rep(log((q0+q1)/2/(1-(q0+q1)/2)),K) ## can consider different settings based on histological data
      ############ Jags model for BHM:
      jags.data <- list("n"=nik.bayes, "Y"=rik.bayes, "K"=K, "q"=q,"g1"=log(q1/(1-q1))-log((q1+q0)/2/(1-(q1+q0)/2)),"g0"=-log((q1+q0)/2/(1-(q1+q0)/2))+log(q0/(1-q0)))
      jags.fit <- jags.model(file = "C:/Users/e0359820/Desktop/Jin/Signal Detection Project/a-ina.txt",data = jags.data,
                             n.adapt=1000,n.chains=1,quiet=T)
      update(jags.fit, 4000)
      bayes.out <- coda.samples(jags.fit,variable.names = c("p","d","pi","delta","tausq","mu1","mu2","tausq2"),n.iter=10000)

      ### Final decision:
      Tstat[sim]=sum(apply(bayes.out[[1]][,sapply(1:K,FUN=function(x){paste("d[",x,"]",sep="")})],1,sum)>0)/nrow(bayes.out[[1]])
      decision.bayes[sim]=ifelse(Tstat[sim]>c.bayes,1,0)
      Decision.bayes[[samplesize.bayes]][scenario-1]=Decision.bayes[[samplesize.bayes]][scenario-1] + ifelse((length(tp)>0)&(decision.bayes[sim]==0),1,0)
      print(sim)
      print(Decision.bayes[[samplesize.bayes]][scenario-1])
    }
  }
}
