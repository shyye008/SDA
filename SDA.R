library(glmnet);library(bootnet);library(arules);library(foreach)
library(doParallel)
#####R function#####
###Single hypothesis testing###
SDA<-function(H,y,X,index,sig.level=0.05,L=1000){
  ts.func<-function(x){
    max(abs(x/sqrt(est.se)))
  }
  ts.func2<-function(x){
    mean(abs(x/sqrt(est.se)))
  }
  
  n<-length(y)
  out<-c()
  for (j in index) {
    X_j<-X[,j];X_c<-X[,-j]
    fit<-cv.glmnet(X_c,X_j)
    pred<-c(predict(fit,X_c,s="lambda.min"))
    KS<-c();CV<-c();CvM<-c();CV2<-c()
    for (i in H) {
      slice<-discretize(y,breaks=i,labels=F)
      H.matrix<-matrix(0,nrow=length(slice),ncol=i)
      H.matrix[cbind(seq_along(slice),slice)]<-1
      m.Z<-kronecker(t(rep(1,i)),(X_j-pred))
      nu<-colMeans(H.matrix*m.Z)
      
      m.nu<-kronecker(rep(1,n),t(nu))
      inf.fun<-H.matrix*m.Z-m.nu
      est.se<-diag(t(inf.fun)%*%inf.fun/n)
      KS<-c(KS,max(abs(sqrt(n)*nu/sqrt(est.se))))
      CvM<-c(CvM,mean(abs(sqrt(n)*nu/sqrt(est.se))))
      
      phi<-c()
      for (l in 1:L) {
        U<-matrix(rnorm(n*i),nrow=n)
        phi_l<-colSums(U*inf.fun)/sqrt(n)
        phi<-cbind(phi,phi_l)
      }
      
      sim_max<-apply(phi,2,ts.func)
      sim_mean<-apply(phi,2,ts.func2)
      CV<-c(CV,quantile(sim_max,probs=1-sig.level))
      CV2<-c(CV2,quantile(sim_mean,probs=1-sig.level))
    }
    
    out<-cbind(out,rbind(KS,CV,CvM,CV2))
  }
  out
}

SDA_chi<-function(H,y,X,index,sig.level=0.05,L=1000){
  n<-length(y)
  out<-c()
  for (j in index) {
    X_j<-X[,j];X_c<-X[,-j]
    fit<-cv.glmnet(X_c,X_j)
    pred<-c(predict(fit,X_c,s="lambda.min"))
    pval<-c()
    for (i in H) {
      slice<-discretize(y,breaks=i,labels=F)
      H.matrix<-matrix(0,nrow=length(slice),ncol=i)
      H.matrix[cbind(seq_along(slice),slice)]<-1
      m.Z<-kronecker(t(rep(1,i)),(X_j-pred))
      nu<-colMeans(H.matrix*m.Z)
      
      m.nu<-kronecker(rep(1,n),t(nu))
      inf.fun<-H.matrix*m.Z-m.nu
      est.se<-t(inf.fun)%*%inf.fun/n
      
      W<-n*t(nu)%*%solve(est.se)%*%nu
      
      pval<-c(pval,pchisq(W,i,lower.tail=F))
    }
    
    out<-c(out,pval)
  }
  out
}

###Mutiple hypothesis testing###
SDA_knock<-function(y,X,H,index=c(1:p)){
  n<-length(y)
  slice<-discretize(y,breaks=H,labels=F)
  H.matrix<-model.matrix(y~as.factor(slice)-1)
  SDA_CvM1<-c();SDA_CvM2<-c()
  for (i in index) {
    X_i<-X[,i]
    X_c<-X[,-i]
    fit<-cv.glmnet(X_c,X_i)
    pred<-c(predict(fit,X_c,s="lambda.min"))
    Z_i<-X_i-pred
    
    s_i<-sd(Z_i)
    Z_Knock<-rnorm(n,0,s_i)
    
    m.Z<-kronecker(t(rep(1,H)),Z_i)
    nu<-colMeans(H.matrix*m.Z)
    
    m.nu<-kronecker(rep(1,n),t(nu))
    inf.fun<-H.matrix*m.Z-m.nu
    est.se<-diag(t(inf.fun)%*%inf.fun/n)
    CvM<-mean(abs(sqrt(n)*nu/sqrt(est.se)))
    
    m.Z.Knock<-kronecker(t(rep(1,H)),Z_Knock)
    nu.Knock<-colMeans(H.matrix*m.Z.Knock)
    
    m.nu.Knock<-kronecker(rep(1,n),t(nu.Knock))
    inf.fun.Knock<-H.matrix*m.Z.Knock-m.nu.Knock
    est.se.Knock<-diag(t(inf.fun.Knock)%*%inf.fun.Knock/n)
    CvM.Knock<-mean(abs(sqrt(n)*nu.Knock/sqrt(est.se.Knock)))
    SDA_CvM1<-c(SDA_CvM1,CvM-CvM.Knock)
    SDA_CvM2<-c(SDA_CvM2,sign(CvM-CvM.Knock)*max(CvM,CvM.Knock))
  }
  list(SDA_CvM1,SDA_CvM2)
}

###Empirical power and FPR###
Output<-function(Stat,beta,q=0.1){
  FDP_hat<-function(Stat,t){
    sum(ifelse(Stat<=-t,1,0))/sum(ifelse(Stat>t,1,0))
  }
  
  FDP<-function(beta,select){
    sum(beta[select]==0)/length(select)
  }
  
  power<-function(beta,select,size){
    sum(beta[select]!=0)/size
  }
  
  p<-length(Stat)
  candidate<-seq(0,max(Stat),length.out=(p+1))[-(p+1)]
  FDP_all<-c()
  for (i in candidate) {
    FDP_all<-c(FDP_all,FDP_hat(Stat,i))
  }
  tau<-candidate[min(which((FDP_all<=q)==T))]
  select<-(1:p)[Stat>tau]
  fdp<-FDP(beta,select)
  power<-power(beta,select,sum(beta))
  data.frame(fdp,power)
}


