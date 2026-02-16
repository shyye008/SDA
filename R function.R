library(glmnet);library(bootnet);library(arules)
#####R function#####
###Single hypothesis testing###
SDA<-function(H,y,X,index,method="Chi",L=1000){
  ##H: number of slices
  ##y: outcome variable
  ##X: covariate matrix
  ##index: index of variable(s) being tested
  ##method: test statistic, can be either "Chi", "CvM", or "KS"
  ##L: number of bootstrap replicates, only for method="CvM" or "KS"
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
      
      if(method=="Chi"){
        est.se<-t(inf.fun)%*%inf.fun/n
        W<-n*t(nu)%*%solve(est.se)%*%nu
        pval<-c(pval,pchisq(W,i,lower.tail=F))
      }else if(method=="CvM"){
        est.se<-diag(t(inf.fun)%*%inf.fun/n)
        ts.func<-function(x){
          mean(abs(x/sqrt(est.se)))
        }
        CvM<-mean(abs(sqrt(n)*nu/sqrt(est.se)))
        
        phi<-c()
        for (l in 1:L) {
          U<-matrix(rnorm(n*i),nrow=n)
          phi_l<-colSums(U*inf.fun)/sqrt(n)
          phi<-cbind(phi,phi_l)
        }
        
        sim_mean<-apply(phi,2,ts.func)
        pval<-c(pval,sum(sim_mean>CvM)/L)
      } else if(method=="KS"){
        est.se<-diag(t(inf.fun)%*%inf.fun/n)
        ts.func<-function(x){
          max(abs(x/sqrt(est.se)))
        }
        KS<-max(abs(sqrt(n)*nu/sqrt(est.se)))
        
        phi<-c()
        for (l in 1:L) {
          U<-matrix(rnorm(n*i),nrow=n)
          phi_l<-colSums(U*inf.fun)/sqrt(n)
          phi<-cbind(phi,phi_l)
        }
        
        sim_max<-apply(phi,2,ts.func)
        pval<-c(pval,sum(sim_max>KS)/L)
      }
    }
    ##out: a vector of p-values
    out<-c(out,pval)
  }
  out
}

###Multiple hypothesis testing###
SDA_knock<-function(y,X,H,method="Chi",feature="Diff",index=c(1:p),q=0.1){
  ##H: number of slices
  ##y: outcome variable
  ##X: covariate matrix
  ##method: test statistic, can be either "Chi", "CvM", or "KS"
  ##feature: feature statistic, can be either "Diff" or "SM"
  ##index: index of variable(s) being tested
  ##q: FDR
  n<-length(y)
  slice<-discretize(y,breaks=H,labels=F)
  H.matrix<-model.matrix(y~as.factor(slice)-1)
  Stat=c()
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
    
    m.Z.Knock<-kronecker(t(rep(1,H)),Z_Knock)
    nu.Knock<-colMeans(H.matrix*m.Z.Knock)
    m.nu.Knock<-kronecker(rep(1,n),t(nu.Knock))
    inf.fun.Knock<-H.matrix*m.Z.Knock-m.nu.Knock
    
    if(method=="CvM"){
      est.se<-diag(t(inf.fun)%*%inf.fun/n)
      est.se.Knock<-diag(t(inf.fun.Knock)%*%inf.fun.Knock/n)
      
      TS<-mean(abs(sqrt(n)*nu/sqrt(est.se)))
      TS.Knock<-mean(abs(sqrt(n)*nu.Knock/sqrt(est.se.Knock)))
    }else if(method=="KS"){
      est.se<-diag(t(inf.fun)%*%inf.fun/n)
      est.se.Knock<-diag(t(inf.fun.Knock)%*%inf.fun.Knock/n)
      
      TS<-max(abs(sqrt(n)*nu/sqrt(est.se)))
      TS.Knock<-max(abs(sqrt(n)*nu.Knock/sqrt(est.se.Knock)))
    }else if(method=="Chi"){
      est.se<-t(inf.fun)%*%inf.fun/n
      est.se.Knock<-t(inf.fun.Knock)%*%inf.fun.Knock/n
      
      TS<-n*t(nu)%*%solve(est.se)%*%nu
      TS.Knock<-n*t(nu.Knock)%*%solve(est.se.Knock)%*%nu.Knock
    }
    
    if(feature=="Diff"){
      Stat<-c(Stat,TS-TS.Knock)
    }else if(feature=="SM"){
      Stat<-c(Stat,sign(TS-TS.Knock)*max(TS,TS.Knock))
    }
  }
  
  FDP_hat<-function(Stat,t){
    sum(ifelse(Stat<=-t,1,0))/sum(ifelse(Stat>t,1,0))
  }
  
  p<-length(Stat)
  candidate<-seq(0,max(Stat),length.out=(p+1))[-(p+1)]
  FDP_all<-c()
  for (i in candidate) {
    FDP_all<-c(FDP_all,FDP_hat(Stat,i))
  }
  tau<-candidate[min(which((FDP_all<=q)==T))]
  select<-(1:p)[Stat>tau]
  
  ##Stat: a vector of feature statistics
  ##select: a vector of indexs for the variables being selected
  list(Stat,select)
}
