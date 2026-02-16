library(glmnet);library(bootnet);library(arules);library(LaplacesDemon)
library(foreach);library(doParallel)
#####R function for multiple hypothesis testing#####
SDA_knock<-function(y,X,H,method="CvM",feature="SM",index=c(1:p),q=0.1){
  ##method can be either "Chi", "CvM", or "KS"
  ##feature can be either "Diff" or "SM"
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
  
  list(Stat,select)
}

#####R function for Empirical power and FPR#####
FDP<-function(beta,select){
  sum(beta[select]==0)/length(select)
}

power<-function(beta,select,size){
  sum(beta[select]!=0)/size
}

#####Data generation#####
n<-400;p<-1000

s1<-1:6;b1<-rep(1,6)
s11<-1:3;s12<-4:6;b11<-rep(1,3);b12<-rep(-1,3)
s2<-1:12;b2<-rep(1,12)
s21<-1:6;s22<-7:12;b21<-rep(1,6);b22<-rep(-1,6)

Omega<-matrix(rep(0.5,5*5),nrow=5)
diag(Omega)<-1
Omega<-kronecker(diag(p/5),Omega)
H<-5

registerDoParallel(35)
###Setting 1###
out<-foreach (i=1:100, .combine=rbind, .packages=c('LaplacesDemon','glmnet','arules')) %dopar% {
  set.seed(i)
  X<-rmvtp(n,rep(0,p),Omega,5)
  error<-rnorm(n)
  y11<-X[,s1]%*%b1+error
  y12<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+0.5*error
  y13<-(3*(X[,s11]%*%b11))/(0.5+(1.5+X[,s12]%*%b12)^2)+0.5*error
  y21<-X[,s2]%*%b2+error
  y22<-sin(X[,s2]%*%b2)*exp(X[,s2]%*%b2)+0.5*error
  y23<-(3*(X[,s21]%*%b21))/(0.5+(1.5+X[,s22]%*%b22)^2)+0.5*error
  
  ###Selection,s=6###
  Betac<-rep(0,p);Betac[s1]<-1
  m11<-SDA_knock(y11,X,ceiling(n^(1/4)))
  est1<-c(FDP(Betac,m11[[2]]),power(Betac,m11[[2]],6))
  
  Betac<-rep(0,p);Betac[s1]<-1
  m12<-SDA_knock(y12,X,ceiling(n^(1/4)))
  est2<-c(FDP(Betac,m12[[2]]),power(Betac,m12[[2]],6))
  
  Betac<-rep(0,p);Betac[s1]<-1
  m13<-SDA_knock(y13,X,ceiling(n^(1/4)))
  est3<-c(FDP(Betac,m13[[2]]),power(Betac,m13[[2]],6))
  
  ###Selection,s=12###
  Betac<-rep(0,p);Betac[s2]<-1
  m21<-SDA_knock(y21,X,ceiling(n^(1/4)))
  est4<-c(FDP(Betac,m21[[2]]),power(Betac,m21[[2]],12))
  
  Betac<-rep(0,p);Betac[s2]<-1
  m22<-SDA_knock(y22,X,ceiling(n^(1/4)))
  est5<-c(FDP(Betac,m22[[2]]),power(Betac,m22[[2]],12))
  
  Betac<-rep(0,p);Betac[s2]<-1
  m23<-SDA_knock(y23,X,ceiling(n^(1/4)))
  est6<-c(FDP(Betac,m12[[2]]),power(Betac,m23[[2]],12))
  
  est<-rbind(est1,est2,est3,est4,est5,est6)
  est
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

###Setting 2###
out<-foreach (i=1:100, .combine=rbind, .packages=c('LaplacesDemon','glmnet','arules')) %dopar% {
  set.seed(i)
  X<-rmvtp(n,rep(0,p),Omega,3)
  error<-rnorm(n)
  y11<-X[,s1]%*%b1+error
  y12<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+0.5*error
  y13<-(3*(X[,s11]%*%b11))/(0.5+(1.5+X[,s12]%*%b12)^2)+0.5*error
  y21<-X[,s2]%*%b2+error
  y22<-sin(X[,s2]%*%b2)*exp(X[,s2]%*%b2)+0.5*error
  y23<-(3*(X[,s21]%*%b21))/(0.5+(1.5+X[,s22]%*%b22)^2)+0.5*error
  
  ###Selection,s=6###
  Betac<-rep(0,p);Betac[s1]<-1
  m11<-SDA_knock(y11,X,ceiling(n^(1/4)))
  est1<-c(FDP(Betac,m11[[2]]),power(Betac,m11[[2]],6))
  
  Betac<-rep(0,p);Betac[s1]<-1
  m12<-SDA_knock(y12,X,ceiling(n^(1/4)))
  est2<-c(FDP(Betac,m12[[2]]),power(Betac,m12[[2]],6))
  
  Betac<-rep(0,p);Betac[s1]<-1
  m13<-SDA_knock(y13,X,ceiling(n^(1/4)))
  est3<-c(FDP(Betac,m13[[2]]),power(Betac,m13[[2]],6))
  
  ###Selection,s=12###
  Betac<-rep(0,p);Betac[s2]<-1
  m21<-SDA_knock(y21,X,ceiling(n^(1/4)))
  est4<-c(FDP(Betac,m21[[2]]),power(Betac,m21[[2]],12))
  
  Betac<-rep(0,p);Betac[s2]<-1
  m22<-SDA_knock(y22,X,ceiling(n^(1/4)))
  est5<-c(FDP(Betac,m22[[2]]),power(Betac,m22[[2]],12))
  
  Betac<-rep(0,p);Betac[s2]<-1
  m23<-SDA_knock(y23,X,ceiling(n^(1/4)))
  est6<-c(FDP(Betac,m12[[2]]),power(Betac,m23[[2]],12))
  
  est<-rbind(est1,est2,est3,est4,est5,est6)
  est
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

###Setting 3###
out<-foreach (i=1:100, .combine=rbind, .packages=c('LaplacesDemon','glmnet','arules')) %dopar% {
  set.seed(i)
  Xn<-rmvnp(n,rep(0,p),Omega)
  X<-c()
  for (l in 1:p) {
    pl<-pnorm(Xn[,l],mean=0,sd=std[l])
    Xl<-qt(pl,df=5)
    X<-cbind(X,Xl)
  }
  error<-rnorm(n)
  y11<-X[,s1]%*%b1+error
  y12<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+0.5*error
  y13<-(3*(X[,s11]%*%b11))/(0.5+(1.5+X[,s12]%*%b12)^2)+0.5*error
  y21<-X[,s2]%*%b2+error
  y22<-sin(X[,s2]%*%b2)*exp(X[,s2]%*%b2)+0.5*error
  y23<-(3*(X[,s21]%*%b21))/(0.5+(1.5+X[,s22]%*%b22)^2)+0.5*error
  
  ###Selection,s=6###
  Betac<-rep(0,p);Betac[s1]<-1
  m11<-SDA_knock(y11,X,ceiling(n^(1/4)))
  est1<-c(FDP(Betac,m11[[2]]),power(Betac,m11[[2]],6))
  
  Betac<-rep(0,p);Betac[s1]<-1
  m12<-SDA_knock(y12,X,ceiling(n^(1/4)))
  est2<-c(FDP(Betac,m12[[2]]),power(Betac,m12[[2]],6))
  
  Betac<-rep(0,p);Betac[s1]<-1
  m13<-SDA_knock(y13,X,ceiling(n^(1/4)))
  est3<-c(FDP(Betac,m13[[2]]),power(Betac,m13[[2]],6))
  
  ###Selection,s=12###
  Betac<-rep(0,p);Betac[s2]<-1
  m21<-SDA_knock(y21,X,ceiling(n^(1/4)))
  est4<-c(FDP(Betac,m21[[2]]),power(Betac,m21[[2]],12))
  
  Betac<-rep(0,p);Betac[s2]<-1
  m22<-SDA_knock(y22,X,ceiling(n^(1/4)))
  est5<-c(FDP(Betac,m22[[2]]),power(Betac,m22[[2]],12))
  
  Betac<-rep(0,p);Betac[s2]<-1
  m23<-SDA_knock(y23,X,ceiling(n^(1/4)))
  est6<-c(FDP(Betac,m12[[2]]),power(Betac,m23[[2]],12))
  
  est<-c(est1,est2,est3,est4,est5,est6)
  est
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

###Setting 4###
out<-foreach (i=1:100, .combine=rbind, .packages=c('LaplacesDemon','glmnet','arules')) %dopar% {
  set.seed(i)
  Xn<-rmvnp(n,rep(0,p),Omega)
  X<-c()
  for (l in 1:p) {
    pl<-pnorm(Xn[,l],mean=0,sd=std[l])
    Xl<-qchisq(pl,df=5)
    X<-cbind(X,Xl)
  }
  error<-rnorm(n)
  y11<-X[,s1]%*%b1+error
  y12<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+0.5*error
  y13<-(3*(X[,s11]%*%b11))/(0.5+(1.5+X[,s12]%*%b12)^2)+0.5*error
  y21<-X[,s2]%*%b2+error
  y22<-sin(X[,s2]%*%b2)*exp(X[,s2]%*%b2)+0.5*error
  y23<-(3*(X[,s21]%*%b21))/(0.5+(1.5+X[,s22]%*%b22)^2)+0.5*error
  
  ###Selection,s=6###
  Betac<-rep(0,p);Betac[s1]<-1
  m11<-SDA_knock(y11,X,ceiling(n^(1/4)))
  est1<-c(FDP(Betac,m11[[2]]),power(Betac,m11[[2]],6))
  
  Betac<-rep(0,p);Betac[s1]<-1
  m12<-SDA_knock(y12,X,ceiling(n^(1/4)))
  est2<-c(FDP(Betac,m12[[2]]),power(Betac,m12[[2]],6))
  
  Betac<-rep(0,p);Betac[s1]<-1
  m13<-SDA_knock(y13,X,ceiling(n^(1/4)))
  est3<-c(FDP(Betac,m13[[2]]),power(Betac,m13[[2]],6))
  
  ###Selection,s=12###
  Betac<-rep(0,p);Betac[s2]<-1
  m21<-SDA_knock(y21,X,ceiling(n^(1/4)))
  est4<-c(FDP(Betac,m21[[2]]),power(Betac,m21[[2]],12))
  
  Betac<-rep(0,p);Betac[s2]<-1
  m22<-SDA_knock(y22,X,ceiling(n^(1/4)))
  est5<-c(FDP(Betac,m22[[2]]),power(Betac,m22[[2]],12))
  
  Betac<-rep(0,p);Betac[s2]<-1
  m23<-SDA_knock(y23,X,ceiling(n^(1/4)))
  est6<-c(FDP(Betac,m12[[2]]),power(Betac,m23[[2]],12))
  
  est<-rbind(est1,est2,est3,est4,est5,est6)
  est
}

stopImplicitCluster()

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)
