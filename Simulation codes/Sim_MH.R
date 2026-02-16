library(arules);library(LaplacesDemon);library(glmnet);library(knockoff)
library(bootnet);library(foreach);library(doParallel)
#####R function for Multiple hypothesis testing#####
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

#####R function for empirical power and FPR#####
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
generator<-ggmGenerator()

registerDoParallel(35)
out<-foreach (i=1:100, .combine=rbind, .packages=c('bootnet','glmnet','arules')) %dopar% {
  set.seed(i)
  trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
  X<-generator(n,trueNet)
  error<-rnorm(n)
  y11<-X[,s1]%*%b1+error
  y12<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+0.5*error
  y13<-(3*(X[,s11]%*%b11))/(0.5+(1.5+X[,s12]%*%b12)^2)+0.5*error
  y21<-X[,s2]%*%b2+error
  y22<-sin(X[,s2]%*%b2)*exp(X[,s2]%*%b2)+0.5*error
  y23<-(3*(X[,s21]%*%b21))/(0.5+(1.5+X[,s22]%*%b22)^2)+0.5*error
  
  ###Selection,s=6###
  Betac<-rep(0,p);Betac[s1]<-1
  m11<-SDA_knock(y11,X,5)
  m11_2<-knockoff.filter(X,y11)
  m11_3<-knockoff.filter(X,y11,offset=0)
  out1_2<-data.frame(FDP(Betac,m11_2$selected),
                     power(Betac,m11_2$selected,sum(Betac)))
  out1_3<-data.frame(FDP(Betac,m11_3$selected),
                     power(Betac,m11_3$selected,sum(Betac)))
  est1<-as.vector(data.frame(Output(m11[[1]],Betac),Output(m11[[2]],Betac),
                             out1_2,out1_3))
  
  Betac<-rep(0,p);Betac[s1]<-1
  m12<-SDA_knock(y12,X,5)
  m12_2<-knockoff.filter(X,y12)
  m12_3<-knockoff.filter(X,y12,offset=0)
  out2_2<-data.frame(FDP(Betac,m12_2$selected),
                     power(Betac,m12_2$selected,sum(Betac)))
  out2_3<-data.frame(FDP(Betac,m12_3$selected),
                     power(Betac,m12_3$selected,sum(Betac)))
  est2<-as.vector(data.frame(Output(m12[[1]],Betac),Output(m12[[2]],Betac),
                             out2_2,out2_3))
  
  Betac<-rep(0,p);Betac[s1]<-1
  m13<-SDA_knock(y13,X,5)
  m13_2<-knockoff.filter(X,y13)
  m13_3<-knockoff.filter(X,y13,offset=0)
  out3_2<-data.frame(FDP(Betac,m13_2$selected),
                     power(Betac,m13_2$selected,sum(Betac)))
  out3_3<-data.frame(FDP(Betac,m13_3$selected),
                     power(Betac,m13_3$selected,sum(Betac)))
  est3<-as.vector(data.frame(Output(m13[[1]],Betac),Output(m13[[2]],Betac),
                             out3_2,out3_3))
  
  ###Selection,s=12###
  Betac<-rep(0,p);Betac[s2]<-1
  m21<-SDA_knock(y21,X,5)
  m21_2<-knockoff.filter(X,y21)
  m21_3<-knockoff.filter(X,y21,offset=0)
  out4_2<-data.frame(FDP(Betac,m21_2$selected),
                     power(Betac,m21_2$selected,sum(Betac)))
  out4_3<-data.frame(FDP(Betac,m21_3$selected),
                     power(Betac,m21_3$selected,sum(Betac)))
  est4<-as.vector(data.frame(Output(m21[[1]],Betac),Output(m21[[2]],Betac),
                             out4_2,out4_3))
  
  Betac<-rep(0,p);Betac[s2]<-1
  m22<-SDA_knock(y22,X,5)
  m22_2<-knockoff.filter(X,y22)
  m22_3<-knockoff.filter(X,y22,offset=0)
  out5_2<-data.frame(FDP(Betac,m22_2$selected),
                     power(Betac,m22_2$selected,sum(Betac)))
  out5_3<-data.frame(FDP(Betac,m22_3$selected),
                     power(Betac,m22_3$selected,sum(Betac)))
  est5<-as.vector(data.frame(Output(m22[[1]],Betac),Output(m22[[2]],Betac),
                             out5_2,out5_3))
  
  Betac<-rep(0,p);Betac[s2]<-1
  m23<-SDA_knock(y23,X,5)
  m23_2<-knockoff.filter(X,y23)
  m23_3<-knockoff.filter(X,y23,offset=0)
  out6_2<-data.frame(FDP(Betac,m23_2$selected),
                     power(Betac,m23_2$selected,sum(Betac)))
  out6_3<-data.frame(FDP(Betac,m23_3$selected),
                     power(Betac,m23_3$selected,sum(Betac)))
  est6<-as.vector(data.frame(Output(m23[[1]],Betac),Output(m23[[2]],Betac),
                             out6_2,out6_3))
  
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

