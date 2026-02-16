library(glmnet);library(bootnet);library(arules);library(foreach)
library(doParallel);library(LaplacesDemon)
#####R function for SDA-KS and SDA-CvM#####
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

s1<-c(1,6,11,12,16,17)
b1<-c(-0.4,0.6,-0.8,-0.8,1,1)
index<-1:100

#####Single index models: fixed precision matrix#####
###Setting 1###
{
  n<-400;p<-1000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out1<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est<-cbind(est1,est2)
    est
  }
}

###Setting 2###
{
  n<-200;p<-1000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out2<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est<-cbind(est1,est2)
    est
  }
}

###Setting 3###
{
  n<-400;p<-2000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out3<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est<-cbind(est1,est2)
    est
  }
}

###Setting 4###
{
  n<-200;p<-2000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out4<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est<-cbind(est1,est2)
    est
  }
}

out<-cbind(out1,out2,out3,out4)

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

#####Single index models: network precision matrix#####
generator<-ggmGenerator()
###Setting 1###
{
  n<-400;p<-1000
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out1<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','arules','glmnet')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est<-cbind(est1,est2)
    est
  }
}

###Setting 2###
{
  n<-200;p<-1000
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out2<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','arules','glmnet')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est<-cbind(est1,est2)
    est
  }
}

###Setting 3###
{
  n<-400;p<-2000
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out3<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','arules','glmnet')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est<-cbind(est1,est2)
    est
  }
}

###Setting 4###
{
  n<-200;p<-2000
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out4<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','arules','glmnet')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est<-cbind(est1,est2)
    est
  }
}

out<-cbind(out1,out2,out3,out4)

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

#####Multiple index models: fixed precision matrix#####
Relu<-function(X){
  ifelse(X>0,X,0)
}

s11<-c(1,6,11);s12<-c(2,7,12)
b11<-c(0.5,-1,0.8);b12<-c(-0.8,-0.5,1)
s21<-1;s22<-6;s23<-11;s24<-12;s25<-c(16,17)
b21<--0.5;b22<-0.8;b23<--1;b24<-1.25;b25<-c(-0.8,1)
index<-1:100
###Setting 1###
{
  n<-400;p<-1000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out1<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<-(X[,s11]%*%b11)^3+(X[,s12]%*%b12)*error
    y3<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est3<-SDA(H,y3,X,index)
    est<-cbind(est1,est2,est3)
    est
  }
}

###Setting 2###
{
  n<-200;p<-1000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out2<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<-(X[,s11]%*%b11)^3+(X[,s12]%*%b12)*error
    y3<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est3<-SDA(H,y3,X,index)
    est<-cbind(est1,est2,est3)
    est
  }
}

###Setting 3###
{
  n<-400;p<-2000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out3<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<-(X[,s11]%*%b11)^3+(X[,s12]%*%b12)*error
    y3<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est3<-SDA(H,y3,X,index)
    est<-cbind(est1,est2,est3)
    est
  }
}

###Setting 4###
{
  n<-200;p<-2000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out4<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<-(X[,s11]%*%b11)^3+(X[,s12]%*%b12)*error
    y3<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est3<-SDA(H,y3,X,index)
    est<-cbind(est1,est2,est3)
    est
  }
}

out<-cbind(out1,out2,out3,out4)

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

#####Multiple index models: network precision matrix#####
generator<-ggmGenerator()
###Setting 1###
{
  n<-400;p<-1000
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out1<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','arules','glmnet')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<-(X[,s11]%*%b11)^3+(X[,s12]%*%b12)*error
    y3<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est3<-SDA(H,y3,X,index)
    est<-cbind(est1,est2,est3)
    est
  }
}

###Setting 2###
{
  n<-200;p<-1000
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out2<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','arules','glmnet')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<-(X[,s11]%*%b11)^3+(X[,s12]%*%b12)*error
    y3<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est3<-SDA(H,y3,X,index)
    est<-cbind(est1,est2,est3)
    est
  }
}

###Setting 3###
{
  n<-400;p<-2000
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out3<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','arules','glmnet')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<-(X[,s11]%*%b11)^3+(X[,s12]%*%b12)*error
    y3<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est3<-SDA(H,y3,X,index)
    est<-cbind(est1,est2,est3)
    est
  }
}

###Setting 4###
{
  n<-200;p<-2000
  H<-ceiling(n^(1/4))
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out4<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','arules','glmnet')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<-(X[,s11]%*%b11)^3+(X[,s12]%*%b12)*error
    y3<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est3<-SDA(H,y3,X,index)
    est<-cbind(est1,est2,est3)
    est
  }
}

out<-cbind(out1,out2,out3,out4)

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)
