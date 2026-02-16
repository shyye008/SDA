library(glmnet);library(bootnet);library(arules);library(foreach)
library(doParallel);library(LaplacesDemon)
#####R function for SDA-Chi-LASSO#####
SDA<-function(H,y,X,index,sig.level=0.05,L=1000){
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

#####Simulation for LASSO#####
###Setting 1###
{
  s1<-c(1,6:7,11:13,16:19,21:25)
  b1<-rep(0.8,15)
  index<-1:25
  n<-400;p<-1000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  H<-5
  
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
    est<-c(est1,est2)
    est
  }
}

###Setting 2###
{
  s1<-c(1,11:12,21:23,31:34,41:45)
  b1<-rep(0.8,15)
  index<-1:50
  n<-400;p<-1000
  Omega<-matrix(rep(0.5,10*10),nrow=10)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/10),Omega)
  H<-5
  
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
    est<-c(est1,est2)
    est
  }
}

###Setting 3###
{
  s1<-c(1,11:12,21:23,31:34,41:45)
  b1<-rep(1,15)
  index<-1:50
  n<-400;p<-1000
  Omega<-matrix(rep(0.5,10*10),nrow=10)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/10),Omega)
  H<-5
  
  out3<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est<-c(est1,est2)
    est
  }
}

out<-cbind(out1,out2,out3)

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

#####R function for SDA-Chi-LASSO+SIS#####
SDA<-function(H,y,X,index,sig.level=0.05,L=1000){
  n<-length(y)
  out<-c()
  for (j in index) {
    X_j<-X[,j];X_c<-X[,-j]
    corr<-cor(X_j,X_c)
    selen<-tail(order(abs(corr)),floor(n/log(n)))
    X_sele<-X_c[,selen]
    fit<-cv.glmnet(X_sele,X_j)
    pred<-c(predict(fit,X_sele,s="lambda.min"))
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

#####Simulation for LASSO+SIS#####
###Setting 1###
{
  s1<-c(1,6:7,11:13,16:19,21:25)
  b1<-rep(0.8,15)
  index<-1:25
  n<-400;p<-1000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  H<-5
  
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
    est<-c(est1,est2)
    est
  }
}

###Setting 2###
{
  s1<-c(1,11:12,21:23,31:34,41:45)
  b1<-rep(0.8,15)
  index<-1:50
  n<-400;p<-1000
  Omega<-matrix(rep(0.5,10*10),nrow=10)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/10),Omega)
  H<-5
  
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
    est<-c(est1,est2)
    est
  }
}

###Setting 3###
{
  s1<-c(1,11:12,21:23,31:34,41:45)
  b1<-rep(1,15)
  index<-1:50
  n<-400;p<-1000
  Omega<-matrix(rep(0.5,10*10),nrow=10)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/10),Omega)
  H<-5
  
  out3<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est<-c(est1,est2)
    est
  }
}

out<-cbind(out1,out2,out3)

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

#####R function for SDA-Chi-LS#####
SDA<-function(H,y,X,index,nei){
  n<-length(y)
  out<-c();k<-1
  for (j in index) {
    X_j<-X[,j];X_c<-X[,nei[[k]]]
    k<-k+1
    datj<-data.frame(X_j,X_c)
    fit<-lm(X_j~.,dat=datj)
    pred<-c(predict(fit))
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

##### Simulation for oracle LS #####
###Setting 1###
{
  s1<-c(1,6:7,11:13,16:19,21:25)
  b1<-rep(0.8,15)
  index<-1:25
  n<-400;p<-1000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  H<-5
  nei<-vector(mode="list",length=25)
  for (i in index) {
    u<-trunc((i-1)/5)
    clusteri<-(u*5+1):(u*5+5)
    ri<-i-u*5
    neii<-clusteri[-ri]
    nei[[i]]<-neii
  }
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out1<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est1<-SDA(H,y1,X,index,nei)
    est2<-SDA(H,y2,X,index,nei)
    est<-c(est1,est2)
    est
  }
}

###Setting 2###
{
  s1<-c(1,11:12,21:23,31:34,41:45)
  b1<-rep(0.8,15)
  index<-1:50
  n<-400;p<-1000
  Omega<-matrix(rep(0.5,10*10),nrow=10)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/10),Omega)
  H<-5
  nei<-vector(mode="list",length=50)
  for (i in index) {
    u<-trunc((i-1)/10)
    clusteri<-(u*10+1):(u*10+10)
    ri<-i-u*10
    neii<-clusteri[-ri]
    nei[[i]]<-neii
  }
  
  library(foreach)
  library(doParallel)
  registerDoParallel(35)
  out2<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est1<-SDA(H,y1,X,index,nei)
    est2<-SDA(H,y2,X,index,nei)
    est<-c(est1,est2)
    est
  }
}

###Setting 3###
{
  s1<-c(1,11:12,21:23,31:34,41:45)
  b1<-rep(1,15)
  index<-1:50
  n<-400;p<-1000
  Omega<-matrix(rep(0.5,10*10),nrow=10)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/10),Omega)
  H<-5
  
  out3<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est1<-SDA(H,y1,X,index,nei)
    est2<-SDA(H,y2,X,index,nei)
    est<-c(est1,est2)
    est
  }
}

out<-cbind(out1,out2,out3)

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)


