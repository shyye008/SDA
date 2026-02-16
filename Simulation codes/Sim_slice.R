library(glmnet);library(bootnet);library(arules);library(foreach)
library(doParallel)
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

Relu<-function(X){
  ifelse(X>0,X,0)
}

###Setting 1###
{
  n<-200;p<-1000;H<-2:20
  b1<-c(1,0.5);b2<-c(1,0.5);b3_1<-c(1,-0.5);b3_2<-c(-1,0.5)
  generator<-ggmGenerator()
  registerDoParallel(35)
  out1<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','arules','glmnet')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-X[,1:2]%*%b1+error
    y2<-sin(X[,1:2]%*%b2)*exp(X[,1:2]%*%b2)+error
    y3<-3*(X[,1:2]%*%b3_1)/(0.5+(1.5+X[,3:4]%*%b3_2)^2)+error
    y4<--2+Relu(-0.5*X[,1])+Relu(0.5*X[,2])+Relu(0.75*X[,3])+
      Relu(X[,4])+Relu(1.25*X[,5])+error
    
    est1<-SDA(H,y1,X,c(1,2,3,20))
    est2<-SDA(H,y2,X,c(1,2,3,20))
    est3<-SDA(H,y3,X,c(1,2,3,4,5,20))
    est4<-SDA(H,y4,X,c(1,2,3,4,5,6,20))
    est<-cbind(est1,est2,est3,est4)
    est
  }
}

###Setting 2###
{
  n<-400;p<-1000;H<-2:20
  b1<-c(1,0.5);b2<-c(1,0.5);b3_1<-c(1,-0.5);b3_2<-c(-1,0.5)
  generator<-ggmGenerator()
  out2<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','arules','glmnet')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-X[,1:2]%*%b1+error
    y2<-sin(X[,1:2]%*%b2)*exp(X[,1:2]%*%b2)+error
    y3<-3*(X[,1:2]%*%b3_1)/(0.5+(1.5+X[,3:4]%*%b3_2)^2)+error
    y4<--2+Relu(-0.5*X[,1])+Relu(0.5*X[,2])+Relu(0.75*X[,3])+
      Relu(X[,4])+Relu(1.25*X[,5])+error
    
    est1<-SDA(H,y1,X,c(1,2,3,20))
    est2<-SDA(H,y2,X,c(1,2,3,20))
    est3<-SDA(H,y3,X,c(1,2,3,4,5,20))
    est4<-SDA(H,y4,X,c(1,2,3,4,5,6,20))
    est<-cbind(est1,est2,est3,est4)
    est
  }
}
stopImplicitCluster()

out<-rbind(out1,out2)
args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

#####R function for SDA-Chi#####
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

###Setting 1###
{
  n<-200;p<-1000;H<-2:20
  b1<-c(1,0.5);b2<-c(1,0.5);b3_1<-c(1,-0.5);b3_2<-c(-1,0.5)
  generator<-ggmGenerator()
  registerDoParallel(35)
  out1<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','arules','glmnet')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-X[,1:2]%*%b1+error
    y2<-sin(X[,1:2]%*%b2)*exp(X[,1:2]%*%b2)+error
    y3<-3*(X[,1:2]%*%b3_1)/(0.5+(1.5+X[,3:4]%*%b3_2)^2)+error
    y4<--2+Relu(-0.5*X[,1])+Relu(0.5*X[,2])+Relu(0.75*X[,3])+
      Relu(X[,4])+Relu(1.25*X[,5])+error
    
    est1<-SDA(H,y1,X,c(1,2,3,20))
    est2<-SDA(H,y2,X,c(1,2,3,20))
    est3<-SDA(H,y3,X,c(1,2,3,4,5,20))
    est4<-SDA(H,y4,X,c(1,2,3,4,5,6,20))
    est<-c(est1,est2,est3,est4)
    est
  }
}

###Setting 2###
{
  n<-400;p<-1000;H<-2:20
  b1<-c(1,0.5);b2<-c(1,0.5);b3_1<-c(1,-0.5);b3_2<-c(-1,0.5)
  generator<-ggmGenerator()
  out2<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','arules','glmnet')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-X[,1:2]%*%b1+error
    y2<-sin(X[,1:2]%*%b2)*exp(X[,1:2]%*%b2)+error
    y3<-3*(X[,1:2]%*%b3_1)/(0.5+(1.5+X[,3:4]%*%b3_2)^2)+error
    y4<--2+Relu(-0.5*X[,1])+Relu(0.5*X[,2])+Relu(0.75*X[,3])+
      Relu(X[,4])+Relu(1.25*X[,5])+error
    
    est1<-SDA(H,y1,X,c(1,2,3,20))
    est2<-SDA(H,y2,X,c(1,2,3,20))
    est3<-SDA(H,y3,X,c(1,2,3,4,5,20))
    est4<-SDA(H,y4,X,c(1,2,3,4,5,6,20))
    est<-c(est1,est2,est3,est4)
    est
  }
}
stopImplicitCluster()

out<-rbind(out1,out2)
args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)




