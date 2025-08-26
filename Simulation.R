#####Simulation#####
Relu<-function(X){
  ifelse(X>0,X,0)
}

#####Number of slices#####
##Setting 1##
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

#####Multiple hypothesis testing#####
###Data generation###
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
  m11<-SDA_knock(y11,X,ceiling(n^(1/4)))
  m11_2<-knockoff.filter(X,y11)
  m11_3<-knockoff.filter(X,y11,offset=0)
  out1_2<-data.frame(FDP(Betac,m11_2$selected),
                     power(Betac,m11_2$selected,sum(Betac)))
  out1_3<-data.frame(FDP(Betac,m11_3$selected),
                     power(Betac,m11_3$selected,sum(Betac)))
  est1<-as.vector(data.frame(Output(m11[[1]],Betac),Output(m11[[2]],Betac),
                             out1_2,out1_3))
  
  Betac<-rep(0,p);Betac[s1]<-1
  m12<-SDA_knock(y12,X,ceiling(n^(1/4)))
  m12_2<-knockoff.filter(X,y12)
  m12_3<-knockoff.filter(X,y12,offset=0)
  out2_2<-data.frame(FDP(Betac,m12_2$selected),
                     power(Betac,m12_2$selected,sum(Betac)))
  out2_3<-data.frame(FDP(Betac,m12_3$selected),
                     power(Betac,m12_3$selected,sum(Betac)))
  est2<-as.vector(data.frame(Output(m12[[1]],Betac),Output(m12[[2]],Betac),
                             out2_2,out2_3))
  
  Betac<-rep(0,p);Betac[s1]<-1
  m13<-SDA_knock(y13,X,ceiling(n^(1/4)))
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
  m21<-SDA_knock(y21,X,ceiling(n^(1/4)))
  m21_2<-knockoff.filter(X,y21)
  m21_3<-knockoff.filter(X,y21,offset=0)
  out4_2<-data.frame(FDP(Betac,m21_2$selected),
                     power(Betac,m21_2$selected,sum(Betac)))
  out4_3<-data.frame(FDP(Betac,m21_3$selected),
                     power(Betac,m21_3$selected,sum(Betac)))
  est4<-as.vector(data.frame(Output(m21[[1]],Betac),Output(m21[[2]],Betac),
                             out4_2,out4_3))
  
  Betac<-rep(0,p);Betac[s2]<-1
  m22<-SDA_knock(y22,X,ceiling(n^(1/4)))
  m22_2<-knockoff.filter(X,y22)
  m22_3<-knockoff.filter(X,y22,offset=0)
  out5_2<-data.frame(FDP(Betac,m22_2$selected),
                     power(Betac,m22_2$selected,sum(Betac)))
  out5_3<-data.frame(FDP(Betac,m22_3$selected),
                     power(Betac,m22_3$selected,sum(Betac)))
  est5<-as.vector(data.frame(Output(m22[[1]],Betac),Output(m22[[2]],Betac),
                             out5_2,out5_3))
  
  Betac<-rep(0,p);Betac[s2]<-1
  m23<-SDA_knock(y23,X,ceiling(n^(1/4)))
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

#####Empirical type I error rate and power#####
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

#####Non-Gaussian covariates#####
s1<-c(1,6,11,12,16,17)
b1<-c(-0.4,0.6,-0.8,-0.8,1,1)
s11<-c(1,6,11);s12<-c(2,7,12)
b11<-c(0.5,-1,0.8);b12<-c(-0.8,-0.5,1)
s21<-1;s22<-6;s23<-11;s24<-12;s25<-c(16,17)
b21<--0.5;b22<-0.8;b23<--1;b24<-1.25;b25<-c(-0.8,1)
index<-1:50
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
    X<-rmvtp(n,rep(0,p),Omega,5)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    y3<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y4<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est3<-SDA(H,y3,X,index)
    est4<-SDA(H,y4,X,index)
    est<-c(est1,est2,est3,est4)
    est
  }
}

###Setting 2###
{
  out2<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    X<-rmvtp(n,rep(0,p),Omega,3)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    y3<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y4<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est3<-SDA(H,y3,X,index)
    est4<-SDA(H,y4,X,index)
    est<-c(est1,est2,est3,est4)
    est
  }
}

###Setting 3###
{
  std<-sqrt(diag(solve(Omega)))
  out3<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    Xn<-rmvnp(n,rep(0,p),Omega)
    X<-c()
    for (l in 1:p) {
      pl<-pnorm(Xn[,l],mean=0,sd=std[l])
      Xl<-qt(pl,df=5)
      X<-cbind(X,Xl)
    }
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    y3<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y4<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est3<-SDA(H,y3,X,index)
    est4<-SDA(H,y4,X,index)
    est<-c(est1,est2,est3,est4)
    est
  }
}

###Setting 4###
{
  out4<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','arules','glmnet')) %dopar% {
    set.seed(i)
    Xn<-rmvnp(n,rep(0,p),Omega)
    X<-c()
    for (l in 1:p) {
      pl<-pnorm(Xn[,l],mean=0,sd=std[l])
      Xl<-qchisq(pl,df=5)
      X<-cbind(X,Xl)
    }
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    y3<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y4<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est1<-SDA(H,y1,X,index)
    est2<-SDA(H,y2,X,index)
    est3<-SDA(H,y3,X,index)
    est4<-SDA(H,y4,X,index)
    est<-c(est1,est2,est3,est4)
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



