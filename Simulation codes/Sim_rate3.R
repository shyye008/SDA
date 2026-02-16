library(foreach);library(doParallel);library(phd);library(LaplacesDemon)
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
  
  registerDoParallel(35)
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','phd')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

###Setting 2###
{
  n<-200;p<-1000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','phd')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

###Setting 3###
{
  n<-400;p<-2000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','phd')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

###Setting 4###
{
  n<-200;p<-2000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','phd')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

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
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','phd')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

###Setting 2###
{
  n<-200;p<-1000
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','phd')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

###Setting 3###
{
  n<-400;p<-2000
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','phd')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

###Setting 4###
{
  n<-200;p<-2000
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','phd')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-X[,s1]%*%b1+error
    y2<-sin(X[,s1]%*%b1)*exp(X[,s1]%*%b1)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

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
  
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','phd')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

###Setting 2###
{
  n<-200;p<-1000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','phd')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

###Setting 3###
{
  n<-400;p<-2000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','phd')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

###Setting 4###
{
  n<-200;p<-2000
  Omega<-matrix(rep(0.5,5*5),nrow=5)
  diag(Omega)<-1
  Omega<-kronecker(diag(p/5),Omega)
  
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('LaplacesDemon','phd')) %dopar% {
    set.seed(i)
    X<-rmvnp(n,rep(0,p),Omega)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

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
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','phd')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

###Setting 2###
{
  n<-200;p<-1000
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','phd')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

###Setting 3###
{
  n<-400;p<-2000
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','phd')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)

###Setting 4###
{
  n<-200;p<-2000
  out<-foreach (i=1:1000, .combine=rbind, .packages=c('bootnet','phd')) %dopar% {
    set.seed(i)
    trueNet<-genGGM(p,nei=5,propPositive=0.8,p=0.25)
    X<-generator(n,trueNet)
    error<-rnorm(n)
    y1<-3*(X[,s11]%*%b11)/(0.5+(1.5+X[,s12]%*%b12)^2)+error
    y2<--2+Relu(b21*X[,s21])+Relu(b22*X[,s22])+Relu(b23*X[,s23])+
      Relu(b24*X[,s24])+Relu(X[,s25]%*%b25)+error
    
    est<-c()
    for (j in index) {
      X_j<-X[,j];X_nj<-X[,-j]
      est1<-doubleres(y1,X_nj,X_j,nperm=1000)
      est2<-doubleres(y2,X_nj,X_j,nperm=1000)
      est<-cbind(est,c(est1,est2))
    }
    est
  }
}

args<-commandArgs(TRUE)
display = function(outfile) {
  write.csv(out, paste(getwd(),"/", outfile, sep =""), row.names = FALSE) 
}  
outfile = args[1]
display(outfile)
