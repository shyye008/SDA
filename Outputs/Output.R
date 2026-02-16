##### Output Figures 1, S1, and S2 #####
library(tidyverse)
library(viridis)
library(here)
library(cowplot)
dat<-read.csv("out_slice.csv")
H<-length(2:20)
index1<-1:(4*H);index2<-(4*H+1):(8*H);index3<-(8*H+1):(14*H);index4<-(14*H+1):(21*H)
dat1<-dat[1:4000,];dat2<-dat[4001:8000,]
dat1_KS<-dat1[seq(1,4000,4),];dat1_CV<-dat1[seq(2,4000,4),]
dat1_CVM<-dat1[seq(3,4000,4),];dat1_CV2<-dat1[seq(4,4000,4),]
dat2_KS<-dat2[seq(1,4000,4),];dat2_CV<-dat2[seq(2,4000,4),]
dat2_CVM<-dat2[seq(3,4000,4),];dat2_CV2<-dat2[seq(4,4000,4),]

f11_KS<-c();f11_CvM<-c();f12_KS<-c();f12_CvM<-c()
f21_KS<-c();f21_CvM<-c();f22_KS<-c();f22_CvM<-c()
f31_KS<-c();f31_CvM<-c();f32_KS<-c();f32_CvM<-c()
f41_KS<-c();f41_CvM<-c();f42_KS<-c();f42_CvM<-c()
for (i in 1:dim(dat)[2]) {
  pw1i_KS<-mean(ifelse(dat1_KS[,i]>dat1_CV[,i],1,0))
  pw1i_CvM<-mean(ifelse(dat1_CVM[,i]>dat1_CV2[,i],1,0))
  pw2i_KS<-mean(ifelse(dat2_KS[,i]>dat2_CV[,i],1,0))
  pw2i_CvM<-mean(ifelse(dat2_CVM[,i]>dat2_CV2[,i],1,0))
  if(i %in% index1){
    f11_KS<-c(f11_KS,pw1i_KS);f11_CvM<-c(f11_CvM,pw1i_CvM)
    f12_KS<-c(f12_KS,pw2i_KS);f12_CvM<-c(f12_CvM,pw2i_CvM)
  }else if(i %in% index2){
    f21_KS<-c(f21_KS,pw1i_KS);f21_CvM<-c(f21_CvM,pw1i_CvM)
    f22_KS<-c(f22_KS,pw2i_KS);f22_CvM<-c(f22_CvM,pw2i_CvM)
  }else if(i %in% index3){
    f31_KS<-c(f31_KS,pw1i_KS);f31_CvM<-c(f31_CvM,pw1i_CvM)
    f32_KS<-c(f32_KS,pw2i_KS);f32_CvM<-c(f32_CvM,pw2i_CvM)
  }else if(i %in% index4){
    f41_KS<-c(f41_KS,pw1i_KS);f41_CvM<-c(f41_CvM,pw1i_CvM)
    f42_KS<-c(f42_KS,pw2i_KS);f42_CvM<-c(f42_CvM,pw2i_CvM)
  }
}

###
dat11<-data.frame(power1=f11_KS,power2=f11_CvM,
                  H=rep(1:H,4),Index=rep(c("1","2","3","20"),each=H))
dat11$Index<-factor(dat11$Index,levels=c("1","2","3","20"))
fit11<-ggplot(dat11,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 1, n=200")+
  theme(legend.position = "bottom")
fit12<-ggplot(dat11,aes(x=H,y=power2,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 1, n=200")+
  theme(legend.position = "bottom")

dat12<-data.frame(power1=f12_KS,power2=f12_CvM,
                  H=rep(1:H,4),Index=rep(c("1","2","3","20"),each=H))
dat12$Index<-factor(dat12$Index,levels=c("1","2","3","20"))
fit13<-ggplot(dat12,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 1, n=400")+
  theme(legend.position = "bottom")
fit14<-ggplot(dat12,aes(x=H,y=power2,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 1, n=400")+
  theme(legend.position = "bottom")

dat21<-data.frame(power1=f21_KS,power2=f21_CvM,
                  H=rep(1:H,4),Index=rep(c("1","2","3","20"),each=H))
dat21$Index<-factor(dat21$Index,levels=c("1","2","3","20"))
fit21<-ggplot(dat21,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 2, n=200")+
  theme(legend.position = "bottom")
fit22<-ggplot(dat21,aes(x=H,y=power2,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 2, n=200")+
  theme(legend.position = "bottom")

dat22<-data.frame(power1=f22_KS,power2=f22_CvM,
                  H=rep(1:H,4),Index=rep(c("1","2","3","20"),each=H))
dat22$Index<-factor(dat22$Index,levels=c("1","2","3","20"))
fit23<-ggplot(dat22,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 2, n=400")+
  theme(legend.position = "bottom")
fit24<-ggplot(dat22,aes(x=H,y=power2,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 2, n=400")+
  theme(legend.position = "bottom")

dat31<-data.frame(power1=f31_KS,power2=f31_CvM,
                  H=rep(1:H,6),Index=rep(c("1","2","3","4","5","20"),each=H))
dat31$Index<-factor(dat31$Index,levels=c("1","2","3","4","5","20"))
fit31<-ggplot(dat31,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 3, n=200")+
  theme(legend.position = "bottom")
fit32<-ggplot(dat31,aes(x=H,y=power2,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 3, n=200")+
  theme(legend.position = "bottom")

dat32<-data.frame(power1=f32_KS,power2=f32_CvM,
                  H=rep(1:H,6),Index=rep(c("1","2","3","4","5","20"),each=H))
dat32$Index<-factor(dat32$Index,levels=c("1","2","3","4","5","20"))
fit33<-ggplot(dat32,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 3, n=400")+
  theme(legend.position = "bottom")
fit34<-ggplot(dat32,aes(x=H,y=power2,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 3, n=400")+
  theme(legend.position = "bottom")

dat41<-data.frame(power1=f41_KS,power2=f41_CvM,
                  H=rep(1:H,7),Index=rep(c("1","2","3","4","5","6","20"),each=H))
dat41$Index<-factor(dat41$Index,levels=c("1","2","3","4","5","6","20"))
fit41<-ggplot(dat41,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 4, n=200")+
  theme(legend.position = "bottom")
fit42<-ggplot(dat41,aes(x=H,y=power2,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 4, n=200")+
  theme(legend.position = "bottom")

dat42<-data.frame(power1=f42_KS,power2=f42_CvM,
                  H=rep(1:H,7),Index=rep(c("1","2","3","4","5","6","20"),each=H))
dat42$Index<-factor(dat42$Index,levels=c("1","2","3","4","5","6","20"))
fit43<-ggplot(dat42,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 4, n=400")+
  theme(legend.position = "bottom")
fit44<-ggplot(dat42,aes(x=H,y=power2,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 4, n=400")+
  theme(legend.position = "bottom")

###CvM###
fig<-plot_grid(fit12+theme(legend.position = "none"), fit22+theme(legend.position = "none"), 
               fit32+theme(legend.position = "none"), fit42+theme(legend.position = "none"),  
               fit14+theme(legend.position = "none"), fit24+theme(legend.position = "none"),
               fit34+theme(legend.position = "none"), fit44+theme(legend.position = "none"),
               #labels = c('A', 'B', 'C', 'D','E','F'), label_size = 12,
               nrow = 2, ncol = 4)

grobs1<-ggplotGrob(fit14)$grobs
leg1<-grobs1[[which(sapply(grobs1, function(x) x$name) == "guide-box")]]
grobs2<-ggplotGrob(fit24)$grobs
leg2<-grobs2[[which(sapply(grobs2, function(x) x$name) == "guide-box")]]
grobs3<-ggplotGrob(fit34)$grobs
leg3<-grobs3[[which(sapply(grobs3, function(x) x$name) == "guide-box")]]
grobs4<-ggplotGrob(fit44)$grobs
leg4<-grobs4[[which(sapply(grobs4, function(x) x$name) == "guide-box")]]
legend<-plot_grid(leg1,leg2,leg3,leg4,nrow=1)

fig_all<-plot_grid(fig,legend,nrow=2,rel_heights=c(1, .1))

ggsave(fig_all, filename = "Fig_H_CvM.pdf", width=14.5,height=7.5)

###KS###
fig<-plot_grid(fit11+theme(legend.position = "none"), fit21+theme(legend.position = "none"), 
               fit31+theme(legend.position = "none"), fit41+theme(legend.position = "none"),  
               fit13+theme(legend.position = "none"), fit23+theme(legend.position = "none"),
               fit33+theme(legend.position = "none"), fit43+theme(legend.position = "none"),
               #labels = c('A', 'B', 'C', 'D','E','F'), label_size = 12,
               nrow = 2, ncol = 4)

grobs1<-ggplotGrob(fit14)$grobs
leg1<-grobs1[[which(sapply(grobs1, function(x) x$name) == "guide-box")]]
grobs2<-ggplotGrob(fit24)$grobs
leg2<-grobs2[[which(sapply(grobs2, function(x) x$name) == "guide-box")]]
grobs3<-ggplotGrob(fit34)$grobs
leg3<-grobs3[[which(sapply(grobs3, function(x) x$name) == "guide-box")]]
grobs4<-ggplotGrob(fit44)$grobs
leg4<-grobs4[[which(sapply(grobs4, function(x) x$name) == "guide-box")]]
legend<-plot_grid(leg1,leg2,leg3,leg4,nrow=1)

fig_all<-plot_grid(fig,legend,nrow=2,rel_heights=c(1, .1))

ggsave(fig_all, filename = "Fig_H_KS.pdf", width=14.5,height=7.5)

###
dat<-read.csv("out_slice2.csv")
H<-length(2:20)
index1<-1:(4*H);index2<-(4*H+1):(8*H);index3<-(8*H+1):(14*H);index4<-(14*H+1):(21*H)
dat1<-dat[1:1000,];dat2<-dat[1001:2000,]

f11<-c();f12<-c();f21<-c();f22<-c();f31<-c();f32<-c()
f41<-c();f42<-c()
for (i in 1:dim(dat)[2]) {
  pw1i<-mean(ifelse(dat1[,i]<0.05,1,0))
  pw2i<-mean(ifelse(dat2[,i]<0.05,1,0))
  if(i %in% index1){
    f11<-c(f11,pw1i);f12<-c(f12,pw2i)
  }else if(i %in% index2){
    f21<-c(f21,pw1i);f22<-c(f22,pw2i)
  }else if(i %in% index3){
    f31<-c(f31,pw1i);f32<-c(f32,pw2i)
  }else if(i %in% index4){
    f41<-c(f41,pw1i);f42<-c(f42,pw2i)
  }
}

dat11<-data.frame(power1=f11,
                  H=rep(1:H,4),Index=rep(c("1","2","3","20"),each=H))
dat11$Index<-factor(dat11$Index,levels=c("1","2","3","20"))
fit11<-ggplot(dat11,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 1, n=200")+
  theme(legend.position = "bottom")

dat12<-data.frame(power1=f12,
                  H=rep(1:H,4),Index=rep(c("1","2","3","20"),each=H))
dat12$Index<-factor(dat12$Index,levels=c("1","2","3","20"))
fit12<-ggplot(dat12,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 1, n=400")+
  theme(legend.position = "bottom")

dat21<-data.frame(power1=f21,
                  H=rep(1:H,4),Index=rep(c("1","2","3","20"),each=H))
dat21$Index<-factor(dat21$Index,levels=c("1","2","3","20"))
fit21<-ggplot(dat21,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 2, n=200")+
  theme(legend.position = "bottom")

dat22<-data.frame(power1=f22,
                  H=rep(1:H,4),Index=rep(c("1","2","3","20"),each=H))
dat22$Index<-factor(dat22$Index,levels=c("1","2","3","20"))
fit22<-ggplot(dat22,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 2, n=400")+
  theme(legend.position = "bottom")

dat31<-data.frame(power1=f31,
                  H=rep(1:H,6),Index=rep(c("1","2","3","4","5","20"),each=H))
dat31$Index<-factor(dat31$Index,levels=c("1","2","3","4","5","20"))
fit31<-ggplot(dat31,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 3, n=200")+
  theme(legend.position = "bottom")

dat32<-data.frame(power1=f32,
                  H=rep(1:H,6),Index=rep(c("1","2","3","4","5","20"),each=H))
dat32$Index<-factor(dat32$Index,levels=c("1","2","3","4","5","20"))
fit32<-ggplot(dat32,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 3, n=400")+
  theme(legend.position = "bottom")

dat41<-data.frame(power1=f41,
                  H=rep(1:H,7),Index=rep(c("1","2","3","4","5","6","20"),each=H))
dat41$Index<-factor(dat41$Index,levels=c("1","2","3","4","5","6","20"))
fit41<-ggplot(dat41,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 4, n=200")+
  theme(legend.position = "bottom")

dat42<-data.frame(power1=f42,
                  H=rep(1:H,7),Index=rep(c("1","2","3","4","5","6","20"),each=H))
dat42$Index<-factor(dat42$Index,levels=c("1","2","3","4","5","6","20"))
fit42<-ggplot(dat42,aes(x=H,y=power1,color=Index))+geom_point()+
  geom_line()+labs(y="Rejection rate",title="Model 4, n=400")+
  theme(legend.position = "bottom")

###Chi###
fig<-plot_grid(fit11+theme(legend.position = "none"), fit21+theme(legend.position = "none"), 
               fit31+theme(legend.position = "none"), fit41+theme(legend.position = "none"),  
               fit12+theme(legend.position = "none"), fit22+theme(legend.position = "none"),
               fit32+theme(legend.position = "none"), fit42+theme(legend.position = "none"),
               #labels = c('A', 'B', 'C', 'D','E','F'), label_size = 12,
               nrow = 2, ncol = 4)

grobs1<-ggplotGrob(fit14)$grobs
leg1<-grobs1[[which(sapply(grobs1, function(x) x$name) == "guide-box")]]
grobs2<-ggplotGrob(fit24)$grobs
leg2<-grobs2[[which(sapply(grobs2, function(x) x$name) == "guide-box")]]
grobs3<-ggplotGrob(fit34)$grobs
leg3<-grobs3[[which(sapply(grobs3, function(x) x$name) == "guide-box")]]
grobs4<-ggplotGrob(fit44)$grobs
leg4<-grobs4[[which(sapply(grobs4, function(x) x$name) == "guide-box")]]
legend<-plot_grid(leg1,leg2,leg3,leg4,nrow=1)

fig_all<-plot_grid(fig,legend,nrow=2,rel_heights=c(1, .1))

ggsave(fig_all, filename = "Fig_H_Chi.pdf", width=14.5,height=7.5)

##### Output results for SDA-KS and SDA-CvM in Tables 1 and S2 #####
dat<-read.csv("out_single1.csv")
dat<-read.csv("out_single2.csv")
dat_KS<-dat[seq(1,4000,4),];dat_CV<-dat[seq(2,4000,4),]
dat_CvM<-dat[seq(3,4000,4),];dat_CV2<-dat[seq(4,4000,4),]

pw_KS11<-c();pw_CvM11<-c();pw_KS21<-c();pw_CvM21<-c()
pw_KS12<-c();pw_CvM12<-c();pw_KS22<-c();pw_CvM22<-c()
pw_KS13<-c();pw_CvM13<-c();pw_KS23<-c();pw_CvM23<-c()
pw_KS14<-c();pw_CvM14<-c();pw_KS24<-c();pw_CvM24<-c()
for (i in 1:400) {
  pwi_KS<-mean(ifelse(dat_KS[,i]>dat_CV[,i],1,0),na.rm=T)
  pwi_CvM<-mean(ifelse(dat_CvM[,i]>dat_CV2[,i],1,0),na.rm=T)
  if(i<=50){
    pw_KS11<-c(pw_KS11,pwi_KS);pw_CvM11<-c(pw_CvM11,pwi_CvM)
  }else if(i<=100){
    pw_KS21<-c(pw_KS21,pwi_KS);pw_CvM21<-c(pw_CvM21,pwi_CvM)
  }else if(i<=150){
    pw_KS12<-c(pw_KS12,pwi_KS);pw_CvM12<-c(pw_CvM12,pwi_CvM)
  }else if(i<=200){
    pw_KS22<-c(pw_KS22,pwi_KS);pw_CvM22<-c(pw_CvM22,pwi_CvM)
  }else if(i<=250){
    pw_KS13<-c(pw_KS13,pwi_KS);pw_CvM13<-c(pw_CvM13,pwi_CvM)
  }else if(i<=300){
    pw_KS23<-c(pw_KS23,pwi_KS);pw_CvM23<-c(pw_CvM23,pwi_CvM)
  }else if(i<=350){
    pw_KS14<-c(pw_KS14,pwi_KS);pw_CvM14<-c(pw_CvM14,pwi_CvM)
  }else if(i<=400){
    pw_KS24<-c(pw_KS24,pwi_KS);pw_CvM24<-c(pw_CvM24,pwi_CvM)
  }
}

s1<-c(1,6,11,12,16,17);s2<-c(1,6,11,12,16,17)
pw_KS11[s1];pw_CvM11[s1];mean(pw_KS11[-s1]);mean(pw_CvM11[-s1])
pw_KS12[s1];pw_CvM12[s1];mean(pw_KS12[-s1]);mean(pw_CvM12[-s1])
pw_KS13[s1];pw_CvM13[s1];mean(pw_KS13[-s1]);mean(pw_CvM13[-s1])
pw_KS14[s1];pw_CvM14[s1];mean(pw_KS14[-s1]);mean(pw_CvM14[-s1])

pw_KS21[s2];pw_CvM21[s2];mean(pw_KS21[-s2]);mean(pw_CvM21[-s2])
pw_KS22[s2];pw_CvM22[s2];mean(pw_KS22[-s2]);mean(pw_CvM22[-s2])
pw_KS23[s2];pw_CvM23[s2];mean(pw_KS23[-s2]);mean(pw_CvM23[-s2])
pw_KS24[s2];pw_CvM24[s2];mean(pw_KS24[-s2]);mean(pw_CvM24[-s2])


###

dat<-read.csv("out_multi1.csv")
dat<-read.csv("out_multi2.csv")
dat_KS<-dat[seq(1,4000,4),];dat_CV<-dat[seq(2,4000,4),]
dat_CvM<-dat[seq(3,4000,4),];dat_CV2<-dat[seq(4,4000,4),]

pw_KS11<-c();pw_CvM11<-c();pw_KS21<-c();pw_CvM21<-c();pw_KS31<-c();pw_CvM31<-c()
pw_KS12<-c();pw_CvM12<-c();pw_KS22<-c();pw_CvM22<-c();pw_KS32<-c();pw_CvM32<-c()
pw_KS13<-c();pw_CvM13<-c();pw_KS23<-c();pw_CvM23<-c();pw_KS33<-c();pw_CvM33<-c()
pw_KS14<-c();pw_CvM14<-c();pw_KS24<-c();pw_CvM24<-c();pw_KS34<-c();pw_CvM34<-c()
for (i in 1:1200) {
  pwi_KS<-mean(ifelse(dat_KS[,i]>dat_CV[,i],1,0))
  pwi_CvM<-mean(ifelse(dat_CvM[,i]>dat_CV2[,i],1,0))
  if(i<=100){
    pw_KS11<-c(pw_KS11,pwi_KS);pw_CvM11<-c(pw_CvM11,pwi_CvM)
  }else if(i<=200){
    pw_KS21<-c(pw_KS21,pwi_KS);pw_CvM21<-c(pw_CvM21,pwi_CvM)
  }else if(i<=300){
    pw_KS31<-c(pw_KS31,pwi_KS);pw_CvM31<-c(pw_CvM31,pwi_CvM)
  }else if(i<=400){
    pw_KS12<-c(pw_KS12,pwi_KS);pw_CvM12<-c(pw_CvM12,pwi_CvM)
  }else if(i<=500){
    pw_KS22<-c(pw_KS22,pwi_KS);pw_CvM22<-c(pw_CvM22,pwi_CvM)
  }else if(i<=600){
    pw_KS32<-c(pw_KS32,pwi_KS);pw_CvM32<-c(pw_CvM32,pwi_CvM)
  }else if(i<=700){
    pw_KS13<-c(pw_KS13,pwi_KS);pw_CvM13<-c(pw_CvM13,pwi_CvM)
  }else if(i<=800){
    pw_KS23<-c(pw_KS23,pwi_KS);pw_CvM23<-c(pw_CvM23,pwi_CvM)
  }else if(i<=900){
    pw_KS33<-c(pw_KS33,pwi_KS);pw_CvM33<-c(pw_CvM33,pwi_CvM)
  }else if(i<=1000){
    pw_KS14<-c(pw_KS14,pwi_KS);pw_CvM14<-c(pw_CvM14,pwi_CvM)
  }else if(i<=1100){
    pw_KS24<-c(pw_KS24,pwi_KS);pw_CvM24<-c(pw_CvM24,pwi_CvM)
  }else if(i<=1200){
    pw_KS34<-c(pw_KS34,pwi_KS);pw_CvM34<-c(pw_CvM34,pwi_CvM)
  }
}

s1<-c(1,2,6,7,11,12);s2<-c(1,6,11,12,16,17)
pw_KS11[s1];pw_CvM11[s1];mean(pw_KS11[-s1]);mean(pw_CvM11[-s1])
pw_KS12[s1];pw_CvM12[s1];mean(pw_KS12[-s1]);mean(pw_CvM12[-s1])
pw_KS13[s1];pw_CvM13[s1];mean(pw_KS13[-s1]);mean(pw_CvM13[-s1])
pw_KS14[s1];pw_CvM14[s1];mean(pw_KS14[-s1]);mean(pw_CvM14[-s1])

pw_KS31[s2];pw_CvM31[s2];mean(pw_KS31[-s2]);mean(pw_CvM31[-s2])
pw_KS32[s2];pw_CvM32[s2];mean(pw_KS32[-s2]);mean(pw_CvM32[-s2])
pw_KS33[s2];pw_CvM33[s2];mean(pw_KS33[-s2]);mean(pw_CvM33[-s2])
pw_KS34[s2];pw_CvM34[s2];mean(pw_KS34[-s2]);mean(pw_CvM34[-s2])

##### Output results for SDA-chi in Tables 1 and S2 #####
dat<-read.csv("out_single1_chi.csv")
dat<-read.csv("out_single2_chi.csv")
pw_11<-c();pw_21<-c()
pw_12<-c();pw_22<-c()
pw_13<-c();pw_23<-c()
pw_14<-c();pw_24<-c()
for (i in 1:400) {
  pwi<-mean(dat[,i]<0.05)
  if(i<=50){
    pw_11<-c(pw_11,pwi)
  }else if(i<=100){
    pw_21<-c(pw_21,pwi)
  }else if(i<=150){
    pw_12<-c(pw_12,pwi)
  }else if(i<=200){
    pw_22<-c(pw_22,pwi)
  }else if(i<=250){
    pw_13<-c(pw_13,pwi)
  }else if(i<=300){
    pw_23<-c(pw_23,pwi)
  }else if(i<=350){
    pw_14<-c(pw_14,pwi)
  }else if(i<=400){
    pw_24<-c(pw_24,pwi)
  }  
}

s1<-c(1,6,11,12,16,17);s2<-c(1,6,11,12,16,17)
pw_11[s1];mean(pw_11[-s1])
pw_12[s1];mean(pw_12[-s1])
pw_13[s1];mean(pw_13[-s1])
pw_14[s1];mean(pw_14[-s1])

pw_21[s2];mean(pw_21[-s2])
pw_22[s2];mean(pw_22[-s2])
pw_23[s2];mean(pw_23[-s2])
pw_24[s2];mean(pw_24[-s2])

###

dat<-read.csv("out_multi1_chi.csv")
dat<-read.csv("out_multi2_chi.csv")
pw_11<-c();pw_21<-c();pw_31<-c()
pw_12<-c();pw_22<-c();pw_32<-c()
pw_13<-c();pw_23<-c();pw_33<-c()
pw_14<-c();pw_24<-c();pw_34<-c()
for (i in 1:1200) {
  pwi<-mean(dat[,i]<0.05)
  if(i<=100){
    pw_11<-c(pw_11,pwi)
  }else if(i<=200){
    pw_21<-c(pw_21,pwi)
  }else if(i<=300){
    pw_31<-c(pw_31,pwi)
  }else if(i<=400){
    pw_12<-c(pw_12,pwi)
  }else if(i<=500){
    pw_22<-c(pw_22,pwi)
  }else if(i<=600){
    pw_32<-c(pw_32,pwi)
  }else if(i<=700){
    pw_13<-c(pw_13,pwi)
  }else if(i<=800){
    pw_23<-c(pw_23,pwi)
  }else if(i<=900){
    pw_33<-c(pw_33,pwi)
  }else if(i<=1000){
    pw_14<-c(pw_14,pwi)
  }else if(i<=1100){
    pw_24<-c(pw_24,pwi)
  }else if(i<=1200){
    pw_34<-c(pw_34,pwi)
  }
}

s1<-c(1,2,6,7,11,12);s2<-c(1,6,11,12,16,17)
pw_11[s1];mean(pw_11[-s1])
pw_12[s1];mean(pw_12[-s1])
pw_13[s1];mean(pw_13[-s1])
pw_14[s1];mean(pw_14[-s1])

pw_31[s2];mean(pw_31[-s2])
pw_32[s2];mean(pw_32[-s2])
pw_33[s2];mean(pw_33[-s2])
pw_34[s2];mean(pw_34[-s2])

##### Output results for HP in Tables 1 and S2 #####
dat<-read.csv("out_single_part_1.csv")
dat<-read.csv("out_single_part_2.csv")
dat<-read.csv("out_single_part_3.csv")
dat<-read.csv("out_single_part_4.csv")
dat<-read.csv("out_single_part_5.csv")
dat<-read.csv("out_single_part_6.csv")
dat<-read.csv("out_single_part_7.csv")
dat<-read.csv("out_single_part_8.csv")
pw1<-c();pw2<-c()
for (i in 1:50) {
  pwi1<-mean(ifelse(dat[1:1000,i]<0.05,1,0))
  pwi2<-mean(ifelse(dat[1001:2000,i]<0.05,1,0))
  pw1<-c(pw1,pwi1);pw2<-c(pw2,pwi2)
}
s1<-c(1,6,11,12,16,17);s2<-c(1,6,11,12,16,17)
pw1[s1];mean(pw1[-s1])
pw2[s2];mean(pw2[-s2])

###

dat<-read.csv("out_multi1_part_1.csv")
dat<-read.csv("out_multi1_part_2.csv")
dat<-read.csv("out_multi1_part_3.csv")
dat<-read.csv("out_multi1_part_4.csv")
dat<-read.csv("out_multi1_part_5.csv")
dat<-read.csv("out_multi1_part_6.csv")
dat<-read.csv("out_multi1_part_7.csv")
dat<-read.csv("out_multi1_part_8.csv")
pw1<-c();pw2<-c()
for (i in 1:50) {
  pwi1<-mean(ifelse(dat[1:1000,i]<0.05,1,0))
  pwi2<-mean(ifelse(dat[1001:2000,i]<0.05,1,0))
  pw1<-c(pw1,pwi1);pw2<-c(pw2,pwi2)
}
s1<-c(1,2,6,7,11,12);s2<-c(1,6,11,12,16,17)
pw1[s1];mean(pw1[-s1])
pw2[s2];mean(pw2[-s2])

##### Output Figure S3 #####
out<-read.csv("out_FDR.csv")
names(out)<-c("FDP1","Power1","FDP2","Power2","FDP3","Power3","FDP4","Power4")
out11<-out[seq(1,600,6),];out12<-out[seq(2,600,6),];out13<-out[seq(3,600,6),]
out21<-out[seq(4,600,6),];out22<-out[seq(5,600,6),];out23<-out[seq(6,600,6),]

fdp11<-out11[,c(1,3,5,7)];pw11<-out11[,c(2,4,6,8)]
fdp12<-out12[,c(1,3,5,7)];pw12<-out12[,c(2,4,6,8)]
fdp13<-out13[,c(1,3,5,7)];pw13<-out13[,c(2,4,6,8)]
fdp21<-out21[,c(1,3,5,7)];pw21<-out21[,c(2,4,6,8)]
fdp22<-out22[,c(1,3,5,7)];pw22<-out22[,c(2,4,6,8)]
fdp23<-out23[,c(1,3,5,7)];pw23<-out23[,c(2,4,6,8)]

colMeans(fdp11,na.rm=T);colMeans(fdp12,na.rm=T);colMeans(fdp13,na.rm=T)
colMeans(pw11,na.rm=T);colMeans(pw12,na.rm=T);colMeans(pw13,na.rm=T)

colMeans(fdp21,na.rm=T);colMeans(fdp22,na.rm=T);colMeans(fdp23,na.rm=T)
colMeans(pw21,na.rm=T);colMeans(pw22,na.rm=T);colMeans(pw23,na.rm=T)

###Setting 1###
dat_f11<-data.frame(Rate=c(as.vector(as.matrix(pw11[,-3])),
                           as.vector(as.matrix(fdp11[,-3]))),
                    Method=rep(c("CvMCD-SDA","CvMSM-SDA","LCD-Knockoff",
                                 "CvMCD-SDA","CvMSM-SDA","LCD-Knockoff"),each=100),
                    Stat=c(rep(c("Power","FDP"),each=300)))
f11<-ggplot(dat_f11,aes(x=Stat,y=Rate,color=Method))+
  geom_boxplot()+ylim(0,1)+labs(x="",title="Model 1, |S|=6")
dat_f12<-data.frame(Rate=c(as.vector(as.matrix(pw12[,-3])),
                           as.vector(as.matrix(fdp12[,-3]))),
                    Method=rep(c("CvMCD-SDA","CvMSM-SDA","LCD-Knockoff",
                                 "CvMCD-SDA","CvMSM-SDA","LCD-Knockoff"),each=100),
                    Stat=c(rep(c("Power","FDP"),each=300)))
f12<-ggplot(dat_f12,aes(x=Stat,y=Rate,color=Method))+
  geom_boxplot()+ylim(0,1)+labs(x="",title="Model 2, |S|=6")
dat_f13<-data.frame(Rate=c(as.vector(as.matrix(pw13[,-3])),
                           as.vector(as.matrix(fdp13[,-3]))),
                    Method=rep(c("CvMCD-SDA","CvMSM-SDA","LCD-Knockoff",
                                 "CvMCD-SDA","CvMSM-SDA","LCD-Knockoff"),each=100),
                    Stat=c(rep(c("Power","FDP"),each=300)))
f13<-ggplot(dat_f13,aes(x=Stat,y=Rate,color=Method))+
  geom_boxplot()+ylim(0,1)+labs(x="",title="Model 3, |S|=6")
###Setting 2###
dat_f21<-data.frame(Rate=c(as.vector(as.matrix(pw21[,-3])),
                           as.vector(as.matrix(fdp21[,-3]))),
                    Method=rep(c("CvMCD-SDA","CvMSM-SDA","LCD-Knockoff",
                                 "CvMCD-SDA","CvMSM-SDA","LCD-Knockoff"),each=100),
                    Stat=c(rep(c("Power","FDP"),each=300)))
f21<-ggplot(dat_f21,aes(x=Stat,y=Rate,color=Method))+
  geom_boxplot()+ylim(0,1)+labs(x="",title="Model 1, |S|=12")
dat_f22<-data.frame(Rate=c(as.vector(as.matrix(pw22[,-3])),
                           as.vector(as.matrix(fdp22[,-3]))),
                    Method=rep(c("CvMCD-SDA","CvMSM-SDA","LCD-Knockoff",
                                 "CvMCD-SDA","CvMSM-SDA","LCD-Knockoff"),each=100),
                    Stat=c(rep(c("Power","FDP"),each=300)))
f22<-ggplot(dat_f22,aes(x=Stat,y=Rate,color=Method))+
  geom_boxplot()+ylim(0,1)+labs(x="",title="Model 2, |S|=12")
dat_f23<-data.frame(Rate=c(as.vector(as.matrix(pw23[,-3])),
                           as.vector(as.matrix(fdp23[,-3]))),
                    Method=rep(c("CvMCD-SDA","CvMSM-SDA","LCD-Knockoff",
                                 "CvMCD-SDA","CvMSM-SDA","LCD-Knockoff"),each=100),
                    Stat=c(rep(c("Power","FDP"),each=300)))
f23<-ggplot(dat_f23,aes(x=Stat,y=Rate,color=Method))+
  geom_boxplot()+ylim(0,1)+labs(x="",title="Model 3, |S|=12")

fig<-plot_grid(f11+theme(legend.position = "none"), f21+theme(legend.position = "none"), 
               f12+theme(legend.position = "none"), f22+theme(legend.position = "none"),  
               f13+theme(legend.position = "none"), f23+theme(legend.position = "none"),
               nrow = 3, ncol = 2)

grobs<-ggplotGrob(f11)$grobs
leg<-grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

fig_all<-plot_grid(fig,leg,ncol=2,rel_widths=c(1, .15))

ggsave(fig_all, filename = "MT_all.pdf", width=10.5,height=12.5)

##### Output results for Table 2 #####
dat<-read.csv("out_miss.csv")
pw_11<-c();pw_21<-c();pw_31<-c();pw_41<-c()
pw_12<-c();pw_22<-c();pw_32<-c();pw_42<-c()
pw_13<-c();pw_23<-c();pw_33<-c();pw_43<-c()
pw_14<-c();pw_24<-c();pw_34<-c();pw_44<-c()
for (i in 1:800) {
  pwi<-mean(dat[,i]<0.05,na.rm=T)
  if(i<=50){
    pw_11<-c(pw_11,pwi)
  }else if(i<=100){
    pw_21<-c(pw_21,pwi)
  }else if(i<=150){
    pw_31<-c(pw_31,pwi)
  }else if(i<=200){
    pw_41<-c(pw_41,pwi)
  }else if(i<=250){
    pw_12<-c(pw_12,pwi)
  }else if(i<=300){
    pw_22<-c(pw_22,pwi)
  }else if(i<=350){
    pw_32<-c(pw_32,pwi)
  }else if(i<=400){
    pw_42<-c(pw_42,pwi)
  }else if(i<=450){
    pw_13<-c(pw_13,pwi)
  }else if(i<=500){
    pw_23<-c(pw_23,pwi)
  }else if(i<=550){
    pw_33<-c(pw_33,pwi)
  }else if(i<=600){
    pw_43<-c(pw_43,pwi)
  }else if(i<=650){
    pw_14<-c(pw_14,pwi)
  }else if(i<=700){
    pw_24<-c(pw_24,pwi)
  }else if(i<=750){
    pw_34<-c(pw_34,pwi)
  }else if(i<=800){
    pw_44<-c(pw_44,pwi)
  }      
}

s1<-c(1,6,11,12,16,17);s2<-c(1,6,11,12,16,17);s3<-c(1,2,6,7,11,12)
s4<-c(1,6,11,12,16,17)
pw_11[s1];mean(pw_11[-s1])
pw_12[s1];mean(pw_12[-s1])
pw_13[s1];mean(pw_13[-s1])
pw_14[s1];mean(pw_14[-s1])

pw_21[s2];mean(pw_21[-s2])
pw_22[s2];mean(pw_22[-s2])
pw_23[s2];mean(pw_23[-s2])
pw_24[s2];mean(pw_24[-s2])

pw_31[s3];mean(pw_31[-s3])
pw_32[s3];mean(pw_32[-s3])
pw_33[s3];mean(pw_33[-s3])
pw_34[s3];mean(pw_34[-s3])

pw_41[s4];mean(pw_41[-s4])
pw_42[s4];mean(pw_42[-s4])
pw_43[s4];mean(pw_43[-s4])
pw_44[s4];mean(pw_44[-s4])

##### Output Figure S4 #####
out1<-read.csv("out_miss_MH1.csv")
out11<-out1[1:100,];out12<-out1[101:200,];out13<-out1[201:300,]
out14<-out1[301:400,];out15<-out1[401:500,];out16<-out1[501:600,]

colMeans(out11,na.rm=T);colMeans(out12,na.rm=T);colMeans(out13,na.rm=T)
colMeans(out14,na.rm=T);colMeans(out15,na.rm=T);colMeans(out16,na.rm=T)

out2<-read.csv("out_miss_MH2.csv")
out21<-out2[1:100,];out22<-out2[101:200,];out23<-out2[201:300,]
out24<-out2[301:400,];out25<-out2[401:500,];out26<-out2[501:600,]

out3<-read.csv("out_miss_MH3.csv")
out31<-out3[,1:2];out32<-out3[,3:4];out33<-out3[,5:6]
out34<-out3[,7:8];out35<-out3[,9:10];out36<-out3[,11:12]

out4<-read.csv("out_miss_MH4.csv")
out41<-out4[1:100,];out42<-out4[101:200,];out43<-out4[201:300,]
out44<-out4[301:400,];out45<-out4[401:500,];out46<-out4[501:600,]

####
rate1<-c(out11[,1],out21[,1],out31[,1],out41[,1])
rate2<-c(out11[,2],out21[,2],out31[,2],out41[,2])
dat_s1<-data.frame(Rate=c(rate1,rate2),
                   Dist=rep(c("D1","D2","D3","D4","D1","D2","D3","D4"),each=100),
                   Stat=c(rep(c("FDP","Power"),each=400)))

f1<-ggplot(dat_s1,aes(x=Stat,y=Rate,color=Dist))+
  geom_boxplot()+ylim(0,1)+labs(x="",title="Model 1, |S|=6")
f1


rate1<-c(out12[,1],out22[,1],out32[,1],out42[,1])
rate2<-c(out12[,2],out22[,2],out32[,2],out42[,2])
dat_s2<-data.frame(Rate=c(rate1,rate2),
                   Dist=rep(c("D1","D2","D3","D4","D1","D2","D3","D4"),each=100),
                   Stat=c(rep(c("FDP","Power"),each=400)))

f2<-ggplot(dat_s2,aes(x=Stat,y=Rate,color=Dist))+
  geom_boxplot()+ylim(0,1)+labs(x="",title="Model 2, |S|=6")
f2

rate1<-c(out13[,1],out23[,1],out33[,1],out43[,1])
rate2<-c(out13[,2],out23[,2],out33[,2],out43[,2])
dat_s3<-data.frame(Rate=c(rate1,rate2),
                   Dist=rep(c("D1","D2","D3","D4","D1","D2","D3","D4"),each=100),
                   Stat=c(rep(c("FDP","Power"),each=400)))

f3<-ggplot(dat_s3,aes(x=Stat,y=Rate,color=Dist))+
  geom_boxplot()+ylim(0,1)+labs(x="",title="Model 3, |S|=6")
f3

rate1<-c(out14[,1],out24[,1],out34[,1],out44[,1])
rate2<-c(out14[,2],out24[,2],out34[,2],out44[,2])
dat_s4<-data.frame(Rate=c(rate1,rate2),
                   Dist=rep(c("D1","D2","D3","D4","D1","D2","D3","D4"),each=100),
                   Stat=c(rep(c("FDP","Power"),each=400)))

f4<-ggplot(dat_s4,aes(x=Stat,y=Rate,color=Dist))+
  geom_boxplot()+ylim(0,1)+labs(x="",title="Model 1, |S|=12")
f4

rate1<-c(out15[,1],out25[,1],out35[,1],out45[,1])
rate2<-c(out15[,2],out25[,2],out35[,2],out45[,2])
dat_s5<-data.frame(Rate=c(rate1,rate2),
                   Dist=rep(c("D1","D2","D3","D4","D1","D2","D3","D4"),each=100),
                   Stat=c(rep(c("FDP","Power"),each=400)))

f5<-ggplot(dat_s5,aes(x=Stat,y=Rate,color=Dist))+
  geom_boxplot()+ylim(0,1)+labs(x="",title="Model 2, |S|=12")
f5

rate1<-c(out16[,1],out26[,1],out36[,1],out46[,1])
rate2<-c(out16[,2],out26[,2],out36[,2],out46[,2])
dat_s6<-data.frame(Rate=c(rate1,rate2),
                   Dist=rep(c("D1","D2","D3","D4","D1","D2","D3","D4"),each=100),
                   Stat=c(rep(c("FDP","Power"),each=400)))

f6<-ggplot(dat_s6,aes(x=Stat,y=Rate,color=Dist))+
  geom_boxplot()+ylim(0,1)+labs(x="",title="Model 3, |S|=12")
f6

####
fig<-plot_grid(f1+theme(legend.position = "none"), f4+theme(legend.position = "none"), 
               f2+theme(legend.position = "none"), f5+theme(legend.position = "none"),  
               f3+theme(legend.position = "none"), f6+theme(legend.position = "none"),
               nrow = 3, ncol = 2)

grobs<-ggplotGrob(f1)$grobs
leg<-grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

fig_all<-plot_grid(fig,leg,ncol=2,rel_widths=c(1, .15))

ggsave(fig_all, filename = "MT_miss.pdf", width=10.5,height=12.5)

##### Output Figures 2 and S5 #####
dat<-read.csv("out_miss_LASSO.csv")

pw_11<-c();pw_21<-c()
pw_12<-c();pw_22<-c()
pw_13<-c();pw_23<-c()
for (i in 1:250) {
  pwi<-mean(dat[,i]<0.05)
  if(i<=25){
    pw_11<-c(pw_11,pwi)
  }else if(i<=50){
    pw_21<-c(pw_21,pwi)
  }else if(i<=100){
    pw_12<-c(pw_12,pwi)
  }else if(i<=150){
    pw_22<-c(pw_22,pwi)
  }else if(i<=200){
    pw_13<-c(pw_13,pwi)
  }else if(i<=250){
    pw_23<-c(pw_23,pwi)
  }
}

###Setting 1###
index<-1:25
Variable<-rep("Null",25)
s1<-c(1,6:7,11:13,16:19,21:25)
Variable[s1]<-rep("Signal",length(s1))
d1<-data.frame(pw_11,pw_21,Variable,index)

plot11<-ggplot(data = d1, aes(x = index, y = pw_11)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 25), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate",
    title="Model 1, q=5, LASSO"
  )
plot11

plot21<-ggplot(data = d1, aes(x = index, y = pw_21)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 25), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate",
    title="Model 2, q=5, LASSO"
  )
plot21

###Setting 2###
index<-1:50
Variable<-rep("Null",50)
s1<-c(1,11:12,21:23,31:34,41:45)
Variable[s1]<-rep("Signal",length(s1))
d1<-data.frame(pw_12,pw_22,Variable,index)

plot12<-ggplot(data = d1, aes(x = index, y = pw_12)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 50), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate",
    title="Model 1, q=10, LASSO"
  )
plot12

plot22<-ggplot(data = d1, aes(x = index, y = pw_22)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 50), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate",
    title="Model 2, q=10, LASSO"
  )
plot22

###Setting 3###
index<-1:50
Variable<-rep("Null",50)
s1<-c(1,11:12,21:23,31:34,41:45)
Variable[s1]<-rep("Signal",length(s1))
d1<-data.frame(pw_13,pw_23,Variable,index)

plot13<-ggplot(data = d1, aes(x = index, y = pw_13)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 50), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate"
  )
plot13

plot23<-ggplot(data = d1, aes(x = index, y = pw_23)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 50), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate"
  )
plot23

###
###
dat2<-read.csv("out_miss_SIS.csv")

pw_11_2<-c();pw_21_2<-c()
pw_12_2<-c();pw_22_2<-c()
pw_13_2<-c();pw_23_2<-c()
for (i in 1:250) {
  pwi<-mean(dat2[,i]<0.05)
  if(i<=25){
    pw_11_2<-c(pw_11_2,pwi)
  }else if(i<=50){
    pw_21_2<-c(pw_21_2,pwi)
  }else if(i<=100){
    pw_12_2<-c(pw_12_2,pwi)
  }else if(i<=150){
    pw_22_2<-c(pw_22_2,pwi)
  }else if(i<=200){
    pw_13_2<-c(pw_13_2,pwi)
  }else if(i<=250){
    pw_23_2<-c(pw_23_2,pwi)
  }
}

###Setting 1###
index<-1:25
Variable<-rep("Null",25)
s1<-c(1,6:7,11:13,16:19,21:25)
Variable[s1]<-rep("Signal",length(s1))
d1_2<-data.frame(pw_11_2,pw_21_2,Variable,index)

plot11_2<-ggplot(data = d1_2, aes(x = index, y = pw_11_2)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 25), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate",
    title="Model 1, q=5, LASSO+SIS"
  )
plot11_2

plot21_2<-ggplot(data = d1_2, aes(x = index, y = pw_21_2)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 25), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate",
    title="Model 2, q=5, LASSO+SIS"
  )
plot21_2

###Setting 2###
index<-1:50
Variable<-rep("Null",50)
s1<-c(1,11:12,21:23,31:34,41:45)
Variable[s1]<-rep("Signal",length(s1))
d1_2<-data.frame(pw_12_2,pw_22_2,Variable,index)

plot12_2<-ggplot(data = d1_2, aes(x = index, y = pw_12_2)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 50), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate",
    title="Model 1, q=10, LASSO+SIS"
  )
plot12_2

plot22_2<-ggplot(data = d1_2, aes(x = index, y = pw_22_2)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 50), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate",
    title="Model 2, q=10, LASSO+SIS"
  )
plot22_2

###Setting 3###
index<-1:50
Variable<-rep("Null",50)
s1<-c(1,11:12,21:23,31:34,41:45)
Variable[s1]<-rep("Signal",length(s1))
d1_2<-data.frame(pw_13_2,pw_23_2,Variable,index)

plot13_2<-ggplot(data = d1, aes(x = index, y = pw_13_2)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 50), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate"
  )
plot13_2

plot23_2<-ggplot(data = d1_2, aes(x = index, y = pw_23_2)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 50), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate"
  )
plot23_2

###
###
dat3<-read.csv("out_miss_LS.csv")

pw_11_3<-c();pw_21_3<-c()
pw_12_3<-c();pw_22_3<-c()
pw_13_3<-c();pw_23_3<-c()
for (i in 1:250) {
  pwi<-mean(dat3[,i]<0.05)
  if(i<=25){
    pw_11_3<-c(pw_11_3,pwi)
  }else if(i<=50){
    pw_21_3<-c(pw_21_3,pwi)
  }else if(i<=100){
    pw_12_3<-c(pw_12_3,pwi)
  }else if(i<=150){
    pw_22_3<-c(pw_22_3,pwi)
  }else if(i<=200){
    pw_13_3<-c(pw_13_3,pwi)
  }else if(i<=250){
    pw_23_3<-c(pw_23_3,pwi)
  }
}

###Setting 1###
index<-1:25
Variable<-rep("Null",25)
s1<-c(1,6:7,11:13,16:19,21:25)
Variable[s1]<-rep("Signal",length(s1))
d1_3<-data.frame(pw_11_3,pw_21_3,Variable,index)

plot11_3<-ggplot(data = d1_3, aes(x = index, y = pw_11_3)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 25), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate",
    title="Model 1, q=5, Oracle"
  )
plot11_3

plot21_3<-ggplot(data = d1_3, aes(x = index, y = pw_21_3)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 25), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate",
    title="Model 2, q=5, Oracle"
  )
plot21_2

###Setting 2###
index<-1:50
Variable<-rep("Null",50)
s1<-c(1,11:12,21:23,31:34,41:45)
Variable[s1]<-rep("Signal",length(s1))
d1_3<-data.frame(pw_12_3,pw_22_3,Variable,index)

plot12_3<-ggplot(data = d1_3, aes(x = index, y = pw_12_3)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 50), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate",
    title="Model 1, q=10, Oracle"
  )
plot12_3

plot22_3<-ggplot(data = d1_3, aes(x = index, y = pw_22_3)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 50), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate",
    title="Model 2, q=10, Oracle"
  )
plot22_3

###Setting 3###
index<-1:50
Variable<-rep("Null",50)
s1<-c(1,11:12,21:23,31:34,41:45)
Variable[s1]<-rep("Signal",length(s1))
d1_3<-data.frame(pw_13_3,pw_23_3,Variable,index)

plot13_3<-ggplot(data = d1, aes(x = index, y = pw_13_3)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 50), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate"
  )
plot13_3

plot23_3<-ggplot(data = d1_3, aes(x = index, y = pw_23_3)) +
  geom_point(aes(color = Variable), alpha = 0.5, size=2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") + 
  coord_cartesian(xlim = c(0, 50), ylim = c(0,1)) + 
  labs(#title = "R1 simulations, C1",
    x = "Index",
    y = "Rejection rate"
  )
plot23_3

###
###
fig<-plot_grid(plot11+theme(legend.position = "none"),
               plot11_2+theme(legend.position = "none"),
               plot11_3+theme(legend.position = "none"),
               plot12+theme(legend.position = "none"),
               plot12_2+theme(legend.position = "none"),
               plot12_3+theme(legend.position = "none"),
               plot21+theme(legend.position = "none"),
               plot21_2+theme(legend.position = "none"),
               plot21_3+theme(legend.position = "none"),
               plot22+theme(legend.position = "none"),
               plot22_2+theme(legend.position = "none"),
               plot22_3+theme(legend.position = "none"),
               #labels = c('A', 'B', 'C', 'D','E','F'), label_size = 12,
               nrow = 4, ncol = 3)

grobs1<-ggplotGrob(plot11)$grobs
leg1<-grobs1[[which(sapply(grobs1, function(x) x$name) == "guide-box")]]

fig_all<-plot_grid(fig,leg1,ncol=2,rel_widths=c(1, .1))

ggsave(fig_all, filename = "Fig_miss.pdf", width=12.5,height=12.5)

#####################
#####################
fig<-plot_grid(plot11+theme(legend.position = "none"),
               plot11_2+theme(legend.position = "none"),
               plot11_3+theme(legend.position = "none"),
               plot12+theme(legend.position = "none"),
               plot12_2+theme(legend.position = "none"),
               plot12_3+theme(legend.position = "none"),
               #labels = c('A', 'B', 'C', 'D','E','F'), label_size = 12,
               nrow = 2, ncol = 3)

grobs1<-ggplotGrob(plot11)$grobs
leg1<-grobs1[[which(sapply(grobs1, function(x) x$name) == "guide-box")]]

fig_all<-plot_grid(fig,leg1,ncol=2,rel_widths=c(1, .1))

ggsave(fig_all, filename = "Fig_miss.pdf", width=12.5,height=6.5)

#####################
#####################
fig<-plot_grid(plot21+theme(legend.position = "none"),
               plot21_2+theme(legend.position = "none"),
               plot21_3+theme(legend.position = "none"),
               plot22+theme(legend.position = "none"),
               plot22_2+theme(legend.position = "none"),
               plot22_3+theme(legend.position = "none"),
               #labels = c('A', 'B', 'C', 'D','E','F'), label_size = 12,
               nrow = 2, ncol = 3)

grobs1<-ggplotGrob(plot11)$grobs
leg1<-grobs1[[which(sapply(grobs1, function(x) x$name) == "guide-box")]]

fig_all<-plot_grid(fig,leg1,ncol=2,rel_widths=c(1, .1))

ggsave(fig_all, filename = "Fig_miss2.pdf", width=12.5,height=6.5)


