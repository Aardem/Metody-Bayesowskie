library(stats)
library(ggplot2)
library(invgamma)
library(mvtnorm)
library(readxl)
library(lmtest)

data <- read_excel("data.xlsx")

#preprocessing
data$Absolwenci<-data$Absolwenci/data$Ludnosc
data$Bezrobotni<-data$Bezrobotni/data$Ludnosc
data$Inwestycje<-data$Inwestycje/data$Ludnosc
data$PracownicyFUN<-data$PracownicyFUN/data$Ludnosc
data<-data[,-3]
data_org<-data

#standaryzacja
data$Wynagrodzenie<-scale(data$Wynagrodzenie)
data$Absolwenci<-scale(data$Absolwenci)
data$Bezrobotni<-scale(data$Bezrobotni)
data$Inwestycje<-scale(data$Inwestycje)
data$KobietyPrac<-scale(data$KobietyPrac)
data$PracownicyFUN<-scale(data$PracownicyFUN)

#podzial danych
set.seed(123)
index<-sample(1:nrow(data),5)
data_test<-data[index,]
data_train<-data[-index,]

#(0;1)  A(0,25;0,1) B(-0,15;0,15) I(0,65;0,25) K(-0,15;0,05) F(0,40;0,5)

#MNK
model<-lm(Wynagrodzenie ~ Absolwenci + Bezrobotni + Inwestycje + KobietyPrac + PracownicyFUN, data=data_train)
summary(model)
parametry_mnk<-model$coefficients
reszty_mnk<-model$residuals

#testy
shapiro.test(reszty_mnk)
bptest(model)
dwtest(model)
#testy git sa

#rozklad a posteriori
#Wynagrodzenie=0+0.25*Absolwenci-0.15*Bezrobotni+0.65*Inwestycje-0.15*KobietyPrac+0.4*PracownicyFUN
parametry_apriori<-c(0,0.25,-0.15,0.65,-0.15,0.4)
odchylenia_apriori<-c(1,0.1,0.15,0.25,0.05,0.5)
wynagrodzenie_apriori<-0+0.25*data_train$Absolwenci-0.15*data_train$Bezrobotni+0.65*data_train$Inwestycje-0.15*data_train$KobietyPrac+0.4*data_train$PracownicyFUN
reszty_apriori<-data_train$Wynagrodzenie-wynagrodzenie_apriori
sd_reszty_apriori<-sd(reszty_apriori)
sigma2<-((4.5/6)^2)/2
#wariancja reszt jest o polowe mniejsza niz wariancja danych

#apriori
alpha0<-nrow(data_train)-(ncol(data_train)-1)
delta0<-sum(reszty_mnk^2)

p = seq(0,1, length=1000)
rapriori<-dinvgamma(p, alpha0/2, delta0/2)

B0<-as.matrix(parametry_apriori)
Y<-as.matrix(data_train$Wynagrodzenie)
X<-as.matrix(data.frame("Stala"=rep(1,nrow(data_train)),data_train[,3:ncol(data_train)]))
E<-diag((odchylenia_apriori^2)/sigma2)

#aposteriori
E_1<-solve(t(X)%*%X+solve(E))
B<-E_1%*%(t(X)%*%Y+solve(E)%*%B0)
alpha1<-alpha0 + nrow(data_train)
delta1<-delta0+t(Y)%*%Y-t(B)%*%solve(E_1)%*%B+t(B0)%*%solve(E)%*%B0

#kolory
cvapr<-"tomato2"
cvapo<-"tomato4"
czapr<-c("turquoise1","violetred1","slateblue1","darkorchid1","indianred1","springgreen1")
czapo<-c("turquoise4","violetred4","slateblue4","darkorchid4","indianred4","springgreen4")
chpdi<-"grey60"
cmint<-"grey5"

#rozklady brzegowe
p = seq(0,1, length=1000)
raposteriori<-dinvgamma(p, alpha1/2, delta1/2)
var_df<-data.frame("var_apriori"=rapriori,"var_aposteriori"=raposteriori)

#wariancja
#plot(y=raposteriori,x=p,type="l",col="red",xlab="x",ylab="gêstoœæ")
#lines(y=rapriori,x=p)
ggplot(var_df,aes(x=p,y=raposteriori))+
  geom_line(col=cvapo)+
  geom_area(fill=cvapo,alpha=0.1)+
  geom_line(aes(x=p,y=rapriori),col=cvapr)+
  geom_area(aes(x=p,y=rapriori),fill=cvapr,alpha=0.1)+
  labs(x = "x",y="density",title = "Rozk³ad wariancji a priori oraz a posteriori")+
  theme_bw()

#funkcja tstudent
t_multi<-function(mi,sigma,k,p,n=1){
  res<-c()
  i<-0
  for (ip in p) {
    i<-i+1
    licznik<-gamma((k+n)/2)
    mianownik<-((pi*k)^(n/2))*gamma(k/2)*abs(sigma)^(1/2)
    iloczyn<-(1+(1/k*((ip-mi)^2)*1/sigma))^(-(k+n)/2)
    res[i]<-(licznik/mianownik)*iloczyn
  }
  return(res)
}

SIGMA<-as.matrix((as.numeric(delta1)/alpha1)*E_1)
p = seq(-2,2, length=1000)
for (i in 1:6) {
  yt<-t_multi(B[i],SIGMA[i,i],alpha1,p)
  yn<-dnorm(p,parametry_apriori[i],odchylenia_apriori[i])
  y_df<-data.frame("y_apriori"=yn,"y_aposteriori"=yt)
  print(ggplot(y_df,aes(x=p,y=y_aposteriori))+
    geom_line(col=czapo[i])+
    geom_area(fill=czapo[i],alpha=0.1)+
    geom_line(aes(x=p,y=y_apriori),col=czapr[i])+
    geom_area(aes(x=p,y=y_apriori),fill=czapr[i],alpha=0.1)+ 
    labs(x = "x",y="density",title = paste0("Rozk³ad a priori
                  oraz a posteriori dla ", rownames(B)[i]))+
    theme_bw())
}

#funkcja hpdi
hpdi = function(x, x.density, coverage){
  best = 0
  for (ai in 1 : (length(x) - 1))
  {for (bi in (ai + 1) : length(x))
  {mass = sum(diff(x[ai : bi]) * x.density[(ai + 1) : bi])
  if (mass >= coverage && mass / (x[bi] - x[ai]) > best)
  {best = mass / (x[bi] - x[ai])
  ai.best = ai
  bi.best = bi}}}
  c(x[ai.best], x[bi.best])
}

hpdi_interval<-list()
p = seq(-2,2, length=1000)
for (i in 1:6) {
  hpdi_interval[[i]]<-hpdi(p,t_multi(B[i],SIGMA[i,i],alpha1,p),0.95)
  yt<-t_multi(B[i],SIGMA[i,i],alpha1,p)
  y_df<-data.frame("y_aposteriori"=yt)
  print(ggplot(y_df,aes(x=p,y=y_aposteriori))+
          geom_line(col=czapo[i])+
          geom_area(fill=czapo[i],alpha=0.1)+ 
          geom_vline(aes(xintercept=hpdi_interval[[i]][1]),
                     color=chpdi, linetype="dashed", size=1)+ 
          geom_vline(aes(xintercept=hpdi_interval[[i]][2]),
                     color=chpdi, linetype="dashed", size=1)+
          labs(x = "x",y="density",title = paste0("Rozk³ad a posteriori z przedzia³ami HPDI dla ", rownames(B)[i]))+
          theme_bw())
  #plot(x=p, y=yt,type = "l",xlab="x",ylab="gêstoœæ")
  #abline(v=hpdi_interval[[i]],col="green")
}

#predykcje
data_pred_m<-as.vector(predict(model,data_test))
m_interval<-predict(model,data_test,interval="prediction")[,2:3]
data_pred_c<-as.vector(unlist(as.matrix(data.frame("stala"=rep(1,5),data_test[3:7]))%*%as.matrix(B)))

Xtau<-as.matrix(data.frame("stala"=rep(1,5),data_test[3:7]))
I<-as.matrix((as.numeric(delta1)/alpha1)*diag(rep(1,nrow(Xtau))))

p = seq(-2,2, length=1000)
hpdi_interval<-list()
for (i in 1:nrow(data_test)) {
  r<-t_multi((Xtau%*%B)[i],(I+Xtau%*%SIGMA%*%t(Xtau))[i,i],alpha1,p)
  hpdi_interval[[i]]<-hpdi(p,r,0.95)
  y_df<-data.frame("y_aposteriori"=r)
  print(ggplot(y_df,aes(x=p,y=y_aposteriori))+
          geom_line(col="darkturquoise")+
          geom_area(fill="darkturquoise",alpha=0.1)+ 
          geom_vline(aes(xintercept=hpdi_interval[[i]][1]),
                     color=chpdi, linetype="longdash", size=1)+ 
          geom_vline(aes(xintercept=hpdi_interval[[i]][2]),
                     color=chpdi, linetype="longdash", size=1)+
          geom_vline(aes(xintercept=m_interval[i,][1]),
                     color=cmint, linetype="dotted", size=1)+ 
          geom_vline(aes(xintercept=m_interval[i,][2]),
                     color=cmint, linetype="dotted", size=1)+
          labs(x = "x",y="density",title = paste0("Rozk³ad predykcyjny z przedzia³ami ufnoœci i HPDI
                                                  dla obserwacji:\n", data_test$Powiat[i]))+
          theme_bw())
}

#porownanie modeli
#model z restrykacjami stala = 0, bezrobotni = 0
model2<-lm(Wynagrodzenie ~ 0 +Absolwenci + Inwestycje + KobietyPrac + PracownicyFUN, data=data_train)
reszty_mnk2<-model2$residuals
#Wynagrodzenie=0.25*Absolwenci+0.65*Inwestycje-0.15*KobietyPrac+0.4*PracownicyFUN
parametry_apriori2<-c(0.25,0.65,-0.15,0.4)
odchylenia_apriori2<-c(0.1,0.25,0.05,0.5)
wynagrodzenie_apriori2<-0.25*data_train$Absolwenci+0.65*data_train$Inwestycje-0.15*data_train$KobietyPrac+0.4*data_train$PracownicyFUN
reszty_apriori2<-data_train$Wynagrodzenie-wynagrodzenie_apriori2
sd_reszty_apriori2<-sd(reszty_apriori2)
sigma22<-((4.5/6)^2)/2
#wariancja danych to range zmiennej objasnianej podzielone na 2 razy regula 3 sigm
#wariancja reszt jest o polowe mniejsza niz wariancja danych

alpha02<-nrow(data_train)-(ncol(data_train)-3)
delta02<-sum(reszty_mnk2^2)

p = seq(0,1, length=1000)
rapriori2<-dinvgamma(p, alpha02/2, delta02/2)

B02<-as.matrix(parametry_apriori2)
Y2<-as.matrix(data_train$Wynagrodzenie)
X2<-as.matrix(data.frame(data_train[,c(3,5:ncol(data_train))]))
E2<-diag((odchylenia_apriori2^2)/sigma22)
E_12<-solve(t(X2)%*%X2+solve(E2))
B2<-E_12%*%(t(X2)%*%Y2+solve(E2)%*%B02)
alpha12<-alpha02 + nrow(data_train)
delta12<-delta02+t(Y2)%*%Y2-t(B2)%*%solve(E_12)%*%B2+t(B02)%*%solve(E2)%*%B02

n<-nrow(data_train)
#czynniki bayesa
PyM<-(1/pi^(n/2))*
  (sqrt(abs(det(E_1)/det(E))))*
  (delta0^(alpha0/2)/delta1^(alpha1/2))*
  (gamma(alpha1/2)/gamma(alpha0/2))

PyMr<-(1/pi^(n/2))*
  (sqrt(abs(det(E_12)/det(E2))))*
  (delta02^(alpha02/2)/delta12^(alpha12/2))*
  (gamma(alpha12/2)/gamma(alpha02/2))

PM<-0.4
PMr<-0.6

iloraz_szans<-(PyM/PyMr)*(PM/PMr)
log_iloraz<-log10(iloraz_szans)

summary(model)
summary(model2)

XX<-t(X)%*%X
XY<-t(X)%*%Y
B1_g<-data.frame("Stala"=c(NA),"Absolwenci"=c(NA),"Bezrobotni"=c(NA),
                 "Inwestycje"=c(NA),"KobietyPrac"=c(NA),"PracownicyFUN"=c(NA))
sigma2_g<-c(sigma2)
m<-2000

for (i in 2:m) {
  E_1i<-solve(1/sigma2_g[i-1]*XX+solve(E))
  B1i<-E_1i%*%(1/sigma2_g[i-1]*XY+solve(E)%*%B0)
  delta1i<-delta0+t(Y-X%*%B1i)%*%(Y-X%*%B1i)
  for (j in 1:length(B)) {
    B1_g[i,j]<-rnorm(1,B1i[j],E_1i[j,j])
  }
  sigma2_g[i]<-rinvgamma(1, alpha1/2, delta1i/2)
}
B1_g<-na.omit(B1_g)
sigma2_g<-na.omit(sigma2_g)
B1_g<-B1_g[1000:2000,]
sigma2_g<-B1_g[1000:2000,]

ggplot(B1_g,aes(x=Stala))+
  geom_density(color=czapr[1], fill=czapr[1],alpha=0.1)+ 
  geom_vline(aes(xintercept=mean(Stala)),
             color=czapr[1], linetype="dashed", size=1)+
  geom_line(aes(x=p,y=y_apriori),col=czapo[1])+
  geom_area(aes(x=p,y=y_apriori),fill=czapo[1],alpha=0.1)+
  labs(x = "x",y="density",title = "Rozk³ad a posteriori sprzê¿ony oraz Gibbsa dla Stala")+
  theme_bw()
ggplot(B1_g,aes(x=Absolwenci))+
  geom_density(color=czapr[2], fill=czapr[2],alpha=0.1)+ 
  geom_vline(aes(xintercept=mean(Absolwenci)),
             color=czapr[2], linetype="dashed", size=1)+
  labs(x = "x",y="density",title = "Rozk³ad a posteriori sprzê¿ony oraz Gibbsa dla Absolwenci")+
  theme_bw()
ggplot(B1_g,aes(x=Bezrobotni))+
  geom_density(color=czapr[3], fill=czapr[3],alpha=0.1)+  
  geom_vline(aes(xintercept=mean(Bezrobotni)),
             color=czapr[3], linetype="dashed", size=1)+
  labs(x = "x",y="density",title = "Rozk³ad a posteriori sprzê¿ony oraz Gibbsa dla Bezrobotni")+
  theme_bw()
ggplot(B1_g,aes(x=Inwestycje))+
  geom_density(color=czapr[4], fill=czapr[4],alpha=0.1)+ 
  geom_vline(aes(xintercept=mean(Inwestycje)),
             color=czapr[4], linetype="dashed", size=1)+
  labs(x = "x",y="density",title = "Rozk³ad a posteriori sprzê¿ony oraz Gibbsa dla Inwestycje")+
  theme_bw()
ggplot(B1_g,aes(x=KobietyPrac))+
  geom_density(color=czapr[5], fill=czapr[5],alpha=0.1)+  
  geom_vline(aes(xintercept=mean(KobietyPrac)),
             color=czapr[5], linetype="dashed", size=1)+
  labs(x = "x",y="density",title = "Rozk³ad a posteriori sprzê¿ony oraz Gibbsa dla KobietyPrac")+
  theme_bw()
ggplot(B1_g,aes(x=PracownicyFUN))+
  geom_density(color=czapr[6], fill=czapr[6],alpha=0.1)+ 
  geom_vline(aes(xintercept=mean(PracownicyFUN)),
             color=czapr[6], linetype="dashed", size=1)+
  labs(x = "x",y="density",title = "Rozk³ad a posteriori sprzê¿ony oraz Gibbsa dla PracownicyFUN")+
  theme_bw()

raposteriori1<-rinvgamma(1000, alpha1/2, delta1/2)
ggplot(data.frame("sigma2"=sigma2_g),aes(x=sigma2))+
  geom_density(color=cvapr, fill=cvapr,alpha=0.1)+ 
  geom_vline(aes(xintercept=mean(sigma2)),
             color=cvapr, linetype="dashed", size=1)+
  geom_density(aes(x=raposteriori1),color=cvapo, fill=cvapo,alpha=0.1)+ 
  geom_vline(aes(xintercept=mean(raposteriori1)),
             color=cvapo, linetype="dashed", size=1)+
  labs(x = "x",y="density",title = "Rozk³ad a posteriori sprzê¿ony oraz Gibbsa dla wariancji")+
  theme_bw()



mv_df<-data.frame("mean_sprz"=rep(0,6),"mean_gibbs"=rep(0,6),"sd_sprz"=rep(0,6),"sd_gibbs"=rep(0,6))
p = seq(-2,2, length=1001)
for (i in 1:6) {
  yt<-t_multi(B[i],SIGMA[i,i],alpha1,p)
  yn<-dnorm(p,parametry_apriori[i],odchylenia_apriori[i])
  y_df<-data.frame("y_apriori"=yn,"y_aposteriori"=yt)
  print(ggplot(B1_g,aes(x=eval(parse(text = rownames(B)[i]))))+
          geom_density(color=czapr[i], fill=czapr[i],alpha=0.1)+ 
          geom_line(aes(x=p,y=yt),col=czapo[i])+
          geom_area(aes(x=p,y=yt),fill=czapo[i],alpha=0.1)+ 
          labs(x = "x",y="density",title = paste0("Rozk³ad a posteriori sprzê¿ony oraz Gibbsa dla ", rownames(B)[i]))+
          xlim(-0.35, 0.75)+
          theme_bw())
  mv_df[i,1]<-B[i]
  mv_df[i,2]<-mean(B1_g[,i])
  mv_df[i,3]<-sqrt(SIGMA[i,i])
  mv_df[i,4]<-sd(B1_g[,i])
  #plot(y=yt,col="red",x=p,type="l",xlab="x",ylab="gêstoœæ")
  #lines(y=yn,x=p)
}

#gibbs<-function(mi0,varp,a0,d0,data,m=1000){
#  h<-c()
#  mi<-c()
#  mi[1]<-mi0
#  h0=1/varp
#  hv=1/sigma2
#  m=m+1
#  n<-nrow(data)
#  y<-as.vector(unlist(data[,4]))
#  for (i in 2:m) {
#    h[i]<-rgamma(1,((a0+n)/2),(1/2)*(d0+sum((y-mi[i-1])^2)))
#    mi[i]<-rnorm(1,(h0*mi0+h[i]*n*mean(y))/(h0+h[i]*n),1/(h0+h[i]*n))
#  }
#  mi[1]<-NA
#  h<-na.omit(h)
#  mi<-na.omit(mi)
#  
#  return(list(h,mi))
#}
#h<-gibbs(parametry_apriori[4],odchylenia_apriori[4]^2,alpha0,delta0,data_train)[[1]]
#mi<-gibbs(parametry_apriori[4],odchylenia_apriori[4]^2,alpha0,delta0,data_train)[[2]]
#plot(density(h))
#plot(density(mi))
