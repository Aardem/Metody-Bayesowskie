library(stats)
library(ggplot2)

BetaParameters <- function(EX, VAR) {
  alpha <- EX*(EX*(1-EX)-VAR)/VAR
  beta <- (1-EX)*(EX*(1-EX)-VAR)/VAR
  return(params = list(alpha = alpha, beta = beta))
}

Beta <- function(alpha, beta) {
  EX <- alpha/(alpha+beta)
  VAR <- (alpha*beta)/((alpha+beta)^2+(alpha+beta+1))
  return(params = list(EX = EX, VAR = VAR))
}

#0.5 i 0.25
#0.5 i 0.01
#0.95 i 0.045
#0.95 i 0.005

aposteriori<-function(mi,tau){
  n<-1000
  set.seed(123)
  vector<-rbinom(n,1,0.8)

  alpha<-BetaParameters(mi,tau)[[1]]
  beta<-BetaParameters(mi,tau)[[2]]
  p = seq(0,1, length=n)
  rapriori<-dbeta(p, alpha, beta)
  
  x<-c()
  k<-c()
  theta_ex<-c()
  theta_var<-c()
  raposteriori_list<-list()
  
  for (i in 1:n) {
    x[i]<-mean(vector[1:i])
    k[i]<-sum(vector[1:i])
    theta_ex[i]<-(k[i]+alpha)/(alpha+i+beta)
    theta_var[i]<-((k[i]+alpha)*(i-k[i]+beta))/((alpha+i+beta)^2*(alpha+i+beta+1))
    alpha_post<-BetaParameters(theta_ex[i],theta_var[i])[[1]]
    beta_post<-BetaParameters(theta_ex[i],theta_var[i])[[2]]
    raposteriori<-dbeta(p, alpha_post, beta_post)
    raposteriori_list[[i]]<-raposteriori
  }
  return(list(x=x,theta_ex=theta_ex,theta_var=theta_var,rapriori=rapriori,raposteriori=raposteriori_list))
}

#===EX&VAR===================================================================================================================================

aposteriori1<-aposteriori(0.5,0.05)
aposteriori2<-aposteriori(0.5,0.03)
aposteriori3<-aposteriori(0.5,0.02)
aposteriori4<-aposteriori(0.5,0.01)

aposteriori1234_df<-data.frame("Obserwacje"=rep(1:1000,4),"Etheta"=c(aposteriori1[[2]],aposteriori2[[2]],aposteriori3[[2]],aposteriori4[[2]]),
                             "VARtheta"=c(aposteriori1[[3]],aposteriori2[[3]],aposteriori3[[3]],aposteriori4[[3]]),
                             "rozklad"=c(rep("(0.5,0.05)",1000),rep("(0.5,0.03)",1000),rep("(0.5,0.02)",1000),rep("(0.5,0.01)",1000)),
                             "srednia"=rep(aposteriori1[[1]],4))
#View(aposteriori1234_df)
ggplot(aposteriori1234_df, aes(Obserwacje, Etheta,col=rozklad)) +
  geom_line() +
  geom_line(aes(Obserwacje,srednia),col="black")

ggplot(aposteriori1234_df, aes(Obserwacje, VARtheta,col=rozklad)) +
  geom_line()

aposteriori5<-aposteriori(0.95,0.03)
aposteriori6<-aposteriori(0.95,0.02)
aposteriori7<-aposteriori(0.95,0.01)
aposteriori8<-aposteriori(0.95,0.005)

aposteriori5678_df<-data.frame("Obserwacje"=rep(1:1000,4),"Etheta"=c(aposteriori5[[2]],aposteriori6[[2]],aposteriori7[[2]],aposteriori8[[2]]),
                               "VARtheta"=c(aposteriori5[[3]],aposteriori6[[3]],aposteriori7[[3]],aposteriori8[[3]]),
                               "rozklad"=c(rep("(0.95,0.03)",1000),rep("(0.95,0.02)",1000),rep("(0.95,0.01)",1000),rep("(0.95,0.005)",1000)),
                               "srednia"=rep(aposteriori5[[1]],4))
#View(aposteriori5678_df)
ggplot(aposteriori5678_df, aes(Obserwacje, Etheta,col=rozklad)) +
  geom_line() +
  geom_line(aes(Obserwacje,srednia),col="black")

ggplot(aposteriori5678_df, aes(Obserwacje, VARtheta,col=rozklad)) +
  geom_line()

#===APRIORI&APOSTERIORI========================================================================================================================

dist_df<-data.frame("X"=rep(seq(0.001,1,0.001),6),"Gestosc"=unlist(aposteriori1$raposteriori[c(4,20,40,100,300,1000)]),
                               "rozklad"=c(rep("4",1000),rep("20",1000),rep("40",1000),rep("100",1000),rep("300",1000),rep("1000",1000)),"apriori"=rep(aposteriori1[[4]],6))

ggplot(dist_df, aes(X, Gestosc,col=rozklad)) +
  geom_line() +
  geom_line(aes(X,apriori),col="black") + 
  scale_colour_discrete(breaks=c("4","20","40","100","300","1000"))+
  labs(y ="Gestosc (0.5,0.05)",colour="Rozklad")

dist_df<-data.frame("X"=rep(seq(0.001,1,0.001),6),"Gestosc"=unlist(aposteriori4$raposteriori[c(4,20,40,100,300,1000)]),
                    "rozklad"=c(rep("4",1000),rep("20",1000),rep("40",1000),rep("100",1000),rep("300",1000),rep("1000",1000)),"apriori"=rep(aposteriori4[[4]],6))

ggplot(dist_df, aes(X, Gestosc,col=rozklad)) +
  geom_line() +
  geom_line(aes(X,apriori),col="black")+ 
  scale_colour_discrete(breaks=c("4","20","40","100","300","1000"))+
  labs(y ="Gestosc (0.5,0.01)",colour="Rozklad")

dist_df<-data.frame("X"=rep(seq(0.001,1,0.001),6),"Gestosc"=unlist(aposteriori5$raposteriori[c(4,20,40,100,300,1000)]),
                    "rozklad"=c(rep("4",1000),rep("20",1000),rep("40",1000),rep("100",1000),rep("300",1000),rep("1000",1000)),"apriori"=rep(aposteriori5[[4]],6))

ggplot(dist_df, aes(X, Gestosc,col=rozklad)) +
  geom_line() +
  geom_line(aes(X,apriori),col="black")+ 
  scale_colour_discrete(breaks=c("4","20","40","100","300","1000"))+
  labs(y ="Gestosc (0.95,0.03)",colour="Rozklad")

dist_df<-data.frame("X"=rep(seq(0.001,1,0.001),6),"Gestosc"=unlist(aposteriori8$raposteriori[c(4,20,40,100,300,1000)]),
                    "rozklad"=c(rep("4",1000),rep("20",1000),rep("40",1000),rep("100",1000),rep("300",1000),rep("1000",1000)),"apriori"=rep(aposteriori8[[4]],6))

ggplot(dist_df, aes(X, Gestosc,col=rozklad)) +
  geom_line() +
  geom_line(aes(X,apriori),col="black")+ 
  scale_colour_discrete(breaks=c("4","20","40","100","300","1000"))+
  labs(y ="Gestosc (0.95,0.005)",colour="Rozklad")