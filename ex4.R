############## HOMEWORK 4


#####a
library(metRology)
library(e1071)
tGARCH= function(n,omega,alpha,beta,nu){
  z=rt.scaled(n, nu, mean = 0, sd = 1)
  y=rep(0,n)
  sigma=rep(0,n)
  sigma[1]=omega
  y[1]=sqrt(omega)*z[1]
  for(t in 2:n){
    sigma[t]= omega+alpha*(y[t-1])^2+beta*sigma[t-1]
    y[t]=sqrt(sigma[t])*z[t]
  }
  result=cbind(y,sigma)

}

d<-tGARCH(1000,0.05,0.1,0.85,4)
plot(d[,1],type="l")
plot(d[,2],type="l")

y<-d[,1]
mean(y)
var(y)
skewness(y)
kurtosis(y)

acf(y,200)
pacf(y,200)
acf(abs(y))

acf(y^2)
pacf(y^2)

gGARCH= function(n,omega,alpha,beta){
  z=rnorm(n,mean=0,sd=1)
  y=rep(0,n)
  sigma=rep(0,n)
  sigma[1]=omega
  y[1]=sqrt(omega)*z[1]
  for(t in 2:n){
    sigma[t]= omega+alpha*(y[t-1])^2+beta*sigma[t-1]
    y[t]=sqrt(sigma[t])*z[t]
  }
  result=cbind(y,sigma)
  
}

dt<-tGARCH(1000,0.05,0.1,0.85,4)
dg<-gGARCH(1000,0.05,0.1,0.85)

plot(dt[,2],type="l")
plot(dg[,2],type="l")


plot(dt[,1],type="l",main = "Estimates for yt t_Garch")
plot(dg[,1],type="l",main = "Estimates for yt g_Garch")



gARCH= function(n,omega,alpha){
  z=rnorm(n,mean=0,sd=1)
  y=rep(0,n)
  sigma=rep(0,n)
  sigma[1]=omega
  y[1]=sqrt(omega)*z[1]
  for(t in 2:n){
    sigma[t]= omega+alpha*(y[t-1])^2
    y[t]=sqrt(sigma[t])*z[t]
  }
  result=cbind(y,sigma)
  
}
dg<-gGARCH(1000,0.05,0.1,0.85)
dga<-gARCH(1000,0.05,0.85)
plot(d[,1],type="l")
plot(d[,2],type="l")

plot(dga[,1],type="l",main="Time varying volatility g_ARCH")
plot(dg[,1],type="l",main="Time varying volatility g_GARCH")

y<-d[,1]
mean(y)
var(y)
skewness(y)
kurtosis(y)

########## HOMEWORK 5
k=function(v){
  ku=((v+3)*(v+1)*(v+2)*3)/(v*(v+5)*(v+7))
  
  return(ku)
  
}

######2
##a

tDCS=function(d,chi,k,nu,lamda,n){
  mu=rep(0,n)
  u=rep(0,n)
  y=rt(n,nu)
  v=rep(0,n)
  mu[1]=d
  
  for(t in 1:n){
    
    v[t]=y[t]-mu[t]
    u[t]=(v[t]/(1+(v[t]^2)/(nu*exp(2*lamda))))
    mu[t+1]=d+ chi*mu[t]+ k*u[t]
  }
  result=list("y"=y,"u"=u,"mean"=mu,"v"=v)
  return(result)
}

sim<-tDCS(0.05,0.2,0.01,200,0,10000)

plot(sim$y,sim$u,pch=20,ylim=c(-5,5),cex=0.5,col="black",main="DCS-t-location model")
sim<-tDCS(0.05,0.2,0.01,10,0,10000)
points(sim$y,sim$u,col="blue",pch=20,cex=0.5)
sim<-tDCS(0.05,0.2,0.01,3,0,10000)
points(sim$y,sim$u,col="red",pch=20,cex=0.5)
legend("topleft", pch=20,legend = c("nu=200","nu=10","nu=3"),col=c("black","blue","red"))


##b
beta_tGARCH=function(d,chi,k,nu,n){
  z=rt.scaled(n, nu, mean = 0, sd = 1)
  sigma=rep(0,n)
  y=rep(0,n)
  u=rep(0,n)
  sigma[1]=d
  for(t in 1:n){
    sigma[t+1]=d+chi*sigma[t]+k*sigma[t]*u[t]
    y[t]=sqrt(sigma[t])*z[t]
    u[t]=(((nu+1)*y[t]^2)/((nu-2)*sigma[t]+y[t]^2))-1
  }
  result= list("y"=y,"u"=u)
  return(result)
  
}
sim<-beta_tGARCH(0.05,0.85,0.01,200,10000)
plot(sim$y,sim$u,pch=20,cex=0.5,col="black",main="Beta-t-GARCH",ylim=c(-2,10),xlim=c(-5,5))
sim<-beta_tGARCH(0.05,0.85,0.01,10,10000)
points(sim$y,sim$u,col="blue",pch=20,cex=0.5)
sim<-beta_tGARCH(0.05,0.85,0.01,3,10000)
points(sim$y,sim$u,col="red",pch=20,cex=0.5)
legend(-5,2, pch=20,legend = c("nu=200","nu=10","nu=3"),col=c("black","blue","red"))


########c
beta_tEGARCH=function(d,chi,k,nu,n){
  z=rt.scaled(n, nu, mean = 0, sd = 1)
  lamda=rep(0,n)
  y=rep(0,n)
  u=rep(0,n)
  lamda[1]=d
  for(t in 1:n){
    lamda[t+1]=d+chi*lamda[t]+k*u[t]
    y[t]=exp(lamda[t])*z[t]
    u[t]=(((nu+1)*y[t]^2)/(nu*exp(2*lamda[t])+y[t]^2))-1
  }
  result= list("y"=y,"u"=u)
  return(result)
  
}
sim<-beta_tEGARCH(0.05,0.85,0.01,200,10000)
plot(sim$y,sim$u,pch=20,cex=0.5,col="black",main="Beta-t-EGARCH",ylim=c(-3,10),xlim=c(-10,10))
sim<-beta_tEGARCH(0.05,0.85,0.01,10,10000)
points(sim$y,sim$u,col="blue",pch=20,cex=0.5)
sim<-beta_tEGARCH(0.05,0.85,0.01,3,10000)
points(sim$y,sim$u,col="red",pch=20,cex=0.5)
legend(-10,1, pch=20,legend = c("nu=200","nu=10","nu=3"),col=c("black","blue","red"))

########3

set.seed(123)
sim<-tDCS(0.05,0.85,0.5,3,0,100)
plot(sim$y,type="l",ylim=c(-15,15),main = "DCS-t location model ")
lines(sim$mean,col="red")
legend(0,-6,legend=c("yt","mean"),col=c("black","red",3),lty=1)

plot(sim$v,type = "l",ylim=c(-15,15))
lines(sim$u,type = "l",col="red")
lines(sim$y,type="l",col=3)
legend(0,-6,legend=c("vt","ut","yt"),col=c("black","red",3),lty=1)
acf(y)
acf(sim$u)
acf(sim$v)

