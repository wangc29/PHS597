#+------------------------------+
#|PHS597 assignment4           	|
#|                              |
#|Author: Chen Wang             |
#|Date:   0922,2020             |
#+------------------------------+
###Implement LDA,QDA; illustration of sphering
library(MASS)
library(mvtnorm)
##simulate data
x1=rmvnorm(100,mean=c(1,1),sigma=0.1*diag(2))
x2=rmvnorm(100,mean=c(2,2),sigma=0.1*diag(2))
x3=rmvnorm(100,mean=c(3,3),sigma=0.1*diag(2))
x=rbind(x1,x2,x3);n=nrow(x)
y=c(rep(0,100),rep(1,100),rep(2,100))
###LDA##
LDA_solver<-function(x,y,QDA=F){
  n<-dim(x)[1]
  p<-dim(x)[2]
  class_label<-unique(y)
  k<-length(class_label)
  mu.mat<-matrix(rep(NA,k*p),ncol = p)
  sigma_k<-matrix(rep(NA,p*p),ncol = p)
  sigma.list<-rep(list(sigma_k),k)
  pi.vec<-rep(NA,k)
  for (i in 1:k) {
    x_i<-x[y==class_label[i],]
    mu.mat[i,]<-colMeans(x_i)
    pi.vec[i]<-nrow(x_i)/n
    sigma_i<-(t(x_i)-mu.mat[i,]) %*% t(t(x_i)-mu.mat[i,])
    sigma.list[[i]]<-sigma_i
  }
  sigma.all<-(Reduce("+",sigma.list))/(n-k)
  if(QDA){
    return(list("pi_hat"=pi.vec,"mu_hat"=mu.mat,sigma_hat=sigma.list))
  }
  return(list("pi_hat"=pi.vec,"mu_hat"=mu.mat,sigma_hat=sigma.all))
}

LDAres<-LDA_solver(x = x,y = y)
QDAres<-LDA_solver(x = x,y = y,QDA = T)
##prediction
LDA_predict<-function(X,LDAres,class_label,QDA=F){
  X<-matrix(X, ncol = 1)
  mu.mat<-LDAres$mu_hat
  pi_hat<-LDAres$pi_hat
  sigma<-LDAres$sigma_hat
  G_X<-rep(NA,nrow(mu.mat))
  if(QDA){
    sigma.list<-LDAres$sigma_hat
    for (i in 1:dim(mu.mat)[1]) {
      mu_i<-mu.mat[i,]
      mu_i<-matrix(mu_i,nrow = 1)
      pi_i<-pi_hat[i]
      sigma_i<-sigma.list[[i]]
      ######
      prob<-(-0.5*log(det(sigma_i))-0.5*(t(X)-mu_i)%*%ginv(sigma_i)%*%t((t(X)-mu_i))+log(pi_i))
      G_X[i]<-prob
    }
    
    return(class_label[which.max(G_X)])
  }
  
  for (i in 1:dim(mu.mat)[1]) {
    mu_i<-mu.mat[i,]
    mu_i<-matrix(mu_i,nrow = 1)
    pi_i<-pi_hat[i]
    prob<-t(X)%*%ginv(sigma)%*%t(mu_i)-0.5*mu_i%*%ginv(sigma)%*%t(mu_i)+log(pi_i)
    G_X[i]<-prob
  }
  return(class_label[which.max(G_X)])
}
##prediction
LDA_predict(X=x[101,],LDAres = LDAres,class_label = c(0,1,2))
LDA_predict(X=x[101,],LDAres = QDAres,QDA = T,class_label = c(0,1,2))



################sphering#####
sigma.mat<-LDAres$sigma_hat
u.mat<-eigen(sigma.mat)$vectors
D.mat<-diag(eigen(sigma.mat)$values)

W.mat<-diag(diag(D.mat)^(-0.5))%*%t(u.mat)
X.star<-t(W.mat%*%t(x))
plot(X.star[1:100,],xlim = c(0,20),ylim = c(-8,3),col='red',xlab = "X1",ylab = "X2");
points(X.star[101:200,],col='blue')
points(X.star[201:300,],col='green')
x1_star<-X.star[1:100,]
x2_star<-X.star[101:200,]
x3_star<-X.star[201:300,]
mu1_star<-colMeans(x1_star)
mu2_star<-colMeans(x2_star)
mu3_star<-colMeans(x3_star)
##
b12<- -1/((mu2_star[2]-mu1_star[2])/(mu2_star[1]-mu1_star[1]))
a12<-0.5*(mu1_star[2]+mu2_star[2])-b12*0.5*(mu1_star[1]+mu2_star[1])
abline(a12,b12)
b13<- -1/((mu3_star[2]-mu1_star[2])/(mu3_star[1]-mu1_star[1]))
a13<-0.5*(mu1_star[2]+mu3_star[2])-b13*0.5*(mu1_star[1]+mu3_star[1])
abline(a13,b13)
b23<- -1/((mu3_star[2]-mu2_star[2])/(mu3_star[1]-mu2_star[1]))
a23<-0.5*(mu2_star[2]+mu3_star[2])-b13*0.5*(mu2_star[1]+mu3_star[1])
abline(a23,b23)

#####transform back#####
p12.1<-as.vector(ginv(W.mat)%*%c(0.5*(mu1_star[1]+mu2_star[1]),0.5*(mu1_star[2]+mu2_star[2])))
p12.2<-as.vector(ginv(W.mat)%*%c(0,a12))
b12.back<-(p12.1[2]-p12.2[2])/(p12.1[1]-p12.2[1])
a12.back<-p12.1[2]-b12.back*p12.1[1]

p13.1<-as.vector(ginv(W.mat)%*%c(0.5*(mu1_star[1]+mu3_star[1]),0.5*(mu1_star[2]+mu3_star[2])))
p13.2<-as.vector(ginv(W.mat)%*%c(0,a13))
b13.back<-(p13.1[2]-p13.2[2])/(p13.1[1]-p13.2[1])
a13.back<-p13.1[2]-b13.back*p13.1[1]

p23.1<-as.vector(ginv(W.mat)%*%c(0.5*(mu2_star[1]+mu3_star[1]),0.5*(mu2_star[2]+mu3_star[2])))
p23.2<-as.vector(ginv(W.mat)%*%c(0,a23))
b23.back<-(p23.1[2]-p23.2[2])/(p23.1[1]-p23.2[1])
a23.back<-p23.1[2]-b23.back*p23.1[1]
##plot 
plot(x1,col='red',xlim=c(0,5),ylim=c(0,5),xlab = "X1",ylab = "X2")
points(x2,col='blue')
points(x3,col='green')
abline(a12.back,b12.back)
abline(a13.back,b13.back)
abline(a23.back,b23.back)




