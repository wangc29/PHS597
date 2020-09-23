#+------------------------------+
#|PHS597 assignment3         	|
#|                              |
#|Author: Chen Wang             |
#|Date:   0922,2020             |
#+------------------------------+
###illustration of classification by regressing indicator matrix

library(MASS)
set.seed(10)
#######simulate three clusters#####
cov.x.mat<-matrix(c(2,0,0,4),ncol = 2)
mu1<-c(3,3)
mu2<-c(10,10)
mu3<-c(15,15)
X1<-mvrnorm(n = 30,mu = mu1,Sigma = cov.x.mat)
X2<-mvrnorm(n = 30,mu = mu2,Sigma = cov.x.mat)
X3<-mvrnorm(n = 30,mu = mu3,Sigma = cov.x.mat)
plot(X1,xlim=c(-5,25),ylim=c(-5,25),col='red',xlab = "X1",ylab = "X2");
points(X2,col='blue');
points(X3,col='green');
X<-rbind(X1,X2,X3)
X<-cbind(rep(1,dim(X)[1]),X)
Y1<-c(rep(1,30),rep(0,60))
Y2<-c(rep(0,30),rep(1,30),rep(0,30))
Y3<-c(rep(0,60),rep(1,30))
Y<-cbind(Y1,Y2,Y3)
####Regression
B<-solve(crossprod(X))%*%crossprod(X,Y)
y.hat<-X%*%B
y.pred<-apply(y.hat,1,which.max);
y.pred
####
# draw decision boundary
b12=B[,1]-B[,2];
b23=B[,2]-B[,3];
b13=B[,1]-B[,3];
abline(-b12[1]/b12[3],-b12[2]/b12[3]);
text(3, 10, "cutoff", col = "red")
abline(-b23[1]/b23[3],-b23[2]/b23[3]);
abline(-b13[1]/b13[3],-b13[2]/b13[3])
