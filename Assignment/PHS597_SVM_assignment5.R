#+------------------------------+
#|PHS597 assignment3         	|
#|                              |
#|Author: Chen Wang             |
#|Date:   1018,2020             |
#+------------------------------+
###implementation of a linear svm solver by using lagrange duality and comparison with known implementation
library(MASS)
library(mvtnorm)
library(quadprog)
library("e1071")
##simulate data
x1=rmvnorm(100,mean=c(1,1),sigma=0.1*diag(2))
x2=rmvnorm(100,mean=c(5,5),sigma=0.1*diag(2))
X=rbind(x1,x2);n=nrow(x)
y=matrix(c(rep(-1,100),rep(1,100)),ncol = 1)
####
svm.solver<-function(X,y){
  err<-5e-10
  N<-nrow(X)
  d.vec<-matrix(1,nrow = N)
  P.mat <- sapply(1:N, function(i) y[i]*t(X)[,i])
  D.mat <- t(P.mat)%*%P.mat
  A.mat<-t(rbind(matrix(y,nrow=1),diag(1,nrow = N)))
  res<-solve.QP(Dmat = D.mat+err*diag(1,N),dvec = d.vec,Amat = A.mat,meq = 1,factorized = F)
  alpha<-res$solution
  beta<-t(X)%*%matrix(alpha*y,nrow =N)
  b0<--0.5*(max(X[y==-1,]%*%beta)+min(X[y==1,]%*%beta))
  return(list(B=beta,B0=b0))
}
###self implemented result
res.own<-svm.solver(X =X,y = y)
B<-res.own$B
B0<-res.own$B0
#### using other package
svm_model<-svm(x = X,y = y,kernel = "linear",type = "C-classification",scale = F)
w <- t(svm_model$coefs) %*% svm_model$SV
slope_1 <- -w[1]/w[2]
intercept_1 <- svm_model$rho/w[2]
###
plot(x1,xlim=c(-1,7),ylim=c(-1,7),col='red',xlab = "X1",ylab = "X2");
points(x2,col='blue');
###plot the separation plane
abline(a = -B0/B[2] ,b = -B[2]/B[1],col="green")
abline(a = intercept_1 ,b = slope_1,col="yellow")
legend("topleft",legend=c("svm e1071", "in house svm"),
       col=c("yellow", "green"), lty=1:1, cex=0.8)




