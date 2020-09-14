#+------------------------------+
#|PHS597 assignment2_2          |
#|                              |
#|Author: Chen Wang             |
#|Date:   0913,2020             |
#+------------------------------+



##implementation of PLS algorithm for dimension reduction
#simulate data
library(MASS)
library(pls)
#simulate data
set.seed(100)
mu.x<-c(2,3,4)
cov.x.mat<-matrix(c(3,1,1,1,4,1,1,1,5),ncol = 3)
X<-scale(mvrnorm(n=30,mu = mu.x,Sigma = cov.x.mat))
y<-rnorm(30,5,2)

##
update.plsX<-function(X.mat,z){
  X<-X.mat
  p<-dim(X)[2]
  for (j in 1:p) {
    X_j<-X[,j]
    X[,j]<- X_j-as.numeric(crossprod(z,X_j)/crossprod(z))*z
  }
  return(X)
}

pls_direction_solver<-function(y,X.mat){
  X<-X.mat
  p<-dim(X.mat)[2]
  N<-dim(X.mat)[1]
  Z<-matrix(rep(0,N*p), ncol = p)
  for ( m in 1:p) {
    phi_m<-t(X)%*%y
    z_m <- X %*% phi_m
    Z[,m]<-z_m
    X<-update.plsX(X.mat = X,z = z_m)
  }
  return(Z)
}

Z.mat<-pls_direction_solver(y = y,X.mat = X)
fit.plsr<-plsr(y~X)
Z.plsr<-fit.plsr$scores

### Check whther directions are the same
angle<-c()
for (i in 1:dim(X)[2]) {
  a<-Z.plsr[,i]
  b<-Z.mat[,i]
  angle<-c(angle,acos(crossprod(a,b)/sqrt(crossprod(a)*crossprod(b))))
}
##All angles are 0
print(angle)          
##[1] 0 0 0

