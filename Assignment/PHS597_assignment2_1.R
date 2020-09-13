#+------------------------------+
#|PHS597 assignment2            |
#|                              |
#|Author: Chen Wang             |
#|Date:   0912,2020             |
#+------------------------------+

##Implentation of a lasso solver with cyclic coordinate descent 
library(glmnet)
#simulate data
library(MASS)
###simulate data
set.seed(100)
mu.x<-c(2,3,4)
cov.x.mat<-matrix(c(3,1,1,1,4,1,1,1,5),ncol = 3)
X<-mvrnorm(n=30,mu = mu.x,Sigma = cov.x.mat)
y<-rnorm(30,5,2)
##define the soft thresholding function
soft_thres_func<-function(lambda,z,var.x){
  b=0
  if (z > lambda) {
    b<- (z-lambda)/var.x
  }
  if (z < -lambda) {
    b<- (z+lambda)/var.x
  }
  return(b)
}

lasso.solver<-function(y,x.mat,intercept=T,lambda,n.iter=10000){
  if (intercept) {
    x.mat<-cbind(rep(1,dim(X)[1]),x.mat)
  }
  N<-dim(x.mat)[1]
  p<-dim(x.mat)[2]
  b.hat<-rep(1,p)
  for (iter in 1:n.iter) {
    for (j in 1:p) {
      X_j <- x.mat[,j]
      var.xj <- crossprod(X_j)/N
      r_j<-(y- x.mat %*% b.hat+b.hat[j]*X_j)
      z.value<- crossprod(X_j,r_j)/N
      if ((j==1) & intercept) {
        b.hat[j]<-z.value/var.xj
        next
      }
      b.hat[j]<- soft_thres_func(lambda = lambda,z = z.value,var.x = var.xj)
    }
  }  
  return(b.hat)
}

lasso.solver(y = y,x.mat = X,lambda = 0.1)
fit.lasso <- glmnet(x = X, y = y,alpha = 1,lambda = 0.1,standardize = F)
c(fit.lasso$a0,as.numeric(fit.lasso$beta))

lasso.solver(y = y,intercept = F,x.mat = X,lambda = 0.5)
fit.lasso2 <- glmnet(x = X, y = y,intercept = F,alpha = 0.5,lambda = 1,standardize = F)
c(fit.lasso2$a0,as.numeric(fit.lasso2$beta))
