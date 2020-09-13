##Implentation of a lasso solver with cyclic coordinate descent 
library(glmnet)
#simulate data
library(MASS)
###simulate data
set.seed(100)
mu.x<-c(2,3,4)
cov.x.mat<-matrix(c(3,1,1,1,4,1,1,1,5),ncol = 3)
X<-mvrnorm(n=30,mu = mu.x,Sigma = cov.x.mat)
x1<-X[,1]
x2<-X[,2]
x3<-X[,3]
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

lasso.solver<-function(y,x.mat,b.hat,lambda,n.iter=1000){
  N<-dim(x.mat)[1]
  p<-dim(x.mat)[2]
  b.hat<-b.hat
  for (iter in 1:n.iter) {
    for (j in 1:p) {
      X_j <- x.mat[,j]
      var.xj <- crossprod(X_j)/N
      r_j<-(y- X %*% b.hat+b.hat[j]*X_j)
      #print(r_j)
      z.value<- crossprod(X_j,r_j)/N
      #print(z.value)
      b.hat[j]<- soft_thres_func(lambda = lambda,z = z.value,var.x = var.xj)
      #print(b.hat)
    }
  }  
  return(b.hat)
}

lasso.solver(y = y,x.mat = X,b.hat = b.start,lambda = 0.5)
fit.lasso <- glmnet(x = X, y = y,alpha = 1,intercept = F,lambda = 0.5,standardize = F)
fit.lasso$beta
lasso.solver(y = y,x.mat = X,b.hat = b.start,lambda = 1)
fit.lasso2 <- glmnet(x = X, y = y,alpha = 1,intercept = F,lambda = 1,standardize = F)
fit.lasso2$beta
