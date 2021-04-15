#+------------------------------+
#|PHS597 Bootstrap          |
#|                              |
#|Author: Chen Wang             |
#|Date:   1109,2020             |
#+------------------------------+

#simulate data
library(MASS)
###simulate data
set.seed(100)
mu.x<-c(2,3)
cov.x.mat<-matrix(c(1,0,0,2),ncol = 2)
X<-mvrnorm(n=30,mu = mu.x,Sigma = cov.x.mat)
x1<-X[,1]
x2<-X[,2]
e<-rnorm(n=30)
###
#beta=0.5
beta_1<-as.vector(c(0.5,0.7))
Y_b1<-X%*%beta_1+e
beta_2<-as.vector(c(0.05,0.07))
Y_b2<-X%*%beta_2+e
############Non-parametric bootstrap functions############
boots.helper<-function(Y,x.mat,size){
  dt.mat<-cbind(Y,x.mat)
  N<-nrow(dt.mat)
  ix<-sample(seq(1:N),size = size,replace = T)
  return(dt.mat[ix,])
}

ols.solver<-function(Y,x.mat){
  return(solve(t(x.mat)%*%x.mat)%*%t(x.mat)%*%Y)
}
boots.solver<-function(Y,x.mat,fold=1000,size){
  n.col<-dim(x.mat)[2]
  res.mat<-matrix(rep(NA,n.col*fold),nrow =  n.col)
  for (i  in 1:fold) {
    boot.mat<-boots.helper(Y = Y,x.mat = x.mat,size = size)
    boot.Y<-boot.mat[,1]
    boot.X<-boot.mat[,-1]
    res.mat[,i]<-ols.solver(Y = boot.Y,x.mat = boot.X)
  }
 return(res.mat) 
}
###1000 fold np bootstrap for beta_1<-as.vector(c(0.5,0.7))
np.boot.res1<-boots.solver(Y = Y_b1,x.mat = X,fold = 1000,size = 100)
#### 95% CI
#as.data.frame(rbind(quantile(np.boot.res1[1,],probs = c(0.025,0.975)),
#      quantile(np.boot.res1[2,],probs = c(0.025,0.975))),row.names = c("b1=0.5,","b2=0.7,"))

###1000 fold np bootstrap for beta_1<-as.vector(c(0.05,0.07))
np.boot.res2<-boots.solver(Y = Y_b2,x.mat = X,fold = 1000,size = 100)
#### 95% CI
#as.data.frame(rbind(quantile(np.boot.res2[1,],probs = c(0.025,0.975)),
#                    quantile(np.boot.res2[2,],probs = c(0.025,0.975))),row.names = c("b1=0.05,","b2=0.07,"))

##########functions for parametric bootstrap#######
p.boots.helper<-function(Y,x.mat,size){
  fit<-lm(Y~x.mat-1)
  beta_hat<-as.vector(fit$coefficients)
  N<-nrow(x.mat)
  err.var.hat<-sum((fit$residuals)^2)/(N-ncol(x.mat))
  ######
  ix<-sample(seq(1:N),size = size,replace = T)
  X.boot<-x.mat[ix,]
  e_star<-rnorm(size,mean = 0,sd = sqrt(err.var.hat))
  Y_star<-X.boot%*%beta_hat+e_star
  return(cbind(Y_star,X.boot))
}
p.boots_solver<-function(Y,x.mat,fold=1000,size){
  n.col<-dim(x.mat)[2]
  res.mat<-matrix(rep(NA,n.col*fold),nrow =  n.col)
  for (i in 1:fold) {
    boot.mat<-p.boots.helper(Y = Y,x.mat = x.mat,size = size)
    boot.Y<-boot.mat[,1]
    boot.X<-boot.mat[,-1]
    res.mat[,i]<-ols.solver(Y = boot.Y,x.mat = boot.X)
  }
  return(res.mat) 
}

###1000 fold parametric bootstrap for beta_1<-as.vector(c(0.5,0.7))
p.boot.res1<-p.boots_solver(Y = Y_b1,x.mat = X,fold = 1000,size = 100)
#### 95% CI
#as.data.frame(rbind(quantile(p.boot.res1[1,],probs = c(0.025,0.975)),
#                    quantile(p.boot.res1[2,],probs = c(0.025,0.975))),row.names = c("b1=0.5,","b2=0.7,"))

###1000 fold np bootstrap for beta_1<-as.vector(c(0.05,0.07))
p.boot.res2<-p.boots_solver(Y = Y_b2,x.mat = X,fold = 1000,size = 100)
#### 95% CI
#as.data.frame(rbind(quantile(p.boot.res2[1,],probs = c(0.025,0.975)),
#                    quantile(p.boot.res2[2,],probs = c(0.025,0.975))),row.names = c("b1=0.05,","b2=0.07,"))

#######Presentation of results######
###for beta=c(0.5,0.7)
b1.1.np<-np.boot.res1[1,]
CI.b1.1.np<-round(quantile(b1.1.np,probs = c(0.025,0.975)),digits = 3)
b2.1.np<-np.boot.res1[2,]
CI.b2.1.np<-round(quantile(b2.1.np,probs = c(0.025,0.975)),digits = 3)
b1.1.p<-p.boot.res1[1,]
CI.b1.1.p<-round(quantile(b1.1.p,probs = c(0.025,0.975)),digits = 3)
b2.1.p<-p.boot.res1[2,]
CI.b2.1.p<-round(quantile(b2.1.p,probs = c(0.025,0.975)),digits = 3)
par(mfrow=c(2,2))
plot(density(b1.1.np),
     main = paste0("non-parametric b1.hat\n\n","95% CI=",
                   paste(CI.b1.1.np,collapse = "-")))
plot(density(b2.1.np),
     main = paste0("non-parametric b2.hat\n\n","95% CI=",
                   paste(CI.b2.1.np,collapse = "-")))
plot(density(b1.1.p),
     main = paste0("parametric b1.hat\n\n","95% CI=",
                   paste(CI.b1.1.p,collapse = "-")))
plot(density(b2.1.p),
     main = paste0("parametric b2.hat\n\n","95% CI=",
                   paste(CI.b2.1.p,collapse = "-")))












