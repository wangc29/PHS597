#+----------------------------+
#|PHS597 assignment1          |
#|                            |
#|Author:	Chen Wang           |
#|Date:		Sept. 6th,2020      |
#+----------------------------+

library(MASS)
###simulate data
set.seed(100)
mu.x<-c(2,3)
cov.x.mat<-matrix(c(1,0,0,2),ncol = 2)
X<-mvrnorm(n=30,mu = mu.x,Sigma = cov.x.mat)
x1<-X[,1]
x2<-X[,2]
y<-rnorm(30,5,2)
###Q1
###############
##one stage multiple regression
beta.hat.full<-(solve(t(X)%*%X))%*%t(X)%*%y 
beta.hat.full
##multiple stage regression
###regress y on x1 get residual e1
b1<-solve(t(x1)%*%x1)%*%t(x1)%*%y
e1<-y-b1*x1
###regress x2 on x1 get residual e2
z1<-solve(t(x1)%*%x1)%*%t(x1)%*%x2
e2<-x2-z1*x1
### regress e1 on e2
b2<-solve(t(e2)%*%e2)%*%t(e2)%*%e1
print(beta.hat.full[2,])
print(b2)
### Thus for b2, estimates are the same from one stage regression 
### and multi-stage regression

###Q2
##orthogorize x1,x2  to get d1,d2 which are orthogonal
d1<-x1
P.mat<-d1%*%solve(t(d1)%*%d1)%*%t(d1)
d2<-x2-P.mat%*%x2
##regress y on d2
b2.ortho<-solve(t(d2)%*%d2)%*%t(d2)%*%y
print(b2.ortho)
### Thus we verified that 
### by successive orthogorization the estimate for b2 
### is the same as the estimate for b2 in multi-variate regression.

print("Number of exploratory variables: p=2")
print(paste0("Estimate of beta_p.hat in one stage multiple regression: ",round(beta.hat.full[2,],digits = 3)))
print(paste0("Estimate of beta_p.hat two stage regression: ",round(b2,digits = 3)))
print(paste0("Estimate of beta_p.hat in sucessive orthogorization: ",round(b2.ortho,digits = 3)))
