#+-----------------------------------------------+
#| PPR implementation using natural cubic spline |
#| Chen Wang, PHS597 DL                          |
#+-----------------------------------------------+
library(MASS)
library(splines2)
##data simulation
set.seed(100)
mu.x<-c(2,3)
cov.x.mat<-matrix(c(1,0,0,2),ncol = 2)
X<-mvrnorm(n=30,mu = mu.x,Sigma = cov.x.mat)
e<-rnorm(30)
y<-y<-exp(X[,1]+X[,2])/(1+exp(X[,1]+X[,2]))+e
####
ppr.fit<-function(X,y){
  p<-dim(X)[2]
  w.rand<-rep(1,p)
  w.old<-w.rand/sqrt(sum(w.rand^2))
  while (T) {
   ##assuming w, estimate g
   V<-X%*%w.old
   g.mat<-naturalSpline(V,df=10,intercept = TRUE)
   #estimate g.beta
   g.beta<-ginv(t(g.mat)%*%g.mat)%*%t(g.mat)%*%y
   ###given g estimate w,
   g.mat.prime<-naturalSpline(V,df=10,derivs = 1,intercept = TRUE)
   g.prime<-g.mat.prime %*%g.beta
   W.mat<-diag(x = as.vector(g.prime^2))
   b_hat<-V+(y- g.mat%*%g.beta)/g.prime
   w.new<-ginv(t(X)%*%W.mat%*%X)%*%t(X)%*%W.mat%*%b_hat
   diff<-max(abs(w.new-w.old))
   #print(diff)
   if (diff<1e-5) {
     return(list(w=w.new,g.beta=g.beta))
   }
   w.old<-w.new
 } 

}
###
res<-ppr.fit(X = X,y = y)
##make prediction
ppr.predict<-function(fit.res,X){
  Basis<-naturalSpline(X%*%fit.res$w,df = 10,intercept = TRUE)
  return(Basis%*%fit.res$g.beta)
}

y_pred<-ppr.predict(fit.res = res,X = X)
sum((y_pred-y)^2)
### use ppr function
fit<-ppr(x = X,y=y,nterms = 1,sm.method = 'spline',df=10)
fit$gof




