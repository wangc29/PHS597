library(splines)
###simulate data
set.seed(20)
x<-rnorm(100)
e<-rnorm(100)
y<-exp(x)/(1+exp(x))+e
###define a hinge helper function
hinge.func<-function(x,cauchy){
  return(max(0,x-cauchy))
}

###function to get cubic basis
get.cubic.basis<-function(x,knot){
  N<-length(x)
  K<-length(knot)
  mat.r<-matrix(rep(NA,N*K),nrow = N)
  for (i in 1:K) {
    mat.r[,i]<-sapply(x, hinge.func,cauchy=knot[i])^3
  }
  return(cbind(rep(1,N),x,x^2,x^3,mat.r))
}
### function to get natural cubic basis
##d_k function in natural cubic spline
get.nat.basis.helper<-function(x,k,knot){
  knot<-sort(knot)
  K<-length(knot)
  C1<-1/(knot[k]-knot[K])
  C2<-1/(knot[K-1]-knot[K])
  l<-C1*(((max(0,x-knot[k]))^3)-((max(0,x-knot[K]))^3))
  r<-C2*(((max(0,x-knot[K-1]))^3)-((max(0,x-knot[K]))^3))
  return(l-r)
}


get.nat.cubic.basis<-function(x,knot){
  knot<-sort(knot)
  K<-length(knot)
  N<-length(x)
  X.mat<-matrix(rep(NA,N*K),nrow = N)
  X.mat[,1]<-rep(1,N)
  X.mat[,2]<-x
  # populate X.mat 
  for (i in 3:K) {
    for (j in 1:N) {
      X.mat[j, i] <- get.nat.basis.helper(x = x[j],k = i-2,knot = knot)
    }
  }
  return(X.mat)
}

##smooth spline
# calculate integral Nk''Nj'' 
## (Shuang Gao helped me get the analytical
##form of a piece-wise function after
## removing the hinge function)
##knots need to be sorted
get.integrate_jk<-function(k,j,knot){
  K<-length(knot)
  Km1<-K-1
  C_k1 <- 1/(knot[k]-knot[K])
  C_k2 <- 1/(knot[Km1]-knot[K])
  C_j1 <- 1/(knot[j]-knot[K])
  C_j2 <- 1/(knot[Km1]-knot[K])
  ##piece-wise integration
  ##piece 1 (j,K-1)
  integrand_p1<-function(x) {36*C_k1*(x-knot[k])*C_j1*(x-knot[j])}
  
  v1 <- integrate(integrand_p1, lower = knot[j], upper=knot[Km1])$value
  ##piece 2 (K-1,K)
  integrand_p2<- function(x) {
    (6*C_k1*(x-knot[k])- 6*C_k2*(x-knot[Km1]))*(6*C_j1*(x-knot[j])- 6*C_j2*(x-knot[Km1]))}
  v2<-integrate(integrand_p2, lower = knot[Km1], upper=knot[K])$value
  return(v1+v2)
}

##function to get the omega matrix
get.omega.mat<- function(x) {
  N<-length(x)
  Omega <- matrix(rep(0, N*N), ncol = N)
  knot<-sort(x)
  
  ##define the sub-matrix by removing first two columns and rows
  Omega.sub<-matrix(rep(0, (N-2)*(N-2)), ncol = N-2)
  ##populate Omega.sub
  for (k in 1:(N-2)) {
    for (j in k:(N-2)) {
      Omega.sub[k, j] <- get.integrate_jk(k,j,knot = knot)
      ##symetric
      Omega.sub[j, k] <- Omega.sub[k, j]
    }
  }
  ##populate Omega
  Omega[3:N,3:N]<-Omega.sub
  return(Omega)
}

smooth.spline.solver<-function(x,y,lambda=100){
  N.mat<-get.nat.cubic.basis(x,knot = sort(x))
  omega.mat<-get.omega.mat(x)
  theta_hat <- solve(t(N.mat)%*%N.mat+lambda*omega.mat) %*% t(N.mat) %*% y
  return(N.mat%*%theta_hat)
}

###fitted cubic spline curve
X.cubic<-get.cubic.basis(x,knot = c(-1,0,1))
beta_hat<-solve(t(X.cubic)%*%X.cubic)%*%t(X.cubic)%*%y
y_hat<-(X.cubic%*%beta_hat)[,1]
##fitted natural cubic spline curve
X.nat.cubic<-get.nat.cubic.basis(x,knot = c(-1,0,1))
beta_nat_hat<-solve(t(X.nat.cubic)%*%X.nat.cubic)%*%t(X.nat.cubic)%*%y
y_nat_hat<-(X.nat.cubic%*%beta_nat_hat)[,1]
##smooth spline
y_hat_smooth<-smooth.spline.solver(x,y,lambda = 10)

##bspline
X.b.spline<-bs(x,knots = c(-1,0,1))
beta_bs_hat<-solve(t(X.b.spline)%*%X.b.spline)%*%t(X.b.spline)%*%y
y_bs_hat<-X.b.spline %*% beta_bs_hat[,1]

plot(x,y,main = "Fitted spline")
points(sort(x),y_hat[order(x)],col="red",type = "l")
points(sort(x),y_nat_hat[order(x)],col="blue",type = "l")
points(sort(x),y_hat_smooth[order(x)],col="green",type = "l")
points(sort(x),y_bs_hat[order(x)],col="yellow",type = "l")
legend("topleft",legend=c("cubic", "natural cubic","smooth spline","b spline"),
       col=c("red", "blue","green","yellow"), lty=1:1, cex=0.8)


