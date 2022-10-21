SVD.fn <-function(Y, K, eps = 0.0001) {

  N = nrow(Y)
  J = ncol(Y)
  #mean imputation
  for(i in 1:ncol(Y)) {
    Y[ , i][is.na(Y[ , i])] <- mean(Y[ , i], na.rm = TRUE)
  }
  temp = svd(Y)
  sigma = temp$d
  U = temp$u
  V = temp$v

  eta = 0.01
  threshold = (1 + eta) * sqrt(N)
  index = 1 : max(which(sigma >= threshold), K + 1)

  X = U[ ,index] %*% diag(sigma[index], length(index), length(index)) %*% t(V[ ,index])
  X[X > 1 - eps] = 1 - eps
  X[X < eps] = eps

  tilde_M = log(X / (1 - X))
  hat_d = apply(tilde_M, 2, mean)

  hat_M = tilde_M - rep(1, N) %*% t(hat_d)

  temp2 = svd(hat_M)
  sigma2 = temp2$d
  U2 = temp2$u
  V2 = temp2$v

  hat_A = V2[, 1 : K] %*% diag(sigma2[1 : K], K, K) / sqrt(N)
  hat_Theta = sqrt(N) * U2[ , 1 : K]

  return(list("A" = hat_A, "Theta" = hat_Theta))

}


identify<-function(u){
  scale=rowSums(u,na.rm = T) #num of item responsed correctly per examinee
  p=apply(u,2,mean,na.rm=T) #the frequency per item
  p=replace(p,p==0,0.001)
  p=replace(p,p==1,0.999)
  q=1-p
  y=qnorm(p,0,1) #inverse of the standard normal CDF
  y=dnorm(y,0,1) #the density function
  s=sd(scale)
  r=NULL
  #r<-ide(u,scale,p,q,y,s)
  for (i in 1:dim(u)[2]) {
    u1=u[!is.na(u[,i]),i]
    scale1=scale[!is.na(u[,i])]
    x1=scale1[u1==1]
    x2=scale1[u1==0]
    if(identical(x1,numeric(0))){
      x1=0
    }
    if(identical(x2,numeric(0))){
      x2=0
    }
    r[i]=(mean(x1)-mean(x2))/s*p[i]*q[i]/y[i]
  }
  return(r)
}

#initialization
init<-function(u,domain,indic){
  r=identify(u)
  person=dim(u)[1]
  r[r>0.9]=0.9
  r[r<0]=abs(r[r<0][1])
  r[r==0]=0.0001
  a0=t(rep(1,domain)%o%(r/sqrt(1-r^2)))*indic
  a0=replace(a0,a0>4,4)
  b0=-qnorm(colSums(u,na.rm=T)/person,0,1)/r
  b0[b0>4]=4
  b0[b0<(-4)]=-4
  Sigma = diag(domain)
  theta=matrix(rnorm(person*domain,0,1),nrow=person)
  #person*item
  xi=array(1,person)%*%t(b0)-theta%*%t(a0)
  eta0=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
  eta0[is.na(u)]=NA
  eps0=xi
  return(list(a0,b0,eta0,eps0,Sigma))
}

#change the sign of a and sigma
rt<-function(A,Sig){
  domain=dim(A)[2]
  #change the sign
  sign_reversed = (-1)^(colSums(A)<0)
  A=A%*%diag(sign_reversed)
  sign_mat=sign_reversed%*%t(sign_reversed)
  Sig=sign_mat*Sig
  return(list(ra=A,rsigma=Sig))
}
