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


#importance gradient descent
importance_gradient_descent <- function(u, ra, rb, rsigma, mu_i, sig_i, indic, lr, S, M, maxiter, beta_1, beta_2) {
  iter <- 1

  new_a <- ra
  new_b <- rb
  new_Sigma_theta <- rsigma

  eps <- 0.0001

  v_a <- matrix(0, nrow = nrow(ra), ncol = ncol(ra))
  v_b <- rep(0,length(rb))
  v_sigma <- matrix(0, nrow = nrow(rsigma), ncol = ncol(rsigma))
  s_a <- matrix(0, nrow = nrow(ra), ncol = ncol(ra))
  s_b <- rep(0,length(rb))
  s_sigma <- matrix(0, nrow = nrow(rsigma), ncol = ncol(rsigma))

  while (iter <= maxiter) {
    old_a <- new_a
    old_b <- new_b
    old_Sigma_theta <- new_Sigma_theta

    # Sampling
    theta_IS <- sampling(mu_i, sig_i, nrow(u), ncol(ra), S, M)
    w_tilde <- importance_weights(theta_IS, mu_i, sig_i, old_a, old_b, old_Sigma_theta, u, nrow(u), ncol(u), S, M)$w_tilde

    grad <- importance_gradient(u, old_a, old_b, old_Sigma_theta, theta_IS, w_tilde, S, M)
    grad_a <- grad$grad_a
    grad_b <- grad$grad_b
    grad_Sigma_theta <- grad$grad_Sigma_theta

    v_a <- beta_1 * v_a + (1 - beta_1) * grad_a
    v_a <- v_a / (1 - beta_1^iter)
    v_b <- beta_1 * v_b + (1 - beta_1) * grad_b
    v_b <- v_b / (1 - beta_1^iter)
    v_sigma <- beta_1 * v_sigma + (1 - beta_1) * grad_Sigma_theta
    v_sigma <- v_sigma / (1 - beta_1^iter)

    s_a <- beta_2 * s_a + (1 - beta_2) * grad_a^2
    s_a <- s_a / (1 - beta_2^iter)
    s_b <- beta_2 * s_b + (1 - beta_2) * (grad_b)^2
    s_b <- s_b / (1 - beta_2^iter)
    s_sigma <- beta_2 * s_sigma + (1 - beta_2) * grad_Sigma_theta^2
    s_sigma <- s_sigma / (1 - beta_2^iter)

    grad_a_adam <- v_a / (sqrt(s_a) + 0.001)
    grad_b_adam <- v_b / (sqrt(s_b) + 0.001)
    grad_sigma_adam <- v_sigma / (sqrt(s_sigma) + 0.001)

    new_a <- old_a + lr * grad_a_adam
    new_a <- new_a * indic
    new_b <- old_b + lr * grad_b_adam
    new_Sigma_theta <- old_Sigma_theta + 0.1 * lr * grad_sigma_adam

    d_temp <- sqrt(diag(diag(new_Sigma_theta)))  # Standardization matrix
    new_Sigma_theta <- solve(d_temp) %*% new_Sigma_theta %*% solve(d_temp)  # (d_temp\Sigma)/d_temp
    new_Sigma_theta <- (new_Sigma_theta + t(new_Sigma_theta)) / 2

    #rescale new_a
    #new_a=new_a%*%d_temp
    e <- max(c(svd(new_a-old_a)$d[1], norm(new_b-old_b,type="2"),svd(new_Sigma_theta-old_Sigma_theta)$d[1]))
    if (e < eps) {
      break
    }

    iter <- iter + 1
  }

  result <- list(new_a=new_a, new_b=new_b, new_Sigma_theta=new_Sigma_theta)
  return(result)
}


