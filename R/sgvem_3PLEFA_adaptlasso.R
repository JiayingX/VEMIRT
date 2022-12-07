#' Stochastic GVEM with Adaptive Lasso Penalty for Exploratory M3PL Analysis
#'
#' @param u a \eqn{N \times J} \code{matrix} or a \code{data.frame} that
#' consists of binary responses of \eqn{N} individuals to \eqn{J} items. The
#' missing values are coded as \code{NA}
#' @param indic a \eqn{J \times K} \code{matrix} or a \code{data.frame} that
#' describes the factor loading structure of \eqn{J} items to \eqn{K} factors. It
#' consists of binary values where 0 refers to the item is irrelevant with this factor,
#' 1 otherwise. For exploratory factor analysis with lasso penalty, \code{indic} should be
#' imposed certain constraints on the a \eqn{K \times K} sub-matrix to ensure identifiability.
#' The remaining parts do not assume any pre-specified zero structure but instead, the
#' appropriate lasso penalty would recover the true zero structure. Also see \code{constrain}
#' @param samp a subsample for each iteration; default is 50
#' @param forgetrate the forget rate for the stochastic algorithm. The value should be within the
#' range from 0.5 to 1. Default is 0.51
#' @param mu_b the mean parameter for the prior distribution of item difficulty parameters
#' @param sigma2_b the variance parameter for the prior distribution of item difficulty parameters
#' @param Alpha the \eqn{\alpha} parameter for the prior distribution of guessing parameters
#' @param Beta the \eqn{\beta} parameter for the prior distribution of guessing parameters
#' @param max.iter the maximum number of iterations for the EM cycle; default is 5000
#' @param constrain the constraint setting: \code{"C1"} or \code{"C2"}. To ensure
#' identifiablity, \code{"C1"} sets a \eqn{K \times K} sub-matrix of \code{indic} to be an
#' identity matrix.This constraint anchor \eqn{K} factors by designating \eqn{K} items that load solely on each factor respectively.
#' Note that the \eqn{K \times K} matrix does not have to appear at the top of the \code{indic} matrix.
#' \code{"C2"} sets the \eqn{K \times K} sub-matrix to be a lower triangular matrix with the diagonal being ones. That is, there
#' are test items associated with each factor for sure and they may be associated
#' with other factors as well. Nonzero entries (in the lower triangular part) except for the diagonal entries of the
#' sub-matrix are penalized during the estimation procedure. For instance, assume \eqn{K=3}, then the \code{"C2"} constraint will
#' imply the following submatrix: \eqn{C2=\begin{bmatrix} 1 & 0 & 0\\ 1 & 1 & 0\\ 1 & 1 & 1\\\end{bmatrix}}. As shown, item 1 is allowed to only
#' load on the first factor, item 2 will for sure load on the second factor but it may also load on the first factor (hence a penalty is added
#' on the \eqn{(2,1)} element of \code{"C1"}, i.e., \eqn{C2_{2,1}} ). Item 3 will for sure load on the third factor but it may also load on the
#' first two factors. However, note that for all remaining items their loading vector will all be \eqn{(1, 1, 1)} hence indistinguishable from the
#' third anchor item. Therefore, we need to alert the algorithm that this third anchor item will for sure load on the third factor, and
#' and whether or not it loads on the first two factors depends on the regularization results. Therefore, we need to specify
#' \code{"non_pen="} to identify the \eqn{K}th anchor item. Although, \code{"C2"} is much weaker than \code{"C1"}, it still ensures empirical identifiability. Default is \code{"C1"}.
#' During estimation, under both the \code{"C1"} and \code{"C1"} constraints, the population means and variances are constrained to be 0 and 1, respectively.
#' @param non_pen the index of an item which is associated with each factor to satisfy \code{"C2"}.
#' For \code{C1}, the input can be \code{NULL}
#' @param gamma a numerical value of adaptive lasso parameter. Zou (2006) recommended three values, 0.5, 1, and 2.
#' The default value is 2.
#'
#' @return a list containing the following objects:
#'   \item{ra}{item discrimination parameters, a \eqn{J \times K} \code{matrix}}
#'   \item{rb}{item difficulty parameters, vector of length \eqn{J}}
#'   \item{rc}{item guessing parameters, vector of length \eqn{J}}
#'   \item{rs}{variational parameters \eqn{s}, a \eqn{N \times J} matrix}
#'   \item{reta}{variational parameters \eqn{\eta(\xi)}, a \eqn{N \times J} matrix}
#'   \item{reps}{variational parameters \eqn{\xi}, a \eqn{N \times J} matrix}
#'   \item{rsigma}{population variance-covariance matrix, a \eqn{K \times K} matrix}
#'   \item{mu_i}{mean parameter for each person, a \eqn{K \times N} matrix}
#'   \item{sig_i}{covariance matrix for each person, a \eqn{K \times K \times N} array}
#'   \item{n}{the number of iterations for the EM cycle}
#'   \item{Q_mat}{factor loading structure, a \eqn{J \times K} matrix}
#'   \item{GIC}{model fit index}
#'   \item{AIC}{model fit index}
#'   \item{BIC}{model fit index}
#'   \item{lbd}{numerical value of lasso penalty parameter \eqn{\lambda}}
#' @references
#' Cho, A. E., Xiao, J., Wang, C., & Xu, G. (2022). Regularized Variational Estimation for Exploratory Item Factor Analysis. \emph{Psychometrika}. https://doi.org/10.1007/s11336-022-09874-6
#'
#' Zou, H. (2006). The adaptive LASSO and its oracle properties.  \emph{Journal of the American Statistical Association, 7}, 1011418–1429.
#'
#' @seealso \code{\link{sgvem_3PLEFA_rot}}, \code{\link{sgvem_3PLEFA_lasso}}, \code{\link{exampleIndic_efa3pl_c1}}, \code{\link{exampleIndic_efa3pl_c2}}
#' @export
#'
#' @examples
#' sgvem_3PLEFA_adaptlasso(exampleData_3pl, exampleIndic_efa3pl_c1,samp=50,
#' forgetrate=0.51,mu_b=0,sigma2_b=4,Alpha=10,Beta=40,max.iter=5000,
#' constrain="C1",non_pen=NULL,gamma=2)
#' @examples
#' sgvem_3PLEFA_adaptlasso(exampleData_3pl, exampleIndic_efa3pl_c2,samp=50,
#' forgetrate=0.51,mu_b=0,sigma2_b=4,Alpha=10,Beta=40,max.iter=5000,
#' constrain="C2",non_pen=19,gamma=2)
#main function for gvem_2PLEFA_adaptlasso
sgvem_3PLEFA_adaptlasso<-function(u,indic,samp=50,forgetrate=0.51,
                                  mu_b,sigma2_b,Alpha,Beta,max.iter=5000,constrain="C1",non_pen=NULL,gamma=2){
  start=Sys.time()
  u=data.matrix(u)
  indic=data.matrix(indic)
  domain=dim(indic)[2]
  if(constrain=="C1"){
    invisible(capture.output(result<-sgvem_3PLEFA_adaptive_const1_all(u,domain,indic,samp,forgetrate,mu_b,sigma2_b,Alpha,Beta,gamma,max.iter)))
  }else{
    if(is.null(non_pen)){
      stop('non_pen argument is required for the C2 constraint',call.=FALSE)
    }else{
      invisible(capture.output(result<-sgvem_3PLEFA_adaptive_const2_all(u,domain,samp,forgetrate,
                                              mu_b,sigma2_b,Alpha,Beta,indic,non_pen,gamma,
                                              max.iter)))
    }
  }
  if(result$lbd==0.1 || result$lbd==40){
    warning("The optimal penalty parameter may be out of range.",call. = FALSE)
  }
  if(result$n==max.iter){
    warning("The maximum number of EM cycles reached!",call.=FALSE)
  }
  end=Sys.time()
  duration=end-start
  cat(paste("Total Execution Time:", round(duration[[1]], 2),  units(duration)),"\n")
  return(result)
}


#adaptive lasso with constraint 1 function
sgvem_3PLEFA_adapt_const1 <- function(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,forgetrate,domain,lbd,indic,
                                      mu_b,sigma2_b,Alpha,Beta,weights,nopenalty_col,max.iter) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 1
  # new_s=matrix(0,nrow=person,ncol=item)
  # aveg_p=colMeans(u==1)
  # for (i in 1:person) {
  #   for (j in 1:item) {
  #     if(u[i,j]==0){
  #       new_s[i,j]=1;
  #     }else{
  #       new_s[i,j]=runif(1,aveg_p[j],1)
  #     }
  #   }
  # }
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(diag(domain),dim=c(domain,domain,person))
  #xi=eta
  prev_a_num = matrix(0,nrow=domain,ncol=item)
  prev_a_denom = array(0,dim=c(domain,domain,item))
  prev_delta = matrix(0,nrow=item,ncol=domain)
  prev_deriv2 = matrix(0,nrow=item,ncol=domain);
  prev_b_num = matrix(0,nrow=item,ncol=1)
  prev_b_denom = matrix(0,nrow=item,ncol=1)
  prev_c_num = matrix(0,nrow=item,ncol=1)
  prev_c_denom = 0
  n = 1
  par_Sigma = Sigma
  #lb=NULL
  #nv<-NULL
  nlb<-n20<-NULL
  prev_lb = 0
  while(converge==1 && Matrix::rankMatrix(Sigma) == domain && n < max.iter){
    if(n==1){
      dec_st=1
      id_1 = sample(person,inite)
    }else{
      dec_st = (n+1)^(-forgetrate)*0.25
      #dec_st = (n+1)^(-forgetrate)
      id_1 = sample(person,samp)
    }
    s_num = (1-new_c)/new_c
    par_Sigma = Sigma
    #update MU and SIGMA,xi,Spart
    rs1<-ecfa3pl(u, Sigma, domain, id_1, item, eta, new_s, new_a, new_b, SIGMA, MU, xi)
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    xi=rs1$xi
    Spart=rs1$Spart

    #update s
    new_s=s3pl(item,id_1,new_s,eta,new_a,s_num,new_b,u,SIGMA,MU,xi)
    #update eta
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=eta3pl(item,id_1,eta,abs(xi))

    Sigma=Spart/length(id_1)
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)

    #update b
    par_b=new_b
    rs2=b3pl(item,id_1,new_s,eta,u,new_a,dec_st,MU,prev_b_num,prev_b_denom,mu_b,sigma2_b)
    new_b=rs2$new_b
    b_denom= rs2$b_denom
    b_num = rs2$b_num


    par_a=new_a

    #update a
    rs3=nalc13pl(u, indic, nopenalty_col, eta, new_s, new_b, SIGMA, MU, new_a, id_1, prev_a_num, prev_a_denom, dec_st)
    new_a1=rs3$new_a
    new_a1[-nopenalty_col,]=new_a[-nopenalty_col,]
    a_denom=rs3$a_denom
    a_num=rs3$a_num
    #L1-penalty
    sdf=setdiff(1:item,nopenalty_col)
    rs4=paal3pl(u,domain, item, id_1, lbd, eta, new_s, new_a1, new_b, SIGMA, MU, dec_st,
                prev_delta, prev_deriv2,weights,sdf)
    new_a=rs4$new_a
    delta=rs4$delta
    deriv2=rs4$deriv2

    #new_a=replace(new_a,new_a< -1,0)

    id_allzero=which(rowSums(new_a)==0)
    if(sum(rowSums(new_a)==0)!=0){
      new_a[id_allzero,] = par_a[id_allzero,]
      break
    }
    #new_a[id_allzero,] = par_a[id_allzero,]

    #update c
    par_c = new_c;
    c_num = colSums(u[id_1,]*(1-new_s[id_1,]))+Alpha - 1
    c_denom = length(id_1) + Alpha + Beta - 2
    new_c = (dec_st*c_num+(1-dec_st)*prev_c_num)/(dec_st*c_denom+(1-dec_st)*prev_c_denom)


    #save old sums for SAEM updates
    prev_a_num = (1-dec_st)*prev_a_num + dec_st*a_num
    prev_a_denom = (1-dec_st)*prev_a_denom + dec_st*a_denom

    prev_b_num = (1-dec_st)*prev_b_num + dec_st*b_num
    prev_b_denom = (1-dec_st)*prev_b_denom + dec_st*b_denom

    prev_c_num = (1-dec_st)*prev_c_num + dec_st*c_num
    prev_c_denom = (1-dec_st)*prev_c_denom + dec_st*c_denom

    prev_delta = (1-dec_st)*prev_delta + dec_st*delta
    prev_deriv2 = (1-dec_st)*prev_deriv2 + dec_st*deriv2


    prev_lb=(1-dec_st)*prev_lb+ dec_st*lb3pl(u, xi, new_s, id_1,
                                             new_a, new_c, Sigma, new_b, SIGMA, MU,
                                             Alpha, Beta, mu_b, sigma2_b)
    nlb<-append(nlb,prev_lb)

    if(n%%20==0 && n>20){
      n20=append(n20,abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n])))
      if(abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n]))<0.006){
        converge=0
      }
    }
    n=n+1
  }
  if (Matrix::rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
    #is_singular = 1
  }else{
    rsigma = Sigma
  }
  #new_a=replace(new_a,new_a< 0,0)
  new_a=replace(new_a,abs(new_a)< 0.001,0)
  Q_mat = (new_a != 0)*1
  #gic
  lbound=lb3pl(u,xi,new_s,c(1:person),new_a,new_c,Sigma, new_b, SIGMA, MU,
               Alpha, Beta, mu_b, sigma2_b)
  gic=log(log(person))*log(person)*sum(Q_mat) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat) - 2*lbound
  aic = 2*sum(Q_mat) -2*lbound
  return(list(ra=new_a,rb=new_b,rc=new_c,rs = new_s,reta = eta,reps=xi,rsigma = Sigma,rs = new_s,
              mu_i = MU,sig_i = SIGMA,n=n,Q_mat=Q_mat,GIC=gic,AIC=aic,BIC=bic))
}


#main function: choose optimal lambda
sgvem_3PLEFA_adaptive_const1_all<-function(u,domain,indic,samp,forgetrate,mu_b,sigma2_b,Alpha,Beta,gamma,max.iter){
  inite=dim(u)[1]
  person=dim(u)[1]
  item=dim(u)[2]
  lbd=seq(2,20,2)
  nopenalty_col=which(rowSums(indic)==1)
  #initialization
  initial=init(u,domain,indic)
  new_a = initial[[1]]
  new_b = initial[[2]]
  eta= initial[[3]]
  xi=initial[[4]]
  Sigma=initial[[5]]
  new_c=rep(0.15,item)
  new_s=matrix(0,nrow=person,ncol=item)
  aveg_p=colMeans(u==1)
  for (i in 1:person) {
    for (j in 1:item) {
      if(u[i,j]==0){
        new_s[i,j]=1;
      }else{
        new_s[i,j]=runif(1,aveg_p[j],1)
      }
    }
  }
  #create weights: use CFA
  wa=sgvem_3PLCFA(u,indic,samp,forgetrate,mu_b,sigma2_b,Alpha,Beta,max.iter)$ra
  weights = abs(wa)^gamma+1e-05
  rl<-vector("list",length(lbd))
  gic<-NULL
  for(j in 1:length(lbd)){
    rl [[j]]=sgvem_3PLEFA_adapt_const1(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,
                                       inite,samp,forgetrate,domain,lbd[j],indic,
                                       mu_b,sigma2_b,Alpha,Beta,weights,nopenalty_col,max.iter)
    lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                 rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
    gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
    new_a = rl[[j]]$ra
    new_b = rl[[j]]$rb
    new_c = rl[[j]]$rc
    new_s = rl[[j]]$rs
    eta= rl[[j]]$reta
    xi=rl[[j]]$reps
  }
  id=which.min(gic)
  #adjust the range of lambda
  if(id==1){
    lbd=seq(0.1,2,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=sgvem_3PLEFA_adapt_const1(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,
                                         inite,samp,forgetrate,domain,lbd[j],indic,
                                         mu_b,sigma2_b,Alpha,Beta,weights,nopenalty_col,
                                         max.iter)
      lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      new_c = rl[[j]]$rc
      new_s = rl[[j]]$rs
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
    }}
  else if(id==10){
    lbd=seq(20,40,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=sgvem_3PLEFA_adapt_const1(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,
                                         inite,samp,forgetrate,domain,lbd[j],indic,
                                         mu_b,sigma2_b,Alpha,Beta,weights,nopenalty_col,
                                         max.iter)
      lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      new_c = rl[[j]]$rc
      new_s = rl[[j]]$rs
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
    }
  }
  id=which.min(gic)
  rs=sgvem_3PLCFA(u,rl[[id]]$Q_mat,samp,forgetrate,mu_b,sigma2_b,Alpha,Beta,max.iter)
  rs$lbd=lbd[id]
  #rs$id=id
  return(rs)
}

#adaptive lasso with constraint 2 function
sgvem_3PLEFA_adapt_const2 <- function(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,forgetrate,domain,lbd,
                                      mu_b,sigma2_b,Alpha,Beta,weights,indic,
                                      nopenalty_col,max.iter) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 1
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(diag(domain),dim=c(domain,domain,person))
  #xi=eta
  prev_a_num = matrix(0,nrow=domain,ncol=item)
  prev_a_denom = array(0,dim=c(domain,domain,item))
  prev_delta = matrix(0,nrow=item,ncol=domain)
  prev_deriv2 = matrix(0,nrow=item,ncol=domain);
  prev_b_num = matrix(0,nrow=item,ncol=1)
  prev_b_denom = matrix(0,nrow=item,ncol=1)
  prev_c_num = matrix(0,nrow=item,ncol=1)
  prev_c_denom = 0
  n = 1
  par_Sigma = Sigma
  #lb=NULL
  #nv<-NULL
  nlb<-n20<-NULL
  prev_lb = 0
  while(converge==1 && Matrix::rankMatrix(Sigma) == domain && n < max.iter){
    if(n==1){
      dec_st=1
      id_1 = sample(person,inite)
    }else{
      dec_st = (n+1)^(-forgetrate)*0.25
      #dec_st = (n+1)^(-forgetrate)
      id_1 = sample(person,samp)
    }
    s_num = (1-new_c)/new_c
    par_Sigma = Sigma
    #update MU and SIGMA,xi,Spart
    rs1<-ecfa3pl(u, Sigma, domain, id_1, item, eta, new_s, new_a, new_b, SIGMA, MU, xi)
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    xi=rs1$xi
    Spart=rs1$Spart

    #update s
    new_s=s3pl(item,id_1,new_s,eta,new_a,s_num,new_b,u,SIGMA,MU,xi)
    #update eta
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=eta3pl(item,id_1,eta,abs(xi))

    Sigma=Spart/length(id_1)
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)

    #update b
    par_b=new_b
    rs2=b3pl(item,id_1,new_s,eta,u,new_a,dec_st,MU,prev_b_num,prev_b_denom,mu_b,sigma2_b)
    new_b=rs2$new_b
    b_denom= rs2$b_denom
    b_num = rs2$b_num


    par_a=new_a

    #update a
    #find the last one for each item by using indicator matrix
    lastone=apply(indic[nopenalty_col,], 1, function(x) tail(which(x!=0),1))
    rs3=nalc23pl(u, domain,item,new_a,nopenalty_col,lastone,new_s,eta,
                 new_b, SIGMA, MU, id_1, prev_a_num, prev_a_denom, dec_st)
    new_a=rs3$new_a
    a_denom=rs3$a_denom
    a_num=rs3$a_num
    #L1-penalty: off-diagnoal
    rs4=paalc23pl(u, new_a, nopenalty_col, new_s, lastone, lbd, id_1, eta, new_b,
                  SIGMA, MU, dec_st, prev_delta, prev_deriv2,weights)
    new_a=rs4$new_a
    delta=rs4$delta
    deriv2=rs4$deriv2

    #upper-tiangular should be zero
    new_a=replace(new_a,indic==0,0)
    #domain+1:item
    #find penaly columns
    pc=setdiff(1:item,nopenalty_col)
    rs5=paalc23pl1(u, new_a, new_s, lbd, id_1, eta, new_b, SIGMA, MU, dec_st,
                   prev_delta, prev_deriv2, pc, delta, deriv2,weights)
    new_a=rs5$new_a
    delta=rs5$delta
    deriv2=rs5$deriv2

    #new_a=replace(new_a,new_a< -1,0)

    id_allzero=which(rowSums(new_a)==0)
    if(sum(rowSums(new_a)==0)!=0){
      new_a[id_allzero,] = par_a[id_allzero,]
      break
    }

    #update c
    par_c = new_c;
    c_num = colSums(u[id_1,]*(1-new_s[id_1,]))+Alpha - 1
    c_denom = length(id_1) + Alpha + Beta - 2
    new_c = (dec_st*c_num+(1-dec_st)*prev_c_num)/(dec_st*c_denom+(1-dec_st)*prev_c_denom)


    #save old sums for SAEM updates
    prev_a_num = (1-dec_st)*prev_a_num + dec_st*a_num
    prev_a_denom = (1-dec_st)*prev_a_denom + dec_st*a_denom

    prev_b_num = (1-dec_st)*prev_b_num + dec_st*b_num
    prev_b_denom = (1-dec_st)*prev_b_denom + dec_st*b_denom

    prev_c_num = (1-dec_st)*prev_c_num + dec_st*c_num
    prev_c_denom = (1-dec_st)*prev_c_denom + dec_st*c_denom

    prev_delta = (1-dec_st)*prev_delta + dec_st*delta
    prev_deriv2 = (1-dec_st)*prev_deriv2 + dec_st*deriv2


    prev_lb=(1-dec_st)*prev_lb+ dec_st*lb3pl(u, xi, new_s, id_1,
                                             new_a, new_c, Sigma, new_b, SIGMA, MU,
                                             Alpha, Beta, mu_b, sigma2_b)
    nlb<-append(nlb,prev_lb)

    if(n%%20==0 && n>20){
      n20=append(n20,abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n])))
      if(abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n]))<0.006){
        converge=0
      }
    }
    n=n+1
  }
  if (Matrix::rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
    #is_singular = 1
  }else{
    rsigma = Sigma
  }
  #new_a=replace(new_a,new_a< 0,0)
  new_a=replace(new_a,abs(new_a)< 0.001,0)
  Q_mat = (new_a != 0)*1
  lbound=lb3pl(u,xi,new_s,c(1:person),new_a,new_c,Sigma, new_b, SIGMA, MU,
               Alpha, Beta, mu_b, sigma2_b)
  gic=log(log(person))*log(person)*sum(Q_mat) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat) - 2*lbound
  aic = 2*sum(Q_mat) -2*lbound
  return(list(ra=new_a,rb=new_b,rc=new_c,rs = new_s,reta = eta,reps=xi,rsigma = Sigma,rs = new_s,
              mu_i = MU,sig_i = SIGMA,n=n,Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic))
}


#main function: choose optimal lambda
sgvem_3PLEFA_adaptive_const2_all<-function(u,domain,samp,forgetrate,
                                           mu_b,sigma2_b,Alpha,Beta,indic,non_pen,gamma,
                                           max.iter){
  lbd=seq(2,20,2)
  inite=dim(u)[1]
  person=dim(u)[1]
  item=dim(u)[2]
  nopenalty_col=c(which(rowSums(indic)<domain),non_pen)
  #initialization
  initial=init(u,domain,indic)
  new_a = initial[[1]]
  new_b = initial[[2]]
  eta= initial[[3]]
  xi=initial[[4]]
  Sigma=initial[[5]]
  new_c=rep(0.15,item)
  new_s=matrix(0,nrow=person,ncol=item)
  aveg_p=colMeans(u==1)
  for (i in 1:person) {
    for (j in 1:item) {
      if(u[i,j]==0){
        new_s[i,j]=1;
      }else{
        new_s[i,j]=runif(1,aveg_p[j],1)
      }
    }
  }
  #create weights: use CFA
  wa=sgvem_3PLCFA(u,indic,samp,forgetrate,mu_b,sigma2_b,Alpha,Beta,max.iter)$ra
  weights = abs(wa)^gamma+1e-05
  rl<-vector("list",length(lbd))
  gic<-NULL
  for(j in 1:length(lbd)){
    rl [[j]]=sgvem_3PLEFA_adapt_const2(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,
                                       forgetrate,domain,lbd[j],mu_b,sigma2_b,Alpha,Beta,
                                       weights,indic,nopenalty_col,max.iter)
    lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                 rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
    gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
    new_a = rl[[j]]$ra
    new_b = rl[[j]]$rb
    new_c = rl[[j]]$rc
    new_s = rl[[j]]$rs
    eta= rl[[j]]$reta
    xi=rl[[j]]$reps
  }
  id=which.min(gic)
  #adjust the range of lambda
  if(id==1){
    lbd=seq(0.1,2,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=sgvem_3PLEFA_adapt_const2(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,
                                         forgetrate,domain,lbd[j],mu_b,sigma2_b,Alpha,Beta,
                                         weights,indic,nopenalty_col,max.iter)
      lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      new_c = rl[[j]]$rc
      new_s = rl[[j]]$rs
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
    }}
  else if(id==10){
    lbd=seq(20,40,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=sgvem_3PLEFA_adapt_const2(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,
                                         forgetrate,domain,lbd[j],mu_b,sigma2_b,Alpha,Beta,
                                         weights,indic,nopenalty_col,max.iter)
      lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      new_c = rl[[j]]$rc
      new_s = rl[[j]]$rs
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
    }
  }
  id=which.min(gic)
  rs=sgvem_3PLCFA(u, rl[[id]]$Q_mat,samp,forgetrate,mu_b,sigma2_b,Alpha,Beta,max.iter)
  rs$lbd=lbd[id]
  #rs$id=id
  return(rs)
}
