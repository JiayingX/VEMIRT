#' Exploratory M2PL Analysis with Lasso Penalty
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
#' @param max.iter the maximum number of iterations for the EM cycle; default is 5000
#' @param constrain the constraint setting: \code{"C1"} or \code{"C2"}. To ensure
#' identifiability, \code{"C1"} sets a \eqn{K \times K} sub-matrix of \code{indic} to be an
#' identity matrix, and the variance of the population variance-covariance matrix to be 1.
#' Instead of setting all off-diagonals of a sub-matrix to be zero, \code{"C2"} keeps the
#' sub-matrix to be a triangular matrix with the diagonal being ones. That is, there
#' are test items associated with each factor for sure and they may be associated
#' with other factors as well. Nonzero entries except for the diagonal entries of the
#' sub-matrix are penalized during the estimation procedure. Although, \code{"C2"}
#' is much weaker than \code{"C1"}, it still ensures empirical identifiability. Default is \code{"C1"}
#' @param non_pen the index of an item which is associated with each factor to satisfy \code{"C2"}.
#' For \code{C1}, the input can be \code{NULL}
#' @return a list containing the following objects:
#'   \item{ra}{item discrimination parameters, a \eqn{J \times K} \code{matrix}}
#'   \item{rb}{item difficulty parameters, vector of length \eqn{J}}
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
#' @seealso \code{\link{gvem_2PLEFA_rot}}, \code{\link{gvem_2PLEFA_adaptlasso}}, \code{\link{exampleIndic_efa2pl_c1}}, \code{\link{exampleIndic_efa2pl_c2}}
#' @export
#'
#' @examples
#' gvem_2PLEFA_lasso(exampleData_2pl, exampleIndic_efa2pl_c1,constrain="C1")
#' @examples
#' gvem_2PLEFA_lasso(exampleData_2pl, exampleIndic_efa2pl_c2,constrain="C2",non_pen=61)

#main function for gvem_2PLEFA_lasso
gvem_2PLEFA_lasso<-function(u,indic,max.iter=5000,constrain="C1",non_pen=NULL){
  u=data.matrix(u)
  indic=data.matrix(indic)
  domain=dim(indic)[2]
  if(constrain=="C1"){
    result=vem_2PLEFA_L1_const1_all(u,domain,indic,max.iter)
  }else{
    result=vem_2PLEFA_L1_const2_all(u,domain,indic,non_pen,max.iter)
  }
  if(result$lbd==0.1 || result$lbd==40){
    warning("The optimal penalty parameter may be out of range.",call. = FALSE)
  }
  return(result)
}



###l1
vem_2PLEFA_L1_const1 <- function(u,new_a,new_b,eta,xi,Sigma, domain,lbd,indic,nopenalty_col,max.iter) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 0
  is_singular = 0

  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(0,dim=c(domain,domain,person))
  n = 0
  par_Sigma = Sigma
  #nv<-NULL
  while(converge==0 && Matrix::rankMatrix(Sigma) == domain && n < max.iter){
    par_Sigma = Sigma
    #update MU, SIGMA, Sigma
    rs1<-ecfa2(u,Sigma,domain, person, item, eta, new_a, new_b)
    xi=rs1$xi
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    Spart=rs1$Spart
    Sigma=Spart/person
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(xi)< 0.01,0.125)


    #update b
    par_b=new_b
    b_part = t(new_a%*%MU)
    new_b=colSums(0.5-u+2*eta*b_part)/colSums(2*eta)

    par_a=new_a
    #update a
    new_a1=nalc12pl(u, indic, nopenalty_col, person, eta, new_b, SIGMA, MU)
    new_a1[-nopenalty_col,]=new_a[-nopenalty_col,]
    #L1-penalty
    sdf=setdiff(1:item,nopenalty_col)
    new_a=palc12pl(u, domain, person, lbd, sdf, eta, new_a1, new_b, SIGMA, MU)
    #par_a=new_a2
    #nv<-append(nv,norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")+
    #             norm(as.vector(Sigma)-as.vector(par_Sigma),type="2"))
    if (norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")+
        norm(as.vector(Sigma)-as.vector(par_Sigma),type="2")<0.001){
      converge=1
    }
    n=n+1
  }
  if (Matrix::rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
    is_singular = 1
  }else{
    rsigma = Sigma
  }
  #new_a=replace(new_a,new_a< 0,0)
  new_a=replace(new_a,abs(new_a)< 0.001,0)
  Q_mat = (new_a != 0)*1
  #gic
  lbound=lb2pl(u,xi,Sigma,new_a,new_b,SIGMA,MU)
  gic=log(log(person))*log(person)*sum(Q_mat) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat) - 2*lbound
  aic = 2*sum(Q_mat) -2*lbound
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,sig_i = SIGMA,n=n,
              Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic))
}


####main function for l1
vem_2PLEFA_L1_const1_all<-function(u,domain,indic,max.iter){
  lbd=seq(2,20,2)
  person=dim(u)[1]
  item=dim(u)[2]
  nopenalty_col=which(rowSums(indic)==1)
  #initialization
  initial=init(u,domain,indic)
  new_a = initial[[1]]
  new_b = initial[[2]]
  eta= initial[[3]]
  xi=initial[[4]]
  Sigma=initial[[5]]
  rl<-vector("list",length(lbd))
  gic<-NULL
  for(j in 1:length(lbd)){
    r0=vem_2PLEFA_L1_const1(u,new_a,new_b,eta,xi,Sigma,domain,lbd[j],indic,nopenalty_col,max.iter)
    rl [[j]]=gvem_2PLCFA(u,r0$Q_mat,max.iter)
    lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                 rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
    gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
    new_a = r0$ra
    new_b = r0$rb
    eta= r0$reta
    xi=r0$reps
    Sigma=r0$rsigma
  }
  id=which.min(gic)
  #adjust the range of lambda
  if(id==1){
    lbd=seq(0.1,2,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      r0=vem_2PLEFA_L1_const1(u,new_a,new_b,eta,xi,Sigma,domain,lbd[j],indic,nopenalty_col,max.iter)
      rl [[j]]=gvem_2PLCFA(u,r0$Q_mat,max.iter)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
      new_a = r0$ra
      new_b = r0$rb
      eta= r0$reta
      xi=r0$reps
      Sigma=r0$rsigma
    }}
  else if(id==10){
    lbd=seq(20,40,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      r0=vem_2PLEFA_L1_const1(u,new_a,new_b,eta,xi,Sigma,domain,lbd[j],indic,nopenalty_col,max.iter)
      rl [[j]]=gvem_2PLCFA(u,r0$Q_mat,max.iter)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
      new_a = r0$ra
      new_b = r0$rb
      eta= r0$reta
      xi=r0$reps
      Sigma=r0$rsigma
    }
  }
  id=which.min(gic)
  rs=rl[[id]]
  rs$lbd=lbd[id]
  #rs$id=id
  return(rs)}


#####l2
vem_2PLEFA_L1_const2 <- function(u,new_a,new_b,eta,xi,Sigma,domain,lbd,
                                 indic,nopenalty_col,max.iter) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 0
  is_singular = 0

  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(0,dim=c(domain,domain,person))
  n = 0

  par_Sigma = Sigma
  while(converge==0 && Matrix::rankMatrix(Sigma) == domain && n < max.iter){
    par_Sigma = Sigma
    #update Sigma, MU, SIGMA
    rs1<-ecfa2(u,Sigma,domain, person, item, eta, new_a, new_b)
    xi=rs1$xi
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    Spart=rs1$Spart
    Sigma=Spart/person
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(xi)< 0.01,0.125)


    #update b
    par_b=new_b
    b_part = t(new_a%*%MU)
    new_b=colSums(0.5-u+2*eta*b_part)/colSums(2*eta)

    par_a=new_a
    #update a
    #find the last one for each item by using indicator matrix
    lastone=apply(indic[nopenalty_col,], 1, function(x) tail(which(x!=0),1))
    new_a=nalc22pl(u, domain,new_a,nopenalty_col,lastone, person, eta, new_b, SIGMA, MU)
    #L1-penalty: off-diagnoal
    new_a=palc22pl(u, new_a,nopenalty_col,lastone, lbd, person, eta, new_b, SIGMA, MU)
    #upper-tiangular should be zero
    new_a=replace(new_a,indic==0,0)
    #domain+1:item
    #find penaly columns
    pc=setdiff(1:item,nopenalty_col)
    new_a=palc22pl1(u, domain, item, person, lbd, eta, new_a, new_b, SIGMA, MU,pc)
    #new_a=replace(new_a,new_a< 0,0)
    #par_a=new_a2
    if (norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")+
        norm(as.vector(Sigma)-as.vector(par_Sigma),type="2")<0.001){
      converge=1
    }
    n=n+1
  }
  if (Matrix::rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
    is_singular = 1
  }else{
    rsigma = Sigma
  }
  new_a=replace(new_a,abs(new_a)< 0.001,0)
  Q_mat = (new_a != 0)*1
  #gic
  lbound=lb2pl(u,xi,Sigma,new_a,new_b,SIGMA,MU)
  gic=log(log(person))*log(person)*sum(Q_mat) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat) - 2*lbound
  aic = 2*sum(Q_mat) -2*lbound
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,sig_i = SIGMA,n=n,
              Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic))
}


####main function for l2
vem_2PLEFA_L1_const2_all<-function(u,domain,indic,non_pen,max.iter){
  lbd=seq(2,20,2)
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
  rl<-vector("list",length(lbd))
  gic<-NULL
  for(j in 1:length(lbd)){
    r0=vem_2PLEFA_L1_const2(u,new_a,new_b,eta,xi,Sigma,  domain,lbd[j],indic,nopenalty_col,max.iter)
    rl [[j]]=gvem_2PLCFA(u,r0$Q_mat,max.iter)
    lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                 rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
    gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
    new_a = r0$ra
    new_b = r0$rb
    eta= r0$reta
    xi=r0$reps
    Sigma=r0$rsigma
  }
  id=which.min(gic)
  #adjust the range of lambda
  if(id==1){
    lbd=seq(0.1,2,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      r0=vem_2PLEFA_L1_const2(u,new_a,new_b,eta,xi,Sigma,  domain,lbd[j],indic,nopenalty_col,max.iter)
      rl [[j]]=gvem_2PLCFA(u, r0$Q_mat,max.iter)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
      new_a = r0$ra
      new_b = r0$rb
      eta= r0$reta
      xi=r0$reps
      Sigma=r0$rsigma
    }}
  else if(id==10){
    lbd=seq(20,40,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      r0=vem_2PLEFA_L1_const2(u,new_a,new_b,eta,xi,Sigma,  domain,lbd[j],indic,nopenalty_col,max.iter)
      rl [[j]]=gvem_2PLCFA(u, r0$Q_mat,max.iter)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
      new_a = r0$ra
      new_b = r0$rb
      eta= r0$reta
      xi=r0$reps
      Sigma=r0$rsigma
    }
  }
  id=which.min(gic)
  rs=rl[[id]]
  rs$lbd=lbd[id]
  #rs$id=id
  return(rs)
}
