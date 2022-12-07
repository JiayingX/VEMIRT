#' Exploratory M2PL Analysis with Post-hoc Rotation
#'
#' @param u a \eqn{N \times J} \code{matrix} or a \code{data.frame} that
#' consists of binary responses of \eqn{N} individuals to \eqn{J} items. The
#' missing values are coded as \code{NA}
#' @param domain the number of factors
#' @param max.iter the maximum number of iterations for the EM cycle; default is 5000
#' @param rot the post-hoc rotation method: Promax or CF-Quartimax; default is \code{"Promax"},
#' but may also be \code{"cfQ"} for conducting the CF-Quartimax rotation
#' @return a list containing the following objects:
#'   \item{ra}{item discrimination parameters, a \eqn{J \times K} \code{matrix}}
#'   \item{rb}{item difficulty parameters, vector of length \eqn{J}}
#'   \item{reta}{variational parameters \eqn{\eta(\xi)}, a \eqn{N \times J} matrix}
#'   \item{reps}{variational parameters \eqn{\xi}, a \eqn{N \times J} matrix}
#'   \item{rsigma}{population variance-covariance matrix, a \eqn{K \times K} matrix}
#'   \item{mu_i}{mean parameter for each person, a \eqn{K \times N} matrix}
#'   \item{sig_i}{covariance matrix for each person, a \eqn{K \times K \times N} array}
#'   \item{n}{the number of iterations for the EM cycle}
#'   \item{rk}{factor loadings, a \eqn{J \times K} \code{matrix}}
#'   \item{Q_mat}{factor loading structure, a \eqn{J \times K} matrix}
#'   \item{GIC}{model fit index}
#'   \item{AIC}{model fit index}
#'   \item{BIC}{model fit index}
#'   \item{ur_a}{item discrimination parameters before conducting the rotation, a \eqn{J \times K} \code{matrix}}
#' @seealso \code{\link{gvem_2PLEFA_lasso}}, \code{\link{gvem_2PLEFA_adaptlasso}}
#' @export
#'
#' @examples
#' gvem_2PLEFA_rot(exampleData_2pl, domain=5,max.iter=3000)
#' @examples
#' gvem_2PLEFA_rot(exampleData_2pl, domain=5,rot="cfQ")
gvem_2PLEFA_rot <- function(u, domain,max.iter=5000,rot="Promax") {
  start=Sys.time()
  u=data.matrix(u)
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 1
  #initialization
  r=identify(u)
  person=dim(u)[1]
  r[r>0.9]=0.9
  r[r<0]=abs(r[r<0][1])
  new_a=t(rep(1,domain)%o%(r/sqrt(1-r^2)))
  new_a=replace(new_a,new_a>4,4)
  theta=matrix(rnorm(person*domain,0,1),nrow=person)
  if(is.na(sum(psych::Promax(new_a)$loadings)) || is.na(sum(GPArotation::cfQ(new_a)$loadings))){
    sda=SVD.fn(u,domain)
    new_a=sda$A
    theta=sda$Theta
  }
  new_b=-qnorm(colSums(u,na.rm=T)/person,0,1)/r
  new_b[new_b>4]=4
  new_b[new_b<(-4)]=-4
  Sigma = diag(domain)

  #person*item
  xi=array(1,person)%*%t(new_b)-theta%*%t(new_a)
  eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
  eta[is.na(u)]=NA

  # initial=init(u,domain,indic)
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(0,dim=c(domain,domain,person))
  n = 0
  # new_a = initial[[1]] * indic
  # new_b = initial[[2]]
  # eta= initial[[3]]
  # xi=initial[[4]]
  # Sigma=initial[[5]]
  # par_Sigma = Sigma
  while(converge==1 && n < max.iter){
    #update mu and sigma for each person
    rs1<-eefa2(u,domain, person, item, eta, new_a, new_b)
    xi=rs1$xi
    xi[is.na(u)]=NA
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(xi)< 0.01,0.125)
    #Sigma=Spart/person
    #d_temp=sqrt(diag(diag(Sigma)))
    #Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)

    #update b
    par_b=new_b
    b_part = t(new_a%*%MU)
    new_b=colSums(0.5-u+2*eta*b_part,na.rm = T)/colSums(2*eta,na.rm = T)
    par_a=new_a
    #update a
    new_a=aefa2(u,domain,person,item,eta,new_b,SIGMA,MU)
    if (norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2") <0.0001){
      converge=0
    }
    n=n+1
  }
  #rotate factors: promax
  ur_a=new_a
  if(rot=="Promax"){rotated=psych::Promax(new_a)}
  else{rotated=GPArotation::cfQ(new_a)}
  new_a=rotated$loadings[1:item,]
  #cut-off points:0.2
  #new_a=replace(new_a,abs(new_a)< 0.2,0)
  Sigma=rotated$Phi
  #change the sign
  rt1<-rt(new_a,Sigma)
  new_a=rt1$ra
  Sigma=rt1$rsigma
  #calculate the factor loadings
  rk<-matrix(NA,item,domain)
  for (j in 1:item) {
    ra1=new_a[j,]
    B=diag(domain)+ra1%*%t(ra1)*Sigma
    rk[j,]=t(solve(chol(B))%*%ra1)
  }
  Q_mat<-(abs(rk)>0.3)*1
  #Q_mat1<-(new_a!=0)*1
  #gic
  lbound=lb2pl(u,xi,Sigma,new_a,new_b,SIGMA,MU)
  gic=log(log(person))*log(person)*sum(Q_mat) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat) - 2*lbound
  aic = 2*sum(Q_mat) -2*lbound
  if(n==max.iter){
    warning("The maximum number of EM cycles reached!",call.=FALSE)
  }
  end=Sys.time()
  duration=end-start
  cat(paste("Total Execution Time:", round(duration[[1]], 2),  units(duration)),"\n")
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,
              sig_i = SIGMA,n=n,rk=rk,Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic,ur_a=ur_a))
}
