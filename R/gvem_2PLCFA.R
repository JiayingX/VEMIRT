#' Confirmatory M2PL Analysis
#'
#' @param u a \eqn{N \times J} \code{matrix} or a \code{data.frame} that
#' consists of binary responses of \eqn{N} individuals to \eqn{J} items. The
#' missing values are coded as \code{NA}
#' @param indic a \eqn{J \times K} \code{matrix} or a \code{data.frame} that
#' describes the factor loading structure of \eqn{J} items to \eqn{K} factors. It
#' consists of binary values where 0 refers to the item is irrelevant with this factor,
#' 1 otherwise
#' @param max.iter the maximum number of iterations for the EM cycle; default is 5000
#'
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
#' @seealso \code{\link{sgvem_3PLCFA}}
#' @export
#'
#' @examples
#' gvem_2PLCFA(exampleData_2pl, exampleIndic_cfa2pl)
gvem_2PLCFA <- function(u,indic,max.iter=5000) {
  u=data.matrix(u)
  indic=data.matrix(indic)
  domain=dim(indic)[2]
  person=dim(u)[1]
  item=dim(u)[2]
  #initialization
  initial=init(u,domain,indic)
  new_a = initial[[1]] * indic
  new_b = initial[[2]]
  eta= initial[[3]]
  xi=initial[[4]]
  Sigma=initial[[5]]
  converge = 0
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(0,dim=c(domain,domain,person))
  n = 0
  par_Sigma = Sigma
  while(converge==0 && Matrix::rankMatrix(Sigma) == domain && n < max.iter){
    par_Sigma = Sigma
    #update MU, SIGMA, sigma, eta
    rs1<-ecfa2(u,Sigma,domain, person, item, eta, new_a, new_b)
    xi=rs1$xi
    xi[is.na(u)]=NA
    Spart=rs1$Spart
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(eta)< 0.01,0.125)
    Sigma=Spart/person
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)

    #update b
    par_b=new_b
    b_part = t(new_a%*%MU)
    new_b=colSums(0.5-u+2*eta*b_part,na.rm = T)/colSums(2*eta,na.rm = T)
    par_a=new_a
    #update a
    new_a=acfa2(u, indic, person, item, domain, eta, new_b, SIGMA, MU,new_a)
    #2022/01/19:rescale a
    new_a=new_a%*%d_temp
    if (norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")+
        norm(as.vector(Sigma)-as.vector(par_Sigma),type="2") <0.0001){
      converge=1
    }
    n=n+1
  }
  if (Matrix::rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
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
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,
              mu_i = MU,sig_i = SIGMA,n=n,Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic))
}
