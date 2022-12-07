#' Stochastic GVEM for Confirmatory M3PL Analysis
#'
#' @param u a \eqn{N \times J} \code{matrix} or a \code{data.frame} that
#' consists of binary responses of \eqn{N} individuals to \eqn{J} items. The
#' missing values are coded as \code{NA}
#' @param indic a \eqn{J \times K} \code{matrix} or a \code{data.frame} that
#' describes the factor loading structure of \eqn{J} items to \eqn{K} factors. It
#' consists of binary values where 0 refers to the item is irrelevant with this factor,
#' 1 otherwise
#' @param samp a subsample for each iteration; default is 50
#' @param forgetrate the forget rate for the stochastic algorithm. The value should be within the
#' range from 0.5 to 1. Default is 0.51
#' @param mu_b the mean parameter for the prior distribution of item difficulty parameters
#' @param sigma2_b the variance parameter for the prior distribution of item difficulty parameters
#' @param Alpha the \eqn{\alpha} parameter for the prior distribution of guessing parameters
#' @param Beta the \eqn{\beta} parameter for the prior distribution of guessing parameters
#' @param max.iter the maximum number of iterations for the EM cycle; default is 5000
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
#' @references
#' Cho, A. E., Wang, C., Zhang, X., & Xu, G. (2021). Gaussian variational estimation for multidimensional item response theory. \emph{British Journal of Mathematical and Statistical Psychology, 74}, 52-85.
#'
#' Cho, A. E., Xiao, J., Wang, C., & Xu, G. (2022). Regularized Variational Estimation for Exploratory Item Factor Analysis. \emph{Psychometrika}. https://doi.org/10.1007/s11336-022-09874-6
#' @seealso \code{\link{gvem_2PLCFA}}
#' @export
#'
#' @examples sgvem_3PLCFA(exampleData_3pl, exampleIndic_cfa3pl,samp=50,forgetrate=0.51,
#' mu_b=0,sigma2_b=4,Alpha=10,Beta=40)
sgvem_3PLCFA <- function(u,indic,samp=50,forgetrate=0.51,mu_b,sigma2_b,Alpha,Beta,max.iter=5000) {
  start=Sys.time()
  u=data.matrix(u)
  indic=data.matrix(indic)
  domain=dim(indic)[2]
  person=dim(u)[1]
  item=dim(u)[2]
  SIGMA=array(0,dim=c(domain,domain,person))
  for (i in 1:person) {
    SIGMA[,,i]=diag(domain)
  }
  MU=matrix(0,nrow=domain,ncol=person)
  prev_a_num = matrix(0,nrow=domain,ncol=item)
  prev_a_denom = array(0,dim=c(domain,domain,item))
  prev_b_num = matrix(0,nrow=item,ncol=1)
  prev_b_denom = matrix(0,nrow=item,ncol=1)
  prev_c_num = matrix(0,nrow=item,ncol=1)
  prev_c_denom = 0
  prev_lb = 0
  converge = 1
  n = 1
  #initialization
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
  initial=init(u,domain,indic)
  new_a = initial[[1]]
  new_b = initial[[2]]
  eta= initial[[3]]
  Sigma=initial[[5]]
  xi=eta
  inite=person
  prev_lb = 0
  nlb<-n20<-NULL
  while(converge==1 && Matrix::rankMatrix(Sigma) == domain && n<max.iter){
    if(n==1){
      #if(n<21){
      dec_st=1
      id_1 = sample(person,inite)
    }else{
      dec_st = (n+1)^(-forgetrate)*0.25
      id_1 = sample(person,samp)
    }
    s_num = (1-new_c)/new_c
    par_Sigma=Sigma
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
    #update Sigma
    Sigma=Spart/length(id_1)
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)

    #update b
    par_b=new_b
    rs2=b3pl(item,id_1,new_s,eta,u,new_a,dec_st,MU,prev_b_num,prev_b_denom,mu_b,sigma2_b)
    new_b=rs2$new_b
    b_denom= rs2$b_denom
    b_num = rs2$b_num

    #update a
    par_a=new_a
    rs3=acfa3(u, indic, id_1,item, domain, eta, new_a,
              new_b, new_s, SIGMA, MU, prev_a_num, prev_a_denom, dec_st)
    new_a=rs3$new_a
    #2022/01/19:rescale a
    new_a=new_a%*%d_temp
    a_denom= rs3$a_denom
    a_num = rs3$a_num

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

    #check the lower bound
    old=prev_lb
    prev_lb=(1-dec_st)*prev_lb+ dec_st*lb3pl(u, xi, new_s, id_1,
                                             new_a, new_c, Sigma, new_b, SIGMA, MU,
                                             Alpha, Beta, mu_b, sigma2_b)
    nlb<-append(nlb,prev_lb)

    #calculate the mean lower bound difference for every 20 iterations
    if(n%%20==0 && n>20){
      n20=append(n20,abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n])))
      if(abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n]))<0.006){
        #if (nrm1<0.32){
        converge=0
      }
    }
    n=n+1
  }
  if (Matrix::rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
  }else{
    rsigma = Sigma
  }
  new_a=replace(new_a,abs(new_a)< 0.001,0)
  Q_mat = (new_a != 0)*1
  #gic
  lbound=lb3pl(u,xi,new_s,person,new_a, new_c, Sigma, new_b, SIGMA, MU,Alpha, Beta, mu_b, sigma2_b)
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
  return(list(ra=new_a,rb=new_b,rc=new_c,rs = new_s,
              reta = eta,reps=xi,rsigma = rsigma,mu_i = MU,sig_i = SIGMA,n=n,
              Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic))
}
