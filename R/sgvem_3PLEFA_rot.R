#' Stochastic GVEM for Exploratory M3PL Analysis
#'
#' @param u a \eqn{N \times J} \code{matrix} or a \code{data.frame} that
#' consists of binary responses of \eqn{N} individuals to \eqn{J} items. The
#' missing values are coded as \code{NA}
#' @param domain the number of factors
#' @param samp a subsample for each iteration; default is 50
#' @param forgetrate the forget rate for the stochastic algorithm. The value should be within the
#' range from 0.5 to 1. Default is 0.51
#' @param mu_b the mean parameter for the prior distribution of item difficulty parameters
#' @param sigma2_b the variance parameter for the prior distribution of item difficulty parameters
#' @param Alpha the \eqn{\alpha} parameter for the prior distribution of guessing parameters
#' @param Beta the \eqn{\beta} parameter for the prior distribution of guessing parameters
#' @param max.iter the maximum number of iterations for the EM cycle; default is 5000
#' @param rot the post-hoc rotation method: Promax or CF-Quartimax; default is \code{"Promax"},
#' but may also be \code{"cfQ"} for conducting the CF-Quartimax rotation
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
#'   \item{rk}{factor loadings, a \eqn{J \times K} \code{matrix}}
#'   \item{GIC}{model fit index}
#'   \item{AIC}{model fit index}
#'   \item{BIC}{model fit index}
#'   \item{ur_a}{item discrimination parameters before conducting the rotation, a \eqn{J \times K} \code{matrix}}
#' @seealso \code{\link{sgvem_3PLEFA_lasso}}, \code{\link{sgvem_3PLEFA_adaptlasso}}
#' @export
#'
#' @examples sgvem_3PLEFA_rot(exampleData_3pl, 3,samp=50,forgetrate=0.51,
#' mu_b=0,sigma2_b=4,Alpha=10,Beta=40,max.iter=5000,rot="Promax")
sgvem_3PLEFA_rot <- function(u,domain, samp=50,forgetrate=0.51,mu_b,sigma2_b,Alpha,Beta,max.iter=5000,rot="Promax") {
  start=Sys.time()
  u=data.matrix(u)
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 1
  aveg_p=colMeans(u==1)
  inite=person
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
  new_c=rep(0.1,item)
  #person*item
  xi=array(1,person)%*%t(new_b)-theta%*%t(new_a)
  eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
  eta[is.na(u)]=NA
  new_s=matrix(NA,nrow=person,ncol=item)
  for (i in 1:person) {
    for (j in 1:item) {
      new_s[i,j]<-ifelse(u[i,j]==0,1,runif(1,aveg_p[j],1))
    }
  }

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

  n = 1
  #lb=NULL
  prev_lb = 0
  nlb<-n20<-NULL
  if(rot=="Promax"){
    while(converge==1 &&  n<max.iter){
      if(n==1){
        dec_st=1
        id_1 = sample(person,inite)
      }else{
        #dec_st = (n+1)^(-forgetrate)*0.25
        dec_st = (n+1)^(-forgetrate)
        id_1 = sample(person,samp)
      }
      s_num = (1-new_c)/new_c
      #par_Sigma=Sigma
      #update MU and SIGMA,xi,Spart
      rs1<-eefa3(u, domain, id_1, item, eta, new_s, new_a, new_b, SIGMA, MU, xi)

      SIGMA=rs1$SIGMA
      MU=rs1$MU
      xi=rs1$xi

      #update s
      new_s=s3pl(item, id_1, new_s, eta, new_a, s_num, new_b, u, SIGMA, MU, xi)

      #update eta
      eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
      eta=eta3pl(item,id_1,eta,abs(xi))


      #update b
      par_b=new_b

      rs2=b3pl(item,id_1,new_s,eta,u,new_a,dec_st,MU,prev_b_num,prev_b_denom,mu_b,sigma2_b)
      new_b=rs2$new_b
      b_denom= rs2$b_denom
      b_num = rs2$b_num
      #update a
      par_a=new_a
      rs3=aefa3(u, id_1, item, domain, eta, new_a, new_b, new_s, SIGMA, MU, prev_a_num, prev_a_denom, dec_st)
      new_a=rs3$new_a
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

      prev_lb=(1-dec_st)*prev_lb+ dec_st*lb3pl(u, xi, new_s, id_1,new_a, new_c,diag(domain), new_b, SIGMA, MU,
                                               Alpha, Beta, mu_b, sigma2_b)
      nlb<-append(nlb,prev_lb)

      if(n%%20==0 && n>20 && testit::has_error(psych::Promax(new_a),silent=T)==0){
        n20=append(n20,abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n])))
        if(abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n]))<0.006 && sum(is.na(psych::Promax(new_a)$loadings[1:item,]))==0){
          converge=0
        }
      }
      n=n+1
    }
    #rotate factors: promax
    ur_a=new_a
    rotated=psych::Promax(new_a)
  }
  else{
    while(converge==1 &&  n<max.iter){
      if(n==1){
        dec_st=1
        id_1 = sample(person,inite)
      }else{
        #dec_st = (n+1)^(-forgetrate)*0.25
        dec_st = (n+1)^(-forgetrate)
        id_1 = sample(person,samp)
      }
      s_num = (1-new_c)/new_c
      #par_Sigma=Sigma
      #update MU and SIGMA,xi,Spart
      rs1<-eefa3(u, domain, id_1, item, eta, new_s, new_a, new_b, SIGMA, MU, xi)

      SIGMA=rs1$SIGMA
      MU=rs1$MU
      xi=rs1$xi

      #update s
      new_s=s3pl(item, id_1, new_s, eta, new_a, s_num, new_b, u, SIGMA, MU, xi)

      #update eta
      eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
      eta=eta3pl(item,id_1,eta,abs(xi))


      #update b
      par_b=new_b

      rs2=b3pl(item,id_1,new_s,eta,u,new_a,dec_st,MU,prev_b_num,prev_b_denom,mu_b,sigma2_b)
      new_b=rs2$new_b
      b_denom= rs2$b_denom
      b_num = rs2$b_num
      #update a
      par_a=new_a
      rs3=aefa3(u, id_1, item, domain, eta, new_a, new_b, new_s, SIGMA, MU, prev_a_num, prev_a_denom, dec_st)
      new_a=rs3$new_a
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

      prev_lb=(1-dec_st)*prev_lb+ dec_st*lb3pl(u, xi, new_s, id_1,new_a, new_c,diag(domain), new_b, SIGMA, MU,
                                               Alpha, Beta, mu_b, sigma2_b)
      nlb<-append(nlb,prev_lb)

      if(n%%20==0 && n>20 && testit::has_error(GPArotation::cfQ(new_a),silent=T)==0){
        n20=append(n20,abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n])))
        if(abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n]))<0.006 && sum(is.na(GPArotation::cfQ(new_a)$loadings[1:item,]))==0){
          converge=0
        }
      }
      n=n+1
    }
    ur_a=new_a
    rotated=GPArotation::cfQ(new_a)
  }

  new_a=rotated$loadings[1:item,]
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
              reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,
              sig_i = SIGMA,n=n,Q_mat=Q_mat,rk=rk,GIC=gic,AIC=aic,
              BIC=bic,ur_a=ur_a))
}
