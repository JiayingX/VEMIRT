#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List eefa3(const arma::mat&u,const int& domain, const arma::vec& id_1,
           const int& item, const arma::mat& eta, const arma::mat& new_s,const arma::mat& new_a,const arma::vec& new_b,
           const arma::cube& SIGMA1,const arma::mat& MU1,const arma::mat& xi1){
  //e step
  arma::mat MU=MU1;
  arma::cube SIGMA=SIGMA1;
  arma::mat xi=xi1;
  for(int i=0; i<id_1.n_elem; ++i){
    int k=id_1(i)-1;
    arma::mat sigma_part=arma::mat(domain,domain,arma::fill::zeros);
    //mu_part.zeros();
    arma::vec mu_part= arma::zeros(domain);
    for(int j=0; j<item; ++j){
      sigma_part=sigma_part+eta(k,j)*(1-u(k,j)+new_s(k,j)*u(k,j))*trans(new_a.row(j))*new_a.row(j);
      mu_part=mu_part+trans((2*eta(k,j)*new_b(j)+u(k,j)-0.5)*(1-u(k,j)+new_s(k,j)*u(k,j))*new_a.row(j));
    }
    arma::mat Sigma= eye(domain,domain);
    arma::mat sigmahat=solve((Sigma+2*sigma_part),eye(domain,domain));
    arma::vec muhat=sigmahat*mu_part;
    SIGMA.slice(k)=sigmahat;
    MU.col(k)=muhat;
    arma::mat mukk=sigmahat+muhat*trans(muhat);
    arma::mat muk=new_a*mukk*trans(new_a);
    xi.row(k)=trans(sqrt(square(new_b)-2*new_b%(new_a*muhat)+muk.diag()));
  }
  return List::create(Named("SIGMA") = SIGMA,
                      Named("MU") = MU,Named("xi") = xi);
}


// [[Rcpp::export]]

List aefa3(const arma::mat&u,const arma::vec& id_1,const int& item,const int& domain, const arma::mat& eta,
           const arma::mat& a,const arma::vec& new_b,const arma::mat& new_s,const arma::cube& SIGMA, const arma::mat& MU,const arma::mat& prev_a_num,
           const arma::cube& prev_a_denom, const double& dec_st){
  //update a
  arma::mat new_a=a;
  arma::mat a_num=zeros(domain,item);
  arma::cube a_denom=zeros(domain,domain,item);
  for(int j=0; j<item; ++j){
    arma::mat a_denom_sub=a_denom.slice(j);
    for(int k=0; k<id_1.n_elem; ++k){
      int i=id_1(k)-1;
      arma::mat sigma=SIGMA.slice(i);
      arma::vec mu=MU.col(i);
      double a1=1-u(i,j)+new_s(i,j)*u(i,j);
      a_denom_sub=a_denom_sub+a1*eta(i,j)*(sigma+mu*trans(mu));
      a_num.col(j)=a_num.col(j)+a1*(u(i,j)-0.5+2*new_b(j)*eta(i,j))*mu;
    }
    if(rank(a_denom_sub) < domain){
      break;
    }
    a_denom.slice(j)=a_denom_sub;
    arma::mat prev_a_denom_sub=prev_a_denom.slice(j);
    arma::mat a2=solve((dec_st*a_denom_sub+(1-dec_st)*prev_a_denom_sub),eye(domain,domain));
    new_a.row(j)=trans(a2*(dec_st*a_num.col(j)+(1-dec_st)*prev_a_num.col(j))/2);
  }
  return List::create(Named("new_a") = new_a,Named("a_denom") = a_denom,Named("a_num") = a_num);
}
