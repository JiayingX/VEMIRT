#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List eefa2(const arma::mat&u, const int& domain, const int& person,
           const int& item, const arma::mat& eta, const arma::mat& new_a,const arma::vec& new_b){
  //estep
  arma::mat MU=mat(domain,person,fill::zeros);
  arma::cube SIGMA=zeros(domain,domain,person);
  arma::mat xi=mat(person,item,fill::zeros);
  for(int i=0; i<person; ++i){
    arma::mat sigma_part=arma::mat(domain,domain,arma::fill::zeros);
    //mu_part.zeros();
    arma::vec mu_part= arma::zeros(domain);
    for(int j=0; j<item; ++j){
      if(NumericVector::is_na(eta(i,j))){
        continue;}
      sigma_part=sigma_part+eta(i,j)*trans(new_a.row(j))*new_a.row(j);
      mu_part=mu_part+trans((2*eta(i,j)*new_b(j)+u(i,j)-0.5)*new_a.row(j));
    }
    arma::mat Sigma= eye(domain,domain);
    arma::mat sigmahat=solve((solve(Sigma,eye(domain,domain))+2*sigma_part),eye(domain,domain));
    arma::vec muhat=sigmahat*mu_part;
    SIGMA.slice(i)=sigmahat;
    MU.col(i)=muhat;
    arma::mat apro=new_a*(sigmahat+muhat*trans(muhat))*trans(new_a);
    xi.row(i)=trans(sqrt(square(new_b)-2*new_b%(new_a*muhat)+apro.diag()));
  }
  return List::create(Named("SIGMA") = SIGMA,
                      Named("MU") = MU,Named("xi") = xi);
}

// [[Rcpp::export]]

arma::mat aefa2(const arma::mat&u,const int& domain,const int& person, const int& item, const arma::mat& eta,
                const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU){
  //update a
  arma::mat new_a=mat(item,domain,fill::zeros);
  for(int j=0; j<item; ++j){
    arma::vec a_nu= zeros(domain);
    arma::mat a_de= zeros(domain,domain);
    for(int i=0; i<person; ++i){
      arma::mat sigma=SIGMA.slice(i);
      arma::vec mu=MU.col(i);
      if(NumericVector::is_na(eta(i,j))){
        continue;}
      a_de=a_de+eta(i,j)*sigma+eta(i,j)*(mu*trans(mu));
      a_nu=a_nu+(u(i,j)-0.5+2*new_b(j)*eta(i,j))*mu;
    }
    new_a.row(j)=trans(solve(a_de,eye(domain,domain))*a_nu/2);
  }
  return new_a;
}
