#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

List ecfa2(const arma::mat&u,const arma::mat& Sigma, const int& domain, const int& person,
           const int& item, const arma::mat& eta, const arma::mat& new_a,const arma::vec& new_b){
  //e step
  arma::mat MU=mat(domain,person,fill::zeros);
  arma::cube SIGMA=zeros(domain,domain,person);
  arma::mat Spart=mat(domain,domain,fill::zeros);
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
    arma::mat sigmahat=solve((solve(Sigma,eye(domain,domain))+2*sigma_part),eye(domain,domain));
    arma::vec muhat=sigmahat*mu_part;
    SIGMA.slice(i)=sigmahat;
    MU.col(i)=muhat;
    arma::mat apro=new_a*(sigmahat+muhat*trans(muhat))*trans(new_a);
    xi.row(i)=trans(sqrt(square(new_b)-2*new_b%(new_a*muhat)+apro.diag()));
    Spart=Spart+sigmahat+muhat*trans(muhat);
  }
  return List::create(Named("Spart") = Spart,Named("SIGMA") = SIGMA,
                      Named("MU") = MU,Named("xi") = xi);
}


// [[Rcpp::export]]

arma::mat acfa2(const arma::mat&u,const arma::mat&indic,const int& person,
                const int& item, const int& domain, const arma::mat& eta,
                const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,const arma::mat& new_a1){
  //update alpha parameters
  arma::mat new_a=new_a1;
  for(int j=0; j<item; ++j){
    arma::rowvec Ind=indic.row(j);
    arma::uvec iind=find(Ind==1);
    int di=std::count(Ind.begin(),Ind.end(), 1);
    arma::vec a_nu= zeros(di);
    arma::mat a_de= zeros(di,di);
    for(int i=0; i<person; ++i){
      arma::mat sigma=SIGMA.slice(i);
      sigma=sigma.submat(iind,iind);
      arma::vec mu=MU.col(i);
      mu=mu.elem(iind);
      if(NumericVector::is_na(eta(i,j))){
        continue;}
      a_de=a_de+eta(i,j)*sigma+eta(i,j)*(mu*trans(mu));
      a_nu=a_nu+(u(i,j)-0.5+2*new_b(j)*eta(i,j))*mu;
    }
    arma::uvec id(1);
    id.at(0)=j;
    //new_a.submat(id,iind)=trans(inv(a_de)*a_nu/2);
    new_a.submat(id,iind)=trans(solve(a_de,eye(di,di))*a_nu/2);
  }
  return new_a;
}


// [[Rcpp::export]]

double lb2pl(const arma::mat&u,const arma::mat& xi,
             const arma::mat& sig,const arma::mat& new_a,
             const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU){
  //lower bound for 2pl
  int person=u.n_rows;
  int item=u.n_cols;
  int domain=sig.n_cols;
  arma::mat eta =(exp(xi)/(1+exp(xi))-0.5)/(2*xi);
  double sum = 0;
  double sum2 = 0;
  for(int i=0; i<person; ++i){
    for(int j=0; j<item; ++j){
      if(NumericVector::is_na(eta(i,j))){
        continue;}
      sum = sum + log(exp(xi(i,j))/(1+exp(xi(i,j))))+(0.5-u(i,j))*new_b(j)+(u(i,j)-0.5)*as_scalar(new_a.row(j)*MU.col(i));
      sum = sum - 0.5*xi(i,j);
      double apro=as_scalar(new_a.row(j)*(SIGMA.slice(i)+MU.col(i)*trans(MU.col(i)))*trans(new_a.row(j)));
      sum = sum - eta(i,j)*(new_b(j)*new_b(j)-2*new_b(j)*as_scalar(new_a.row(j)*MU.col(i))+apro - xi(i,j)*xi(i,j));
    }
    sum2=sum2+trace(solve(sig,eye(domain,domain))*(SIGMA.slice(i)+MU.col(i)*trans(MU.col(i))));
  }
  double lb=sum - 0.5*sum2 + person/2*log(det(solve(sig,eye(domain,domain))));
  return lb;
}
