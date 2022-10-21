#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

arma::mat paal2pl(const arma::mat&u,const int& domain,const int& person,const double& lbd, const arma::vec& sdf, const arma::mat& eta,
                  const arma::mat& a,
                  const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,const arma::mat& weights){
  //adaptive lasso l1: update penalty a
  arma::mat new_a=a;
  for(int j=0; j<sdf.n_elem; ++j){
    int m=sdf(j)-1;
    for(int k=0; k<domain; ++k){
      double delta = 0;
      double deriv2 = 0;
      for(int i=0; i<person; ++i){
        arma::mat sigma=SIGMA.slice(i);
        arma::vec mu=MU.col(i);
        arma::mat sigmumu=sigma+mu*trans(mu);
        arma::uvec A=regspace<uvec>(0, 1,k-1);
        arma::uvec B=regspace<uvec>(k+1, 1, domain-1);
        arma::uvec ex_k=join_vert(A,B);
        arma::uvec id(1);
        id.at(0)=m;
        arma::uvec id2(1);
        id2.at(0)=k;
        delta = delta + (u(i,m)-0.5)*mu(k)+2*new_b(m)*eta(i,m)*mu(k)-2*eta(i,m)*dot(trans(new_a.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
        deriv2 = deriv2 + 2*eta(i,m)*sigmumu(k,k);
      }
      if(abs(delta)>lbd*(1/weights(m,k))){
        double S=sign(delta)*(abs(delta) - lbd*(1/weights(m,k)));
        new_a(m,k) = S/deriv2;
      }else{
        new_a(m,k) = 0;
      }
    }
  }
  return new_a;
}


// [[Rcpp::export]]

arma::mat paalc22pl(const arma::mat&u,const arma::mat& a,
                    const arma::mat&indic,const arma::vec& nopenalty_col,
                    const arma::vec& lastone,const double& lbd,
                    const int& person, const arma::mat& eta,const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,
                    const arma::mat&weights){
  //adaptive lasso l2:update a with penalty: 1:domain, off-diagonal
  arma::mat new_a1=a;
  int domain=new_a1.n_cols;
  for(int j=0; j<nopenalty_col.n_elem; ++j){
    int n=nopenalty_col(j)-1;
    int l=lastone(j)-1;
    arma::uvec A=regspace<uvec>(0, 1,l-1);
    arma::uvec B=regspace<uvec>(l+1, 1, domain-1);
    arma::uvec sdf=join_vert(A,B);
    for(int k=0; k<sdf.n_elem; ++k){
      double delta = 0;
      double deriv2 = 0;
      int m=sdf(k);
      for(int i=0; i<person; ++i){
        arma::mat sigma=SIGMA.slice(i);
        arma::vec mu=MU.col(i);
        arma::mat sigmumu=sigma+mu*trans(mu);
        arma::uvec A1=regspace<uvec>(0, 1,m-1);
        arma::uvec B1=regspace<uvec>(m+1, 1, domain-1);
        arma::uvec ex_k=join_vert(A1,B1);
        arma::uvec id(1);
        id.at(0)=n;
        arma::uvec id2(1);
        id2.at(0)=m;
        delta = delta + (u(i,n)-0.5)*mu(m)+2*new_b(n)*eta(i,n)*mu(m)-2*eta(i,n)*dot(trans(new_a1.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
        deriv2 = deriv2 + 2*eta(i,n)*sigmumu(m,m);
      }
      if(abs(delta)>lbd*(1/weights(n,m))){
        double S=sign(delta)*(abs(delta) - lbd*(1/weights(n,m)));
        new_a1(n,m) = S/deriv2;
      }else{
        new_a1(n,m) = 0;
      }
    }
  }
  return new_a1;
}


// [[Rcpp::export]]

arma::mat paalc22pl1(const arma::mat&u,const int& domain,const int& item,const int& person,const double& lbd, const arma::mat& eta, arma::mat& new_a,
                     const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,const arma::vec& pc,const arma::mat&weights){
  //adaptive lasso l2: update a with penalty: domain+1:item
  for(int n=0; n<pc.n_elem; ++n){
    int j=pc(n)-1;
    for(int k=0; k<domain; ++k){
      double delta = 0;
      double deriv2 = 0;
      for(int i=0; i<person; ++i){
        arma::mat sigma=SIGMA.slice(i);
        arma::vec mu=MU.col(i);
        arma::mat sigmumu=sigma+mu*trans(mu);
        arma::uvec A=regspace<uvec>(0, 1,k-1);
        arma::uvec B=regspace<uvec>(k+1, 1, domain-1);
        arma::uvec ex_k=join_vert(A,B);
        arma::uvec id(1);
        id.at(0)=j;
        arma::uvec id2(1);
        id2.at(0)=k;
        delta = delta + (u(i,j)-0.5)*mu(k)+2*new_b(j)*eta(i,j)*mu(k)-2*eta(i,j)*dot(trans(new_a.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
        deriv2 = deriv2 + 2*eta(i,j)*sigmumu(k,k);
      }
      if(abs(delta)>lbd*(1/weights(j,k))){
        double S=sign(delta)*(abs(delta) - lbd*(1/weights(j,k)));
        new_a(j,k) = S/deriv2;
      }else{
        new_a(j,k) = 0;
      }
    }
  }
  return new_a;
}
