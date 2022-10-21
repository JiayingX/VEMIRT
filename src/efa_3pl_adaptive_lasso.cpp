#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

List paal3pl(const arma::mat&u,const int& domain,const int& item,const arma::vec& id_1,const double& lbd,
             const arma::mat& eta,const arma::mat& s,const arma::mat& a,const arma::vec& new_b,const arma::cube& SIGMA,
             const arma::mat& MU,const double& dec_st,const arma::mat& prev_delta,const arma::mat& prev_deriv2,
             const arma::mat& weights,const arma::vec& sdf){
  //adaptive l1:update penalty a
  arma::mat new_a=a;
  arma::mat delta=mat(item,domain,fill::zeros);
  arma::mat deriv2=mat(item,domain,fill::zeros);
  for(int n=0; n<sdf.n_elem; ++n){
    int j=sdf(n)-1;
    for(int k=0; k<domain; ++k){
      for(int i=0; i<id_1.n_elem; ++i){
        int l=id_1(i)-1;
        arma::mat sigma=SIGMA.slice(l);
        arma::vec mu=MU.col(l);
        arma::mat sigmumu=sigma+mu*trans(mu);
        arma::uvec A=regspace<uvec>(0, 1,k-1);
        arma::uvec B=regspace<uvec>(k+1, 1, domain-1);
        arma::uvec ex_k=join_vert(A,B);
        double a1=1-u(l,j)+s(l,j)*u(l,j);
        arma::uvec id(1);
        id.at(0)=j;
        arma::uvec id2(1);
        id2.at(0)=k;
        delta(j,k) = delta(j,k) + a1*((u(l,j)-0.5)*mu(k)+2*new_b(j)*eta(l,j)*mu(k)-2*eta(l,j)*dot(trans(new_a.submat(id,ex_k)),sigmumu.submat(ex_k,id2)));
        deriv2(j,k) = deriv2(j,k) + a1*2*eta(l,j)*sigmumu(k,k);
      }
      double avg_delta = dec_st*delta(j,k) + (1-dec_st)*prev_delta(j,k);
      if(abs(avg_delta)>lbd*(1/weights(j,k))){
        double S=sign(avg_delta)*(abs(avg_delta) - lbd*(1/weights(j,k)));
        new_a(j,k) = S/(dec_st*deriv2(j,k) + (1-dec_st)*prev_deriv2(j,k));
      }else{
        new_a(j,k) = 0;
      }
    }
  }
  return List::create(Named("new_a") = new_a,Named("delta") = delta,Named("deriv2") = deriv2);
}


// [[Rcpp::export]]

List paalc23pl(const arma::mat&u,const arma::mat& a,
               const arma::vec& nopenalty_col,const arma::mat& s,
               const arma::vec& lastone,const double& lbd,
               const arma::vec& id_1, const arma::mat& eta,const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,
               const double& dec_st,const arma::mat& prev_delta,const arma::mat& prev_deriv2,const arma::mat&weights){
  //adapt l2: update a with penalty: 1:domain, off-diagonal
  arma::mat new_a1=a;
  int domain=new_a1.n_cols;
  int item=u.n_cols;
  arma::mat delta=mat(item,domain,fill::zeros);
  arma::mat deriv2=mat(item,domain,fill::zeros);
  for(int j=0; j<nopenalty_col.n_elem; ++j){
    int n=nopenalty_col(j)-1;
    int l=lastone(j)-1;
    arma::uvec A=regspace<uvec>(0, 1,l-1);
    arma::uvec B=regspace<uvec>(l+1, 1, domain-1);
    arma::uvec sdf=join_vert(A,B);
    for(int k=0; k<sdf.n_elem; ++k){
      int m=sdf(k);
      for(int p=0; p<id_1.n_elem; ++p){
        int i=id_1(p)-1;
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
        double a2=1-u(i,n)+s(i,n)*u(i,n);
        delta(n,m) = delta(n,m) + a2*(u(i,n)-0.5)*mu(m)+2*new_b(n)*eta(i,n)*mu(m)-2*eta(i,n)*dot(trans(new_a1.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
        deriv2(n,m) = deriv2(n,m) + 2*a2*eta(i,n)*sigmumu(m,m);
      }
      double avg_delta = dec_st*delta(n,m) + (1-dec_st)*prev_delta(n,m);
      if(avg_delta>lbd*(1/weights(n,m))){
        double S=sign(avg_delta)*(abs(avg_delta) - lbd*(1/weights(n,m)));
        new_a1(n,m) = S/(dec_st*deriv2(n,m) + (1-dec_st)*prev_deriv2(n,m));
      }else{
        new_a1(n,m) = 0;
      }
    }
  }
  return List::create(Named("new_a") = new_a1,Named("delta") = delta,Named("deriv2") = deriv2);
}


// [[Rcpp::export]]

List paalc23pl1(const arma::mat&u,const arma::mat& a,const arma::mat& s,const double& lbd,
                const arma::vec& id_1, const arma::mat& eta,const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,
                const double& dec_st,const arma::mat& prev_delta,const arma::mat& prev_deriv2,const arma::vec& pc,
                const arma::mat& delta1, const arma::mat& deriv22,const arma::mat&weights){
  //adapt l2: update a with penalty: domain+1:item
  arma::mat new_a1=a;
  int domain=new_a1.n_cols;
  arma::mat delta=delta1;
  arma::mat deriv2=deriv22;
  for(int j=0; j<pc.n_elem; ++j){
    int n=pc(j)-1;
    for(int m=0; m<domain; ++m){
      for(int p=0; p<id_1.n_elem; ++p){
        int i=id_1(p)-1;
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
        double a2=1-u(i,n)+s(i,n)*u(i,n);
        delta(n,m) = delta(n,m) + a2*(u(i,n)-0.5)*mu(m)+2*new_b(n)*eta(i,n)*mu(m)-2*eta(i,n)*dot(trans(new_a1.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
        deriv2(n,m) = deriv2(n,m) + 2*a2*eta(i,n)*sigmumu(m,m);
      }
      double avg_delta = dec_st*delta(n,m) + (1-dec_st)*prev_delta(n,m);
      if(avg_delta>lbd*(1/weights(n,m))){
        double S=sign(avg_delta)*(abs(avg_delta) - lbd*(1/weights(n,m)));
        new_a1(n,m) = S/(dec_st*deriv2(n,m) + (1-dec_st)*prev_deriv2(n,m));
      }else{
        new_a1(n,m) = 0;
      }
    }
  }
  return List::create(Named("new_a") = new_a1,Named("delta") = delta,Named("deriv2") = deriv2);
}

