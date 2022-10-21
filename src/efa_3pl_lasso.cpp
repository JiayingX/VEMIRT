#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

List nalc13pl(const arma::mat&u,const arma::mat&indic,const arma::vec& nopenalty_col,const arma::mat& eta,const arma::mat& s,
              const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,
              const arma::mat& a,const arma::vec& id_1,const arma::mat& prev_a_num,const arma::cube& prev_a_denom,
              const double& dec_st){
  //lasso l1:update a without penalty
  int domain=indic.n_cols;
  int item=u.n_cols;
  arma::mat new_a=a;
  arma::mat a_num=mat(domain,item,fill::zeros);
  arma::cube a_denom=zeros(domain,domain,item);
  for(int j=0; j<nopenalty_col.n_elem; ++j){
    int k=nopenalty_col(j)-1;
    arma::rowvec Ind=indic.row(k);
    arma::uvec iind=find(Ind==1);
    arma::mat a_denom_sub=a_denom.slice(j);
    arma::uvec id2(1);
    id2(0)=k;
    for(int i=0; i<id_1.n_elem; ++i){
      int l=id_1(i)-1;
      double a1=1-u(l,k)+s(l,j)*u(l,k);
      arma::mat sigma=SIGMA.slice(l);
      arma::uvec id(1);
      id(0)=l;
      arma::mat a2=sigma.submat(iind,iind)+MU.submat(iind,id)*trans(MU.submat(iind,id));
      a_denom_sub.submat(iind,iind)=a_denom_sub.submat(iind,iind)+a1*eta(l,k)*a2;
      a_num.submat(iind,id2)=a_num.submat(iind,id2)+a1*(u(l,k)-0.5+2*new_b(k)*eta(l,k))*MU(iind,id);
    }
    arma::mat prev_a_denom_sub=prev_a_denom.slice(k);
    arma::mat a3=dec_st*a_denom_sub.submat(iind,iind) + (1-dec_st)*prev_a_denom_sub.submat(iind,iind);
    arma::mat a4=(dec_st*a_num.submat(iind,id2)+(1-dec_st)*prev_a_num.submat(iind,id2))/2;
    new_a.submat(id2,iind)=trans(solve(a3,eye(iind.n_elem,iind.n_elem))*a4);
    a_denom.slice(k)=a_denom_sub;
  }
  return List::create(Named("new_a") = new_a,Named("a_denom") = a_denom,Named("a_num") = a_num);
}


// [[Rcpp::export]]

List palc13pl(const arma::mat&u,const int& domain,const int& item,const arma::vec& id_1,const double& lbd,
              const arma::mat& eta,const arma::mat& s,const arma::mat& a,const arma::vec& new_b,const arma::cube& SIGMA,
              const arma::mat& MU,const double& dec_st,const arma::mat& prev_delta,const arma::mat& prev_deriv2,const arma::vec& sdf){
  //lasso l1:update a with penalty
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
      if(abs(avg_delta)>lbd){
        double S=sign(avg_delta)*(abs(avg_delta) - lbd);
        new_a(j,k) = S/(dec_st*deriv2(j,k) + (1-dec_st)*prev_deriv2(j,k));
      }else{
        new_a(j,k) = 0;
      }
    }
  }
  return List::create(Named("new_a") = new_a,Named("delta") = delta,Named("deriv2") = deriv2);
}


// [[Rcpp::export]]

List nalc23pl(const arma::mat&u,const int& domain,const int& item,
              const arma::mat& a,const arma::vec& nopenalty_col,
              const arma::vec& lastone,const arma::mat& s, const arma::mat& eta,const arma::vec& new_b,
              const arma::cube& SIGMA, const arma::mat& MU,const arma::vec& id_1,const arma::mat& prev_a_num,const arma::cube& prev_a_denom,
              const double& dec_st){
  //laso l2: non penalty a
  arma::mat new_a=a;
  arma::mat a_num=mat(domain,item,fill::zeros);
  arma::cube a_denom=zeros(domain,domain,item);
  for(int j=0; j<nopenalty_col.n_elem; ++j){
    int k=nopenalty_col(j)-1;
    int l=lastone(j)-1;
    for(int i=0; i<id_1.n_elem; ++i){
      int n=id_1(i)-1;
      double sigma=SIGMA(l,l,n);
      double mu=MU(l,n);
      double a1=1-u(n,k)+s(n,k)*u(n,k);
      a_denom(l,l,k)=a_denom(l,l,k)+a1*eta(n,k)*sigma+eta(n,k)*(mu*mu);
      a_num(l,k)=a_num(l,k)+a1*(u(n,k)-0.5+2*new_b(k)*eta(n,k))*mu;
    }
    new_a(k,l)=1/(dec_st*a_denom(l,l,k)+(1-dec_st)*prev_a_denom(l,l,k))*(dec_st*a_num(l,k)+(1-dec_st)*prev_a_num(l,k))/2;

  }
  return List::create(Named("new_a") = new_a,Named("a_denom") = a_denom,Named("a_num") = a_num);
}


// [[Rcpp::export]]

List palc23pl(const arma::mat&u,const arma::mat& a,
              const arma::vec& nopenalty_col,const arma::mat& s,
              const arma::vec& lastone,const double& lbd,
              const arma::vec& id_1, const arma::mat& eta,const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,
              const double& dec_st,const arma::mat& prev_delta,const arma::mat& prev_deriv2){
  //lasso l2: update a with penalty: 1:domain, off-diagonal
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
      if(abs(avg_delta)>lbd){
        double S=sign(avg_delta)*(abs(avg_delta) - lbd);
        new_a1(n,m) = S/(dec_st*deriv2(n,m) + (1-dec_st)*prev_deriv2(n,m));
      }else{
        new_a1(n,m) = 0;
      }
    }
  }
  return List::create(Named("new_a") = new_a1,Named("delta") = delta,Named("deriv2") = deriv2);
}


// [[Rcpp::export]]

List palc23pl1(const arma::mat&u,const arma::mat& a,const arma::mat& s,const double& lbd,
               const arma::vec& id_1, const arma::mat& eta,const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,
               const double& dec_st,const arma::mat& prev_delta,const arma::mat& prev_deriv2,const arma::vec& pc,
               const arma::mat& delta1, const arma::mat& deriv22){
  //lasso l2: update a with penalty: domain+1:item
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
      if(abs(avg_delta)>lbd){
        double S=sign(avg_delta)*(abs(avg_delta) - lbd);
        new_a1(n,m) = S/(dec_st*deriv2(n,m) + (1-dec_st)*prev_deriv2(n,m));
      }else{
        new_a1(n,m) = 0;
      }
    }
  }
  return List::create(Named("new_a") = new_a1,Named("delta") = delta,Named("deriv2") = deriv2);
}
