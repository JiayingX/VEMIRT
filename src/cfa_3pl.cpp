#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

List ecfa3pl(const arma::mat&u,const arma::mat& Sigma, const int& domain, const arma::vec& id_1,
             const int& item, const arma::mat& eta, const arma::mat& new_s,const arma::mat& new_a,const arma::vec& new_b,
             const arma::cube& SIGMA1,const arma::mat& MU1,const arma::mat& xi1){
  //estep
  arma::mat MU=MU1;
  arma::cube SIGMA=SIGMA1;
  arma::mat Spart=mat(domain,domain,fill::zeros);
  arma::mat xi=xi1;
  for(int i=0; i<id_1.n_elem; ++i){
    int k=id_1(i)-1;
    arma::mat sigma_part=arma::mat(domain,domain,arma::fill::zeros);
    arma::vec mu_part= arma::zeros(domain);
    for(int j=0; j<item; ++j){
      sigma_part=sigma_part+eta(k,j)*(1-u(k,j)+new_s(k,j)*u(k,j))*trans(new_a.row(j))*new_a.row(j);
      mu_part=mu_part+trans((2*eta(k,j)*new_b(j)+u(k,j)-0.5)*(1-u(k,j)+new_s(k,j)*u(k,j))*new_a.row(j));
    }
    arma::mat sigmahat=solve((solve(Sigma,eye(domain,domain))+2*sigma_part),eye(domain,domain));
    arma::vec muhat=sigmahat*mu_part;
    SIGMA.slice(k)=sigmahat;
    MU.col(k)=muhat;
    arma::mat mukk=sigmahat+muhat*trans(muhat);
    arma::mat muk=new_a*mukk*trans(new_a);
    xi.row(k)=trans(sqrt(square(new_b)-2*new_b%(new_a*muhat)+muk.diag()));
    Spart=Spart+sigmahat+muhat*trans(muhat);
  }
  return List::create(Named("SIGMA") = SIGMA,Named("MU") = MU,Named("xi") = xi,Named("Spart") = Spart);
}


// [[Rcpp::export]]

arma::mat s3pl(const int& item, const arma::vec& id_1,const arma::mat& s,const arma::mat& eta,
               const arma::mat& new_a,const arma::vec& s_num,const arma::vec& new_b,const arma::mat&u,const arma::cube& SIGMA,
               const arma::mat& MU,const arma::mat& xi){
  //update s
  arma::mat new_s=s;
  for(int i=0; i<id_1.n_elem; ++i){
    int k=id_1(i)-1;
    for(int j=0; j<item; ++j){
      double s_denom = 0;
      if (u(k,j) == 0){
        new_s(k,j) = 1;
      }else{
        s_denom = s_denom + log(exp(xi(k,j))/(1+exp(xi(k,j))))+(0.5-u(k,j))*new_b(j)+(u(k,j)-0.5)*as_scalar(new_a.row(j)*MU.col(k))-0.5*xi(k,j);
        arma::mat mid=SIGMA.slice(k)+MU.col(k)*trans(MU.col(k));
        double aa=as_scalar(new_a.row(j)*mid*trans(new_a.row(j)));
        s_denom = s_denom - eta(k,j)*as_scalar(new_b(j)*new_b(j)-2*new_b(j)*new_a.row(j)*MU.col(k)+aa-xi(k,j)*xi(k,j));
        s_denom = exp(-s_denom);
        new_s(k,j) = s_num(j)/(s_num(j)+s_denom);
      }
    }
  }
  return new_s;
}

// [[Rcpp::export]]

arma::mat eta3pl(const int& item, const arma::vec& id_1,const arma::mat& ea,const arma::mat& xi){
  //update eta
  arma::mat eta=ea;
  for(int i=0; i<id_1.n_elem; ++i){
    int k=id_1(i)-1;
    for(int j=0; j<item; ++j){
      if(xi(k,j)<0.01){
        eta(k,j)=0.125;
      }
    }
  }
  return eta;
}


// [[Rcpp::export]]

List b3pl(const int& item, const arma::vec& id_1,const arma::mat& new_s,const arma::mat& eta,
          const arma::mat& u, const arma::mat& new_a,const double& dec_st,const arma::mat& MU,
          const arma::vec& prev_b_num,const arma::vec& prev_b_denom,const double& mu_b,const double& sigma2_b){
  arma::vec b_denom= arma::zeros(item);
  arma::vec b_num= arma::zeros(item);
  //update b
  for(int j=0; j<item; ++j){
    for(int i=0; i<id_1.n_elem; ++i){
      int k=id_1(i)-1;
      b_denom(j) = b_denom(j) + 2*(1-u(k,j)+new_s(k,j)*u(k,j))*eta(k,j);
      b_num(j) = b_num(j) + (1-u(k,j)+new_s(k,j)*u(k,j))*(0.5-u(k,j)+2*eta(k,j)*as_scalar(new_a.row(j)*MU.col(k)));
    }
  }
  b_denom = b_denom + 1/sigma2_b;
  b_num = b_num + mu_b/sigma2_b;
  arma::vec new_b=(dec_st*b_num+(1-dec_st)*prev_b_num)/(dec_st*b_denom+(1-dec_st)*prev_b_denom);
  return List::create(Named("new_b") = new_b,Named("b_denom") = b_denom,Named("b_num") = b_num);
}


// [[Rcpp::export]]

List acfa3(const arma::mat&u,const arma::mat&indic,const arma::vec& id_1,const int& item,const int& domain, const arma::mat& eta,
           const arma::mat& a,const arma::vec& new_b,const arma::mat& new_s,const arma::cube& SIGMA, const arma::mat& MU,const arma::mat& prev_a_num,
           const arma::cube& prev_a_denom, const double& dec_st){
  //update a
  arma::mat new_a=a;
  arma::mat a_num=zeros(domain,item);
  arma::cube a_denom=zeros(domain,domain,item);
  for(int j=0; j<item; ++j){
    arma::rowvec Ind=indic.row(j);
    arma::uvec iind=find(Ind==1);
    arma::mat a_denom_sub=a_denom.slice(j);
    arma::vec a_num_sub=a_num.col(j);
    arma::vec a_num_sub1=a_num_sub.elem(iind);
    for(int i=0; i<id_1.n_elem; ++i){
      int k=id_1(i)-1;
      arma::mat sigma=SIGMA.slice(k);
      sigma=sigma.submat(iind,iind);
      arma::vec mu=MU.col(k);
      mu=mu.elem(iind);
      double a1=1-u(k,j)+new_s(k,j)*u(k,j);
      a_denom_sub.submat(iind,iind)=a_denom_sub.submat(iind,iind)+a1*eta(k,j)*(sigma+mu*trans(mu));
      a_num_sub1=a_num_sub1+a1*(u(k,j)-0.5+2*new_b(j)*eta(k,j))*mu;
    }
    a_denom.slice(j)=a_denom_sub;
    a_num_sub.elem(iind)=a_num_sub1;
    a_num.col(j)=a_num_sub;
    arma::mat prev_a_denom_sub=prev_a_denom.slice(j);
    prev_a_denom_sub=prev_a_denom_sub.submat(iind,iind);
    arma::mat a2=solve((dec_st*a_denom_sub.submat(iind,iind)+(1-dec_st)*prev_a_denom_sub),eye(iind.n_elem,iind.n_elem));
    arma::vec prev_a_num_sub=prev_a_num.col(j);
    prev_a_num_sub=prev_a_num_sub.elem(iind);
    arma::uvec id(1);
    id.at(0)=j;
    new_a.submat(id,iind)=trans(a2*(dec_st*a_num_sub1+(1-dec_st)*prev_a_num_sub)/2);
  }
  return List::create(Named("new_a") = new_a,Named("a_denom") = a_denom,Named("a_num") = a_num);
}

// [[Rcpp::export]]

double lb3pl(const arma::mat&u,const arma::mat& xi, const arma::mat& s, const arma::vec& id_1,
             const arma::mat& new_a,const arma::vec& new_c,const arma::mat& sig,
             const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,const double& Alpha,
             const double& Beta,const double& mu_b,const double& sigma2_b){
  //compute the lower bound
  int item=u.n_cols;
  int domain=new_a.n_cols;
  arma::mat eta =(exp(xi)/(1+exp(xi))-0.5)/(2*xi);
  double sum1 = 0;
  double sum2 = 0;
  double sum3 = 0;
  for(int k=0; k<id_1.n_elem; ++k){
    int i=id_1(k)-1;
    for(int j=0; j<item; ++j){
      double a1=1-u(i,j)+s(i,j)*u(i,j);
      double a2=log(exp(xi(i,j))/(1+exp(xi(i,j))))+(0.5-u(i,j))*new_b(j);
      double a3=(u(i,j)-0.5)*as_scalar(new_a.row(j)*MU.col(i));
      double apro=as_scalar(new_a.row(j)*(SIGMA.slice(i)+MU.col(i)*trans(MU.col(i)))*trans(new_a.row(j)));
      double a4=eta(i,j)*(new_b(j)*new_b(j)-2*new_b(j)*as_scalar(new_a.row(j)*MU.col(i))+apro - xi(i,j)*xi(i,j));
      sum1 = sum1 + a1*(a2+a3-0.5*xi(i,j)-a4);
      sum3=sum3 + a1*log(1-new_c(j)) + u(i,j)*(1-s(i,j))*log(new_c(j));
    }
    sum2=sum2+trace(solve(sig,eye(domain,domain))*(SIGMA.slice(i)+MU.col(i)*trans(MU.col(i))));
  }
  double sum_b = as_scalar(-1/(2*sigma2_b)*trans(new_b-mu_b)*(new_b-mu_b));
  double sum_c = sum((Alpha-1)*log(new_c)+(Beta-1)*log(1-new_c));
  double lb=sum1 - 0.5*sum2 + id_1.n_elem/2*log(det(eye(domain,domain)))+sum3+ sum_b + sum_c;
  return lb;
}
