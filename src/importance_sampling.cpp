#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

arma::cube sampling(const arma::mat& mu_i, const arma::cube& sig_i, const int& person,
                    const int& domain, const int& S, const int& M) {
  arma::cube theta_IS(S * M, domain,person);
  for (int i = 0; i < person; i++){
    arma::vec mu_i_slice = mu_i.col(i);
    arma::mat sig_i_slice = sig_i.slice(i);
    arma::mat theta_ism = arma::randn(S*M, domain);
    arma::mat theta_i = arma::repmat(mu_i_slice, 1, S*M).t() + theta_ism * arma::chol(sig_i_slice);
    theta_IS.slice(i) = theta_i;
  }
  return theta_IS;
}


// [[Rcpp::export]]

List importance_weights(const arma::cube theta_IS, const arma::mat& mu_i, const arma::cube sig_i,
                        const arma::mat& a, const arma::vec& b, const arma::mat& Sigma_theta, const arma::mat& y,
                        const int& person, const int& item, const int& S, const int& M) {
  // Calculate item probability [s*m, item, person]
  arma::mat part1=zeros(S*M,person);
  arma::mat part2=zeros(S*M,person);
  arma::mat part3=zeros(S*M,person);
  arma::mat q=zeros(S*M,person);
  // Calculate part3
  double log_det_sigma_theta = log(det(Sigma_theta));
  arma::mat inv_sigma_theta = inv(Sigma_theta);
  arma::mat theta_mul = arma::mat(S*M,person,arma::fill::zeros);
  for (int i = 0; i < person; i++) {
    arma::mat p_item_part=arma::mat(S*M,item,arma::fill::zeros);
    arma::mat b_rep = arma::mat(S*M,item,arma::fill::zeros);
    b_rep.each_row() = trans(b);
    // Calculate q
    arma::mat sig_subi = sig_i.slice(i);
    double log_det_sigma_i = log(det(sig_subi));
    arma::mat inv_sigma_i = inv(sig_subi);
    // Calculate theta_part
    arma::mat theta_part = theta_IS.slice(i);
    for (int s = 0; s < S*M; s++) {
      part3(s,i) = -(as_scalar(theta_part.row(s) * inv_sigma_theta * trans(theta_part.row(s))))/2 - log_det_sigma_theta/2;
      q(s,i) = -(as_scalar((theta_part.row(s)-trans(mu_i.col(i))) * inv_sigma_i * trans(theta_part.row(s)-trans(mu_i.col(i)))))/2 - log_det_sigma_i/2;
    }
    p_item_part=1/(1+exp(-theta_part * trans(a) + b_rep));
    part1.col(i) = log(p_item_part) * trans(y.row(i));
    part2.col(i) = log(1-p_item_part) * trans(1-y.row(i));
  }
  arma::mat p = part1 + part2 + part3;
  arma::mat w = exp(p - q);
  arma::mat w_tildesub = zeros(S, person);
  for (int i = 0; i < person; i++) {
    for (int j = 0; j < S; j++) {
      int k = j * M;
      arma::mat sm = w(arma::span(k, k + M-1), i);
      w_tildesub(j, i) = arma::accu(sm);
    }
  }
  arma::mat w_tilde = w / repelem(w_tildesub, M, 1);
  return List::create(Named("w") = w,Named("w_tilde") = w_tilde,Named("w_tildesub") = w_tildesub);
}


// [[Rcpp::export]]

List importance_gradient(const arma::mat& y, const arma::mat& a, const arma::vec& b,
                         const arma::mat&  Sigma_theta, const arma::cube theta_IS,
                         const arma::mat& w_tilde,const int& S, const int& M) {
  int person = y.n_rows;
  int item = y.n_cols;
  int domain = a.n_cols;
  arma::vec grad_b= zeros(item);
  arma::mat grad_a = arma::mat(item, domain, arma::fill::zeros);
  arma::mat grad_Sigma_theta = arma::mat(domain, domain, arma::fill::zeros);
  arma::mat b_rep = arma::mat(S*M,item,arma::fill::zeros);
  b_rep.each_row() = trans(b);
  // Calculate
  for (int i = 0; i < person; i++) {
    arma::mat theta_part = theta_IS.slice(i);
    arma::mat temp = 1/(1+exp(theta_part * trans(a) - b_rep));
    for (int s = 0; s < S; s++) {
      for (int m = 0; m < M; m++){
        int k = s * M + m;
        grad_a += w_tilde(k,i)*(y.row(i) - 1 + temp.row(k)).t()*theta_part.row(k);
        grad_b += w_tilde(k,i)*(1 - y.row(i) - temp.row(k)).t();
        grad_Sigma_theta += w_tilde(k,i)*(Sigma_theta - (theta_part.row(k)).t() * theta_part.row(k))/2;
      }
    }
  }
  grad_a = grad_a / (S*person);
  grad_b = grad_b / (S*person);
  grad_Sigma_theta = grad_Sigma_theta / (S*person);
  return List::create(Named("grad_a") = grad_a,Named("grad_b") = grad_b,Named("grad_Sigma_theta") = grad_Sigma_theta);
}
