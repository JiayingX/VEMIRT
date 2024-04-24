#' Importance Weighted Version of GVEM Analysis for M2PL Models
#' @description An importance weighted version of GVEM (i.e., IW-GVEM) can be implemented
#' to correct the bias on item parameters under M2PL models
#'
#' @param u a \eqn{N \times J} \code{matrix} or a \code{data.frame} that
#' consists of binary responses of \eqn{N} individuals to \eqn{J} items. The
#' missing values are coded as \code{NA}
#' @param gvem_result a list that includes exploratory or confirmatory GVEM
#' results for M2PL models.
#' @param S the number of times to draw samples;default is 10
#' @param M the number of samples drawn from the variational distributions;default is 10
#' @param max.iter the maximum number of iterations for the EM cycle; default is 10
#'
#' @return a list containing the following objects:
#'   \item{ra}{item discrimination parameters estimated by GVEM, a \eqn{J \times K} \code{matrix}}
#'   \item{rb}{item difficulty parameters estimated by GVEM, vector of length \eqn{J}}
#'   \item{reta}{variational parameters \eqn{\eta(\xi)}, a \eqn{N \times J} matrix}
#'   \item{reps}{variational parameters \eqn{\xi}, a \eqn{N \times J} matrix}
#'   \item{rsigma}{population variance-covariance matrix estimated by GVEM, a \eqn{K \times K} matrix}
#'   \item{mu_i}{mean parameter for each person, a \eqn{K \times N} matrix}
#'   \item{sig_i}{covariance matrix for each person, a \eqn{K \times K \times N} array}
#'   \item{n}{the number of iterations for the EM cycle}
#'   \item{rk}{factor loadings, a \eqn{J \times K} \code{matrix}, for exploratory analysis only}
#'   \item{Q_mat}{factor loading structure, a \eqn{J \times K} matrix}
#'   \item{GIC}{model fit index}
#'   \item{AIC}{model fit index}
#'   \item{BIC}{model fit index}
#'   \item{SE}{Standard errors of item parameters, a \eqn{J \times (K+1)} matrix where the last column includes SE estimates for item difficulty parameters, for confirmatory analysis only}
#'   \item{ur_a}{item discrimination parameters before conducting the rotation, a \eqn{J \times K} \code{matrix}, for exploratory analysis only}
#'   \item{new_a}{item discrimination parameters estimated by IW-GVEM, a \eqn{J \times K} \code{matrix}}
#'   \item{new_b}{item difficulty parameters estimated by IW-GVEM, vector of length \eqn{J}}
#'   \item{new_Sigma_theta}{population variance-covariance matrix estimated by IV-GVEM, a \eqn{K \times K} matrix}
#'   \item{best_lr}{The learning rate used for importance sampling}
#'   \item{best_lb}{The lower bound value for importance sampling}
#' @seealso \code{\link{gvem_2PLCFA}}, \code{\link{gvem_2PLEFA_rot}},\code{\link{bs_2PLCFA}}
#' @export
#'
#' @examples
#' \dontrun{
#' CFA_result <- gvem_2PLCFA(exampleData_2pl, exampleIndic_cfa2pl)
#' importanceSampling(exampleData_2pl,CFA_result)
#' }
importanceSampling<-function(u,gvem_result,S=10,M=10,max.iter=10){
  u=data.matrix(u)
  person<-dim(u)[1]
  indic<-gvem_result$Q_mat
  domain=dim(indic)[2]
  beta_1=0.9
  beta_2=0.99
  theta_IS <- sampling(gvem_result$mu_i,gvem_result$sig_i,person,domain,S,M)
  lr_list <- c(0.5, 0.1, 0.05, 0.01)
  best_lr = 0
  best_lb = -Inf
  best_result=NULL
  for (i in 1:length(lr_list)) {
    lr =lr_list[i]
    new_results<-importance_gradient_descent(u, gvem_result$ra, gvem_result$rb, gvem_result$rsigma,
                                             gvem_result$mu_i, gvem_result$sig_i, indic, lr, S, M,
                                             max.iter, beta_1, beta_2)
    w_tildesub<-importance_weights(theta_IS, gvem_result$mu_i, gvem_result$sig_i, new_results$new_a, new_results$new_b,
                                   new_results$new_Sigma_theta, u, person,dim(u)[2], S, M)$w_tildesub
    #lower bound
    lb_val <- sum(colMeans(log(w_tildesub)))
    if(lb_val>best_lb){
      best_lr = lr
      best_lb = lb_val
      best_result = new_results
    }
  }
  best_result<-c(gvem_result,best_result)
  best_result$best_lr = best_lr
  best_result$best_lb = best_lb
  return(best_result)
}
