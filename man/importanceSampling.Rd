% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importanceSampling.R
\name{importanceSampling}
\alias{importanceSampling}
\title{Importance Weighted Version of GVEM Analysis for M2PL Models}
\usage{
importanceSampling(u, gvem_result, S = 10, M = 10, max.iter = 10)
}
\arguments{
\item{u}{a \eqn{N \times J} \code{matrix} or a \code{data.frame} that
consists of binary responses of \eqn{N} individuals to \eqn{J} items. The
missing values are coded as \code{NA}}

\item{gvem_result}{a list that includes exploratory or confirmatory GVEM
results for M2PL models.}

\item{S}{the number of times to draw samples;default is 10}

\item{M}{the number of samples drawn from the variational distributions;default is 10}

\item{max.iter}{the maximum number of iterations for the EM cycle; default is 10}
}
\value{
a list containing the following objects:
  \item{ra}{item discrimination parameters estimated by GVEM, a \eqn{J \times K} \code{matrix}}
  \item{rb}{item difficulty parameters estimated by GVEM, vector of length \eqn{J}}
  \item{reta}{variational parameters \eqn{\eta(\xi)}, a \eqn{N \times J} matrix}
  \item{reps}{variational parameters \eqn{\xi}, a \eqn{N \times J} matrix}
  \item{rsigma}{population variance-covariance matrix estimated by GVEM, a \eqn{K \times K} matrix}
  \item{mu_i}{mean parameter for each person, a \eqn{K \times N} matrix}
  \item{sig_i}{covariance matrix for each person, a \eqn{K \times K \times N} array}
  \item{n}{the number of iterations for the EM cycle}
  \item{rk}{factor loadings, a \eqn{J \times K} \code{matrix}, for exploratory analysis only}
  \item{Q_mat}{factor loading structure, a \eqn{J \times K} matrix}
  \item{GIC}{model fit index}
  \item{AIC}{model fit index}
  \item{BIC}{model fit index}
  \item{SE}{Standard errors of item parameters, a \eqn{J \times (K+1)} matrix where the last column includes SE estimates for item difficulty parameters, for confirmatory analysis only}
  \item{ur_a}{item discrimination parameters before conducting the rotation, a \eqn{J \times K} \code{matrix}, for exploratory analysis only}
  \item{new_a}{item discrimination parameters estimated by IW-GVEM, a \eqn{J \times K} \code{matrix}}
  \item{new_b}{item difficulty parameters estimated by IW-GVEM, vector of length \eqn{J}}
  \item{new_Sigma_theta}{population variance-covariance matrix estimated by IV-GVEM, a \eqn{K \times K} matrix}
  \item{best_lr}{The learning rate used for importance sampling}
  \item{best_lb}{The lower bound value for importance sampling}
}
\description{
An importance weighted version of GVEM (i.e., IW-GVEM) can be implemented
to correct the bias on item parameters under M2PL models
}
\examples{
\dontrun{
CFA_result <- gvem_2PLCFA(exampleData_2pl, exampleIndic_cfa2pl)
importanceSampling(exampleData_2pl,CFA_result)
}
}
\seealso{
\code{\link{gvem_2PLCFA}}, \code{\link{gvem_2PLEFA_rot}},\code{\link{bs_2PLCFA}}
}
