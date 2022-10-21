\name{VEMIRT-package}
\alias{VEMIRT-package}
\alias{VEMIRT}
\docType{package}
\title{VEMIRT: Variational Expectation Maximization for High-dimensional IRT Models}
\description{
\packageDescription{VEMIRT}
}
\details{
The package includes three modules: parallel analysis, exploratory and confirmatory analysis for M2PL and M3PL models. The number of factors can be identified via parallel analysis using the \code{\link{pa_poly}} function. To conduct the exploratory analysis, the  Gaussian Variational EM (GVEM) algorithms with post-hoc rotation (Promax & CF-Quartimax), Lasso, or adaptive Lasso are provided. The stochastic GVEM is implemented to futher imporve the computational efficacy when analyzing M3PL models. The package also supports the confirmatory analysis by using the  \code{\link{gvem_2PLCFA}} and  \code{\link{sgvem_3PLCFA}} functions.
}
\author{
\itemize{
  \item Jiaying Xiao \email{jxiao6@uw.edu}
  \item Gongjun Xu \email{gongjun@umich.edu}
  \item Chun Wang \email{wang4066@uw.edu}
}

}
\references{
Cho, A. E., Wang, C., Zhang, X., & Xu, G. (2021). Gaussian variational estimation for multidimensional item response theory. \emph{British Journal of Mathematical and Statistical Psychology, 74}, 52-85.

Cho, A. E., Xiao, J., Wang, C., & Xu, G. (2022). Regularized Variational Estimation for Exploratory Item Factor Analysis. \emph{Psychometrika}. https://doi.org/10.1007/s11336-022-09874-6

}

\keyword{ package }

