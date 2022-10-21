#' Parallel analysis using polychoric correlation
#' @description Identify the number of factors
#' @param data a \eqn{N \times J} \code{matrix} or a \code{data.frame} that
#' consists of the responses of \eqn{N} individuals to \eqn{J} items
#' without any missing values. The responses are binary or polytomous.
#' @param n.iter Number of simulated analyses to perform
#' @param figure By default, \code{pa_poly} draws an eigenvalue plot.
#' If FALSE, it suppresses the graphic output
#' @return \code{pa_poly} returns a \code{data.frame} with the eigenvalues
#' for the real data and the simulated data.
#' @export
#'
#' @examples
#' pa_poly(exampleData_2pl, n.iter=20)
pa_poly=function(data,n.iter=10,figure=TRUE){
  J<-dim(data)[2]
  N<-dim(data)[1]


  tetrachoric<-matrix(0,nrow=J,ncol=J)

  for (i in 1:J){
    for (j in 1:J){
      if(i==j){tetrachoric[i,j]=1.0}
      if(i<j){
        tetrachoric[i,j]=polycor::polychor(data[,i],data[,j])}
      else{
        tetrachoric[i,j]=tetrachoric[j,i]
      }
    }
  }

  real_eigen<-eigen(tetrachoric)$value


  p_correct<-matrix(0,J,2)
  for (i in 1:J){
    p_correct[i,1]=length(which(data[,i]==0))/N
    p_correct[i,2]=length(which(data[,i]==1))/N
  }
  rep_eigen<-matrix(0,nrow=J, ncol=n.iter)
  for (ss in 1:n.iter){

    new_response<-matrix(0,nrow=N,ncol=J)
    for (i in 1:N){

      for (j in 1:J){
        p_temp=matrix(0,2,1)
        p_temp[1]=p_correct[j,1]
        temp<-runif(1,0,1)
        if(temp<p_temp[1]){
          new_response[i,j]=1 }
      }

    }
    tetrachoric<-matrix(0,nrow=J,ncol=J)

    data<-new_response

    for (i in 1:J){
      for (j in 1:J){
        if(i==j){tetrachoric[i,j]=1.0}
        if(i<j){
          tetrachoric[i,j]=polycor::polychor(data[,i],data[,j])}
        else{
          tetrachoric[i,j]=tetrachoric[j,i]
        }
      }
    }
    rep_eigen[, ss]<-eigen(tetrachoric)$value

  }
  sim_eigen<-apply(rep_eigen,1,mean)
  eigen_parallel<-data.frame(real_eigen,sim_eigen)
  colnames(eigen_parallel)<-c("Actual Data",
                              "Simulated Data")
  com=eigen_parallel[,1]-eigen_parallel[,2]
  n=Position(function(x) x<0, com)-1
  cat("Parallel analysis suggests that ")
  cat("the number of factors = ",n,"\n")
  if(figure == TRUE){
    xrange<-range(1,dim(eigen_parallel)[1])
    yrange <- range(c(eigen_parallel[,1],eigen_parallel[,2]))
    plot(xrange, yrange, type="n", xlab="Factor Number",
         ylab="Eigenvalues" )
    colors <- rainbow(2)
    linetype <- c(1:2)
    plotchar <- seq(18,18+2,1)
    lines(1:dim(eigen_parallel)[1], eigen_parallel[,1], type="b", lwd=1.5,
          lty=linetype[1], col=colors[1], pch=plotchar[1])
    lines(1:dim(eigen_parallel)[1], eigen_parallel[,2], type="b", lwd=1.5,
          lty=linetype[2], col=colors[2], pch=plotchar[2])
    abline(h=1)
    title("Parallel Analysis Scree Plots")
    legend("topright",legend=c("Actual Data", "Simulated Data"),cex=0.8, col=colors,
           pch=plotchar, lty=linetype)
  }
  return(eigen_parallel)
}
