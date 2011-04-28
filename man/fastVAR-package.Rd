\name{fastVAR-package}
\alias{fastVAR-package}
\alias{fastVAR-package}
\docType{package}
\title{
Compute large VAR and VARX models
}
\description{
This package is designed for time series data.  Uses fast implementations to 
estimate Vector Autoregressive models and Vector Autoregressive models with 
Exogenous Inputs.  For speedup, fastVAR can use multiple cpu cores to
calculate the estimates.  For very large systems, fastVAR uses Lasso penalty to 
return very sparse coefficient matrices.  Regression diagnostics can be used 
to compare models, and prediction functions can be used to calculate the 
n-step ahead prediction.

Map-Reduce functions are in the works (Beta) for estimating large VAR models
on a compute cluster.

Includes Canada data set from package vars to compare results and performance.
}
\details{
\tabular{ll}{
Package: \tab fastVAR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-05-06\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
Depends: glmnet, multicore
}

Very few functions:
\itemize{
  \item fastVAR - estimate a standard VAR model
  \item fastVARX - estimate a standard VARX model with exogenous inputs
  \item VARlasso - use lasso penalty to estimate a VAR model
  \item VARXlasso - use lasso penalty to estimate a VARX model
  \item VARpredict - get the n-step ahead prediction for a VAR model
  \item VARXpredict - get the n-step ahead prediction for a VARX model
  \item VARdiag - get regression diagnostics for a VAR or VARX model
}

Cluster functions:
\itemize{
  \item map.VARXlasso - map a large VARXlasso function to multiple compute nodes
  \item reduce.VARXlasso - compute partitions of the coefficient matrix on compute nodes
  \item gather.VARXlasso - retrieve the partitions and combine them to return the final coefficient matrix
}
}
\author{
Jeffrey Wong <jeffwong@stanford.edu>
}
\references{
Robert Tibshirani <http://www-stat.stanford.edu/~tibs/lasso.html>
}
\keyword{ package }
\seealso{
glmnet, vars
}
\examples{

}
