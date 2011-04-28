\name{map.VARXlasso}
\alias{map.VARXlasso}
\title{
Map step for Map-Reduce
}
\description{
Break down the original input into components that are distributed to computenodes
}
\usage{
map.VARXlasso(y, x, p, b, y.spec = matrix(1, nrow = nrow(y), ncol = nrow(y)), x.spec = matrix(1, nrow = nrow(y), ncol = nrow(x)), computenodes = 4, mapFolder)
}
\arguments{
  \item{y}{
Response matrix where each column is a separate time series
}
  \item{x}{
Matrix of exogenous inputs, where each column is a separate exogenous variable
}
  \item{p}{
The autoregressive lag of the response variable
}
  \item{b}{
The autoregressive lag of the exogenous variables
}
  \item{y.spec}{
A binary matrix that can constrain the number of lagged predictor variables.  
If y.spec[i][j] = 0, the ith time series in y will not be regressed on the jth
time series of y, or any of its lags.
}
  \item{x.spec}{
A binary matrix that can constrain the number of lagged predictor variables.  
If x.spec[i][j] = 0, the ith time series in y will not be regressed on the jth
time series of x, or any of its lags.
}
  \item{computenodes}{
Number of compute nodes in the cluster
}
  \item{mapFolder}{
The path to the folder which contains the independent components of the input data
}
}
\details{
This function is the first step of a Map-Reduce procedure to calculate the
coefficient matrix of the VARXlasso.

The first step maps the reduced response matrix and the design matrix by writing
out individual files representing components of the problem, which are then read
into individual compute nodes.

The second step reduces the problem where each compute node calculates a certain set
of columns of the final coefficient matrix.  These compute nodes receive only components
of the original input data.  The partitions of the coefficient matrix are written out to files.

The final step is done at the head, where the components of the coefficient matrix
are read in and combined to form the final coefficient matrix.
}
\value{

}
\references{

}
\author{
Jeffrey Wong
}
\note{

}

\seealso{
reduce.VARXlasso, gather.VARXlasso
}
\examples{

}
