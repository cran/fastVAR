\name{reduce.VARXlasso}
\alias{reduce.VARXlasso}
\title{
Reduce step for Map-Reduce
}
\description{
Run this function on individal compute nodes to calculate
partitions of the coefficient matrix
}
\usage{
reduce.VARXlasso(y, x, Z, p, b, y.spec = matrix(1, nrow = nrow(y), ncol = nrow(y)), x.spec = matrix(1, nrow = nrow(y), ncol = nrow(x)), lengthY, lengthX, id, reduceFolder, numcore = 1, ...)
}
\arguments{
  \item{y}{
A partition of the response matrix where each column is a separate time series
}
  \item{x}{
A partition of the matrix of exogenous inputs, where each column is a separate exogenous variable
}
  \item{Z}{
The original design matrix
}
  \item{p}{
The autoregressive lag of the response variable
}
  \item{b}{
The autoregressive lag of the exogenous variables
}
  \item{y.spec}{
The original y.spec matrix.
A binary matrix that can constrain the number of lagged predictor variables.  
If y.spec[i][j] = 0, the ith time series in y will not be regressed on the jth
time series of y, or any of its lags.
}
  \item{x.spec}{
The original x.spec matrix.
A binary matrix that can constrain the number of lagged predictor variables.  
If x.spec[i][j] = 0, the ith time series in y will not be regressed on the jth
time series of x, or any of its lags.
}
  \item{lengthY}{
The number of endogenous time series being measured in the original problem
}
  \item{lengthX}{
The number of exogenous variables being used for prediction in the original problem
}
  \item{id}{
The id number of the partition that this compute node should reduce.  
This number must lie in 1:computenodes (see map.VARXlasso)
}
  \item{reduceFolder}{
The path to the folder where compute nodes can write out their output
}
  \item{numcore}{
The number of cores that can be used to run VARXlasso
}
  \item{\dots}{
Additional arguments to pass to mclapply of the 'multicore' package, which is used in VARXlasso
}
}
\details{
This function is the second step of a Map-Reduce procedure to calculate the
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
map.VARXlasso, gather.VARXlasso
}
\examples{

}
