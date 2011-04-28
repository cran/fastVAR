\name{gather.VARXlasso}
\alias{gather.VARXlasso}
\title{
Gather step for Map-Reduce
}
\description{
Gather the components of the coefficient matrix that were computed in compute nodes
and produce the final coefficient matrix
}
\usage{
gather.VARXlasso(foldername)
}
\arguments{
  \item{foldername}{
The path to the folder that contains the components of the coefficient matrix
}
}
\details{
This function is the third step of a Map-Reduce procedure to calculate the
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
The sparse coefficient matrix returned by VARXlasso
}
\references{

}
\author{
Jeffrey Wong
}
\note{

}

\seealso{
map.VARXlasso, reduce.VARXlasso
}
\examples{

}
