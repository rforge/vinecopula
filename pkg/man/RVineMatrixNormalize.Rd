\name{RVineMatrixNormalize}
\Rdversion{1.1}
\alias{RVineMatrixNormalize}

\title{Permute the variables to achieve a natural ordering}
\description{
A RVineMatrix is permuted to achieve a natural ordering (i.e. \code{diag(RVM$Matrix)==d:1})
}

\usage{
RVineMatrixNormalize(RVM)
}
\arguments{
  \item{RVM}{\code{\link{RVineMatrix}} defining the R-vine structure}

}
\value{
  \item{RVM}{A RVineMatrix in natural ordering with entries in \code{RVM$names} keeping track of the reordering}
}

\examples{
corMat <- matrix(c(1.00, 0.17, 0.15, 0.14, 0.13,
                   0.17, 1.00, 0.30, 0.28, 0.05,
                   0.15, 0.30, 1.00, 0.17, 0.05,
                   0.14, 0.28, 0.17, 1.00, 0.04,
                   0.13, 0.05, 0.05, 0.04, 1.00),5,5)

Matrix = matrix(c(5,2,3,1,4,
                  0,2,3,4,1,
                  0,0,3,4,1,
                  0,0,0,4,1,
                  0,0,0,0,1),5,5)
family = matrix(1,5,5)

par = matrix(c(0,0.2,0.9,0.5,0.8,
               0,  0,0.1,0.6,0.9,
               0,  0,  0,0.7,0.5,
               0,  0,  0,  0,0.8,
               0,  0,  0,  0,  0),5,5)

# define RVineMatrix object
RVM = RVineMatrix(Matrix,family,par)

# normalise the RVine
RVineMatrixNormalize(RVM)
}

\keyword{vine}