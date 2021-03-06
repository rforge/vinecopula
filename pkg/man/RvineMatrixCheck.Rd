\name{RVineMatrixCheck}
\Rdversion{1.1}
\alias{RVineMatrixCheck}

\title{
R-Vine Matrix Check
}
\description{
The given matrix is tested to be a valid R-vine matrix.
}
\usage{RVineMatrixCheck(M)}

\arguments{
  \item{M}{A dxd vine matrix: only lower triangle is used;
For the check, M is assumed to be in natural order, i.e. d:1 on diagonal.
Further M[j+1,j]=d-j and M[j,j]=d-j
}
}

\value{
\item{code}{
	\code{1} for OK; \cr
	\code{-3} diagonal can not be put in order d:1; \cr
	\code{-2} for not permutation of j:d in column d-j; \cr
	\code{-1} if cannot find proper binary array from array in natural order. 
}
}

\note{
The matrix M do not have to be given in natural order or the diagonal in order d:1. The test checks if it can be done in order to be a valid R-vine matrix. \cr
If a function in this package needs the natural order the \code{RVineMatrix} object is automatically "normalized". \cr
The function \code{\link{RVineMatrix}} automatically checks if the given R-vine matrix is valid.
}

\author{Harry Joe}

\references{
Joe H, Cooke RM and Kurowicka D (2011). Regular vines: generation algorithm 
and number of equivalence classes.
In Dependence Modeling: Vine Copula Handbook, pp 219--231.
World Scientific, Singapore.
}

\seealso{\code{\link{RVineMatrix}}}

\examples{
A1 <- matrix(c(6, 0, 0, 0, 0, 0, 
			         5, 5, 0, 0, 0, 0,
			         3, 4, 4, 0, 0, 0,
			         4, 3, 3, 3, 0, 0,
			         1, 1, 2, 2, 2, 0,
			         2, 2, 1, 1, 1, 1), 6, 6, byrow = TRUE)
b1 <- RVineMatrixCheck(A1)
print(b1)
# improper vine matrix, code=-1
A2 <- matrix(c(6, 0, 0, 0, 0, 0, 
			         5, 5, 0, 0, 0, 0,
			         4, 4, 4, 0, 0, 0,
			         1, 3, 3, 3, 0, 0,
			         3, 1, 2, 2, 2, 0,
			         2, 2, 1, 1, 1,1 ), 6, 6, byrow = TRUE)
b2 <- RVineMatrixCheck(A2)
print(b2)
# improper vine matrix, code=-2
A3 <- matrix(c(6, 0, 0, 0, 0, 0, 
			         3, 5, 0, 0, 0, 0,
			         3, 4, 4, 0, 0, 0,
			         4, 3, 3, 3, 0, 0,
			         1, 1, 2, 2, 2, 0,
			         2, 2, 1, 1, 1, 1), 6, 6, byrow = TRUE)
b3 <- RVineMatrixCheck(A3)
print(b3)
}
