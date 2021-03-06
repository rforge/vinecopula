\name{RVineCopSelect}       
\alias{RVineCopSelect}

\title{Sequential Pair-Copula Selection and Estimation for R-Vine Copula Models}

\description{
This function fits a R-vine copula model to a d-dimensional copula data set.
Pair-copula families are selected using \code{\link{BiCopSelect}} and estimated sequentially.
}

\usage{
RVineCopSelect(data, familyset = NA, Matrix, selectioncrit = "AIC", indeptest = FALSE,
               level = 0.05, trunclevel = NA, rotations = TRUE)
}

\arguments{
  \item{data}{An N x d data matrix (with uniform margins).}
  \item{familyset}{An integer vector of pair-copula families to select from (the independence copula MUST NOT be specified in this vector unless one wants to fit an independence vine!).
    The vector has to include at least one pair-copula family that allows for positive and one that allows for negative dependence. Not listed copula families might be included to better handle limit cases. 
    If \code{familyset = NA} (default), selection among all possible families is performed.
    The coding of pair-copula families is shown below.}
  \item{Matrix}{Lower or upper triangular d x d matrix that defines the R-vine tree structure.}  
  \item{selectioncrit}{Character indicating the criterion for pair-copula selection. Possible choices: \code{selectioncrit = "AIC"} (default) or \code{"BIC"} (see \code{\link{BiCopSelect}}).}
  \item{indeptest}{Logical; whether a hypothesis test for the independence of \code{u1} and \code{u2} is performed before bivariate copula selection
    (default: \code{indeptest = FALSE}; see \code{\link{BiCopIndTest}}).
    The independence copula is chosen for a (conditional) pair if the null hypothesis of independence cannot be rejected.}
  \item{level}{Numeric; significance level of the independence test (default: \code{level = 0.05}).}
  \item{trunclevel}{Integer; level of truncation.}
  \item{rotations}{If \code{TRUE}, all rotations of the families in \code{familyset} are included.}
}

\details{
R-vine copula models with unknown structure can be specified using \code{\link{RVineStructureSelect}}.
}

\value{
	An \code{\link{RVineMatrix}} object with the following matrix components
	\item{Matrix}{R-vine tree structure matrix as given by the argument \code{Matrix}.}
	\item{family}{Selected pair-copula family matrix with values corresponding to\cr
		\code{0} = independence copula \cr
	        \code{1} = Gaussian copula \cr
	        \code{2} = Student t copula (t-copula) \cr
	        \code{3} = Clayton copula \cr
	        \code{4} = Gumbel copula \cr
	        \code{5} = Frank copula \cr
	        \code{6} = Joe copula \cr 
		\code{7} = BB1 copula \cr
		\code{8} = BB6 copula \cr
		\code{9} = BB7 copula \cr
		\code{10} = BB8 copula \cr
		\code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
		\code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
		\code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr 
		\code{17} = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
		\code{18} = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
		\code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
		\code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
		\code{23} = rotated Clayton copula (90 degrees) \cr
		\code{24} = rotated Gumbel copula (90 degrees) \cr
		\code{26} = rotated Joe copula (90 degrees) \cr
		\code{27} = rotated BB1 copula (90 degrees) \cr
		\code{28} = rotated BB6 copula (90 degrees) \cr
		\code{29} = rotated BB7 copula (90 degrees) \cr
		\code{30} = rotated BB8 copula (90 degrees) \cr
		\code{33} = rotated Clayton copula (270 degrees) \cr
		\code{34} = rotated Gumbel copula (270 degrees) \cr
		\code{36} = rotated Joe copula (270 degrees) \cr
		\code{37} = rotated BB1 copula (270 degrees) \cr
		\code{38} = rotated BB6 copula (270 degrees) \cr
		\code{39} = rotated BB7 copula (270 degrees) \cr
		\code{40} = rotated BB8 copula (270 degrees) \cr
    \code{104} = Tawn type 1 copula \cr
    \code{114} = rotated Tawn type 1 copula (180 degrees) \cr
    \code{124} = rotated Tawn type 1 copula (90 degrees)  \cr
    \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
    \code{204} = Tawn type 2 copula  \cr
    \code{214} = rotated Tawn type 2 copula (180 degrees) \cr
    \code{224} = rotated Tawn type 2 copula (90 degrees)  \cr
    \code{234} = rotated Tawn type 2 copula (270 degrees) \cr
		}
	\item{par}{Estimated pair-copula parameter matrix.}
	\item{par2}{Estimated second pair-copula parameter matrix with parameters of pair-copula families with two parameters.}
}

\author{Eike Brechmann}

\seealso{\code{\link{RVineStructureSelect}}, \code{\link{BiCopSelect}}, \code{\link{RVineSeqEst}}}

\references{
Brechmann, E. C., C. Czado, and K. Aas (2012).
Truncated regular vines in high dimensions with applications to financial data.
Canadian Journal of Statistics 40 (1), 68-85.

Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
Selecting and estimating regular vine copulae and application to financial returns.
Computational Statistics & Data Analysis, 59 (1), 52-69.
}

\examples{
# define 5-dimensional R-vine tree structure matrix
Matrix <- c(5, 2, 3, 1, 4,
            0, 2, 3, 4, 1,
            0, 0, 3, 4, 1,
            0, 0, 0, 4, 1,
            0, 0, 0, 0, 1)
Matrix <- matrix(Matrix, 5, 5)

# define R-vine pair-copula family matrix
family <- c(0, 1, 3, 4, 4,
            0, 0, 3, 4, 1,
            0, 0, 0, 4, 1,
            0, 0, 0, 0, 3,
            0, 0, 0, 0, 0)
family <- matrix(family, 5, 5)

# define R-vine pair-copula parameter matrix
par <- c(0, 0.2, 0.9, 1.5, 3.9,
         0, 0, 1.1, 1.6, 0.9,
         0, 0, 0, 1.9, 0.5,
         0, 0, 0, 0, 4.8,
         0, 0, 0, 0, 0)
par <- matrix(par, 5, 5)

# define second R-vine pair-copula parameter matrix
par2 <- matrix(0, 5, 5)

# define RVineMatrix object
RVM <- RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("V1", "V2", "V3", "V4", "V5"))
                  
# simulate a sample of size 1000 from the R-vine copula model
set.seed(123)
simdata <- RVineSim(1000, RVM)

# determine the pair-copula families and parameters
RVM1 <- RVineCopSelect(simdata, familyset = c(1, 3, 4, 5 ,6), Matrix)
RVM1$family
round(RVM1$par, 2)
}

