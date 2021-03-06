\name{BiCopGofTest}
\alias{BiCopGofTest}

\title{Goodness-of-Fit Test for Bivariate Copulas}

\description{
This function performs a goodness-of-fit test for bivariate copulas,
either based on White's information matrix equality (White 1982) as introduced by Huang and Prokhorov (2011)
or based on Kendall's process.
It computes the test statistics and p-values.
}

\usage{
BiCopGofTest(u1, u2, family, par = 0, par2 = 0, method = "white",
             max.df = 30, B = 100, obj = NULL)
}

\arguments{
  \item{u1, u2}{Numeric vectors of equal length with values in [0,1].}
  \item{family}{An integer defining the bivariate copula family: \cr
		\code{0} = independence copula \cr
	    \code{1} = Gaussian copula \cr
	    \code{2} = Student t copula (t-copula) (only for \code{method = "white"}; see details)\cr
	    \code{3} = Clayton copula \cr
	    \code{4} = Gumbel copula \cr
	    \code{5} = Frank copula \cr
	    \code{6} = Joe copula (only for \code{method = "kendall"}) \cr
		\code{7} = BB1 copula (only for \code{method = "kendall"})\cr
		\code{8} = BB6 copula (only for \code{method = "kendall"})\cr
		\code{9} = BB7 copula (only for \code{method = "kendall"})\cr
		\code{10} = BB8 copula (only for \code{method ="kendall"})\cr
		\code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
		\code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
		\code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr 
		\code{17} = rotated BB1 copula (180 degrees; ``survival BB1''; only for \code{method = "kendall"})\cr
		\code{18} = rotated BB6 copula (180 degrees; ``survival BB6''; only for \code{method = "kendall"})\cr
		\code{19} = rotated BB7 copula (180 degrees; ``survival BB7''; only for \code{method = "kendall"})\cr
		\code{20} = rotated BB8 copula (180 degrees; ``survival BB8''; only for \code{method = "kendall"})\cr
		\code{23} = rotated Clayton copula (90 degrees) \cr
		\code{24} = rotated Gumbel copula (90 degrees) \cr
		\code{26} = rotated Joe copula (90 degrees) \cr
		\code{27} = rotated BB1 copula (90 degrees; only for \code{method = "kendall"}) \cr
		\code{28} = rotated BB6 copula (90 degrees; only for \code{method = "kendall"}) \cr
		\code{29} = rotated BB7 copula (90 degrees; only for \code{method = "kendall"}) \cr
		\code{30} = rotated BB8 copula (90 degrees; only for \code{method = "kendall"}) \cr
		\code{33} = rotated Clayton copula (270 degrees) \cr
		\code{34} = rotated Gumbel copula (270 degrees) \cr
		\code{36} = rotated Joe copula (270 degrees) \cr
		\code{37} = rotated BB1 copula (270 degrees; only for \code{method = "kendall"}) \cr
		\code{38} = rotated BB6 copula (270 degrees; only for \code{method = "kendall"}) \cr
		\code{39} = rotated BB7 copula (270 degrees; only for \code{method = "kendall"}) \cr
		\code{40} = rotated BB8 copula (270 degrees; only for \code{method = "kendall"})
		}
  \item{par}{Copula parameter (optional).}
  \item{par2}{Second parameter for bivariate t-copula (optional); default: \code{par2 = 0}.}
  \item{max.df}{Numeric; upper bound for the estimation of the degrees of freedom parameter of the t-copula (default: \code{max.df = 30}).}
  \item{method}{A string indicating the goodness-of-fit method:\cr
		\code{"white"} = goodness-of-fit test based on White's information matrix equality (default) \cr
		\code{"kendall"} = goodness-of-fit test based on Kendall's process }
  \item{B}{Integer; number of bootstrap samples (default: \code{B = 100}).
    For \code{B = 0} only the the test statistics are returned.\cr
    WARNING: If \code{B} is chosen too large, computations will take very long.} 
  \item{obj}{\code{BiCop} object containing the family and parameter specification.}
}

\value{
For \code{method = "white"}:
	\item{p.value}{Asymptotic p-value.}
	\item{statistic}{The observed test statistic.}\cr
For \code{method = "kendall"}
	\item{p.value.CvM}{Bootstrapped p-value of the goodness-of-fit test using the Cramer-von Mises statistic (if \code{B > 0}).}
  \item{p.value.KS}{Bootstrapped p-value of the goodness-of-fit test using the Kolmogorov-Smirnov statistic (if \code{B > 0}).}
  \item{statistic.CvM}{The observed Cramer-von Mises test statistic.}
  \item{statistic.KS}{The observed Kolmogorov-Smirnov test statistic.}
	
}

\details{
\code{method = "white"}:\cr
This goodness-of fit test uses the information matrix equality of White (1982) and was investigated by Huang and Prokhorov (2011).
The main contribution is that under correct model specification the Fisher Information can be equivalently calculated as minus the expected Hessian matrix or as the expected outer product of the score function.
The null hypothesis is 
\deqn{
	H_0: \boldsymbol{H}(\theta) + \boldsymbol{C}(\theta) = 0 
}
against the alternative
\deqn{
	H_0: \boldsymbol{H}(\theta) + \boldsymbol{C}(\theta) \neq 0 ,
}
where \eqn{\boldsymbol{H}(\theta)} is the expected Hessian matrix and \eqn{\boldsymbol{C}(\theta)} is the expected outer product of the score function.
For the calculation of the test statistic we use the consistent maximum likelihood estimator \eqn{\hat{\theta}} and the sample counter parts of \eqn{\boldsymbol{H}(\theta)} and \eqn{\boldsymbol{C}(\theta)}.
The correction of the covariance-matrix in the test statistic for the uncertainty in the margins is skipped. The implemented tests assumes that where is no uncertainty in the margins.
The correction can be found in Huang and Prokhorov (2011). It involves two-dimensional integrals.\cr
WARNING: For the t-copula the test may be instable. The results for the t-copula therefore have to be treated carefully.\cr
\cr
\code{method = "kendall"}:\cr
This copula goodness-of-fit test is based on Kendall's process as investigated by Genest and Rivest (1993) and Wang and Wells (2000).
For rotated copulas the input arguments are transformed and the goodness-of-fit procedure for the corresponding non-rotated copula is used.
}

\author{Ulf Schepsmeier, Wanling Huang, Jiying Luo, Eike Brechmann}

\references{
Genest, C. and L.-P. Rivest (1993).
Statistical inference procedures for bivariate Archimedean copulas.
Journal of the American Statistical Association, 88 (423), 1034-1043.

Huang, w. and A. Prokhorov (2011).
A goodness-of-fit test for copulas.
to appear in Econometric Reviews

Luo J. (2011).
Stepwise estimation of D-vines with arbitrary specified copula pairs and EDA tools.
Diploma thesis, Technische Universitaet Muenchen.\cr
\url{http://mediatum.ub.tum.de/?id=1079291}.

Wang, W. and M. T. Wells (2000).
Model selection and semiparametric inference for bivariate failure-time data.
Journal of the American Statistical Association, 95 (449), 62-72.

White, H. (1982)
Maximum likelihood estimation of misspecified models,
Econometrica, 50, 1-26.
}

\seealso{\code{\link{BiCopDeriv2}}, \code{\link{BiCopDeriv}}, \code{\link{BiCopIndTest}}, \code{\link{BiCopVuongClarke}}}

\examples{
# simulate from a bivariate Clayton copula
set.seed(123)
simdata <- BiCopSim(300, 3, 2)
u1 <- simdata[,1]
u2 <- simdata[,2]

# perform White's goodness-of-fit test for the true copula
BiCopGofTest(u1, u2, family = 3)

# perform White's goodness-of-fit test for the Frank copula
BiCopGofTest(u1, u2, family = 5)

# perform Kendall's goodness-of-fit test for the true copula
gof <- BiCopGofTest(u1, u2, family = 3, method = "kendall", B=50)
gof$p.value.CvM
gof$p.value.KS

# perform Kendall's goodness-of-fit test for the Frank copula
gof <- BiCopGofTest(u1, u2, family = 5, method = "kendall", B=50)
gof$p.value.CvM
gof$p.value.KS
}

