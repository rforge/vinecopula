\name{RVineStructureSelect}
\alias{RVineStructureSelect}

\title{Sequential Specification of R- and C-Vine Copula Models}

\description{
This function fits either an R- or a C-vine copula model to a d-dimensional copula data set.
Tree structures are determined and appropriate pair-copula families are selected using \code{\link{BiCopSelect}} and estimated sequentially (forward selection of trees).
}

\usage{
RVineStructureSelect(data, familyset = NA, type = 0, selectioncrit = "AIC",
                      indeptest = FALSE, level = 0.05, trunclevel = NA,
                      progress = FALSE, weights = NA, rotations = TRUE)
}

\arguments{
  \item{data}{An N x d data matrix (with uniform margins).}
  \item{familyset}{An integer vector of pair-copula families to select from (the independence copula MUST NOT be specified in this vector unless one wants to fit an independence vine!).
    The vector has to include at least one pair-copula family that allows for positive and one that allows for negative dependence. Not listed copula families might be included to better handle limit cases.
    If \code{familyset = NA} (default), selection among all possible families is performed.
    Coding of pair-copula families: \cr
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
  \item{type}{Type of the vine model to be specified:\cr
    \code{0} or \code{"RVine"} = R-vine (default)\cr
    \code{1} or \code{"CVine"} = C-vine\cr
    C- and D-vine copula models with pre-specified order can be specified using \code{CDVineCopSelect} of the package CDVine.
    Similarly, R-vine copula models with pre-specified tree structure can be specified using \code{\link{RVineCopSelect}}.}
  \item{selectioncrit}{Character indicating the criterion for pair-copula selection. Possible choices: \code{selectioncrit = "AIC"} (default) or \code{"BIC"} (see \code{\link{BiCopSelect}}).}
  \item{indeptest}{Logical; whether a hypothesis test for the independence of \code{u1} and \code{u2} is performed before bivariate copula selection
    (default: \code{indeptest = FALSE}; see \code{\link{BiCopIndTest}}).
    The independence copula is chosen for a (conditional) pair if the null hypothesis of independence cannot be rejected.}
  \item{level}{Numerical; significance level of the independence test (default: \code{level = 0.05}).}
  \item{trunclevel}{Integer; level of truncation.}
  \item{progress}{Logical; whether the tree-wise specification progress is printed (default: \code{progress = FALSE}).}
  \item{weights}{Numerical; weights for each observation (opitional).}
  \item{rotations}{If \code{TRUE}, all rotations of the families in \code{familyset} are included.}
}

\details{
R-vine trees are selected using maximum spanning trees with absolute values of pairwise Kendall's taus as weights, i.e.,
the following optimization problem is solved for each tree:
\deqn{
\max \sum_{edges\ e_{ij}\ in\ spanning\ tree} |\hat{\tau}_{ij}|,
}{
\max \sum_{edges e_{ij} in spanning tree} |\hat{\tau}_{ij}|,
}
where \eqn{\hat{\tau}_{ij}} denote the pairwise empirical Kendall's taus and a spanning tree is a tree on all nodes.
The setting of the first tree selection step is always a complete graph.
For subsequent trees, the setting depends on the R-vine construction principles, in particular on the proximity condition.

The root nodes of C-vine trees are determined similarly by identifying the node with strongest dependencies to all other nodes.
That is we take the node with maximum column sum in the empirical Kendall's tau matrix.

Note that a possible way to determine the order of the nodes in the D-vine is to identify a shortest Hamiltonian path in terms
of weights \eqn{1-|\tau_{ij}|}.
This can be established for example using the package TSP.
Example code is shown below.
}

\value{
  An \code{\link{RVineMatrix}} object with the selected structure (\code{RVM$Matrix}) and families (\code{RVM$family})
  as well as sequentially estimated parameters stored in \code{RVM$par} and \code{RVM$par2}.
}

\references{
Brechmann, E. C., C. Czado, and K. Aas (2012).
Truncated regular vines in high dimensions with applications to financial data.
Canadian Journal of Statistics 40 (1), 68-85.

Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
Selecting and estimating regular vine copulae and application to financial returns.
Computational Statistics & Data Analysis, 59 (1), 52-69.
}

\author{Jeffrey Dissmann, Eike Brechmann, Ulf Schepsmeier, Thomas Nagler}

\seealso{\code{\link{RVineTreePlot}}, \code{\link{RVineCopSelect}}}

\examples{
# load data set
data(daxreturns)

# select the R-vine structure, families and parameters
# using only the first 4 variables and the first 750 observations
# we allow for the copula families: Gauss, t, Clayton, Gumbel, Frank and Joe
RVM <- RVineStructureSelect(daxreturns[1:750,1:4], c(1:6), progress = TRUE)

\dontrun{
# specify a C-vine copula model with only Clayton, Gumbel and Frank copulas (time-consuming)
CVM <- RVineStructureSelect2(daxreturns, c(3,4,5), "CVine")
}

\dontrun{
# determine the order of the nodes in a D-vine using the package TSP (time-consuming)
library(TSP)
d <- dim(daxreturns)[2]
M <- 1 - abs(TauMatrix(daxreturns))
hamilton <- insert_dummy(TSP(M), label = "cut")
sol <- solve_TSP(hamilton, method = "repetitive_nn")
order <- cut_tour(sol, "cut")
DVM <- D2RVine(order, family = rep(0,d*(d-1)/2), par = rep(0, d*(d-1)/2))
RVineCopSelect(daxreturns, c(1:6), DVM$Matrix)
}
}

